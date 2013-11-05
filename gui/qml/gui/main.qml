import QtQuick 2.0
import QtGraphicalEffects 1.0
import Qt3D 2.0
import Qt3D.Shapes 2.0
import MolecularDynamics 1.0
import QtQuick.Controls 1.0
import QtQuick.Layouts 1.0
import StereoViewport 1.0
import FileIO 1.0
import OculusReader 1.0

Rectangle {
    id: rectRoot
    property point lensOffsetFromCenter: Qt.point(0,0)
    property rect distortion: Qt.rect(1, 0.22, 0.24, 0.0)
    property real aspectRatio: width / height;
    property real fillScale: 1.8;
    width: 1280
    height: 800
    color: "black"

    OculusReader {
        id: oculusNavigator
        camera: viewportRoot.camera
    }

    StereoViewport {
        id: viewportRoot
        blending: true
        anchors.fill: parent
//        navigation: false

        light: Light {
            ambientColor: Qt.rgba(1,1,1,1)
            position: myCamera.eye.plus(myCamera.center.minus(myCamera.eye).times(0.2))
            quadraticAttenuation: 0.000005
        }

        camera: Camera {
            id: myCamera
            eye: Qt.vector3d(10,0,0)
            center: Qt.vector3d(0,0,0)
            nearPlane: 0.1
            farPlane: 50
            fieldOfView: 120
            onFieldOfViewChanged: {
                console.log(fieldOfView)
                if(fieldOfView > 160) {
                    fieldOfView = 160
                } else if(fieldOfView < 20) {
                    fieldOfView = 20
                }
            }
            eyeSeparation: 0.1
        }

        MolecularDynamics {
            id: moleculeSystem
            x: 2
            y: 2
            z: 2
            effect: Effect {
                blending: true
                texture: "particle.png"
            }
            targetTemperature: targetTemperatureSlider.value
            useThermostat: thermostatCheckBox.checked
        }

        // TODO Fix bug in MolecularDynamics class and remove this Sphere
        Sphere {
            x: 100
            effect: Effect {
                color: "blue"
            }
        }

        Timer {
            id: timer
            interval: 50
            repeat: true
            onTriggered: {
                moleculeSystem.stepForward()
            }
        }
        MouseArea {
            propagateComposedEvents: false
            anchors.fill: parent
            onPressed: {
                oculusNavigator.enabled = false
                menuRect.state = ""
                mouse.accepted = false
            }
            onReleased: {
                oculusNavigator.enabled = true
            }
        }
        PinchArea {
            property real startingFieldOfView: 90
            anchors.fill: parent
            enabled: true
            onPinchStarted: {
                startingFieldOfView = myCamera.fieldOfView
            }

            onPinchUpdated: {
                var pinchDiff = pinch.scale;
//                console.log("Scale: " + pinch.scale + " previousScale: " + pinch.previousScale + " diff " + pinchDiff)
                console.log("Starting " + startingFieldOfView + " scale " + pinch.scale)
                if(pinchDiff > 1) {
                    pinchDiff *= 10
                } else {
                    pinchDiff = -1/pinchDiff * 10
                }
                console.log("pinchDiff " + pinchDiff)
                myCamera.fieldOfView = startingFieldOfView - pinchDiff
//                myCamera.fieldOfView = startingFieldOfView / (Math.sqrt(Math.sqrt(pinch.scale)))
            }
        }
    }
    Rectangle {
        id: menuRect
        color: Qt.rgba(0.85, 0.85, 0.95, 1.0)
        width: parent.width * 0.2
        radius: parent.width * 0.01
        anchors {
            left: parent.left
            top: parent.top
            bottom: parent.bottom
            margins: width * 0.05
            leftMargin: -(width - radius)
        }
        Behavior on anchors.leftMargin {
            NumberAnimation { duration: 400; easing.type: Easing.InOutCubic }
        }

        states: [
            State {
                name: "hovered"
                PropertyChanges {
                    target: menuRect
                    anchors.leftMargin: -menuRect.radius
                }
            }

        ]
        MouseArea {
            anchors.fill: parent
            hoverEnabled: true
            onPressed: {
                menuRect.state = "hovered"
            }
        }
        ColumnLayout {
            anchors.fill: parent
            anchors.margins: menuRect.radius
            anchors.leftMargin: 2 * menuRect.radius
            Label {
                text: "Thermostat:"
            }
            CheckBox {
                id: thermostatCheckBox
                checked: false
                text: "Enabled"
            }
            Label {
                text: "Target temperature:"
//                enabled: thermostatCheckBox.checked
            }
            Slider {
                id: targetTemperatureSlider
                Layout.fillWidth: true
                minimumValue: 0.0001
                maximumValue: 500
                value: 1.0
//                enabled: thermostatCheckBox.checked
            }
            Label {
                text: "Simulation:"
            }
            Button {
                text: timer.running ? "Pause" : "Play"
                onClicked: {
                    timer.running = !timer.running
                }
            }

            Item {
                Layout.fillHeight: true
            }
        }
    }
    FileIO {
        id: vertexShaderFile
        source: "oculus.vert"
        onError: console.log(msg)
    }

    FileIO {
        id: fragmentShaderFile
        source: "oculus.frag"
        onError: console.log(msg)
    }

    ShaderEffectSource {
        id: shaderEffectSourceLeft
        width: rectRoot.width / 2
        anchors {
            left: rectRoot.left
            top: rectRoot.top
            bottom: rectRoot.bottom
        }
        visible: false

        hideSource: true
        sourceItem: viewportRoot
        sourceRect: Qt.rect(0, 0, viewportRoot.width / 2, viewportRoot.height)
    }

    ShaderEffectSource {
        id: shaderEffectSourceRight
        width: rectRoot.width / 2
        anchors {
            right: rectRoot.right
            top: rectRoot.top
            bottom: rectRoot.bottom
        }
        visible: false

        hideSource: true
        sourceItem: viewportRoot
        sourceRect: Qt.rect(viewportRoot.width / 2, 0, viewportRoot.width / 2, viewportRoot.height)
    }

    Item {
        width: rectRoot.width / 2
        anchors {
            left: rectRoot.left
            top: rectRoot.top
            bottom: rectRoot.bottom
        }
        clip: true
        ShaderEffect {
            width: parent.width + 100
            height: parent.height
            x: 0

            property variant qt_Texture0: shaderEffectSourceLeft
            property point lensOffsetFromCenter: rectRoot.lensOffsetFromCenter
            property rect distortion: rectRoot.distortion
            property real aspectRatio: rectRoot.aspectRatio
            property real fillScale: rectRoot.fillScale
            vertexShader: vertexShaderFile.read()
            fragmentShader: fragmentShaderFile.read()
        }
    }

    Item {
        width: rectRoot.width / 2
        anchors {
            right: rectRoot.right
            top: rectRoot.top
            bottom: rectRoot.bottom
        }
        clip: true
        ShaderEffect {
            width: parent.width + 100
            height: parent.height
            x: -100

            property variant qt_Texture0: shaderEffectSourceRight
            property point lensOffsetFromCenter: Qt.point(-rectRoot.lensOffsetFromCenter.x, rectRoot.lensOffsetFromCenter.y)
            property rect distortion: rectRoot.distortion
            property real aspectRatio: rectRoot.aspectRatio
            property real fillScale: rectRoot.fillScale
            vertexShader: vertexShaderFile.read()
            fragmentShader: fragmentShaderFile.read()
        }
    }
    NavigationPad {
        anchors {
            right: parent.right
            bottom: parent.bottom
            margins: parent.width * 0.01
        }

        width: parent.width * 0.1
        height: parent.width * 0.1
    }
}
