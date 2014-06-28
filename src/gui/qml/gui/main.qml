import QtQuick 2.0
import QtGraphicalEffects 1.0
import Qt3D 2.0
import Qt3D.Shapes 2.0
import MolecularDynamics 1.0
import QtQuick.Controls 1.0
import QtQuick.Layouts 1.0
import StereoViewport 1.0
import FileIO 1.0
//import OculusReader 1.0

Rectangle {
    id: rectRoot
    property point lensOffsetFromCenter: Qt.point(0,0)
    property rect distortion: Qt.rect(1, 0.22, 0.24, 0.0)
    property real aspectRatio: width / height;
    property real fillScale: 1.3;
    width: 1280
    height: 800
    color: "black"

    StereoViewport {
        id: viewportRoot
        blending: true
        anchors.fill: parent
        navigation: true

        light: Light {
            ambientColor: Qt.rgba(1,1,1,1)
            position: Qt.vector3d(myCamera.eye.x + 0.3 * (myCamera.center.x - myCamera.eye.x),
                                  myCamera.eye.y + 0.3 * (myCamera.center.y - myCamera.eye.y),
                                  myCamera.eye.z + 0.3 * (myCamera.center.z - myCamera.eye.z))
            quadraticAttenuation: 0.005
        }

        camera: Camera {
            id: myCamera
            eye: Qt.vector3d(23,0,0)
            center: Qt.vector3d(0,0,0)
            nearPlane: 0.1
            farPlane: 5000
            fieldOfView: 60
//            eyeSeparation: 0.1
        }

        MolecularDynamics {
            id: moleculeSystem
            x: -2
            y: -2
            z: -2
            sortPoints: MolecularDynamics.BackToFront
            effect: Effect {
                blending: true
                texture: "sphere3-purple.png"
            }
            targetTemperature: tools.thermostat.targetTemperature
            useThermostat: tools.thermostat.thermostatEnabled
        }

        // TODO Fix bug in MolecularDynamics class and remove this Sphere
//        Sphere {
//            x: 100
//            effect: Effect {
//                color: "blue"
//            }
//        }

        Text {
            anchors.left: parent.left
            anchors.top: parent.top
            text: moleculeSystem.fps.toFixed(1) + " FPS"
            color: "white"
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
            anchors.fill: parent
            onPressed: {
                myCamera.center = Qt.vector3d(0,0,0)
                mouse.accepted = false
            }
        }

//        MouseArea {
//            propagateComposedEvents: false
//            anchors.fill: parent
//            onPressed: {
////                menuRect.state = ""
//                mouse.accepted = false
//            }
//        }
//        PinchArea {
//            property real startingFieldOfView: 90
//            anchors.fill: parent
//            enabled: true
//            onPinchStarted: {
//                startingFieldOfView = myCamera.fieldOfView
//            }

//            onPinchUpdated: {
//                var pinchDiff = pinch.scale;
////                console.log("Scale: " + pinch.scale + " previousScale: " + pinch.previousScale + " diff " + pinchDiff)
//                console.log("Starting " + startingFieldOfView + " scale " + pinch.scale)
//                if(pinchDiff > 1) {
//                    pinchDiff *= 10
//                } else {
//                    pinchDiff = -1/pinchDiff * 10
//                }
//                console.log("pinchDiff " + pinchDiff)
//                myCamera.fieldOfView = startingFieldOfView - pinchDiff
////                myCamera.fieldOfView = startingFieldOfView / (Math.sqrt(Math.sqrt(pinch.scale)))
//            }
//        }
    }

    NavigationPad {
        anchors {
            right: parent.right
            verticalCenter: parent.verticalCenter
            margins: parent.width * 0.01
        }

        width: parent.width * 0.1
        height: parent.width * 0.1
    }

    RotationPad {
        anchors {
            left: parent.left
            verticalCenter: parent.verticalCenter
            margins: parent.width * 0.01
        }

        width: parent.width * 0.1
        height: parent.width * 0.1
    }

    PlaybackControls {

    }

    Tools {
        id: tools
    }

    StatisticsView {
        pressure: moleculeSystem.pressure
        temperature: moleculeSystem.temperature
    }

    //    OculusReader {
    //        id: oculusNavigator
    //        camera: viewportRoot.camera
    //    }
//    FileIO {
//        id: vertexShaderFile
//        source: "oculus.vert"
//        onError: console.log(msg)
//    }

//    FileIO {
//        id: fragmentShaderFile
//        source: "oculus.frag"
//        onError: console.log(msg)
//    }

//    ShaderEffectSource {
//        id: shaderEffectSourceLeft
//        width: rectRoot.width / 2
//        anchors {
//            left: rectRoot.left
//            top: rectRoot.top
//            bottom: rectRoot.bottom
//        }
//        visible: false

//        hideSource: true
//        sourceItem: viewportRoot
//        sourceRect: Qt.rect(0, 0, viewportRoot.width / 2, viewportRoot.height)
//    }

//    ShaderEffectSource {
//        id: shaderEffectSourceRight
//        width: rectRoot.width / 2
//        anchors {
//            right: rectRoot.right
//            top: rectRoot.top
//            bottom: rectRoot.bottom
//        }
//        visible: false

//        hideSource: true
//        sourceItem: viewportRoot
//        sourceRect: Qt.rect(viewportRoot.width / 2, 0, viewportRoot.width / 2, viewportRoot.height)
//    }

//    Item {
//        width: rectRoot.width / 2
//        anchors {
//            left: rectRoot.left
//            top: rectRoot.top
//            bottom: rectRoot.bottom
//        }
//        clip: true
//        ShaderEffect {
//            width: parent.width + 100
//            height: parent.height
//            x: 0

//            property variant qt_Texture0: shaderEffectSourceLeft
//            property point lensOffsetFromCenter: rectRoot.lensOffsetFromCenter
//            property rect distortion: rectRoot.distortion
//            property real aspectRatio: rectRoot.aspectRatio
//            property real fillScale: rectRoot.fillScale
//            vertexShader: vertexShaderFile.read()
//            fragmentShader: fragmentShaderFile.read()
//        }
//    }

//    Item {
//        width: rectRoot.width / 2
//        anchors {
//            right: rectRoot.right
//            top: rectRoot.top
//            bottom: rectRoot.bottom
//        }
//        clip: true
//        ShaderEffect {
//            width: parent.width + 100
//            height: parent.height
//            x: -100

//            property variant qt_Texture0: shaderEffectSourceRight
//            property point lensOffsetFromCenter: Qt.point(-rectRoot.lensOffsetFromCenter.x, rectRoot.lensOffsetFromCenter.y)
//            property rect distortion: rectRoot.distortion
//            property real aspectRatio: rectRoot.aspectRatio
//            property real fillScale: rectRoot.fillScale
//            vertexShader: vertexShaderFile.read()
//            fragmentShader: fragmentShaderFile.read()
//        }
//    }
}
