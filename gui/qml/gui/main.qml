import QtQuick 2.0
import Qt3D 2.0
import Qt3D.Shapes 2.0
import MolecularDynamics 1.0
import QtQuick.Controls 1.0
import QtQuick.Layouts 1.0

Rectangle {
    width: 1280
    height: 800
    color: "black"
    Viewport {
        blending: true
        anchors.fill: parent

        light: Light {
            ambientColor: Qt.rgba(1,1,1,1)
            position.x: myCamera.eye.x / 2
            position.y: myCamera.eye.y / 2
            position.z: myCamera.eye.z / 2
            quadraticAttenuation: 0.005
        }

        camera: Camera {
            id: myCamera
            eye: Qt.vector3d(21,7,19)
            nearPlane: 5
            farPlane: 50

        }

        MolecularDynamics {
            id: moleculeSystem
            effect: Effect {
                blending: true
                texture: "particle.png"
            }
            targetTemperature: targetTemperatureSlider.value
            useThermostat: thermostatCheckBox.checked
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
                menuRect.state = ""
                mouse.accepted = false
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
                enabled: thermostatCheckBox.checked
            }
            Slider {
                id: targetTemperatureSlider
                Layout.fillWidth: true
                minimumValue: 0.0001
                maximumValue: 100
                value: 1.0
                enabled: thermostatCheckBox.checked
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
}
