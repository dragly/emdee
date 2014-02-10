import QtQuick 2.0
import QtQuick.Layouts 1.0

import "tools"

Rectangle {
    id: toolsRect
    property alias thermostat: temperatureControl
    color: "black"
    anchors {
        right: parent.right
        bottom: parent.bottom
        rightMargin: -parent.width * 0.05 * 2
        bottomMargin: -parent.width * 0.05 * 1
    }
    width: parent.width * 0.05 * 3
    height: parent.width * 0.05 * 2

    states: [
        State {
            name: "toggled"
            PropertyChanges {
                target: toolsRect
                anchors.rightMargin: 0
                anchors.bottomMargin: 0
            }
        }

    ]

    transitions: [
        Transition {
            from: ""
            to: "toggled"
            SequentialAnimation {
                NumberAnimation {
                    properties: "anchors.bottomMargin"
                    duration: 200
                    easing.type: Easing.InQuad
                }
                NumberAnimation {
                    properties: "anchors.rightMargin"
                    duration: 400
                    easing.type: Easing.OutQuad
                }
            }
        },

        Transition {
            from: "toggled"
            to: ""
            SequentialAnimation {
                NumberAnimation {
                    properties: "anchors.rightMargin"
                    duration: 400
                    easing.type: Easing.InQuad
                }
                NumberAnimation {
                    properties: "anchors.bottomMargin"
                    duration: 200
                    easing.type: Easing.OutQuad
                }
            }
        }
    ]

    Thermostat {
        id: temperatureControl
        height: parent.height * 2
    }

    GridLayout {
        anchors.fill: parent
        columns: 3
        rows: 2
        rowSpacing: 0
        columnSpacing: 0
        Rectangle {
            Layout.fillWidth: true
            Layout.fillHeight: true
            color: "black"
            MouseArea {
                anchors.fill: parent
                hoverEnabled: true
                onPressed: {
                    if(toolsRect.state === "toggled") {
                        toolsRect.state = ""
                        temperatureControl.state = ""
                    } else {
                        toolsRect.state = "toggled"
                    }
                }
            }
            Text {
                anchors.centerIn: parent
                text: "Tools"
                color: "white"
            }
        }
        Rectangle {
            Layout.fillWidth: true
            Layout.fillHeight: true
            color: "black"
            Text {
                anchors.centerIn: parent
                text: "Select"
                color: "white"
            }
        }
        Rectangle {
            Layout.fillWidth: true
            Layout.fillHeight: true
            color: "black"
            Text {
                anchors.centerIn: parent
                text: "Create"
                color: "white"
            }
        }
        Rectangle {
            Layout.fillWidth: true
            Layout.fillHeight: true
            color: "black"
            MouseArea {
                anchors.fill: parent
                onPressed: {
                    if(temperatureControl.state === "") {
                        temperatureControl.state = "toggled"
                    } else {
                        temperatureControl.state = ""
                    }
                }
            }

            Text {
                anchors.centerIn: parent
                text: "Temp"
                color: "white"
            }
        }
        Rectangle {
            Layout.fillWidth: true
            Layout.fillHeight: true
            color: "black"
            Text {
                anchors.centerIn: parent
                text: "Force"
                color: "white"
            }
        }
        Rectangle {
            Layout.fillWidth: true
            Layout.fillHeight: true
            color: "black"
            Text {
                anchors.centerIn: parent
                text: "Delete"
                color: "white"
            }
        }
    }
}
