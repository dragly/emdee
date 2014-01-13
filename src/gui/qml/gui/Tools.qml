import QtQuick 2.0
import QtQuick.Layouts 1.0

import "tools"

Rectangle {
    id: toolsRect
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

    TemperatureControl {
        id: temperatureControl
        anchors {
            bottom: parent.top
            bottomMargin: -temperatureControl.height
            left: parent.left
            right: parent.right
        }
        height: parent.height * 2
        states: [
            State {
                name: "toggled"
                PropertyChanges {
                    target: temperatureControl
                    anchors.bottomMargin: 0
                }
            }
        ]
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
            MouseArea {
                anchors.fill: parent
                hoverEnabled: true
                onPressed: {
                    if(toolsRect.state === "toggled") {
                        toolsRect.state = ""
                    } else {
                        toolsRect.state = "toggled"
                    }
                }
            }
            Text {
                anchors.centerIn: parent
                text: "Tools"
            }
        }
        Rectangle {
            Layout.fillWidth: true
            Layout.fillHeight: true
            Text {
                anchors.centerIn: parent
                text: "Select"
            }
        }
        Rectangle {
            Layout.fillWidth: true
            Layout.fillHeight: true
            Text {
                anchors.centerIn: parent
                text: "Create"
            }
        }
        Rectangle {
            Layout.fillWidth: true
            Layout.fillHeight: true
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
            }
        }
        Rectangle {
            Layout.fillWidth: true
            Layout.fillHeight: true
            Text {
                anchors.centerIn: parent
                text: "Force"
            }
        }
        Rectangle {
            Layout.fillWidth: true
            Layout.fillHeight: true
            Text {
                anchors.centerIn: parent
                text: "Delete"
            }
        }
    }
}
