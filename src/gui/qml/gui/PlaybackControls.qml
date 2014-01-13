import QtQuick 2.0
import QtQuick.Layouts 1.0
import QtQuick.Controls 1.0

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
