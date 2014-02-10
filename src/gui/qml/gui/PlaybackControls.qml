import QtQuick 2.0
import QtQuick.Layouts 1.0
import QtQuick.Controls 1.0

Rectangle {
    id: menuRect
    //    color: Qt.rgba(0.85, 0.85, 0.95, 1.0)
    color: "black"
    width: parent.width * 0.1
    height: parent.width * 0.1
    anchors {
        left: parent.left
        bottom: parent.bottom
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
    //    RowLayout {
    //        anchors.fill: parent
    //        anchors.margins: menuRect.radius
    //        anchors.leftMargin: 2 * menuRect.radius
    Text {
        anchors.centerIn: parent
        text: timer.running ? "Pause" : "Play"
        color: "white"
        MouseArea {
            anchors.fill: parent

            onClicked: {
                timer.running = !timer.running
            }
        }
    }

    //        Item {
    //            Layout.fillHeight: true
    //        }
    //    }
}
