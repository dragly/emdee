import QtQuick 2.0
import QtQuick.Controls 1.0
import QtQuick.Layouts 1.0

Rectangle {
    id: toolPaneRoot
    default property alias content: theContent.data
    property alias title: titleText.text
    width: 200
    height: 100
    color: "black"

    anchors {
        bottom: parent.top
        bottomMargin: -temperatureControl.height
        left: parent.left
        right: parent.right
    }

    Rectangle {
        id: header
        color: "darkblue"
        anchors {
            top: parent.top
            left: parent.left
            right: parent.right
        }
        height: 30

        Label {
            id: titleText
            text: "Tool Name"
            horizontalAlignment: Text.AlignHCenter
            verticalAlignment: Text.AlignVCenter
            anchors.fill: parent
            color: "white"
        }
    }
    Rectangle {
        id: theContent
        color: "black"
        anchors {
            left: parent.left
            right: parent.right
            bottom: parent.bottom
            top: header.bottom
        }
    }

    // Hide show states
    states: [
        State {
            name: "toggled"
            PropertyChanges {
                target: toolPaneRoot
                anchors.bottomMargin: 0
            }
        }
    ]

    transitions: [
        Transition {
            NumberAnimation {
                properties: "anchors.bottomMargin"
                duration: 250
                easing.type: Easing.InOutQuad
            }
        }
    ]
}
