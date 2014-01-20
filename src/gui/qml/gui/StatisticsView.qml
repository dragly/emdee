import QtQuick 2.0

Rectangle {
    property double temperature
    property double pressure
    width: parent.width * 0.2
    height: parent.width * 0.025
    anchors {
        right: parent.right
        top: parent.top
    }

    function kelvinToCelcius(kelvin) {
        return kelvin - 273.15
    }

    Text {
        anchors {
            left: parent.left
        }
        width: parent.width / 2
        height: parent.height
        horizontalAlignment: Text.AlignHCenter
        verticalAlignment: Text.AlignVCenter

        text: kelvinToCelcius(temperature).toFixed(1) + " °C(a.u.)"
    }

    Text {
        anchors {
            right: parent.right
        }
        width: parent.width / 2
        height: parent.height
        horizontalAlignment: Text.AlignHCenter
        verticalAlignment: Text.AlignVCenter

        text: pressure.toFixed(1) + " Pa(a.u.)"
    }
}
