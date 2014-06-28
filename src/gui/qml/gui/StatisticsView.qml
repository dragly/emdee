import QtQuick 2.0

Rectangle {
    property double temperature
    property double pressure
    color: "black"
    width: parent.width * 0.2
    height: parent.width * 0.025
    anchors {
        right: parent.right
        top: parent.top
    }

    function kelvinToCelcius(kelvin) {
        return kelvin - 273.15
    }

    function temperatureToKelvin(temperature) {
        return temperature*119.57
    }

    function temperatureToCelsius(temperature) {
        return kelvinToCelcius(temperatureToKelvin(temperature))
    }

    function pressureToPascal(pressure) {
        return pressure * 41818086.5802
    }

    Text {
        anchors {
            left: parent.left
        }
        width: parent.width / 2
        height: parent.height
        horizontalAlignment: Text.AlignHCenter
        verticalAlignment: Text.AlignVCenter

        text: temperatureToCelsius(temperature).toFixed(1) + " °C"
        color: "white"
    }

//    Text {
//        anchors {
//            right: parent.right
//        }
//        width: parent.width / 2
//        height: parent.height
//        horizontalAlignment: Text.AlignHCenter
//        verticalAlignment: Text.AlignVCenter

//        text: (pressureToPascal(pressure)/1000000 * 3.3).toFixed(1) + " kPa"
//        color: "white"
//    }
}
