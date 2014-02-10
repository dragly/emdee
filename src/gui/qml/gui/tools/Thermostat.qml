import QtQuick 2.0
import QtQuick.Controls 1.0
import QtQuick.Layouts 1.0

ToolPane {
    property alias targetTemperature: targetTemperatureSlider.value
    property alias thermostatEnabled: thermostatEnabledCheckbox.checked
    height: layout.height
    title: "Thermostat"

    function kelvinToCelcius(kelvin) {
        return kelvin - 273.15
    }

    function temperatureToKelvin(temperature) {
        return temperature*119.57
    }

    function temperatureToCelsius(temperature) {
        return kelvinToCelcius(temperatureToKelvin(temperature))
    }

    ColumnLayout {
        id: layout
        anchors.fill: parent
        anchors.margins: 10
        Label {
            text: "Target temperature:"
            horizontalAlignment: Text.AlignHCenter
            Layout.fillWidth: true
            color: "white"
        }
        Label {
            text: temperatureToCelsius(targetTemperatureSlider.value).toFixed(1) + " °C"
            horizontalAlignment: Text.AlignHCenter
            Layout.fillWidth: true
            color: "white"
        }
        Slider {
            id: targetTemperatureSlider
            Layout.fillWidth: true
            minimumValue: 0.0001
            maximumValue: 5.0
            value: 1.0
        }
        Label {
            text: "Enabled:"
            horizontalAlignment: Text.AlignHCenter
            Layout.fillWidth: true
            color: "white"
        }
        CheckBox{
            id: thermostatEnabledCheckbox
            anchors.horizontalCenter: parent.horizontalCenter
        }
    }
}
