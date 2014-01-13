import QtQuick 2.0
import QtQuick.Controls 1.0
import QtQuick.Layouts 1.0

ToolPane {
    height: layout.height
    title: "Thermostat"

    function kelvinToCelcius(kelvin) {
        return kelvin - 273.15
    }

    ColumnLayout {
        id: layout
        anchors.fill: parent
        anchors.margins: 10
        Label {
            text: "Target temperature:"
            horizontalAlignment: Text.AlignHCenter
            Layout.fillWidth: true
        }
        Label {
            text: kelvinToCelcius(targetTemperatureSlider.value).toFixed(1) + " °C"
            horizontalAlignment: Text.AlignHCenter
            Layout.fillWidth: true
        }
        Slider {
            id: targetTemperatureSlider
            Layout.fillWidth: true
            minimumValue: 0.0001
            maximumValue: 500
            value: 1.0
        }
        Label {
            text: "Enabled:"
            horizontalAlignment: Text.AlignHCenter
            Layout.fillWidth: true
        }
        CheckBox{
            anchors.horizontalCenter: parent.horizontalCenter
        }
    }
}
