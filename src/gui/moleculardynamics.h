#ifndef MOLECULARDYNAMICS_H
#define MOLECULARDYNAMICS_H

#include <QObject>
#include <QQuickItem3D>
#include <QElapsedTimer>
#include <vector>

using std::vector;

class MoleculeSystem;
class BerendsenThermostat;

class MolecularDynamics : public QQuickItem3D
{
    Q_OBJECT
    Q_PROPERTY(double targetTemperature READ targetTemperature WRITE setTargetTemperature NOTIFY targetTemperatureChanged)
    Q_PROPERTY(double temperature READ temperature NOTIFY temperatureChanged)
    Q_PROPERTY(double pressure READ pressure NOTIFY pressureChanged)
    Q_PROPERTY(bool useThermostat READ useThermostat WRITE setUseThermostat NOTIFY useThermostatChanged)
    Q_PROPERTY(SortMode sortPoints READ sortPoints WRITE setSortPoints NOTIFY sortPointsChanged)
    Q_PROPERTY(double fps READ fps NOTIFY fpsChanged)
public:
    explicit MolecularDynamics(QQuickItem *parent = 0);

    void drawItem(QGLPainter *painter);

    Q_INVOKABLE void stepForward();
    double targetTemperature() const;

    bool useThermostat() const;

    double fps() const
    {
        return m_fps;
    }

    SortMode sortPoints() const
    {
        return m_sortPoints;
    }

    double temperature() const
    {
        return m_temperature;
    }

    double pressure() const
    {
        return m_pressure;
    }

signals:

    void targetTemperatureChanged(double arg);

    void useThermostatChanged(bool arg);

    void fpsChanged(double arg);

    void sortPointsChanged(SortMode arg);

    void temperatureChanged(double arg);

    void pressureChanged(double arg);

public slots:

void setTargetTemperature(double arg);

void setUseThermostat(bool arg);

void setSortPoints(SortMode arg)
{
    if (m_sortPoints != arg) {
        m_sortPoints = arg;
        emit sortPointsChanged(arg);
    }
}

private:
    QList<QVector3D> m_points;
    MoleculeSystem *m_moleculeSystem;
    QGLSceneNode* m_geometry;
    BerendsenThermostat* m_thermostat;
    double m_targetTemperature;
    bool m_useThermostat;

    QArray<QVector3D> vertices;
    QArray<QVector3D> normals;
    QArray<QVector2D> texCoords;
    QArray<uint> indexes;
    double m_fps;
    QElapsedTimer fpsTimer;
    SortMode m_sortPoints;
    double m_pressure;
    double m_temperature;
};

#endif // MOLECULARDYNAMICS_H
