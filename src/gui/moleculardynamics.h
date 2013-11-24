#ifndef MOLECULARDYNAMICS_H
#define MOLECULARDYNAMICS_H

#include <QObject>
#include <QQuickItem3D>
#include <vector>

using std::vector;

class MoleculeSystem;
class BerendsenThermostat;

class MolecularDynamics : public QQuickItem3D
{
    Q_OBJECT
    Q_PROPERTY(double targetTemperature READ targetTemperature WRITE setTargetTemperature NOTIFY targetTemperatureChanged)
    Q_PROPERTY(bool useThermostat READ useThermostat WRITE setUseThermostat NOTIFY useThermostatChanged)
public:
    explicit MolecularDynamics(QQuickItem *parent = 0);

    void drawItem(QGLPainter *painter);

    Q_INVOKABLE void stepForward();
    double targetTemperature() const;

    bool useThermostat() const;

signals:

    void targetTemperatureChanged(double arg);

    void useThermostatChanged(bool arg);

public slots:

void setTargetTemperature(double arg);

void setUseThermostat(bool arg);

private:
    QList<QVector3D> m_points;
    MoleculeSystem *m_moleculeSystem;
    QGLSceneNode* m_geometry;
    BerendsenThermostat* m_thermostat;
    double m_temperature;
    bool m_useThermostat;

    QArray<QVector3D> vertices;
    QArray<QVector3D> normals;
    QArray<QVector2D> texCoords;
    QArray<uint> indexes;
};

#endif // MOLECULARDYNAMICS_H
