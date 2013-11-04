#include "moleculardynamics.h"

// Qt includes
#include <QQuickItem>
#include <QQuickItem3D>
#include <QGLBuilder>

// Molecular dynamics includes
#include <src/moleculesystem.h>
#include <src/generator.h>
#include <src/force/lennardjonesforce.h>
#include <src/integrator/velocityverletintegrator.h>
#include <src/atom.h>
#include <src/moleculesystemcell.h>
#include <src/modifier/berendsenthermostat.h>

MolecularDynamics::MolecularDynamics(QQuickItem *parent) :
    QQuickItem3D(parent),
    m_temperature(1.0),
    m_useThermostat(false),
    m_thermostat(0)
{
    m_moleculeSystem = new MoleculeSystem();
    m_moleculeSystem->setOutputEnabled(false);
    m_moleculeSystem->setSaveEnabled(false);

    double potentialConstant = 1;
    double bUnit = 5.620 / 3.405;

    Generator generator;
    //    generator.setUnitLength(unitLength);
    LennardJonesForce *force = new LennardJonesForce();
    force->setPotentialConstant(potentialConstant);
    vector<Atom*> atoms = generator.generateFcc(bUnit, 4, AtomType::argon());
    generator.boltzmannDistributeVelocities(3, atoms);

    VelocityVerletIntegrator *integrator = new VelocityVerletIntegrator(m_moleculeSystem);
    integrator->setTimeStep(0.001);
    m_moleculeSystem->setIntegrator(integrator);
    m_moleculeSystem->addTwoParticleForce(force);
    //    system.setPotentialConstant(potentialConstant);
    mat lastBoundaries = generator.lastBoundaries();
    lastBoundaries(1,0) += 5;
    lastBoundaries(1,1) += 5;
    lastBoundaries(1,2) += 5;
    m_moleculeSystem->setBoundaries(lastBoundaries);
    m_moleculeSystem->addAtoms(atoms);
    m_moleculeSystem->setupCells(potentialConstant * 3);
}

void MolecularDynamics::drawItem(QGLPainter *painter)
{
    if(!painter) {
        return;
    }
    //    qDebug() << "Painting...";
    QGLBuilder builder;
    builder.newSection(QGL::NoSmoothing);
    const QMatrix4x4 &modelViewMatrix = painter->modelViewMatrix();
    QVector3D right;
    right.setX(modelViewMatrix(0,0));
    right.setY(modelViewMatrix(0,1));
    right.setZ(modelViewMatrix(0,2));
    QVector3D up;
    up.setX(modelViewMatrix(1,0));
    up.setY(modelViewMatrix(1,1));
    up.setZ(modelViewMatrix(1,2));
    QGeometryData quad;

    //    if(m_sortPoints == BackToFront) {
    QMultiMap<double, QVector3D> sortedPoints;
    for(MoleculeSystemCell* cell : m_moleculeSystem->allCells()) {
        for(Atom* atom : cell->atoms()) {
            QVector3D center = QVector3D(atom->position().x(), atom->position().y(), atom->position().z());
            const QVector4D &depthVector = painter->modelViewMatrix() * center;
            double depth = depthVector.z();
            sortedPoints.insert(depth, center);
        }
    }
    m_points.clear();
    QMapIterator<double, QVector3D> i(sortedPoints);
    while(i.hasNext()) {
        m_points.push_back(i.next().value());
    }
    sortedPoints.clear();
    //    }

    QVector3D a;
    QVector3D b;
    QVector3D c;
    QVector3D d;
    QVector2D ta(0,0);
    QVector2D tb(0,1);
    QVector2D tc(1,1);
    QVector2D td(1,0);
    //    for(MoleculeSystemCell* cell : m_moleculeSystem->allCells()) {
    //        for(Atom* atom : cell->atoms()) {
    for(const QVector3D& centerIn : m_points) {
        QVector3D center = centerIn - QVector3D(5,5,5);
        if(painter->isCullable(center)) {
            continue;
        }
        double size = 0.2;
        a = center - right * (size * 0.5);
        b = center + right * size * 0.5;
        c = center + right * size * 0.5 + up * size;
        d = center - right * size * 0.5 + up * size;
        quad.appendVertex(a,b,c,d);
        quad.appendTexCoord(ta, tb, tc, td);
    }
    //        }
    //    }

    // }
    builder.addQuads(quad);
    QGLSceneNode* geometry = builder.finalizedSceneNode();
    if(m_geometry) {
        delete m_geometry;
    }
    m_geometry = geometry;
    if(m_geometry) {
        m_geometry->draw(painter);
    }
}

void MolecularDynamics::stepForward()
{
    m_moleculeSystem->setNSimulationSteps(2);
    m_moleculeSystem->simulate();
    update();
}

double MolecularDynamics::targetTemperature() const
{
    return m_temperature;
}

bool MolecularDynamics::useThermostat() const
{
    return m_useThermostat;
}

void MolecularDynamics::setTargetTemperature(double arg)
{
    if (m_temperature != arg) {
        if(m_thermostat) {
            m_thermostat->setTargetTemperature(arg);
        }
        m_temperature = arg;
        emit targetTemperatureChanged(arg);
    }
}

void MolecularDynamics::setUseThermostat(bool arg)
{
    if (m_useThermostat != arg) {
        m_useThermostat = arg;
        if(arg) {
            m_thermostat = new BerendsenThermostat(m_moleculeSystem);
            m_thermostat->setTargetTemperature(m_temperature);
            m_thermostat->setRelaxationTime(0.01);
            m_moleculeSystem->addModifier(m_thermostat);
        } else {
            m_moleculeSystem->removeModifier(m_thermostat);
            delete m_thermostat;
        }
        emit useThermostatChanged(arg);
    }
}
