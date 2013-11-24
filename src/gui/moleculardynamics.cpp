#include "moleculardynamics.h"

// Qt includes
#include <QQuickItem>
#include <QQuickItem3D>
#include <QGLBuilder>

// Molecular dynamics includes
#include <moleculesystem.h>
#include <generator.h>
#include <force/lennardjonesforce.h>
#include <force/vashishtatwoparticleforce.h>
#include <force/vashishtathreeparticleforce.h>
#include <integrator/velocityverletintegrator.h>
#include <atom.h>
#include <moleculesystemcell.h>
#include <modifier/berendsenthermostat.h>

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
//    LennardJonesForce *force = new LennardJonesForce();
//    force->setPotentialConstant(potentialConstant);
    VashishtaTwoParticleForce *force = new VashishtaTwoParticleForce();
    VashishtaThreeParticleForce *threeParticleForce = new VashishtaThreeParticleForce();
    m_moleculeSystem->setBoundaries(0,10,0,10,0,10);
    vector<Atom*> atoms;
    int idCounter = 1;
    int type = 0;
    AtomType silicon(0);
    silicon.setEffectiveCharge(1.6);
    silicon.setElectronicPolarizability(0.0);
    silicon.setMass(28.0855);
    AtomType oxygen(1);
    oxygen.setEffectiveCharge(-0.8);
    oxygen.setElectronicPolarizability(2.4);
    oxygen.setMass(15.9994);
    for(int i = 0; i < 5; i++) {
        for(int j = 0; j < 5; j++) {
            type = ((i+j)%2);
            for(int k = 0; k < 5; k++) {
                Atom* atom;
                if(type % 2) {
                    atom = new Atom(silicon);
                } else {
                    atom = new Atom(oxygen);
                }
                atom->setID(idCounter);
                Vector3 position(i + sin(i / 5.0),j + cos(i / 5.0),k);
                position *= 3.0;
                atom->setPosition(position);
                atoms.push_back(atom);
                idCounter++;
                type++;
            }
        }
    }
    m_moleculeSystem->addAtoms(atoms);
    generator.boltzmannDistributeVelocities(1000, atoms);

    VelocityVerletIntegrator *integrator = new VelocityVerletIntegrator(m_moleculeSystem);
    integrator->setTimeStep(0.001);
    m_moleculeSystem->setIntegrator(integrator);
    m_moleculeSystem->addTwoParticleForce(force);
    m_moleculeSystem->addThreeParticleForce(threeParticleForce);
    m_moleculeSystem->setupCells(3);
}

void MolecularDynamics::drawItem(QGLPainter *painter)
{
    if(!painter) {
        return;
    }
    //    qDebug() << "Painting...";
//    QGLBuilder builder;
//    builder.newSection(QGL::NoSmoothing);
    const QMatrix4x4 &modelViewMatrix = painter->modelViewMatrix();
    QVector3D right;
    right.setX(modelViewMatrix(0,0));
    right.setY(modelViewMatrix(0,1));
    right.setZ(modelViewMatrix(0,2));
    QVector3D up;
    up.setX(modelViewMatrix(1,0));
    up.setY(modelViewMatrix(1,1));
    up.setZ(modelViewMatrix(1,2));


    QGLVertexBundle vertexBundle;
    QGLIndexBuffer indexBuffer;

    vertices.clear();
//    normals.clear();
    texCoords.clear();
    indexes.clear();

    vertices.reserve(4*m_points.length());
//    normals.reserve(4*m_points.length());
    texCoords.reserve(4*m_points.length());
    indexes.reserve(6*m_points.length());

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
    QVector3D aOffset =  - right * 0.5 - up * 0.5;
    QVector3D bOffset =  right * 0.5 - up * 0.5;
    QVector3D cOffset =  right * 0.5 + up * 0.5;
    QVector3D dOffset =  - right * 0.5 + up * 0.5;
    QVector2D ta(0,0);
    QVector2D tb(0,1);
    QVector2D tc(1,1);
    QVector2D td(1,0);
    //    for(MoleculeSystemCell* cell : m_moleculeSystem->allCells()) {
    //        for(Atom* atom : cell->atoms()) {
    QVector3D normal(QVector3D::crossProduct(right, up));
    int count = 0;
    for(const QVector3D& centerIn : m_points) {
        QVector3D center = centerIn - QVector3D(5,5,5);
//        if(painter->isCullable(center)) {
//            continue;
//        }
        double size = 0.2;
        a = center + size * aOffset;
        b = center + size * bOffset;
        c = center + size * cOffset;
        d = center + size * dOffset;
        vertices.append(a, b, c, d);
        texCoords.append(ta, tb, tc, td);
//        normals.append(normal, normal, normal, normal);
        indexes.append(count*4 + 0, count*4 + 1, count*4 + 2);
        indexes.append(count*4 + 2, count*4 + 3, count*4 + 0);
        count++;
    }
    vertexBundle.addAttribute(QGL::Position, vertices);
    vertexBundle.addAttribute(QGL::TextureCoord0, texCoords);
//    vertexBundle.addAttribute(QGL::Normal, normals);
    indexBuffer.setIndexes(indexes);

    painter->clearAttributes();
    // Set up normal attributes to use only one element
    painter->glDisableVertexAttribArray(GLuint(QGL::Normal));
    painter->glVertexAttrib3f(GLuint(QGL::Normal), normal.x(), normal.y(), normal.z());

    // Set the rest of the vertex bundle (basically only positions)
    painter->setVertexBundle(vertexBundle);
    painter->draw(QGL::DrawingMode(QGL::Triangles), indexBuffer, 0, indexBuffer.indexCount());
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
