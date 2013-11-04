#ifndef SUPERCAMERA_H
#define SUPERCAMERA_H

#include <QGLCamera>

class SuperCamera : public QGLCamera
{
    Q_OBJECT
public:
    explicit SuperCamera(QObject *parent = 0);

signals:

public slots:

};

#endif // SUPERCAMERA_H
