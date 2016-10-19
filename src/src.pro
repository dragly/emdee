TEMPLATE=subdirs
CONFIG+=ordered
SUBDIRS+=libs
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

!noapp {
    SUBDIRS += app
    app.depends = libs
}
#!nogui {
#    SUBDIRS += gui
#    gui.depends = libs
#}
mpi:!nogui {
    message(Cannot compile gui with mpi enabled)
    SUBDIRS -= gui
}
