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
!nomdgui {
    SUBDIRS += gui
    gui.depends = libs
}
mpi:!nomdgui {
    message(Cannot compile gui with mpi enabled)
    SUBDIRS -= gui
}
