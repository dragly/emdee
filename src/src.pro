TEMPLATE=subdirs
CONFIG+=ordered
SUBDIRS+=libs
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += mdgui

app {
    SUBDIRS += app
}
!mpi:mdgui {
    SUBDIRS += gui
}


mpi:mdgui {
    message(Cannot compile gui with mpi enabled)
}
