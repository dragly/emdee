TEMPLATE=subdirs
CONFIG+=ordered
SUBDIRS+=libs

CONFIG += gui

app {
    SUBDIRS += app
}
!mpi:gui {
    SUBDIRS += gui
}


mpi:gui {
    message(Cannot compile gui with mpi enabled)
}
