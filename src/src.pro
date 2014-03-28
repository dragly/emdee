TEMPLATE=subdirs
CONFIG+=ordered
SUBDIRS+=libs
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

!noapp {
    SUBDIRS += app
}
!nomdgui {
    SUBDIRS += gui
}
mpi:!nomdgui {
    message(Cannot compile gui with mpi enabled)
    SUBDIRS -= gui
}
