TEMPLATE=subdirs
CONFIG=ordered
SUBDIRS+=libs

CONFIG += gui

app {
    SUBDIRS += app
}
gui {
    SUBDIRS += gui
}
