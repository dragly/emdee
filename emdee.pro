TEMPLATE=subdirs
SUBDIRS+=src \
    tools
CONFIG+=ordered

tools.depends = src

!notests {
    SUBDIRS += tests
    tests.depends = src
}
