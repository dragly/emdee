TEMPLATE=subdirs
SUBDIRS+=src
CONFIG+=ordered

!notools {
    SUBDIRS += tools
    tools.depends = src
}
!notests {
    SUBDIRS += tests
    tests.depends = src
}
!noexamples {
    SUBDIRS += examples
    examples.depends = src
}
