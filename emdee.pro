TEMPLATE=subdirs
SUBDIRS+=src \
    tools

!notests {
    SUBDIRS += tests
}
