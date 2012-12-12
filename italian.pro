#-------------------------------------------------
#
# Project created by QtCreator 2012-10-13T13:38:59
#
#-------------------------------------------------

QT          += core gui
CONFIG      += qwt

TARGET      = italian
TEMPLATE    = app


SOURCES     += main.cpp\
            mainwindow.cpp

HEADERS     += mainwindow.h \
            trimatrix.h

FORMS       += mainwindow.ui

LIBS += -lqwt
INCLUDEPATH += /usr/include/qwt
DEPENDPATH += /usr/include/qwt

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp
