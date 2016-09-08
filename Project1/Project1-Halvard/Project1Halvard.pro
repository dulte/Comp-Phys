TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    armadillo.cpp

HEADERS += \
    balle.h


LIBS += -llapack -lblas -larmadillo
