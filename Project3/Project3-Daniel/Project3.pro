TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    odesolver.cpp \
    vec3.cpp \
    particle.cpp \
    system.cpp

HEADERS += \
    odesolver.h \
    vec3.h \
    particle.h \
    system.h
