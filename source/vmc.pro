TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo -lunittest++

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp
QMAKE_CFLAGS_DEBUG += -fopenmp

SOURCES += main.cpp \
    vmcsolver.cpp \
    trialFunctions/trialfunction.cpp \
    trialFunctions/heliumsimpleanalytical.cpp \
    trialFunctions/heliumjastrownumerical.cpp \
    trialFunctions/heliumjastrowanalytical.cpp \
    trialFunctions/heliumsimplenumerical.cpp \
    trialFunctions/beryllium.cpp \
    trialFunctions/hydrogen.cpp \
    trialFunctions/neon.cpp \
    libm.cpp \
    lib.cpp

HEADERS += \
    vmcsolver.h \
    trialFunctions/trialfunction.h \
    trialFunctions/heliumsimpleanalytical.h \
    trialFunctions/heliumjastrownumerical.h \
    trialFunctions/heliumjastrowanalytical.h \
    trialFunctions/heliumsimplenumerical.h \
    trialFunctions/beryllium.h \
    trialFunctions/hydrogen.h \
    trialFunctions/neon.h \
    libm.h \
    lib.h
