TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo

SOURCES += main.cpp \
    vmcsolver.cpp \
    lib.cpp \
    trialFunctions/trialfunction.cpp \
    trialFunctions/heliumsimpleanalytical.cpp \
    trialFunctions/heliumjastrownumerical.cpp \
    trialFunctions/heliumjastrowanalytical.cpp \
    trialFunctions/heliumsimplenumerical.cpp \
    trialFunctions/beryllium.cpp \
    trialFunctions/hydrogen.cpp

HEADERS += \
    vmcsolver.h \
    lib.h \
    trialFunctions/trialfunction.h \
    trialFunctions/heliumsimpleanalytical.h \
    trialFunctions/heliumjastrownumerical.h \
    trialFunctions/heliumjastrowanalytical.h \
    trialFunctions/heliumsimplenumerical.h \
    trialFunctions/beryllium.h \
    trialFunctions/hydrogen.h
