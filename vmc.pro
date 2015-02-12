TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo

SOURCES += main.cpp \
    vmcsolver.cpp \
    lib.cpp \
    trialFunctions/trialfunction.cpp \
    trialFunctions/heliumsimplenumerically.cpp

HEADERS += \
    vmcsolver.h \
    lib.h \
    trialFunctions/trialfunction.h \
    trialFunctions/heliumsimplenumerically.h
