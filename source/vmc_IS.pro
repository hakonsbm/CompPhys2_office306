TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo

SOURCES += \
    lib.cpp \
    trialFunctions/trialfunction.cpp \
    trialFunctions/heliumsimpleanalytical.cpp \
    trialFunctions/heliumjastrownumerical.cpp \
    trialFunctions/heliumjastrowanalytical.cpp \
    trialFunctions/heliumsimplenumerical.cpp \
    vmcsolver_IS.cpp \
    main_IS.cpp

HEADERS += \
    lib.h \
    trialFunctions/trialfunction.h \
    trialFunctions/heliumsimpleanalytical.h \
    trialFunctions/heliumjastrownumerical.h \
    trialFunctions/heliumjastrowanalytical.h \
    trialFunctions/heliumsimplenumerical.h \
    vmcsolver_IS.h
