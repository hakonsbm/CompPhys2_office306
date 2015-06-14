TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo -lunittest++

#QMAKE_CXXFLAGS+= -fopenmp
#QMAKE_LFLAGS +=  -fopenmp
#QMAKE_CFLAGS_DEBUG += -fopenmp

#MPI settings
# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

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
#    GTO/basisbank.cpp \
    lib.cpp \
    tests.cpp \
    derivatives.cpp \
    slaterdeterminant.cpp \
    trialFunctions/helium.cpp \
    trialFunctions/hydrogentwo.cpp \
    GTO/basisbank.cpp \
    GTO/contracted.cpp \
    GTO/gto.cpp \
    trialFunctions/berylliumtwo.cpp

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
#    GTO/basisbank.h \
    lib.h \
    derivatives.h \
    slaterdeterminant.h \
    trialFunctions/helium.h \
    trialFunctions/hydrogentwo.h \
    GTO/basisbank.h \
    GTO/contracted.h \
    GTO/gto.h \
    trialFunctions/berylliumtwo.h
