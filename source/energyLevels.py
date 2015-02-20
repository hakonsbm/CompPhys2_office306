import numpy as np
import pylab as pl
import scipy.stats as stats

def LowestEnergy(data):
	
	bestRun = np.argmin(data[:,0])

	alpha = data[bestRun, 2]
	beta = data[bestRun, 3]
	energy = data[bestRun, 0]

	return alpha, beta, energy


HeliumSimpleAnalytical = np.genfromtxt("outfiles/HeliumSimpleAnalytical")
HeliumSimpleNumerical = np.genfromtxt("outfiles/HeliumSimpleNumerical")
HeliumJastrowAnalytical = np.genfromtxt("outfiles/HeliumJastrowAnalytical")
HeliumJastrowNumerical = np.genfromtxt("outfiles/HeliumJastrowNumerical")



alpha, beta, energy = LowestEnergy(HeliumSimpleAnalytical)

print "With trialfunction HeliumSimpleAnalytical"
print "alpha = ", alpha , ", beta = ", beta ,  " and the lowest energy was : " , energy

alpha, beta, energy = LowestEnergy(HeliumSimpleNumerical)

print "With trialfunction HeliumSimpleNumerical"
print "alpha = ", alpha , ", beta = ", beta ,  " and the lowest energy was : " , energy

alpha, beta, energy = LowestEnergy(HeliumJastrowAnalytical)

print "With trialfunction HeliumJastrowAnalytical"
print "alpha = ", alpha , ", beta = ", beta ,  " and the lowest energy was : " , energy

alpha, beta, energy = LowestEnergy(HeliumJastrowNumerical)

print "With trialfunction HeliumJastrowNumerical"
print "alpha = ", alpha , ", beta = ", beta ,  " and the lowest energy was : " , energy
