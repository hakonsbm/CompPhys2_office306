import numpy as np
import pylab as pl
import scipy.stats as stats


def findLowestEnergy():

	def LowestEnergy(data):
		bestRun = np.argmin(data[:,0])

		alpha = data[bestRun, 2]
		beta = data[bestRun, 3]
		energy = data[bestRun, 0]

		return alpha, beta, energy


	HeliumSimpleAnalytical = np.genfromtxt("outfiles/HeliumSimpleAnalytical_alpha_beta")
	alpha, beta, energy = LowestEnergy(HeliumSimpleAnalytical)
	print "With trialfunction HeliumSimpleAnalytical"
	print "alpha = ", alpha , ", beta = ", beta ,  " and the lowest energy was : " , energy

	# HeliumSimpleNumerical = np.genfromtxt("outfiles/HeliumSimpleNumerical_alpha_beta")
	# alpha, beta, energy = LowestEnergy(HeliumSimpleNumerical)
	# print "With trialfunction HeliumSimpleNumerical"
	# print "alpha = ", alpha , ", beta = ", beta ,  " and the lowest energy was : " , energy

	HeliumJastrowAnalytical = np.genfromtxt("outfiles/HeliumJastrowAnalytical_alpha_beta")
	alpha, beta, energy = LowestEnergy(HeliumJastrowAnalytical)
	print "With trialfunction HeliumJastrowAnalytical"
	print "alpha = ", alpha , ", beta = ", beta ,  " and the lowest energy was : " , energy

	# HeliumJastrowNumerical = np.genfromtxt("outfiles/HeliumJastrowNumerical_alpha_beta")
	# alpha, beta, energy = LowestEnergy(HeliumJastrowNumerical)
	# print "With trialfunction HeliumJastrowNumerical"
	# print "alpha = ", alpha , ", beta = ", beta ,  " and the lowest energy was : " , energy


#Making the plots timesteps and energy
def plotResultsVSTimestep(data):

	# print data
	pl.figure()
	pl.title('Energy vs timestep, alpha = %.2f, beta = %.2f' % (data[0,2] , data[0,3]))

	pl.plot(data[:,5] , data[:,0])
	pl.show()

	return


#Decide what we want to plot this time

HeliumSimpleAnalyticalTime = np.genfromtxt("outfiles/HeliumSimpleAnalytical_timeStep")
HeliumJastrowAnalyticalTime = np.genfromtxt("outfiles/HeliumJastrowAnalytical_timeStep")



# findLowestEnergy()
plotResultsVSTimestep(HeliumJastrowAnalyticalTime)



