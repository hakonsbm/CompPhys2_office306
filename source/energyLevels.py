import numpy as np
import pylab as pl
import scipy.stats as stats
import mpl_toolkits.mplot3d.axes3d as p3


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

def plotEnergyVsAlphaBeta(data):

	alpha = data[:,2]
	beta = data[:,3]
	energy = data[:,0]
	variance = data[:,1] - data[:,0]**2


	fig=pl.figure()
	ax = p3.Axes3D(fig)
	ax.title("  'variance' vs alpha and beta ")
	ax.plot_trisurf(alpha,beta,energy)
	ax.set_xlabel('alpha')
	ax.set_ylabel('beta')
	ax.set_zlabel('Energy')
	pl.show()


	fig=pl.figure()
	ax = p3.Axes3D(fig)
	ax.title("  'variance' vs alpha and beta ")
	ax.plot_trisurf(alpha,beta,variance)
	ax.set_xlabel('alpha')
	ax.set_ylabel('beta')
	ax.set_zlabel('variance')
	pl.show()

	return

def plotEnergyVsAlpha(data):
	pl.figure()
	pl.title("Energy vs alpha")

	alpha = data[:,2]
	energy = data[:,0]
	variance = data[:,1] - data[:,0]**2 



	pl.plot(alpha,energy)

	pl.xlabel("alpha")
	pl.ylabel("Energy")

	pl.figure()
	pl.title(" 'Variance vs alpha'")
	pl.plot (alpha,variance)
	pl.show()



#Decide what we want to plot this time

HeliumSimpleAnalyticalTime = np.genfromtxt("outfiles/HeliumSimpleAnalytical_timeStep")
HeliumJastrowAnalyticalTime = np.genfromtxt("outfiles/HeliumJastrowAnalytical_timeStep")

HeliumSimpleNumericalAlphaBeta = np.genfromtxt("outfiles/HeliumSimpleNumerical_alpha_beta")
HeliumSimpleAnalyticalAlphaBeta = np.genfromtxt("outfiles/HeliumSimpleAnalytical_alpha_beta")

HeliumJastrowAnalyticalAlphaBeta = np.genfromtxt("outfiles/HeliumJastrowAnalytical_alpha_beta")
HeliumJastrowNumericalAlphaBeta = np.genfromtxt("outfiles/HeliumJastrowNumerical_alpha_beta")

# findLowestEnergy()
# plotResultsVSTimestep(HeliumJastrowAnalyticalTime)
# plotEnergyVsAlphaBeta(HeliumJastrowAnalyticalAlphaBeta)
plotEnergyVsAlpha(HeliumSimpleAnalyticalAlphaBeta)



