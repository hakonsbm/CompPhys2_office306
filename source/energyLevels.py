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
def plotResultsVSTimestep(data, name):

	energy = data[:,0]
	variance = data[:,1] - data[:,0]**2
	timeStep = data[:,5]


	timeFig = pl.figure()
	pl.title('Energy vs timestep, alpha = %.2f, beta = %.2f'  % (data[0,2] , data[0,3]) + " for " + name)

	pl.plot(timeStep , energy)
	pl.xlabel("Energy")
	pl.ylabel("Timestep")

	timeFig.savefig("../Report/figures/" + name + "TimeEnergy")

	timeFig = pl.figure()
	pl.title('Energy vs timestep, alpha = %.2f, beta = %.2f'  % (data[0,2] , data[0,3]) + " for " + name)

	pl.plot(timeStep , variance)
	pl.xlabel("Variance")
	pl.ylabel("Timestep")

	timeFig.savefig("../Report/figures/" + name + "TimeVariance")




	return

def plotEnergyVsAlphaBeta(data, name):

	alpha = data[:,2]
	beta = data[:,3]
	energy = data[:,0]
	variance = data[:,1] - data[:,0]**2



	fig = pl.figure()
	ax = p3.Axes3D(fig)
	pl.title("Energy vs alpha and beta, " + name)
	ax.plot_trisurf(alpha,beta,energy)
	ax.set_xlabel('alpha')
	ax.set_ylabel('beta')
	ax.set_zlabel('Energy')

	fig.savefig("../Report/figures/" + name + "AlphaBetaEnergy")



	fig=pl.figure()
	ax = p3.Axes3D(fig2)
	pl.title("Variance vs Alpha and beta, " + name)
	ax.plot_trisurf(alpha,beta,variance)
	ax.set_xlabel('alpha')
	ax.set_ylabel('beta')
	ax.set_zlabel('variance')

	fig.savefig("../Report/figures/" + name + "AlphaBetaVariance")


	return

def plotEnergyVsAlpha(data, name):

	alpha = data[:,2]
	energy = data[:,0]
	variance = data[:,1] - data[:,0]**2 

	energyFig = pl.figure()
	pl.title("Energy vs alpha, " + name)

	pl.plot(alpha,energy)

	pl.xlabel("alpha")
	pl.ylabel("Energy")

	energyFig.savefig("../Report/figures/EnergyVsAlpha" + name)

	varianceFig = pl.figure()
	pl.title(" Variance vs alpha, " + name)
	pl.plot (alpha,variance)

	varianceFig.savefig("../Report/figures/VarianceVsAlpha" + name)



#Decide what we want to plot this time

name = "HeliumJastrowAnalytical"

data = np.genfromtxt("outfiles/" + name + "_alpha_beta")
datatime = np.genfromtxt("outfiles/" + name +"_timeStep")



# findLowestEnergy()
plotResultsVSTimestep(datatime , name)
# plotEnergyVsAlphaBeta(data, name)
# plotEnergyVsAlpha(data, name)

pl.show()

