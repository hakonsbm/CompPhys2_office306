import numpy as np
import pylab as pl
import scipy.stats as stats
import mpl_toolkits.mplot3d.axes3d as p3


def findLowestEnergy(data, name):

	def LowestEnergy(data):
		bestRun = np.argmin(data[:,0])

		alpha = data[bestRun, 3]
		beta = data[bestRun, 4]
		energy = data[bestRun, 0]
		variance = data[bestRun, 2] 

		return alpha, beta, energy, variance


	alpha, beta, energy, variance = LowestEnergy(data)
	print "With trialfunction " + name
	print "alpha = ", alpha , ", beta = ", beta ,  ", the lowest energy was : " , energy, " with sample variance " , variance




#Making the plots timesteps and energy
def plotResultsVSTimestep(data, name):

	energy = data[:,0]
	variance = data[:,1] - data[:,0]**2
	timeStep = data[:,6]


	timeFig = pl.figure()
	pl.title('Energy vs timestep, alpha = %.2f, beta = %.2f'  % (data[0,3] , data[0,4]) + " for " + name)

	pl.plot(timeStep , energy)
	pl.xlabel("Energy")
	pl.ylabel("Timestep")

	timeFig.savefig("../Report/figures/" + name + "TimeEnergy")

	timeFig = pl.figure()
	pl.title('Energy vs timestep, alpha = %.2f, beta = %.2f'  % (data[0,3] , data[0,4]) + " for " + name)

	pl.plot(timeStep , variance)
	pl.xlabel("Variance")
	pl.ylabel("Timestep")

	timeFig.savefig("../Report/figures/" + name + "TimeVariance")




	return

def plotEnergyVsAlphaBeta(data, name):

	alpha = data[:,3]
	beta = data[:,4]
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
	ax = p3.Axes3D(fig)
	pl.title("Variance vs Alpha and beta, " + name)
	ax.plot_trisurf(alpha,beta,variance)
	ax.set_xlabel('alpha')
	ax.set_ylabel('beta')
	ax.set_zlabel('variance')

	fig.savefig("../Report/figures/" + name + "AlphaBetaVariance")


	return

def plotEnergyVsAlpha(data, name):

	alpha = data[:,3]
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

def plotChargeDensity(data, name):

	x1 = data[:,2]
	y1 = data[:,3]
	z1 = data[:,4]

	x2 = data[:,5]
	y2 = data[:,6]
	z2 = data[:,7]

	fig = pl.figure()
	ax = p3.Axes3D(fig)
	pl.title("Energy vs alpha and beta, " + name)
	ax.scatter(x1,y1,z1)
	ax.set_xlabel('alpha')
	ax.set_ylabel('beta')
	ax.set_zlabel('Energy')

	fig.savefig("../Report/figures/VarianceVsAlpha" + name)



	return




#Decide what we want to plot this time

# name = "HeliumSimpleAnalytical"
name = "Beryllium";

data = np.genfromtxt("outfiles/" + name + "_alpha_beta_10M")
# datatime = np.genfromtxt("outfiles/" + name +"_timeStep_10M")
# dataSample = np.genfromtxt("outfiles/" + name +"_blockingSamples_10M")


findLowestEnergy(data, name)
# plotResultsVSTimestep(datatime , name)
# plotEnergyVsAlphaBeta(data, name)
# plotEnergyVsAlpha(data, name)
# plotResultsVSTimestep(datatime , name)
# plotChargeDensity(dataSample[0 : : 100, :], name)

pl.show()

