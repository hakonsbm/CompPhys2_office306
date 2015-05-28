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
	variance = np.sqrt(data[:,1] - data[:,0]**2)/np.sqrt(1000000)
	timeStep = data[:,6]
	
	print energy

	timeFig = pl.figure()
	pl.title('Energy vs timestep, alpha = %.2f, beta = %.2f'  % (data[0,3] , data[0,4]) + " for " + name)

	pl.plot(timeStep , energy)
	pl.ylabel("Energy")
	pl.xlabel("Timestep")

	timeFig.savefig("../../Report/figures/" + name + "TimeEnergy")

	timeFig = pl.figure()
	pl.title('Energy vs timestep, alpha = %.2f, beta = %.2f'  % (data[0,3] , data[0,4]) + " for " + name)

	pl.plot(timeStep , variance)
	pl.ylabel("Variance")
	pl.xlabel("Timestep")

	timeFig.savefig("../../Report/figures/" + name + "TimeVariance")




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

	fig.savefig("../../Report/figures/" + name + "AlphaBetaEnergy")



	fig=pl.figure()
	ax = p3.Axes3D(fig)
	pl.title("Variance vs Alpha and beta, " + name)
	ax.plot_trisurf(alpha,beta,variance)
	ax.set_xlabel('alpha')
	ax.set_ylabel('beta')
	ax.set_zlabel('variance')

	fig.savefig("../../Report/figures/" + name + "AlphaBetaVariance")

	findLowestEnergy(data, name)


def plotEnergyVsAlpha(data, name):

	alpha = data[:,3]
	energy = data[:,0]
	variance = data[:,2] 

	energyFig = pl.figure()
	pl.title("Energy vs alpha, " + name)

	pl.plot(alpha,energy)

	pl.xlabel("alpha")
	pl.ylabel("Energy")

	energyFig.savefig("../../Report/figures/EnergyVsAlpha" + name)

	varianceFig = pl.figure()
	pl.title(" Variance vs alpha, " + name)
	pl.plot (alpha,variance)

	varianceFig.savefig("../../Report/figures/VarianceVsAlpha" + name)

	findLowestEnergy(data, name)

def plotChargeDensity(data, name):

	if(name == "Beryllium"):
		nElectrons = 4
	else: 
		if (name == "Neon"):
			nElectrons = 10
		else:
			if (name == "HeliumJastrowAnalytical" or "HydrogenTwo"):
				nElectrons = 2



	norm = np.ndarray(shape=(0), dtype=float, order='F')

	datapoints = data.shape[0]
	ylimit = (0, 3)
	xlimit = (0, 4)


	#This sorts out the coordinates of th eelectrons and stores them as r values
	for i in range(0 , nElectrons):
		lower = 4 + 3*i
		upper = 7 + 3*i
		r = data[: , lower : upper ]
		normTemp = r[:, 0]**2 + r[:,1]**2 + r[:,2]**2
		norm = np.append(norm,normTemp)
		normTemp = sorted(normTemp)


####################################################
	#If only one electron is wanted enable this

	# electron = 0
	# lower = 4 + 3*electron
	# upper = 7 + 3*electron
	# r = data[: , lower : upper ]
####################################################



	norm = sorted(norm)
	norm = np.sqrt(norm)
# 
	# print norm

	fig = pl.figure()
	# ax = p3.Axes3D(fig)
	pl.title("Charge Density of " + name)
	pl.xlabel("r")
	pl.xlim([0,4])
	pl.ylim([0,1.2])

	pl.hist(norm, normed=True, bins=100)


	# ax.scatter(x1,y1,z1)
	# ax.set_xlabel('x')
	# ax.set_ylabel('y')
	# ax.set_zlabel('z')

	fig.savefig("../../Report/figures/ChargeDensity" + name)

	xySlice = np.ndarray(shape=(3), dtype=float, order='F')
	xySlice[2] = 0

	planeIgnore = 0

	#Now we want to do a plot of a slice around the x-y axis
	for i in range(0, r.shape[0] ):
		if np.abs(r[i,planeIgnore]) < 0.01: 
			xySlice = np.vstack((xySlice, r[i]))
	
	# print xySlice

	fig = pl.figure()
	pl.plot(xySlice[:,np.abs(planeIgnore - 1)], xySlice[:,np.abs(planeIgnore - 2)], ".")
	pl.title("Charge Density of " + name)
	pl.xlabel("x")
	pl.ylabel("z")
	# pl.ylim([-1,1])
	# pl.xlim([-1,1])

	fig.savefig("../../Report/figures/ChargeDensityXYSlice" + name)


	###################################################
	#And at the last we want to create a 3D plot of it
	###################################################

	# print r.shape
	# print xySlice.shape

	# xySlice = xySlice[1:xySlice.shape[0]:50,:]
	# print xySlice.shape


	# return
	# r = xySlice
	print r.shape
	r = r[10000::10,:]
	print r.shape
	fig=pl.figure()
	ax = p3.Axes3D(fig)
	pl.title("Charge Density, " + name)
	ax.scatter(r[:,0],r[:,1],r[:,2])
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_zlabel('z')
	ax.set_xlim([-3,3])
	ax.set_ylim([-3,3])
	ax.set_zlim([-3,3])



	fig.savefig("../../Report/figures/ChargeDensity3D" + name)


	return

def plotVarVSnCycles(data, name):
	
	variance = np.sqrt(data[:,1] - data[:,0]**2)/np.sqrt(data[:,7]) 
	nCycles = data[:, 7 ]



	energyFig = pl.figure()
	pl.title("Variance vs Cycles, " + name)

	pl.plot(nCycles , variance)

	pl.xlabel("Cycles")
	pl.ylabel("Variance")

	energyFig.savefig("../../Report/figures/VarianceNCycles" + name)


	return





#Decide what we want to plot this time

# name = "HeliumSimpleAnalytical"
# name = "HeliumJastrowAnalytical"
name = "Beryllium"
# name = "Neon"
# name = "HydrogenTwo"
# name = "Helium"




#	Picks the relevant data sample
# data = np.genfromtxt("../outfiles/" + name + "_alpha_beta")
# datatime = np.genfromtxt("../outfiles/" + name +"_timeStep")
dataSample = np.genfromtxt("../outfiles/" + name +"_blockingSamples")
# dataCycles = np.genfromtxt("../outfiles/" + name +"_nCycles")


# findLowestEnergy(data, name)
# plotResultsVSTimestep(datatime , name)
# plotEnergyVsAlphaBeta(data, name)
# plotEnergyVsAlpha(data, name)
# plotResultsVSTimestep(datatime , name)
plotChargeDensity(dataSample, name)
# plotVarVSnCycles(dataCycles[1:,:], name)


# pl.show()

