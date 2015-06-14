import numpy as np
import pylab as pl
import scipy.stats as stats
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib.colors import LogNorm







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
	if (name == "Neon"):
		nElectrons = 10
	if (name == "HeliumJastrowAnalytical" or "HydrogenTwo"):
		nElectrons = 2
	if(name == "BerylliumTwo"):
		nElectrons = 8


	datapoints = data.shape[0]


	norm = np.ndarray(shape=(0), dtype=float, order='F')
	positions = np.ndarray(shape=(datapoints*nElectrons, 3), dtype = float)
	electrons = np.ndarray(shape = (datapoints, 3*nElectrons) , dtype = float)


	print name
	print data.shape
	print positions.shape
	print electrons.shape

	ylimit = (0, 3)
	xlimit = (0, 4)


	
	#This sorts out the coordinates of th eelectrons and stores them as r values
	for i in range(0 , nElectrons):
		lower = 4 + 3*i
		upper = 7 + 3*i
		r = data[: , lower : upper ]
		positions[ i*datapoints : (i + 1)*datapoints , : ] = r
		# positions = np.vstack((positions, r))
		print "Done with atom %d", i

		# normTemp = r[:, 0]**2 + r[:,1]**2 + r[:,2]**2
		# norm = np.append(norm,normTemp)
		# normTemp = sorted(normTemp)

	print positions.shape



	print "This was not the problem"


	#Testing to see how large the array needs to be first to avoid having to resize the array
	#So we check how many entries should be in the 2D slice, then creates the array the right size before putting them there
	counter = 0;

	for i in range(0,positions.shape[0]):
		if np.abs(positions[i,1]) < 0.5:
			counter = counter + 1

	slice2D = np.ndarray(shape=(counter, 3), dtype = float)
	print counter
	counter = 0
	for i in range(0,positions.shape[0]):
		if np.abs(positions[i,1]) < 0.2:
			# slice2D = np.vstack((slice2D,positions[i,:]))
			slice2D[counter,:] = positions[i,:]
			counter = counter + 1

	print "Was there a problem?"

####################################################
	#If only one electron is wanted enable this

	# electron = 0
	# lower = 4 + 3*electron
	# upper = 7 + 3*electron
	# r = data[: , lower : upper ]
####################################################
	norm = positions[:, 0]**2 + positions[:,1]**2 + positions[:,2]**2
	norm = sorted(norm)
	norm = np.sqrt(norm)
# 
	# print norm

	fig = pl.figure()
	# ax = p3.Axes3D(fig)
	pl.title("Radial Charge Density of " + name)
	pl.xlabel("r [a.u.]")
	pl.xlim([0,6])
	pl.ylim([0,1.6])

	pl.hist(norm, normed=True, bins=200)


	fig.savefig("../../Report/figures/ChargeDensity" + name)


	# 	#Creating a better slicething
	# for i in range (0,nElectrons):
	fig = pl.figure()
	pl.hist2d(slice2D[:,0], slice2D[:,2], bins=100, norm=LogNorm())
	pl.title("Charge density " + name)
	pl.xlabel("x axis [a.u.]")
	pl.ylabel("z axis [a.u.]")
	pl.xlim([-6,6])
	pl.ylim([-6,6])
	pl.colorbar()
	fig.savefig("../../Report/figures/ChargeDensity3D" + name)


	if(name == "BerylliumTwo"):
		#Plotting each electron by itself
		# Three subplots sharing both x/y axes
		fig, ((ax1, ax2), (ax3 , ax4), (ax5, ax6), (ax7, ax8)) = pl.subplots(4,2, sharex='col', sharey='row')
		# fig, ((ax1, ax2), (ax3, ax4)) = pl.subplots(2,2)#, sharex=True, sharey=True)
		pl.xlim([-6,6])
		pl.ylim([-6,6])
		pl.xlabel("x axis [a.u.]")
		pl.ylabel("z axis [a.u.]")

		k = 0
		lower = 4 + 3*k
		upper = 7 + 3*k
		r = data[: , lower : upper ]

		ax1.hist2d(r[:,0], r[:,2], bins=100, norm=LogNorm())

		k = 1
		lower = 4 + 3*k
		upper = 7 + 3*k
		r = data[: , lower : upper ]

		ax2.hist2d(r[:,0], r[:,2], bins=100, norm=LogNorm())

		k = 2
		lower = 4 + 3*k
		upper = 7 + 3*k
		r = data[: , lower : upper ]

		ax3.hist2d(r[:,0], r[:,2], bins=100, norm=LogNorm())

		k = 3
		lower = 4 + 3*k
		upper = 7 + 3*k
		r = data[: , lower : upper ]

		ax4.hist2d(r[:,0], r[:,2], bins=100, norm=LogNorm())

		k = 4
		lower = 4 + 3*k
		upper = 7 + 3*k
		r = data[: , lower : upper ]

		ax5.hist2d(r[:,0], r[:,2], bins=100, norm=LogNorm())

		k = 5
		lower = 4 + 3*k
		upper = 7 + 3*k
		r = data[: , lower : upper ]

		ax6.hist2d(r[:,0], r[:,2], bins=100, norm=LogNorm())

		k = 6
		lower = 4 + 3*k
		upper = 7 + 3*k
		r = data[: , lower : upper ]

		ax7.hist2d(r[:,0], r[:,2], bins=100, norm=LogNorm())

		k = 7
		lower = 4 + 3*k
		upper = 7 + 3*k
		r = data[: , lower : upper ]

		ax8.hist2d(r[:,0], r[:,2], bins=100, norm=LogNorm())




		ax1.set_title("                              One Body Density electron " + name)
		fig.subplots_adjust(hspace=0)

		fig.savefig("../../Report/figures/OneBodyDensityElectrons" + name)

	if(name == "HydrogenTwo"):
		#Plotting each electron by itself
		# Three subplots sharing both x/y axes
		fig, (ax1, ax2) = pl.subplots(2,1, sharex='col', sharey='row')
		# fig, ((ax1, ax2), (ax3, ax4)) = pl.subplots(2,2)#, sharex=True, sharey=True)
		

		k = 0
		lower = 4 + 3*k
		upper = 7 + 3*k
		r = data[: , lower : upper ]

		ax1.hist2d(r[:,0], r[:,2], bins=100, norm=LogNorm())
		ax1.set_xlim([-3,3])
		ax1.set_ylim([-3,3])

		k = 1
		lower = 4 + 3*k
		upper = 7 + 3*k
		r = data[: , lower : upper ]

		ax2.hist2d(r[:,0], r[:,2], bins=100, norm=LogNorm())

		pl.xlim([-3,3])
		pl.ylim([-3,3])
		pl.xlabel("x axis [a.u.]")
		pl.ylabel("z axis [a.u.]")


		ax1.set_title("One Body Density electron " + name)
		fig.subplots_adjust(hspace=0)

		fig.savefig("../../Report/figures/OneBodyDensityElectrons" + name)


	# ax1.colorbar()

	# pl.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)


		# pl.xlim([-6,6])
		# pl.ylim([-6,6])


	
	# ax1.plot(x, y)
	# ax1.set_title('Sharing both axes')
	# ax2.scatter(x, y)
	# ax3.scatter(x, 2 * y ** 2 - 1, color='r')
	# # Fine-tune figure; make subplots close to each other and hide x ticks for
	# # all but bottom plot.
	# f.subplots_adjust(hspace=0)
	# plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

	# counter = 0
	# for i in range(k*datapoints,(k+1)*datapoints):
	# 	if np.abs(positions[i,1]) < 0.5:
	# 		counter = counter + 1



	# slice2D = np.ndarray(shape=(counter, 3), dtype = float)
	# counter = 0
	# # for i in range(k*data.shape[0],(k+1)*data.shape[0]):
	# # 	if np.abs(positions[i,1]) < 0.2:
	# # 		# slice2D = np.vstack((slice2D,positions[i,:]))
	# # 		slice2D[counter,:] = positions[i,:]
	# # 		counter = counter + 1
	

	pl.show()



	# xySlice = np.ndarray(shape=(3), dtype=float, order='F')
	# xySlice[2] = 0

	# planeIgnore = 0

	# #Now we want to do a plot of a slice around the x-y axis
	# for i in range(0, r.shape[0] ):
	# 	if np.abs(r[i,planeIgnore]) < 0.05: 
	# 		xySlice = np.vstack((xySlice, r[i]))
	

	# norm = r[:, 0]**2 + r[:,1]**2
	# norm = sorted(norm)
	# norm = np.sqrt(norm)

	# fig = pl.figure()
	# # ax = p3.Axes3D(fig)
	# pl.title("Charge Density of " + name)
	# pl.xlabel("r")
	# pl.xlim([0,3])
	# pl.ylim([0,1.6])

	# pl.hist(norm, normed=True, bins=200)
	# # # print xySlice

	# # fig = pl.figure()
	# # pl.plot(xySlice[:,np.abs(planeIgnore - 1)], xySlice[:,np.abs(planeIgnore - 2)], ".")
	# # pl.title("Charge Density of " + name)
	# # pl.xlabel("x")
	# # pl.ylabel("z")
	# # # pl.ylim([-1,1])
	# # # pl.xlim([-1,1])

	# # fig.savefig("../../Report/figures/ChargeDensityXYSlice" + name)

	# #Creating a better slicething
	# fig = pl.figure()
	# pl.hist2d(xySlice[:,np.abs(planeIgnore - 1)], xySlice[:,np.abs(planeIgnore - 2)], bins=100, norm=LogNorm())
	# pl.colorbar()
	# pl.show()


	# ###################################################
	# #And at the last we want to create a 3D plot of it
	# ###################################################

	# print r.shape
	# print xySlice.shape

	# xySlice = xySlice[1:xySlice.shape[0]:50,:]
	# print xySlice.shape


	# return
	# r = xySlice

	# r = r[10000::10,:]
	# fig=pl.figure()
	# ax = p3.Axes3D(fig)
	# pl.title("Charge Density, " + name)
	# ax.scatter(r[:,0],r[:,1],r[:,2])
	# ax.set_xlabel('x')
	# ax.set_ylabel('y')
	# ax.set_zlabel('z')
	# ax.set_xlim([-3,3])
	# ax.set_ylim([-3,3])
	# ax.set_zlim([-3,3])



	# fig.savefig("../../Report/figures/ChargeDensity3D" + name)


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



pl.close('all')

#Decide what we want to plot this time

# name = "HeliumSimpleAnalytical"
# name = "HeliumJastrowAnalytical"
# name = "Beryllium"
# name = "Neon"
name = "HydrogenTwo"
# name = "Helium"
# name = "BerylliumTwo"




#	Picks the relevant data sample
# data = np.genfromtxt("../outfiles/" + name + "_alpha_beta")
# datatime = np.genfromtxt("../outfiles/" + name +"_timeStep")
dataSample = np.genfromtxt("../outfiles/" + name +"_blockingSamples")
# dataCycles = np.genfromtxt("../outfiles/" + name +"_nCycles")

print "Datareading sure took some time"

# findLowestEnergy(data, name)
# plotResultsVSTimestep(datatime , name)
# plotEnergyVsAlphaBeta(data, name)
# plotEnergyVsAlpha(data, name)
# plotResultsVSTimestep(datatime , name)
plotChargeDensity(dataSample, name)
# plotVarVSnCycles(dataCycles[1:,:], name)


pl.show()

