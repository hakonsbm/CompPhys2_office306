import numpy as np
import matplotlib.pyplot as plt
from sys import argv

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator
from matplotlib import cm



atom = "HeliumJastrowAnalytical"
#atom = "Beryllium"
#atom = "Neon"

for quant in ["energy"]:

	filename = atom + "_alpha_beta"
	infile = open("./../outfiles/"+filename, "r")
	energy = []
	energySq = []
	var = []
	alpha = []
	beta = []
	dt = []
	for line in infile:
		energy.append(float(line.split()[0]))
		energySq.append(float(line.split()[1]))
		var.append(float(line.split()[2]))
		alpha.append(float(line.split()[3]))
		beta.append(line.split()[4])
		dt.append(line.split()[5])
		
	infile.close
	np.array(energy)
	np.array(energySq)
	np.array(var)
	np.array(alpha)
	np.array(beta)
	np.array(dt)

	min_energy = 1000
	min_energy_alpha = 0;
	min_energy_beta = 0;
	
	min_var = 1000
	min_var_alpha = 0;
	min_var_beta = 0;

	for i in range(len(energy)):
		if energy[i] < min_energy:
			min_energy = energy[i]
			min_energy_alpha = alpha[i]
			min_energy_beta = beta[i]
		if var[i] < min_var:
			min_var = var[i]
			min_var_alpha = alpha[i]
			min_var_beta = beta[i]

	print "Min var = ", min_var
	print "Min var alpha = ", min_var_alpha
	print "Min var beta = ", min_var_beta
	print 
	print "Min energy = ", min_energy
	print "Min energy alpha = ", min_energy_alpha
	print "Min energy beta = ", min_energy_beta
			



	fig = plt.figure()
	ax  = fig.gca(projection="3d")

	if quant == "variance":
		ax.plot_trisurf(alpha, beta, var, cmap=cm.jet, linewidth=0.2)
	elif quant == "energy":
		ax.plot_trisurf(alpha, beta, energy, cmap=cm.jet, linewidth=0.2)

	#ax.zaxis.set_major_locator(LinearLocator(10))
	#fig.colorbar(surf, shrink=0.5, aspect=5)
	#plt.plot(alpha, energy, "x-", label="")
	#plt.hold("on")
	plt.title(atom + " " + quant)
	ax.set_xlabel('alpha')
	ax.set_ylabel('beta')
	ax.set_zlabel(quant+" [au]")
	#plt.axis([dt[0], dt[-1], 0, 10])
	#plt.legend()
	save_path = "../../Report/figures/" + filename + "_" + quant + ".png"
	print "saving to ", save_path
	plt.savefig(save_path) #"../Report/figures/" + 

	plt.show()
