import numpy as np
import matplotlib.pyplot as plt
from sys import argv

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator
from matplotlib import cm




filename = "./outfiles/Neon_alpha_beta"
infile = open(filename, "r")
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

def test(x, y):
	return x*x + y


fig = plt.figure()
ax  = fig.gca(projection="3d")
ax.plot_trisurf(alpha, beta, var, cmap=cm.jet, linewidth=0.2)

#ax.zaxis.set_major_locator(LinearLocator(10))
#fig.colorbar(surf, shrink=0.5, aspect=5)
#plt.plot(alpha, energy, "x-", label="")
#plt.hold("on")
#plt.title("")
ax.set_xlabel('alpha')
ax.set_ylabel('beta')
ax.set_zlabel('variance')
#plt.axis([dt[0], dt[-1], 0, 10])
#plt.legend()
plt.savefig(filename + "_var.png")

plt.show()
