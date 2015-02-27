import numpy as np
import matplotlib.pyplot as plt
from sys import argv


filename = "./outfiles/HeliumSimpleAnalytical_IS_dt_0.1to1.0"
infile = open(filename, "r")
energy = []
energySq = []
alpha = []
beta = []
dt = []
for line in infile:
	energy.append(float(line.split()[0]))
	energySq.append(float(line.split()[1]))
	alpha.append(float(line.split()[2]))
	beta.append(line.split()[3])
	dt.append(line.split()[4])
	
infile.close
np.array(energy)
np.array(energySq)
np.array(alpha)
np.array(beta)
np.array(dt)

plt.plot(dt, energy, "x-", label="")
#plt.hold("on")
plt.title("")
plt.xlabel("dt")
plt.ylabel("energy")
#plt.axis([dt[0], dt[-1], 0, 10])
#plt.legend()
plt.savefig(filename + ".png")

plt.show()
