import numpy as np
import matplotlib.pyplot as plt
from sys import argv


filename = "../../build-blocking-Desktop-Debug/outputSTD.dat"
infile = open(filename, "r")
blocksize = []
std = []
mean = []
for line in infile:
	blocksize.append(float(line.split()[0]))
	std.append(float(line.split()[1]))
	mean.append(float(line.split()[2]))
	
infile.close
np.array(blocksize)
np.array(std)
np.array(mean)

#ax.zaxis.set_major_locator(LinearLocator(10))
#fig.colorbar(surf, shrink=0.5, aspect=5)

plt.plot(blocksize, std, "o-", label="")
plt.xlabel("Blocksize")
plt.ylabel("Standard deviation")
#plt.hold("on")
plt.title("Beryllium")
#plt.axis([dt[0], dt[-1], 0, 10])
#plt.legend()
plt.savefig("../../Report/figures/" + "Beryllium_blocking" + ".png")

plt.show()
