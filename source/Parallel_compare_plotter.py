import numpy as np
import matplotlib.pyplot as plt
from sys import argv

processors = [1,2,3,4]
time = [26.6016, 13.4549, 9.10255, 7.305]
speedup = np.zeros(4)

for i in xrange(0,4):
	speedup[i] = time[0]/time[i]
	print "Processors: ", processors[i]
	print "Speedup: ", speedup[i], "\n"

np.array(processors)
np.array(time)
np.array(speedup)

plt.plot(processors, time, "x-", label="")
#plt.hold("on")
# plt.title("")
# plt.xlabel("dt")
# plt.ylabel("energy")
#plt.axis([dt[0], dt[-1], 0, 10])
#plt.legend()
# plt.savefig(filename + ".png")

plt.show()
