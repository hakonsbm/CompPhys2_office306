# Python script to run vmc

#	Options for runtype
#
# runWithDiffConstants
# runSIWithDiffTimesteps
# runBlockingSampledRun
# runCompareAnalytical
# runDiffNCycles
# runFindAlphaBeta
# runCompareParallelize
# runTests
#

import os
import subprocess as sp
import re

#-n 4 vmc Beryllium runFindAlphaBeta 5000000

#-n 4 vmc HeliumJastrowAnalytical runFindAlphaBeta 200000000

#atom = "Beryllium" # atom
#cycles = "1000000" # Monte Carlo cycles
path = "../build-vmc-hakonsbm-release" # Change to your QT build directory

for rundata in [["HeliumJastrowAnalytical", "200000000"], ["Beryllium", "5000000"]]:

	n_procs = "4" # processes
	runtype = "runFindAlphaBeta" # runtype
	atom = rundata[0]
	cycles = rundata[1]
	print atom
	os.chdir(path)
	output = sp.check_output(["/usr/bin/mpirun","-n", n_procs, "vmc", atom, runtype, cycles])
	#print "Python output:"
	#print output
	#print 

	# Using regex to find values
	#energy = re.search("(Energy: )(.*.\d+)( Energy)", output)
	#variance = re.search("(Variance: )(.*.\d+)", output)
	#time = re.search("(processes: )(\d+.\d+)", output)
	variance_a_b = re.search("The lowest variance, (.+), is found with alpha  (.+) and beta (.+)", output)
	energy_a_b = re.search("The lowest energy, (.+), is found with alpha  (.+) and beta (.+)", output)

	print "The lowest variance, %s, is found with alpha  %s and beta %s" % (variance_a_b.group(1), variance_a_b.group(2), variance_a_b.group(3))
	print "The lowest energy, %s, is found with alpha  %s and beta %s" % (energy_a_b.group(1), energy_a_b.group(2), energy_a_b.group(3))
	#print energy.group(2)
	#print variance.group(2)
	#print time.group(2)
	print 




# expected output:
"""

HELIUM:
The lowest variance, 0.00485688, is found with alpha  1.9 and beta 0.375
The lowest energy, -3.14961, is found with alpha  1.5 and beta 0.741667

BERYLLIUM:
The lowest variance, 0.046652, is found with alpha  4 and beta 0.925
The lowest energy, -19.1429, is found with alpha  2.4 and beta 1.01667

"""

