# Python script to run vmc

import os
import subprocess as sp
import re

path = "../build-vmc-hakonsbm" # Change to your QT build directory

n_procs = "2" # processes
atom = "Beryllium" # atom
runtype = "runCompareParallelize" # runtype
cycles = "1000" # Monte Carlo cycles

os.chdir(path)
output = sp.check_output(["/usr/bin/mpirun","-n", n_procs, "vmc", atom, runtype, cycles])
print "Python output:"
print output
print 

# Using regex to find values
energy = re.search("(Energy: )(.*.\d+)( Energy)", output)
variance = re.search("(Variance: )(.*.\d+)", output)
time = re.search("(processes: )(\d+.\d+)", output)

print energy.group(2)
print variance.group(2)
print time.group(2)