# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 10:59:55 2013

@author: svenni
"""
from pylab import *
from sys import argv
from glob import glob
from time import time
from fys4460 import dataType, boltzmannConstant, headerType

fileNames = argv[1:]

if len(fileNames) == 1:
    fileNames = glob(fileNames[0])
fileNames.sort()
loadTime = 0
calculateTime = 0
plotTime = 0

#skipFiles = 10
myBins = array([])
normBins = array([])
averageDisplacements = zeros(len(fileNames))
averageSquareDisplacements = zeros(len(fileNames))
times = zeros(len(fileNames))
dt = 0.01
i = 0
for fileName in fileNames:
    print "Loading data for " + fileName
    #velocities = loadtxt(fileName, skiprows=2, usecols=[4,5,6], unpack=True)
    f = open(fileName, "rb")
    header = fromfile(f, dtype=headerType, count=1)
    atoms = fromfile(f, dtype=dataType)
#    f = h5py.File(fileName, "r")
#    atoms = f.get("ArrayOfStructures")
    averageDisplacements[i] = header["averageDisplacement"][0]
    averageSquareDisplacements[i] = header["averageSquareDisplacement"][0]

    print "averageSquareDisplacements[i]: ", averageSquareDisplacements[i]
    times[i] = header["time"][0]
    i += 1
#    f.close()

figure()
plot(times, averageSquareDisplacements, label="Average square displacement")
legend()
grid()

show()
