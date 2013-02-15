# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 09:41:51 2013

@author: svenni
"""
from pylab import *
from sys import argv
from glob import glob
from time import time
from fys4460 import dataType, boltzmannConstant, headerType
import h5py

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
kineticEnergies = zeros(len(fileNames))
#kineticEnergies2 = zeros(len(fileNames))
potentialEnergies = zeros(len(fileNames))
times = zeros(len(fileNames))
cTemperatures = zeros(len(fileNames))
dt = 0.01
i = 0
for fileName in fileNames:
    t1 = time()
    print "Loading data for " + fileName
    #velocities = loadtxt(fileName, skiprows=2, usecols=[4,5,6], unpack=True)
    f = open(fileName, "rb")
    header = fromfile(f, dtype=headerType, count=1)
    atoms = fromfile(f, dtype=dataType)
#    f = h5py.File(fileName, "r")
#    atoms = f.get("ArrayOfStructures")
    atomMass = 6.6353628e-26
    velocityMagnitude = sqrt(atoms["velocity"][:,0]**2 + atoms["velocity"][:,1]**2 + atoms["velocity"][:,2]**2)
    kineticEnergy = 0.5 * atomMass * sum(velocityMagnitude**2)
    potentialEnergy = sum(atoms["potential"])

    print "kineticEnergy: ", kineticEnergy
    print "potentialEnergy: ", potentialEnergy
    kineticEnergies[i] = kineticEnergy
#    kineticEnergies2[i] = header["kineticEnergy"][0]
    potentialEnergies[i] = potentialEnergy
    cTemperatures[i] = header["temperature"][0]
    
    t2 = time()
    calculateTime += t2 - t1
    print calculateTime
    times[i] = dt * i
    i += 1
#    f.close()

figure()
plot(times, kineticEnergies, label="Kinetic")
plot(times, potentialEnergies, label="Potential")
plot(times, kineticEnergies + potentialEnergies, label="Sum")
legend()
grid()

nMovingAverage = 100
figure()
unitTemperature = 119.74
temperatures = kineticEnergies / (3. / 2. * len(atoms) * boltzmannConstant)
plot(times, temperatures, label="Temperature")
plot(times, cTemperatures, label="Temperature")
#plot(times[nMovingAverage - 1:], movavg(temperatures, nMovingAverage), label="Moving average")
legend()
grid()

show()