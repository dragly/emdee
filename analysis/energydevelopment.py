# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 09:41:51 2013

@author: svenni
"""
from pylab import *
from sys import argv
from glob import glob
from time import time
from fys4460 import loadAtoms, boltzmannConstant
#import os
#import h5py

saveDir = argv[1]

fileNames = argv[2:]

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
pressures = zeros(len(fileNames))
times = zeros(len(fileNames))
cTemperatures = zeros(len(fileNames))
dt = 0.01
i = 0
for fileName in fileNames:
    t1 = time()
    print "Loading data for " + fileName
    #velocities = loadtxt(fileName, skiprows=2, usecols=[4,5,6], unpack=True)
    header, atoms = loadAtoms(fileName)
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
    pressures[i] = header["pressure"][0]
    print "pressure: ", pressures[i]
    
    t2 = time()
    calculateTime += t2 - t1
    print calculateTime
    times[i] = header["time"][0]
    i += 1
#    f.close()

# Energy plot
figure()
picoTimes = times / 1e-12
electronVolt = 1.6e-19 # J
plot(picoTimes, kineticEnergies / electronVolt, label="Kinetic")
plot(picoTimes, potentialEnergies / electronVolt, label="Potential")
plot(picoTimes, (kineticEnergies + potentialEnergies) / electronVolt, label="Sum")
xlabel("t [ps]")
ylabel(r"energy [eV]")
legend()
grid()
savefig(saveDir + "/energy.pdf")

# Temperature plot
nMovingAverage = 100
figure()
unitTemperature = 119.74
temperatures = kineticEnergies / (3. / 2. * len(atoms) * boltzmannConstant)
#plot(times, temperatures, label="Temperature")
plot(picoTimes, cTemperatures, label="Temperature")
plot(picoTimes[nMovingAverage - 1:], movavg(temperatures, nMovingAverage), label="Moving average")
xlabel("t [ps]")
ylabel(r"temperature [K]")
legend()
grid()
savefig(saveDir + "/temperature.pdf")

# Pressure plot
figure()
plot(picoTimes, pressures, label="Pressure")
xlabel("t [ps]")
ylabel(r"pressure [Pa]")
legend()
grid()
savefig(saveDir + "/pressure.pdf")

# Pressure temperature plot
figure()
scatter(temperatures, pressures, label="Pressure vs temperature")
xlabel(r"temperature [K]")
ylabel(r"pressure [Pa]")
legend()
grid()
savefig(saveDir + "/pressure-vs-temperature.pdf")

show()