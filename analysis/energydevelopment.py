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
from os.path import expanduser, join, split
from pylibconfig import Config
#matplotlib.rcParams["mathtext.default"] = "regular"
#matplotlib.rcParams["font.size"] = 16
#import os
#import h5py
configFilePaths = argv[1:]
if len(configFilePaths) == 1:
    configFilePaths = glob(configFilePaths[0])
for configFilePath in configFilePaths:
    saveDir, configFileName = split(configFilePath)
    
    config = Config()
    config.readFile(configFilePath)
    fileNames = config.value("simulation.saveFileName")[0]
    fileNames = expanduser(fileNames)
    fileNames = glob(fileNames)
    
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
    
    for percentage in array([0.9, 1.0]):
        cutoff = int(percentage * len(times))
        cutoffString = ("%.2f" % percentage).replace(".", "_")
        # Energy plot
        figure()
        picoTimes = times[-cutoff:] / 1e-12
        kineticEnergiesCutoff = kineticEnergies[-cutoff:]
        potentialEnergiesCutoff = potentialEnergies[-cutoff:]
        cTemperaturesCutoff = cTemperatures[-cutoff:]
        pressuresCutoff = pressures[-cutoff:]
        
        
        electronVolt = 1 #1.6e-19 # J
        plot(picoTimes, kineticEnergiesCutoff / electronVolt, label="Kinetic")
        plot(picoTimes, potentialEnergiesCutoff / electronVolt, label="Potential")
        plot(picoTimes, (kineticEnergiesCutoff + potentialEnergiesCutoff) / electronVolt, label="Sum")
        xlabel("t [ps]")
        ylabel(r"energy [eV]")
        legend()
        grid()
        savefig(saveDir + "/energy" + cutoffString + ".pdf")
        
        # Temperature plot
        nMovingAverage = 100
        figure()
        unitTemperature = 119.74
        #temperatures = kineticEnergies / (3. / 2. * len(atoms) * boltzmannConstant)
        #plot(times, temperatures, label="Temperature")
        plot(picoTimes, cTemperaturesCutoff, label="Temperature")
        plot(picoTimes[nMovingAverage - 1:], movavg(cTemperaturesCutoff, nMovingAverage), label="Moving average")
        xlabel("t [ps]")
        ylabel(r"temperature [K]")
        legend()
        grid()
        savefig(saveDir + "/temperature" + cutoffString + ".pdf")
        
        # Pressure plot
        figure()
        plot(picoTimes, pressuresCutoff, label="Pressure")
        xlabel("t [ps]")
        ylabel(r"pressure [Pa]")
        legend()
        grid()
        savefig(saveDir + "/pressure" + cutoffString + ".pdf")
        
        # Pressure temperature plot
        figure()
        scatter(cTemperaturesCutoff, pressuresCutoff, label="Pressure vs temperature")
        xlabel(r"temperature [K]")
        ylabel(r"pressure [Pa]")
        legend()
        grid()
        savefig(saveDir + "/pressure-vs-temperature" + cutoffString + ".pdf")
    
show()