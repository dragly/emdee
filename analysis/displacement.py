# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 10:59:55 2013

@author: svenni
"""
from pylab import *
from sys import argv
from glob import glob
from time import time
from fys4460 import loadAtoms
from os.path import expanduser, join, split, isdir
from scipy import stats
from pylibconfig import Config
matplotlib.rcParams["mathtext.default"] = "regular"

configFilePaths = argv[1:]
if len(configFilePaths) == 1:
    configFilePaths = glob(configFilePaths[0])
for configFilePath in configFilePaths:
    if isdir(configFilePath):
        print "Info: Input is directory, selecting first config file"
        configFilePath = glob(join(configFilePath, "*.cfg"))[0]
        
    saveDir, configFileName = split(configFilePath)
    
    config = Config()
    config.readFile(configFilePath)
    fileNames = config.value("simulation.saveFileName")[0]
    fileNames = expanduser(fileNames)
    fileNames = glob(fileNames)
    fileNames.sort()
#    fileNames = fileNames[400:1400]
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
        header, lammps, atoms = loadAtoms(fileName)
    #    f = h5py.File(fileName, "r")
    #    atoms = f.get("ArrayOfStructures")
        averageDisplacements[i] = header["averageDisplacement"][0]
        averageSquareDisplacements[i] = header["averageSquareDisplacement"][0]
    
        print "averageSquareDisplacements[i]: ", averageSquareDisplacements[i]
        times[i] = header["time"][0]
        i += 1
    #    f.close()
    
    figure()
    picoTimes = times / 1e-12
    plot(picoTimes, averageSquareDisplacements, label="Average square displacement")
    xlabel("t [ps]")
    ylabel(r"distance [$m^2$]")
    legend()
    grid()
    savefig(saveDir + "/displacement.pdf")
    gradient = stats.linregress(times, averageSquareDisplacements)[0]
    diffConstant = gradient / 6
    print "diffConstant ", diffConstant
    f = open(saveDir + "/displacement-diffusion.tex", "w")
    f.write(str(diffConstant))
    f.close()
    show()