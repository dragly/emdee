# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 09:41:51 2013

@author: svenni
"""
from pylab import *
from sys import argv
from glob import glob
from time import time
from fys4460 import loadAtoms
from os.path import expanduser, join, split
from pylibconfig import Config
matplotlib.rcParams["mathtext.default"] = "regular"
matplotlib.rcParams["font.size"] = 16

configFilePath = argv[1]
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
myBins = array([0,18])
normBins = array([0,25])
for fileName in fileNames[0:100:5]:
    print "Loading data for " + fileName
    
    t1 = time()
    header, atoms = loadAtoms(fileName)
    t2 = time()
    loadTime += t2 - t1
    
    t1 = time()
    velocityMagnitude = sqrt(atoms["velocity"][:,0]**2 + atoms["velocity"][:,1]**2 + atoms["velocity"][:,2]**2)
    t2 = time()
    calculateTime += t2 - t1
    
    #figure()
    clf()

    myBins = hist(atoms["velocity"][:,0], bins=myBins[1], alpha=0.3, label=r"$v_x$")
    hist(atoms["velocity"][:,1], bins=myBins[1], alpha=0.3, label=r"$v_y$")
    hist(atoms["velocity"][:,2], bins=myBins[1], alpha=0.3, label=r"$v_z$")
    legend()
    xlabel("velocity [m/s]")
    ylabel("number of atoms")
    savefig(join(saveDir, split(fileName)[1].replace(".bin", "") + "-comp.pdf"))
    
    t1 = time()
    print "Plotting for " + fileName
    clf()
    normBins = hist(velocityMagnitude, bins=normBins[1])
    xlabel("speed [m/s]")
    ylabel("number of atoms")
    savefig(join(saveDir, split(fileName)[1].replace(".bin", "") + "-norm.pdf"))
    #clf()
    #plot(normBins[1][1:len(normBins[0]) + 1], normBins[0] / (4 * pi * normBins[1][1:len(normBins[0])+1]**2))
    #show()
    #plot(velocityMagnitude)
    
    #show()
    print loadTime
    print calculateTime
    t2 = time()
    plotTime += t2 - t1
    print plotTime
#    show()