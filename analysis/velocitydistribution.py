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
for fileName in fileNames:
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
    if len(myBins) == 0:
        myBins = hist(atoms["velocity"][:,0], bins=25, alpha=0.3)
    else:
        hist(atoms["velocity"][:,0], bins=myBins[1], alpha=0.3)
        
    hist(atoms["velocity"][:,1], bins=myBins[1], alpha=0.3)
    hist(atoms["velocity"][:,2], bins=myBins[1], alpha=0.3)
    savefig(saveDir + "/" + fileName.split("/")[-1] + "-comp.pdf")
    
    t1 = time()
    print "Plotting for " + fileName
    clf()
    if len(normBins) == 0:
        normBins = hist(velocityMagnitude, bins=40)
    else:
        hist(velocityMagnitude, bins=normBins[1])
    #clf()
    #plot(normBins[1][1:len(normBins[0]) + 1], normBins[0] / (4 * pi * normBins[1][1:len(normBins[0])+1]**2))
    #show()
    #plot(velocityMagnitude)
    
    #show()
    print "Saving for " + fileName
    print loadTime
    print calculateTime
    savefig(saveDir + "/" + fileName.split("/")[-1] + "-norm.pdf")
    t2 = time()
    plotTime += t2 - t1
    print plotTime
#    show()