# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 09:41:51 2013

@author: svenni
"""
from pylab import *
from sys import argv
from glob import glob
from time import time
from fys4460 import dataType

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
for fileName in fileNames:
    print "Loading data for " + fileName
    
    t1 = time()
    atoms = fromfile(fileName, dtype=dataType)
    t2 = time()
    loadTime += t2 - t1
    
    t1 = time()
    velocityMagnitude = sqrt(atoms["velocityX"]**2 + atoms["velocityY"]**2 + atoms["velocityZ"]**2)
    t2 = time()
    calculateTime += t2 - t1
    
    #figure()
    clf()
    if len(myBins) == 0:
        myBins = hist(atoms["velocityX"], bins=25, alpha=0.3)
    else:
        hist(atoms["velocityX"], bins=myBins[1], alpha=0.3)
        
    hist(atoms["velocityY"], bins=myBins[1], alpha=0.3)
    hist(atoms["velocityZ"], bins=myBins[1], alpha=0.3)
    savefig(fileName + "-comp.png")
    
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
    savefig(fileName + "-norm.png")
    t2 = time()
    plotTime += t2 - t1
    print plotTime