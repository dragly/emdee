# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 11:15:36 2013

@author: svenni
"""

from pylab import *
from fys4460 import loadAtoms, listSaveFileNames, saveAtoms
from sys import argv
from glob import glob
from os.path import expanduser, join, split, isdir, islink
from os import readlink
from pylibconfig import Config
from time import time

startTime = time()

configFilePaths = argv[1:]
if len(configFilePaths) == 1:
    configFilePaths = glob(configFilePaths[0])
for configFilePath in configFilePaths:
    if isdir(configFilePath):
        print "Info: Input is directory, selecting first config file"
        configFilePath = glob(join(configFilePath, "*.cfg"))[0]
    
    if islink(configFilePath):
        configFilePath = readlink(configFilePath)

    config = Config()
    config.readFile(configFilePath)
    fileNames = config.value("simulation.saveFileName")[0]
    fileNames = expanduser(fileNames)
    fileNames = glob(fileNames)
    
    fileNames.sort()
    
    prevAtoms = []
    
    nBins = 50
    bins = linspace(0.0, 3e-9, nBins + 1)
    xVals = bins[:-1] + (bins[1] - bins[0]) / 2
    yVals = zeros(nBins)
    nAtoms = zeros(nBins)
    
#    fileNames = fileNames

    times =  zeros(len(fileNames))   
    drifts = zeros(len(fileNames))
    
#    fileNames = [fileNames[0]]
    
    counter = 0
    for fileName in fileNames:
        header, lammpsHeader, atoms = loadAtoms(fileName)
        atoms = sort(atoms)
        zLength = header["upperBounds"][0,2] - header["lowerBounds"][0,2]
        zMin = header["lowerBounds"][0,2]
        # pick atoms around center
    #    atoms = atoms[where(atoms["position"][:,2] < zMin + zLength / 2 + zLength * 0.4)]
    #    atoms = atoms[where(atoms["position"][:,2] > zMin + zLength / 2 - zLength * 0.4)]
        atoms = atoms[where(atoms["isPositionFixed"] < 1)]
    #    
    #    if len(prevAtoms) > 0:
    #        positiveAtoms = atoms[where(atoms["position"][:,2] > zMin + zLength / 2)]
    #        negativeAtoms = prevAtoms[where(prevAtoms["position"][:,2] < zMin + zLength / 2)]
    #        
    #        movedThrough = intersect1d(positiveAtoms["id"], negativeAtoms["id"])
    #        
    #        if len(movedThrough) > 0:
    #            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    #            print str(len(movedThrough)) + " atoms moved through!"
    #            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        
        #    prevAtoms = atoms.copy()
        centerPoint = header["lowerBounds"] + (header["upperBounds"] - header["lowerBounds"]) / 2
        diff = atoms["position"] - centerPoint
        radius = sqrt(diff[:,0]**2 + diff[:,1]**2)
        
        nAtomsCurrent = histogram(radius, bins=bins)[0]
        nAtoms += nAtomsCurrent
        nAtomsCurrent[where(nAtomsCurrent==0)] = 1
        velocityBins = histogram(radius, bins=bins, weights=atoms["velocity"][:,2])[0]
        yVals += velocityBins / nAtomsCurrent
        
        times[counter] = header["time"]
        drifts[counter] = sum(atoms["velocity"][:,2]) / len(atoms)
        counter += 1
        
#        clf()
        
#        for i in range(nBins):        
#            minRadius = bins[i]
#            maxRadius = bins[i + 1]
#            indicesAbove = where(radius > minRadius)
#            radiusAbove = radius[indicesAbove]
#            atomsAbove = atoms[indicesAbove]
#            indicesInBin = where(radiusAbove <= maxRadius)
#            radiusInBin = radiusAbove[indicesInBin]
#            atomsInBin = atomsAbove[indicesInBin]
#            nAtoms[i] += len(atomsInBin)
#            if len(atomsInBin) > 0:
#                yVals[i] += mean(atomsInBin["velocity"][:,2])
    
    V = pi * (bins[1:]**2 - bins[:-1]**2) * (zLength)
    yVals = yVals / len(fileNames)            
    nAtoms = nAtoms / len(fileNames)            
    figure()            
    plot(xVals, yVals)
    xlabel("distance from cylinder center")
    ylabel("velocity")
    
    figure()
    F = 4.84848e-11 
    Fz = F * (nAtoms / V)
    mu = Fz * (2e-9**2 - xVals**2) / (4 * yVals)
    plot(xVals[:15], mu[:15])
    xlabel("distance from cylinder center")
    ylabel("viscosity")
    
    figure()
    plot(xVals, nAtoms / V)
    xlabel("distance from cylinder center")
    ylabel("atom density")
    
    figure()
    plot(times, drifts)

    endTime = time()
    
    print "Total time taken: ", endTime - startTime    
    
    
    show()
