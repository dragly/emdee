# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 11:15:36 2013

@author: svenni
"""

from pylab import *
from fys4460 import loadAtoms, listSaveFileNames, saveAtoms
from sys import argv
from glob import glob
fileNames = argv[1:]
if len(fileNames) == 1:
    fileNames = glob(fileNames[0])
fileNames = sort(fileNames)
prevAtoms = []

nBins = 25
bins = linspace(0.0, 3e-9, nBins + 1)
xVals = bins[:-1] + (bins[1] - bins[0]) / 2
yVals = zeros(nBins)
nAtoms = zeros(nBins)

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
    
    for i in range(nBins):        
        minRadius = bins[i]
        maxRadius = bins[i + 1]
        indicesAbove = where(radius > minRadius)
        radiusAbove = radius[indicesAbove]
        atomsAbove = atoms[indicesAbove]
        indicesInBin = where(radiusAbove < maxRadius)
        radiusInBin = radiusAbove[indicesInBin]
        atomsInBin = atomsAbove[indicesInBin]
        nAtoms[i] += len(atomsInBin)
        if len(atomsInBin) > 0:
            yVals[i] += mean(atomsInBin["velocity"][:,2])

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


show()