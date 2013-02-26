# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 08:54:28 2013

@author: svenni
"""

from pylab import *
from sys import argv
from glob import glob
from time import time
from fys4460 import loadAtoms, boltzmannConstant
from os.path import expanduser, join, split
from pylibconfig import Config
#import os
#import h5py
configFilePath = argv[1]
saveDir, configFileName = split(configFilePath)

config = Config()
config.readFile(configFilePath)
fileNames = config.value("simulation.saveFileName")[0]
fileNames = expanduser(fileNames)
fileNames = glob(fileNames)

fileNames.sort()
distancesList = []
for fileName in fileNames[400:1400:100]:
    print "Loading file\n" + fileName
    header, atoms = loadAtoms(fileName)

    lowerBounds = header["lowerBounds"][0]
    upperBounds = header["upperBounds"][0]
    
    sideLengths = upperBounds - lowerBounds
    
    positions = atoms["position"]
    distances = zeros(int(len(positions) * (len(positions) - 1) / 2.0))
    counter = 0
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            distance = positions[j] - positions[i]
            for k in range(3):
                if abs(distance[k]) < abs(distance[k] + sideLengths[k]) and abs(distance[k]) < abs(distance[k] - sideLengths[k]):
                    continue
                elif abs(distance[k] + sideLengths[k]) <  abs(distance[k] - sideLengths[k]):
                    distance[k] = distance[k] + sideLengths[k]
                else:
                    distance[k] = distance[k] - sideLengths[k]
            distances[counter] = linalg.norm(distance, 2)
            counter += 1
        print i
    distancesList.append(distances)

distBins = [0,20]
totalBins = zeros(distBins[1])
figure()
for distances in distancesList:
    distBins = hist(distances, bins=distBins[1])
    V = 4./3. * pi * distBins[1]**3
    Vdiff = V[1:] - V[:-1]
    newBins = distBins[0] / Vdiff
    totalBins = totalBins + newBins

figure()
plot(distBins[1][1:], newBins)
savefig(join(saveDir, "radial-distribution-latest.pdf"))
figure()
plot(distBins[1][1:], totalBins)
savefig(join(saveDir, "radial-distribution.pdf"))
show()