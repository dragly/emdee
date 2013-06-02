# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 11:15:36 2013

@author: svenni
"""

from pylab import *
from fys4460 import loadAtoms, listSaveFileNames, saveAtoms
from sys import argv
from glob import glob
import os
from os.path import join

import re

def tryint(s):
    try:
        return int(s)
    except:
        return s
    
def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)

baseDir = argv[1]
baseStartFileName = argv[2]
baseEndFileName = argv[3]

porosities = []
discharges = []
dirs = os.listdir(baseDir)
sort_nicely(dirs)
for subdir in dirs:
    print subdir
    startFileName = join(baseDir, subdir, baseStartFileName)
    endFileName = join(baseDir, subdir, baseEndFileName)
    
    startHeader, startLammpsHeader, startAtoms = loadAtoms(startFileName)
    startAtoms = sort(startAtoms)
    
    endHeader, endLammpsHeader, endAtoms = loadAtoms(endFileName)
    endAtoms = sort(endAtoms)
    
    nAtoms = len(startAtoms)
    nFreeAtoms = len(where(startAtoms["isPositionFixed"] == 0)[0])
    porosity = float(nFreeAtoms) / float(nAtoms)
    
    systemVector = startHeader["upperBounds"][0] - startHeader["lowerBounds"][0]
    systemVolume = linalg.norm(systemVector)
    atomVolume = systemVolume / nAtoms
    
    nAtomsThroughSurface = len(where(startAtoms["position"][:,2] + endAtoms["displacement"][:,2] > startHeader["upperBounds"][0][2])[0])
    nAtomsThroughSurface -= len(where(startAtoms["position"][:,2] + endAtoms["displacement"][:,2] < startHeader["lowerBounds"][0][2])[0])
    print "nAtomsThroughSurface: " + str(nAtomsThroughSurface)
    
    volumeFlux = nAtomsThroughSurface * atomVolume
    
    timeDiff = endHeader["time"] - startHeader["time"]
    
    porosities.append(porosity)
    print "porosity: " + str(porosity)
    discharges.append( volumeFlux / timeDiff )
    print "discharge: " + str(volumeFlux / timeDiff) 
figure()
plot(porosities, discharges, "-o")
plot(porosities[-6:-3], discharges[-6:-3], "-o", color="red")

xlabel("porosity (%)")
ylabel(u"discharge (m³/s)")
print "Incline: " + str((discharges[1] - discharges[3]) / (porosities[1] - porosities[3]))
show()