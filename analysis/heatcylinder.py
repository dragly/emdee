# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 16:19:22 2013

@author: svenni
"""

from pylab import *
from fys4460 import loadAtoms, listSaveFileNames, saveAtoms
from sys import argv
fileName = argv[1]
header, lammpsHeader, atoms = loadAtoms(fileName)

cylinderRadius = 2e-9

systemLengths = header["upperBounds"][0] - header["lowerBounds"][0]

nFixed = 0
for atom in atoms:
    diff = (atom["position"] - systemLengths / 2)
    if(sqrt(diff[0]**2 + diff[1]**2) < cylinderRadius):
        atom["velocity"][2] = 10*abs(atom["velocity"][2])
#        atom["isPositionFixed"] = 0
    else:
        atom["velocity"] *= 0.001
#        atom["isPositionFixed"] = 1
        nFixed += 1
    
print "Fixed", nFixed, "of", len(atoms), "atoms"
saveAtoms(header, lammpsHeader, atoms, argv[2])
