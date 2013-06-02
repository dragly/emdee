# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 16:19:22 2013

@author: svenni
"""

from pylab import *
from fys4460 import loadAtoms, listSaveFileNames, saveAtoms, dataType
from sys import argv
fileName = argv[1]
header, lammpsHeader, atoms = loadAtoms(fileName)

cylinderRadius = 2e-9

systemLengths = header["upperBounds"][0] - header["lowerBounds"][0]

newAtoms = []

nFixed = 0
for atom in atoms:
    diff = (atom["position"] - systemLengths / 2)
    if(sqrt(diff[0]**2 + diff[1]**2) < cylinderRadius):
        atom["isPositionFixed"] = 0
        if random() < 0.5:
            newAtoms.append(atom)
    else:
        newAtoms.append(atom)
        atom["isPositionFixed"] = 1
        nFixed += 1

newAtoms = asarray(newAtoms, dtype=dataType)
    
print "Fixed", nFixed, "of", len(atoms), "atoms"
saveAtoms(header, lammpsHeader, newAtoms, argv[2])
