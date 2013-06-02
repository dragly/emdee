# -*- coding: utf-8 -*-
"""
Created on Wed May 22 11:54:14 2013

@author: svenni
"""

from fys4460 import *
from pylab import *
from sys import argv

#sideLengths = ones(3) * 7.166e-10 * (1/2.)**(1/3.)
sideLengths = ones(3) * 7.166e-10 * (1/2.)**(1/3.)
n = ones(3,int) * int(argv[1])
boxLengths = n * sideLengths

#test1, test2, test3 = loadAtoms("/home/svenni/scratch/fys4460/tmp/banana/data000000.bin")

header = zeros(1, dtype=headerType)
lammpsHeader = zeros(1, dtype=lammpsHeaderType)

positions1 = array([
                    [0,0,0],
                   [0.5 , 0.5 , 0],
                    [0.5 , 0 , 0.5], 
                    [0 , 0.5 , 0.5]#,
                    ])
                    
positions2 = array([
                    [0.25 , 0.25 , 0.25],
                    [0.75 , 0.75 , 0.25],
                    [0.75 , 0.25 , 0.75],
                    [0.25 , 0.75 , 0.75],
                [0.5,0.5,0.5],
               [0, 0, 0.5],
                [0,0.5,0],
                [0.5,0,0]
#                [0.75,0.75,0.75],
#                [0.25,0.25,0.75],
#                [0.25,0.75,0.25],
#                [0.75,0.25,0.25]
            ])
#positions1 = array([
#                    [0,0,0],
#                   [0.5 , 0.5 , 0],
#                    [0.5 , 0 , 0.5], 
#                    [0 , 0.5 , 0.5],
##                    [0.25 , 0.25 , 0.25],
##                    [0.75 , 0.75 , 0.25],
##                    [0.75 , 0.25 , 0.75],
##                    [0.25 , 0.75 , 0.75]
#                    ])
#                    
#positions2 = array([
#                    [0.25 , 0.25 , 0.25],
#                    [0.75 , 0.75 , 0.25],
#                    [0.75 , 0.25 , 0.75],
#                    [0.25 , 0.75 , 0.75]
#            ])

atoms = []
idCounter = 1
totalMass = 0
for i in range(n[0]):
    for j in range(n[1]):
        for k in range(n[2]):
            offset = array([i,j,k]) * sideLengths
#            newAtoms[0]["position"][:] = array([0 , 0 , 0]) * sideLength
#            newAtoms[1]["position"][:] = array([i*3 - 1,  j*3    , k*3 - 1/sqrt(2)]) * sideLength / 3
#            newAtoms[2]["position"][:] = array([i*3    ,  j*3 + 1, k*3 + 1/sqrt(2)]) * sideLength / 3
#            newAtoms[3]["position"][:] = array([i*3    ,  j*3 - 1, k*3 + 1/sqrt(2)]) * sideLength / 3
            newAtoms = zeros(len(positions1), dtype=dataType)
            for a in range(len(positions1)):
                newAtoms[a]["position"][:] = positions1[a] * sideLengths + offset
                newAtoms[a]["type"] = 14
                newAtoms[a]["id"] = idCounter
                atoms.append(newAtoms[a])
                totalMass += 4.66370658657455e-26
                idCounter += 1
            newAtoms = zeros(len(positions2), dtype=dataType)
            for a in range(len(positions2)):
                newAtoms[a]["position"][:] = positions2[a] * sideLengths + offset
                newAtoms[a]["type"] = 8
                newAtoms[a]["id"] = idCounter
                atoms.append(newAtoms[a])
                totalMass += 2.65676264126474e-26
                idCounter += 1

atoms = asarray(atoms, dtype=dataType)
numberDensity = len(atoms) / prod(boxLengths)
print "Mass density (g/cm^3): " + str(totalMass * 1e3 / (prod(boxLengths) * 1e6))
print "Number density (10^22/cm^3): " + str(numberDensity / (1e22 * 1e6))
header["upperBounds"][:] = boxLengths
#obeyBoundaries(header, atoms)
saveFileName = argv[2]
saveFileDir = os.path.dirname(saveFileName)
if not os.path.exists(saveFileDir):
    os.makedirs(saveFileDir)
saveAtoms(header, lammpsHeader, atoms, saveFileName)