# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 17:02:54 2013

@author: svenni
"""

from pylab import *
import time
from sys import argv
from fys4460 import dataType, headerType, loadAtoms, loadHeader
from glob import glob
fileNames = argv[1:]
if len(fileNames) == 1:
    fileNames = glob(fileNames[0])
    
for fileName in fileNames:
    print "Converting " + fileName
    if not fileName[-4:] == ".bin":
        print "Filename should end with .bin"
        exit()
    t1 = time.time()
    
    header, lammpsHeader, atoms = loadAtoms(fileName)
    
    f = open(fileName.replace(".bin", ".xyz"), "w")
    f.write(str(len(atoms)) + "\n")
    f.write("Some nice comment\n")
    for atom in atoms:
        atomString = ("%s %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %d\n" % (atom["type"], 
                                                                                     atom["position"][0] / 1e-10,atom["position"][1] / 1e-10,atom["position"][2] / 1e-10,
#                                                                                     atom["positionX"],atom["positionY"],atom["positionZ"],
                                                                                atom["velocity"][0], atom["velocity"][1], atom["velocity"][2],
                                                                                atom["force"][0],atom["force"][1],atom["force"][2],
                                                                                atom["potential"], atom["cellID"]))
        f.write(atomString)
            
        
    f.close()
    
    t2 = time.time()
    print "Conversion took " + str(t2 - t1)