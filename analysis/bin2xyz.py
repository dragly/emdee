# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 17:02:54 2013

@author: svenni
"""

from pylab import *
import time
from sys import argv
from fys4460 import dataType, headerType
fileNames = argv[1:]
for fileName in fileNames:
    print "Converting " + fileName
    if not fileName[-4:] == ".bin":
        print "Filename should end with .bin"
        exit()
    t1 = time.time()
    f1 = open(fileName, "rb")
    header = fromfile(f1, dtype=headerType, count=1)
    atoms = fromfile(f1, dtype=dataType)
    
    f = open(fileName.replace(".bin", ".xyz"), "w")
    f.write(str(len(atoms)) + "\n")
    f.write("Some nice comment\n")
    for atom in atoms:
        atomString = ("%s %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %d\n" % (atom["type"], 
                                                                                     atom["positionX"] / 1e-10,atom["positionY"] / 1e-10,atom["positionZ"] / 1e-10,
#                                                                                     atom["positionX"],atom["positionY"],atom["positionZ"],
                                                                                atom["velocityX"],atom["velocityY"],atom["velocityZ"],
                                                                                atom["forceX"],atom["forceY"],atom["forceZ"],
                                                                                atom["potential"], atom["cellID"]))
        f.write(atomString)
            
        
    f.close()
    
    t2 = time.time()
    print "Conversion took " + str(t2 - t1)