# -*- coding: utf-8 -*-
"""
Created on Tue May 21 09:36:03 2013

@author: svenni
"""

from fys4460 import *
from os.path import join
from sys import argv

dataDir1 = argv[1]
dataDir2 = argv[2]

for i in range(1,10):
#    print i
    header1,lammpsHeader1,atoms1 = loadAtoms(join(dataDir1,"data%06d.bin" % i))
    header2,lammpsHeader2,atoms2 = loadAtoms(join(dataDir2,"data%06d.bin" % i))
    
    atoms1 = sort(atoms1, order=("id"))
    atoms2 = sort(atoms2, order=("id"))
    
    forceerror = abs((atoms1["force"] - atoms2["force"]).max() / (atoms1["force"].max() + atoms2["force"].max()))
    if forceerror > 1e-10:
        print "Force:"
        print (atoms1["force"] - atoms2["force"]).max()
        print forceerror
#        raise Exception
    
    poserror = abs((atoms1["position"] - atoms2["position"]).max() / (atoms1["position"].max() + atoms2["position"].max()))
    if poserror > 1e-10:
        print "Position:"
        print (atoms1["position"] - atoms2["position"]).max()
        print poserror
#        raise Exception