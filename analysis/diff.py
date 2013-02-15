# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 19:18:44 2013

@author: svenni
"""

from fys4460 import *
from pylab import *
from sys import argv
from glob import glob

fileNames = glob("/tmp/mdout/*.bin")
fileNames2 = glob("/tmp/mdoutmirror/*.bin")
fileNames.sort()
fileNames2.sort()

for i in range(len(fileNames)):
    fileName = fileNames[i]
    fileName2 = fileNames2[i]
    print fileName
    atoms = fromfile(fileName, dtype=dataType)
    atoms2 = fromfile(fileName2, dtype=dataType)
    
    diff = atoms["positionX"] - atoms2["positionX"]
    print diff