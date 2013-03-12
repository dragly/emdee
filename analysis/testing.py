# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 11:16:34 2013

@author: svenni
"""

from pylab import *
from fys4460 import loadAtoms, listSaveFileNames, saveAtoms
from sys import argv
fileName = argv[1]
header, atoms = loadAtoms(fileName)