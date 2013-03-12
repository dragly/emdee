# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:31:33 2013

@author: svenni
"""
from sys import argv
from fys4460 import loadHeader
from pylab import dtype, fromfile, concatenate
from pylibconfig import Config
from os.path import expanduser, join, split
from glob import glob
import subprocess
import os
import errno

fileNames = argv[1:]
if len(fileNames) < 2:
    fileNames = glob(expanduser(argv[1]))
fileNames.sort()
print fileNames
for fileName in fileNames:
    header = loadHeader(fileName)
    
    nProcessors = header['nProcessors'][0]
    print "Has", nProcessors, "processor(s)"
    
    lammpsFileName = fileName.replace(".bin", ".lmp")
    runList = "/bin/cat " + lammpsFileName + " "
    if nProcessors > 1:
        for i in range(1, nProcessors):
            runList += (lammpsFileName + (".%04d " % i))
    runList += " > "
    runList += lammpsFileName.replace(".lmp", ".lmp.merged")
    print runList
    subprocess.check_output(runList, shell=True)
    print "done"