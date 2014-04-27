# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:31:33 2013

@author: svenni
"""
from sys import argv
from fys4460 import loadHeader
from pylab import dtype, fromfile, concatenate
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
    
    lammpsFileName = fileName.replace(".bin", ".lmp")
    runList = "/bin/cat " + lammpsFileName + " "
    nAtomsTotal = header["nAtoms"]
    if nProcessors < 2:
        print "Skipped " + fileName
        continue
    
    for i in range(1, nProcessors):
        runList += (lammpsFileName + (".%04d " % i))
        header2 = loadHeader(fileName + (".%04d" % i))
        nAtomsTotal += header2["nAtoms"]
            
    runList += " > "
    mergedFileName = lammpsFileName.replace(".lmp", ".lmp.merged")
    runList += mergedFileName
    print runList
    subprocess.check_output(runList, shell=True)
    
    header["nProcessors"] = 1
    header["nAtoms"] = nAtomsTotal
    
    for i in range(1, nProcessors):
        os.remove(fileName  + (".%04d" % i))
        os.remove(lammpsFileName  + (".%04d" % i))
    os.rename(mergedFileName, lammpsFileName)
    
    headerFile = open(fileName, "wb")
    header.tofile(headerFile)
    headerFile.close()
    
print "Done"
