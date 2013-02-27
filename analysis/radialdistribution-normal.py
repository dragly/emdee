# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 08:54:28 2013

@author: svenni
"""

from pylab import *
from sys import argv
from glob import glob
from time import time
from fys4460 import loadAtoms, boltzmannConstant
from os.path import expanduser, join, split
from pylibconfig import Config
import os
import subprocess
#import os
#import h5py
configFilePaths = argv[1:]
for configFilePath in configFilePaths:
    print "Loading config file"
    configFilePath
    saveDir, configFileName = split(configFilePath)
    
    config = Config()
    config.readFile(configFilePath)
    fileNames = config.value("simulation.saveFileName")[0]
    fileNames = expanduser(fileNames)
    fileNames = glob(fileNames)
    
    fileNames.sort()
    distBins = [0,1000]
    totalBins = zeros(distBins[1])
    iFiles = 0
    temperature = 0
    for fileName in fileNames[400:1400:100]:
        header, atoms = loadAtoms(fileName)
        temperature += header["temperature"]
        outFileName = fileName.replace("data", "distances")
        process = subprocess.call(["../tools/radial-distribution-build-Desktop_Qt_5_0_1_GCC_64bit-Release/radial-distribution", fileName, outFileName])
        distances = fromfile(outFileName)
        os.remove(outFileName)
        
        figure("histFigure")
        distBins = hist(distances, bins=distBins[1])
        V = 4./3. * pi * distBins[1]**3
        Vdiff = V[1:] - V[:-1]
        newBins = distBins[0] / (Vdiff * len(distances))
        totalBins = totalBins + newBins
        
        iFiles += 1
    temperature /= iFiles
    temperatureLabel = "%.0f K" % temperature
    totalBins /= iFiles
    latestFigure = figure("latestFigure")
    plot(distBins[1][1:] * 1e9, newBins, label=temperatureLabel)
    xlabel("r [nm]")
    ylabel("g(r)")
    legend()
    distributionFigure = figure("distributionFigure")
    plot(distBins[1][1:] * 1e9, totalBins, label=temperatureLabel)
    xlabel("r [nm]")
    ylabel("g(r)")
    legend()
    
    
latestFigure.savefig(join(saveDir, "..", "..", "radial-distribution-latest.pdf"))
distributionFigure.savefig(join(saveDir, "..", "..", "radial-distribution.pdf"))
show()