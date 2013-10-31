# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 08:54:28 2013

@author: svenni
"""

from pylab import *
from sys import argv
from glob import glob
from fys4460 import loadHeader
from os.path import expanduser, join, split, isdir, islink
from os import readlink
from pylibconfig import Config
from itertools import chain
import os
import subprocess
#import os
#import h5py
combinations = [["8", "8", "8", "3.0", "3.0"],
                ["8", "8", "14", "3.0", "2.0"],
                ["8", "14", "8", "2.0", "2.0"],
                ["8", "14", "14", "2.0", "3.3"],
                ["14", "8", "14", "2.0", "2.0"],
                ["14", "14", "14", "3.3", "3.3"]]
configFilePaths = argv[1:]
for configFilePath in configFilePaths:
    print "Loading config file"
    configFilePath
    
    if islink(configFilePath):
        configFilePath = readlink(configFilePath)
        
    saveDir, configFileName = split(configFilePath)
    
    config = Config()
    config.readFile(configFilePath)
    fileNames = config.value("simulation.saveFileName")[0]
    fileNames = expanduser(fileNames)
    fileNames = glob(fileNames)
    
    fileNames.sort()
    
    for combination in combinations:
        figure()
        distBins = [0,360]
        totalBins = zeros(distBins[1])
        iFiles = 0
        temperature = 0
        for fileName in fileNames[-20:-1]:
            header = loadHeader(fileName)
            temperature += header["temperature"]
            outFileName = fileName.replace("data", "distances")
            arguments = ["../tools/build-angular-distribution-Qt_5_0_1_gcc_64-Release/angular-distribution"]
            arguments += [fileName, outFileName]
            arguments += combination
            arguments += [str(len(totalBins))]
            process = subprocess.call(arguments)
            outFile = open(outFileName, "rb")
            nBins = fromfile(outFile, dtype='int32', count=1)
            binEdges = fromfile(outFile, dtype=float, count=nBins + 1) 
            binContents = fromfile(outFile, dtype=float, count=nBins)
            outFile.close()
    #        distances = fromfile(outFileName)
            os.remove(outFileName)
            
            #figure("histFigure")
            #distBins = hist(distances, bins=distBins[1])
    #        V = 4./3. * pi * binEdges**3
    #        Vdiff = V[1:] - V[:-1]
            newBins = binContents # / (Vdiff * sum(binContents))
            totalBins = totalBins + newBins
                
            iFiles += 1
        temperature /= iFiles
        temperatureLabel = "%.0f K" % temperature
        totalBins /= iFiles
#        latestFigure = figure("latestFigure")
#        plot(binEdges[:-1] * 360 / (2*pi), newBins, label=temperatureLabel)
#        xlabel("r [nm]")
#        ylabel("g(r)")
#        legend()
#        distributionFigure = figure("distributionFigure")
        title("-".join(combination[:3]).replace("8", "O").replace("14", "Si"))
        plot(binEdges[:-1] * 360 / (2*pi), totalBins, label=temperatureLabel)
        xlabel(r"$\theta$ [deg]")
        ylabel("g(r)")
        grid()
        legend()
        
        
#        latestFigure.savefig(join(saveDir, "angular-distribution-latest-" + " ".join(combination) + ".pdf"))
        savefig(join(saveDir, "angular-distribution-" + "-".join(combination).replace(".", "_") + ".pdf"))
    show()