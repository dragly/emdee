# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 15:33:40 2013

@author: svenni
"""

from pylab import *
from pylibconfig import Config
from os import listdir
from os.path import join, split, expanduser, isdir
from sys import argv
from fys4460 import loadAtoms
from glob import glob


runDir = argv[1]
saveDir = split(argv[1])[0]
    
temperatureDirs = listdir(runDir)
temperatureDirs.sort()
saveFileNamesList = []
nFilenames = 0
for temperatureDir in temperatureDirs:
    if not isdir(join(runDir, temperatureDir)):
        continue
    print "Loading dir " + temperatureDir
    config = Config()
    configFileName = join(runDir, temperatureDir, "step1", "step1.cfg")
    print "Loading config " + configFileName
    config.readFile(configFileName)
    saveFileNames = config.value("simulation.saveFileName")[0]
    saveFileNames = expanduser(saveFileNames)
    saveFileNames = glob(saveFileNames)
    saveFileNames.sort()
    saveFileNamesSelection = saveFileNames[300:399]
    saveFileNamesList.append(saveFileNamesSelection)
    nFilenames += len(saveFileNamesSelection)

pressures = zeros(nFilenames)
temperatures = zeros(nFilenames)

pressuresAverage = zeros(len(saveFileNamesList))
temperaturesAverage = zeros(len(saveFileNamesList))
i = 0
avgi = 0
for fileNames in saveFileNamesList:
    starti = i
    for fileName in fileNames:
        print "Loading file\n" + fileName
        header, atoms = loadAtoms(fileName)
        pressures[i] = header["pressure"][0]
        temperatures[i] = header["temperature"][0]
        i += 1
    pressuresAverage[avgi] = mean(pressures[starti:i - 1])
    temperaturesAverage[avgi] = mean(temperatures[starti:i - 1])
    avgi += 1
    
figure()
scatter(temperatures, pressures)
xlabel(r"temperature [K]")
ylabel(r"pressure [Pa]")
grid()
savefig(saveDir + "/pressure-vs-temperature.pdf")
#show()
figure()
plot(temperaturesAverage, pressuresAverage, "-o")
xlabel(r"temperature [K]")
ylabel(r"pressure [Pa]")
grid()
savefig(saveDir + "/pressure-vs-temperature-mean.pdf")
show()