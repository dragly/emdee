# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 15:20:49 2013

@author: svenni
"""

#from runconfig import run
from pylibconfig import Config
from pylab import *
from sys import argv
from os.path import split, join
from fys4460 import makedirsSilent

defaultConfigFilePath = "1k-temperature/default.cfg"
defaultConfigDir, defaultConfigFileName = split(defaultConfigFilePath)

systemSizes = [6, 7, 8, 9, 10, 11, 12]

config = Config()
config.addGroup("", "runConfig")
config.addList("runConfig", "subConfigs")

iRun = 0
for systemSize in systemSizes:
    runName = "systemsize%04d" % systemSize
    subConfigFileName = join(runName, runName + ".cfg")
    subConfigFilePath = join(defaultConfigDir, subConfigFileName)
    print split(subConfigFilePath)[0]
    makedirsSilent(split(subConfigFilePath)[0])
    subConfig = Config()
    subConfig.readFile(defaultConfigFilePath)
    subConfig.setValue("initialization.[0].nCells", systemSize)
    saveFileName = subConfig.value("simulation.saveFileName")[0]
    saveFileName = saveFileName.replace("$DATEDIR", "$DATEDIR/" + runName)
    subConfig.setValue("simulation.saveFileName", saveFileName)
    subConfig.writeFile(subConfigFilePath)
    config.appendToList("runConfig.subConfigs", subConfigFileName)
    iRun += 1
    
config.writeFile(join(defaultConfigDir, "systemsizes.cfg"))