# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 12:15:18 2013

@author: svenni
"""

#from runconfig import run
from pylibconfig import Config
from pylab import *
from sys import argv
from os.path import split, join
from fys4460 import makedirsSilent

defaultConfigFilePath = "1-equilibrated-fast-zoom/default-start.cfg"
defaultEndConfigFilePath = "1-equilibrated-fast-zoom/default-end.cfg"
defaultConfigDir, defaultConfigFileName = split(defaultConfigFilePath)
defaultEndConfigDir, defaultEndConfigFileName = split(defaultEndConfigFilePath)

temperatures = linspace(100, 600, 40)

config = Config()
config.addGroup("", "runConfig")
config.addList("runConfig", "subConfigs")

iRun = 0
for temperature in temperatures:
    runName = "temperature%04d" % iRun
    subConfigFileName = join(runName, "step0.cfg")
    subConfigEndFileName = join(runName, "step1.cfg")
    subConfigFilePath = join(defaultConfigDir, subConfigFileName)
    subConfigEndFilePath = join(defaultConfigDir, subConfigEndFileName)
    makedirsSilent(split(subConfigFilePath)[0])
    subConfigStart = Config()
    subConfigStart.readFile(defaultConfigFilePath)
    oldSaveFile = subConfigStart.value("simulation.saveFileName")[0]
    newSaveFile = oldSaveFile.replace("$DATEDIR", "$DATEDIR/" + runName)
    subConfigStart.setValue("simulation.saveFileName", newSaveFile)
    subConfigStart.setValue("initialization.[1].initialTemperature", temperature)
    subConfigStart.setValue("modifiers.[0].targetTemperature", temperature)
    subConfigStart.writeFile(subConfigFilePath)
    
    subConfigEnd = Config()
    subConfigEnd.readFile(defaultEndConfigFilePath)
    oldSaveFile = subConfigEnd.value("simulation.saveFileName")[0]
    newSaveFile = oldSaveFile.replace("$DATEDIR", "$DATEDIR/" + runName)
    subConfigEnd.setValue("simulation.saveFileName", newSaveFile)
    
    oldLoadFile = subConfigEnd.value("initialization.[0].fileName")[0]
    newLoadFile = oldLoadFile.replace("$DATEDIR", "$DATEDIR/" + runName)
    subConfigEnd.setValue("initialization.[0].fileName", newLoadFile)
    subConfigEnd.writeFile(subConfigEndFilePath)
    
    stepConfig = Config()
    stepConfig.addGroup("", "runConfig")
    stepConfig.addList("runConfig", "subConfigs")
    stepConfigFileName = join(defaultConfigDir, runName, runName + ".cfg")
    stepConfig.appendToList("runConfig.subConfigs", split(subConfigFileName)[1])
    stepConfig.appendToList("runConfig.subConfigs", split(subConfigEndFileName)[1])
    stepConfig.writeFile(stepConfigFileName)
    config.appendToList("runConfig.subConfigs", join(runName, runName + ".cfg"))
    iRun += 1
    
config.writeFile(join(defaultConfigDir, "temperatures.cfg"))