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

defaultConfigFilePath1 = "1-push-default.cfg"
defaultConfigDir, defaultConfigFileName1 = split(defaultConfigFilePath1)
defaultConfigFilePath2 = "2-collect-default.cfg"
defaultConfigDir2, defaultConfigFileName2 = split(defaultConfigFilePath2)

systemSizes = [8, 16, 20, 24, 28, 32, 64]

config = Config()
config.addGroup("", "runConfig")
config.addList("runConfig", "subConfigs")

iRun = 0
for systemSize in systemSizes:
    runName = "frozenpores%d" % systemSize
    subConfigFileName = join(runName, "1-" + runName + "-push.cfg")
    subConfigFilePath = join(defaultConfigDir, subConfigFileName)
    subConfigFileName2 = join(runName, "2-" + runName + "-collect.cfg")
    subConfigFilePath2 = join(defaultConfigDir, subConfigFileName2)
    makedirsSilent(split(subConfigFilePath)[0])
    
    subConfig = Config()
    subConfig.readFile(defaultConfigFilePath1)
    saveFileName = subConfig.value("simulation.saveFileName")[0]
    print saveFileName
    saveFileName = saveFileName.replace("$NPORES", str(systemSize))
    subConfig.setValue("simulation.saveFileName", saveFileName)
    fileName = subConfig.value("initialization.[0].fileName")[0]
    fileName = fileName.replace("$NPORES", str(systemSize))
    subConfig.setValue("initialization.[0].fileName", fileName)
    subConfig.writeFile(subConfigFilePath)
    
    subConfig2 = Config()
    subConfig2.readFile(defaultConfigFilePath2)
    saveFileName2 = subConfig2.value("simulation.saveFileName")[0]
    saveFileName2 = saveFileName2.replace("$NPORES", str(systemSize))
    subConfig2.setValue("simulation.saveFileName", saveFileName2)
    fileName2 = subConfig2.value("initialization.[0].fileName")[0]
    fileName2 = fileName2.replace("$NPORES", str(systemSize))
    subConfig2.setValue("initialization.[0].fileName", fileName2)
    subConfig2.writeFile(subConfigFilePath2)
    
    stepConfig = Config()
    stepConfig.addGroup("", "runConfig")
    stepConfig.addList("runConfig", "subConfigs")
    stepConfigFileName = join(defaultConfigDir, runName, runName + ".cfg")
    stepConfig.appendToList("runConfig.subConfigs", split(subConfigFileName)[1])
    stepConfig.appendToList("runConfig.subConfigs", split(subConfigFileName2)[1])
    stepConfig.writeFile(stepConfigFileName)
    config.appendToList("runConfig.subConfigs", join(runName, runName + ".cfg"))
    iRun += 1
    
config.writeFile(join(defaultConfigDir, "permeability.cfg"))