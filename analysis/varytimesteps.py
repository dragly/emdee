#from runconfig import run
from pylibconfig import Config
from pylab import *
from sys import argv
from os.path import split, join
from fys4460 import makedirsSilent

defaultConfigFilePath = "1j-energy/default.cfg"
defaultConfigDir, defaultConfigFileName = split(defaultConfigFilePath)

timeSteps = linspace(1e-15, 1e-13, 10)

config = Config()
config.addGroup("", "runConfig")
config.addList("runConfig", "subConfigs")

iRun = 0
for timeStep in timeSteps:
    runName = "timestep%04d" % iRun
    subConfigFileName = join(runName, runName + ".cfg")
    subConfigFilePath = join(defaultConfigDir, subConfigFileName)
    print split(subConfigFilePath)[0]
    makedirsSilent(split(subConfigFilePath)[0])
    subConfig = Config()
    subConfig.readFile(defaultConfigFilePath)
    subConfig.setValue("integrator.timeStep", timeStep)
    saveFileName = subConfig.value("simulation.saveFileName")[0]
    saveFileName = saveFileName.replace("$DATEDIR", "$DATEDIR/" + runName)
    subConfig.setValue("simulation.saveFileName", saveFileName)
    subConfig.writeFile(subConfigFilePath)
    config.appendToList("runConfig.subConfigs", subConfigFileName)
    iRun += 1
    
config.writeFile(join(defaultConfigDir, "timesteps.cfg"))