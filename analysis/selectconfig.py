# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 16:54:02 2013

@author: svenni
"""

from sys import argv
from fys4460 import createSymlink
from pylibconfig import Config
from os.path import join, split, expanduser, isdir
from os import listdir

configFilePath = argv[1]

if isdir(configFilePath):
    foundConfig = False    
    for fileName in listdir(configFilePath):
        if fileName[-4:] == ".cfg":
            configFilePath = join(configFilePath, fileName)
            foundConfig = True
            break
    if not foundConfig:
        print "Could not find config..."
        raise FileIOException

print "Loading config " + configFilePath + "\n"

config = Config()
config.readFile(configFilePath)

symlinkNameConfig = "/tmp/selectedconfig.cfg"
print "Target: " + configFilePath
print "Name: " + symlinkNameConfig
createSymlink(configFilePath, symlinkNameConfig)

saveFilePath, saveFileName = split(config.value("simulation.saveFileName")[0])
saveFilePath = expanduser(saveFilePath)
symLinkNameSave = "/tmp/selectedsave"
print "Target: " + saveFilePath
print "Name: " + symLinkNameSave
createSymlink(saveFilePath, symLinkNameSave)