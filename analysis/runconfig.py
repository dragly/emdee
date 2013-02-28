# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 21:52:16 2013

@author: svenni
"""

from pylibconfig import Config
from sys import argv
import subprocess
import sqlite3
import datetime
import os
import errno
from os.path import expanduser, join, split
from git import Repo
from fys4460 import makedirsSilent, createSymlink
import socket

def printChildren(database, lastrowid, config, name):
    for child in config.children(name):
        if len(config.children(child)) > 0:
            printChildren(database, lastrowid, config, child)
        else:
            values = [lastrowid, child, config.value(child)[0]]
            database.execute("INSERT INTO runconfig (runid, name, value) VALUES (?, ?, ?)", values)

def replaceDateDir(config, dateDir, parent = ""):
    for child in config.children(parent):
        if len(config.children(child)) > 0:
            replaceDateDir(config, dateDir, child)
        else:
            oldValue = config.value(child)[0]
            if type(oldValue) == type("abc") and "$DATEDIR" in oldValue:
                newValue = oldValue.replace("$DATEDIR", dateDir)
                config.setValue(child, newValue)

def addRunInformation(config, name, value):
    config.addString("runInformation", name)
    config.setValue("runInformation." + name, value)
    
def parseSaveFile(saveFile, dateDir):
    saveFile = saveFile.replace("$DATEDIR", dateDir)
    saveFile = expanduser(saveFile)
    return saveFile
            
dateDir = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")

def parseAndRun(executable, configFile, runDir):
    # Build run dir name
    
    print "Parsing configuration...\nConfigFile: " + configFile + "\nRunDir: " + runDir
    
    configFileDir, configFileName = split(configFile)
    if runDir == "":
        runDir = join(configFileDir, "runs", dateDir)
    makedirsSilent(runDir)  
    
    config = Config()
    config.readFile(configFile)
    
    if config.exists("runConfig.subConfigs"):
        print "This is a parallel config. Launching subconfigs..."
        for subConfigValue in config.children("runConfig.subConfigs"):
            subConfigFile = join(configFileDir, config.value(subConfigValue)[0])
            subConfigFileDir, subConfigFileName = split(subConfigFile)
            subConfigRunDir = join(runDir, subConfigFileName).replace(".cfg", "")
            parseAndRun(executable, subConfigFile, subConfigRunDir)
#            prefix = subConfigFileName.replace(".cfg", "")
#            subRunDir = join(runDir, prefix)
#            subDateDir = join(dateDir, prefix)
#            makedirsSilent(subRunDir)  
#            subRunConfigFile = join(subRunDir, subConfigFileName)
#            subConfig = Config()
#            print "Loading " + subConfigFile
#            subConfig.readFile(subConfigFile)
#            subConfig.writeFile(subRunConfigFile)
#            
#            # Create symlink to root datedir
#            saveFile = subConfig.value("simulation.saveFileName")[0]
#            saveFile = parseSaveFile(saveFile, dateDir)
#            saveDir = os.path.dirname(saveFile)
#            createSymlink(saveDir, "/tmp/latestsaveroot")            
#            
#            run(executable, subRunConfigFile, subDateDir)
        
        config.writeFile(join(configFileDir, configFileName))
    else:
        print "This is a singular config. Launching it alone..."
        runConfigFile = join(runDir, configFileName)
        config.writeFile(runConfigFile)
        run(executable, runConfigFile, dateDir, runDir)
    

def run(executable, configFile, dateDir, runDir):
    executable = os.path.realpath(executable)
    temp,configFileName = os.path.split(configFile)    
    
    config = Config()
    config.readFile(configFile)
    
    replaceDateDir(config, dateDir)
    
    if(config.exists("runInformation")):
        print "Config cannot contain section runInformation already. It would be overwritten."
        raise Exception
    
    ### Git info ###
    repo = Repo(".")
    head = repo.heads[0]
    commitSHA = getattr(head.commit, "id")
    
    ### Run information ###
    config.addGroup("", "runInformation")
    addRunInformation(config, "hostname", socket.gethostname())
    addRunInformation(config, "timestamp", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    addRunInformation(config, "gitCommitSHA", commitSHA)
    addRunInformation(config, "dateDir", dateDir)
    addRunInformation(config, "runDir", runDir)
    addRunInformation(config, "configFile", configFile)
    
    config.writeFile(configFile)

    # Input stuff to database
    database = sqlite3.connect("test.db")
    values = [executable, runDir]
    insertID = database.execute("INSERT INTO run (executable, dir) VALUES (?, ?)", values)
    printChildren(database, insertID.lastrowid, config, "")
    database.commit()
    database.close()
    
    # Create symlinks
    createSymlink(os.getcwd() + "/" + runDir, "/tmp/latest")
    saveFile = config.value("simulation.saveFileName")[0]
    saveFile = expanduser(saveFile)
    saveDir = os.path.dirname(saveFile)
    createSymlink(saveDir, "/tmp/latestsavedir")
    
    logFilePath = join(runDir, "run_" + configFileName.replace(".cfg", "") + ".log")
    f = open(logFilePath, "w")
    print "Executable: " + executable + "\nConfig: " + configFileName + "\nRunDir: " + runDir + "\nLogFile: " + logFilePath
    print "Starting..."
    process = subprocess.Popen([executable, configFileName], cwd=runDir, stdout=f, stderr=subprocess.STDOUT)
    process.wait()
    f.close()
        
if __name__ == "__main__":
    executable = argv[1]
    configFiles = argv[2]
    
    parseAndRun(executable, configFiles, "")
