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
from git import Repo

executable = argv[1]
configFile = argv[2]

rootRunDir = "runs"

def printChildren(database, lastrowid, config, name):
    for child in config.children(name):
        if len(config.children(child)) > 0:
            printChildren(database, lastrowid, config, child)
        else:
            values = [lastrowid, child, config.value(child)[0]]
            database.execute("INSERT INTO runconfig (runid, name, value) VALUES (?, ?, ?)", values)

def addRunInformation(config, name, value):
    config.addString("runInformation", name)
    config.setValue("runInformation." + name, value)
        
# Build run dir name
dateDir = datetime.datetime.now().strftime("%Y/%m/%d/%H%M%S_%s/")
runDir = os.getcwd() + "/" + rootRunDir + "/" + dateDir

try:
    os.makedirs(runDir)
except OSError as exc: # Python >2.5
    if exc.errno == errno.EEXIST and os.path.isdir(runDir):
        pass
    else: 
        raise    

config = Config()
config.readFile(configFile)
if(config.exists("runInformation")):
    print "Config cannot contain section runInformation already. It would be overwritten."
    raise Exception


### Git info ###
repo = Repo(".")
head = repo.heads[0]
commitSHA = getattr(head.commit, "id")

### Run information ###
config.addGroup("", "runInformation")
addRunInformation(config, "timestamp", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
addRunInformation(config, "gitCommitSHA", commitSHA)
addRunInformation(config, "dateDir", dateDir)
addRunInformation(config, "runDir", runDir)

database = sqlite3.connect("test.db")

values = [executable, runDir]

insertID = database.execute("INSERT INTO run (executable, dir) VALUES (?, ?)", values)

printChildren(database, insertID.lastrowid, config, "")

database.commit()
database.close()

output = subprocess.check_output([executable, configFile], cwd=runDir)

f = open(runDir + "run.log", "w")
f.write(output)
f.close()

config.writeFile(runDir + "config.cfg")

print runDir