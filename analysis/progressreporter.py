# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 00:20:37 2013

@author: svenni
"""

from pylab import *
from sys import argv
import json
import urllib2
import time
import os
import ConfigParser
import uuid

runName = argv[1]
progressFileName = argv[2]

runID = str(uuid.uuid1())

# UPDATE THESE VARIABLES
baseUrl = "http://compphys.dragly.org"
# END UPDATE VARIABLES

def loadUserNameFromFile():
    global username
    configParser = ConfigParser.ConfigParser()
    configParser.read(configPath)
    if configParser.has_section("General"):
        if configParser.has_option("General", "username"):
            username = configParser.get("General", "username")
            print("Updating username to", username)
    

configPath = "/etc/cpu-usage-reporter.conf"
username = "unnamed" + ("%.0f" % (random() * 100))


# If a username is supplied by argument
if len(argv) > 1:
    username = argv[1]
# If a config file has been set up
elif os.path.exists(configPath):
    loadUserNameFromFile()

print("Using username " + username)

samples = 1
usageSum = 0
if os.path.exists(configPath):
    loadUserNameFromFile()
progress = 0
attempts = 0
while(progress < 1 and attempts < 6):
    try:
        if os.path.exists(progressFileName):
            progressFile = open(progressFileName, "r")
            progress = float(progressFile.readline())
            progressFile.close()
            print "Pushing progress to server", progress
            calledURL = baseUrl + "/wp-content/plugins/run-reporter/submit.php?user=dragly&progress=" + str(progress) + "&runName=" + urllib2.quote(runName) + "&runID=" + urllib2.quote(runID)
            response = urllib2.urlopen(calledURL)
            runData = json.load(response)
        else:
            attempts+=1
    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except:
        print "Something wrong happened... Whatever..."
    if progress < 1:
        time.sleep(2)

print "Progressreporter done!"