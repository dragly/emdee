#!/usr/bin/python
from os import system
from sys import argv
sourceDir = argv[1]

system("cp -r " + sourceDir + "/systems .")
system("cp -r " + sourceDir + "/testconfig.cfg .")