#!/usr/bin/python
from os import system, makedirs
from os.path import isdir
from sys import argv
sourceDir = argv[1]

system("cp -r " + sourceDir + "/systems .")
system("cp -r " + sourceDir + "/testconfig.cfg .")
