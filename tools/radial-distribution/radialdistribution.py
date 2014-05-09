# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 08:54:28 2013

@author: svenni
"""

sys.path.append("../../analysis")
from pylab import *
from sys import argv
from glob import glob
import sys
from time import time
from fys4460 import loadAtoms, boltzmannConstant
from os.path import expanduser, join, split, isdir, islink, exists
from os import readlink
import shutil
from argparse import ArgumentParser

import os
import subprocess

parser = ArgumentParser()
parser.add_argument("--lammps_files", required=True)
parser.add_argument("--id", default="tmp")
args = parser.parse_args()

combination = ["1","1"]
fileNames = args.lammps_files
fileNames = glob(fileNames)

fileNames.sort()

current_path = os.path.dirname(os.path.realpath(__file__))
build_path = os.path.abspath(os.path.join(current_path, "..", "..", "..", "build-radial-distribution"))
project_path = os.path.abspath(os.path.join(current_path, "..", ".."))

print "Building in:\n", build_path

if os.path.exists(build_path):
    shutil.rmtree(build_path)

os.makedirs(build_path)

subprocess.call(["qmake", project_path, "CONFIG+=nogui", "CONFIG+=mpi", "CONFIG+=notests", "CONFIG+=noapp"], cwd=build_path)
subprocess.call(["make", "-j", "8"], cwd=build_path)

staterunner_path = os.path.join(build_path, "tools", "radial-distribution")
lib_path = os.path.abspath(os.path.join(build_path, "src", "libs"))

env = dict(os.environ)
env['LD_LIBRARY_PATH'] = lib_path
binary = os.path.join(staterunner_path, "radial-distribution")

distBins = [0,500]
totalBins = zeros(distBins[1])
iFiles = 0
temperature = 0

for fileName in fileNames:
    header, lammps, atoms = loadAtoms(fileName)
    temperature += header["temperature"]
    outFileName = fileName + ".out"
    arguments = [binary, fileName, outFileName] + combination + [str(distBins[1])]
    process = subprocess.call(arguments, env=env)
    outFile = open(outFileName, "rb")
    nBins = fromfile(outFile, dtype='int32', count=1)
    binEdges = fromfile(outFile, dtype=float, count=nBins + 1) 
    binContents = fromfile(outFile, dtype=float, count=nBins)
    outFile.close()
#        distances = fromfile(outFileName)
    os.remove(outFileName)
    
    binEdges *= 1e10
    sideLengths = 1e10*(header["upperBounds"] - header["lowerBounds"])[0]
    atoms_count = len(atoms)
    pair_count = atoms_count * (atoms_count-1) / 2.0
    system_number_density = pair_count / prod(sideLengths)
    
    #figure("histFigure")
    #distBins = hist(distances, bins=distBins[1])
    V = 4./3. * pi * binEdges**3
    Vdiff = V[1:] - V[:-1]

    slice_number_densities = binContents / Vdiff
    newBins = slice_number_densities / system_number_density
    
    #newBins = binContents / Vdiff
    totalBins = totalBins + newBins
        
    iFiles += 1
temperature /= iFiles
temperatureLabel = "%.0f K" % temperature
totalBins /= iFiles
#        latestFigure = figure("latestFigure")
#        plot(binEdges[:-1] * 1e10, newBins, label=temperatureLabel)
#        xlabel(u"r [Å]")
#        ylabel("g(r)")
#        legend()
figure()
title("-".join(combination))

# Change of units
#    binEdges /= 0.458

plot(binEdges[:-1], totalBins, label=temperatureLabel)
xlabel(u"r")
ylabel("g(r)")
xlim(0, min(sideLengths) / 2.0)
ylim(0,5.0)
legend()
numberDensity = len(atoms) / prod(header["upperBounds"] - header["lowerBounds"])
grid()

if args.id != "tmp":
    from sumatra.projects import load_project
    from os import makedirs
    output_dir = join(abspath(load_project().data_store.root), args.id)
    if not exists(output_dir):
        makedirs(output_dir)
    savefig(join(output_dir, "radial-distribution.pdf"))
    savefig(join(output_dir, "radial-distribution.png"))
    savetxt(join(output_dir, "radial-distribution.dat"), array([binEdges[:-1], totalBins]).transpose())
else:
    show()
