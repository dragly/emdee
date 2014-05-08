# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 08:54:28 2013

@author: svenni
"""

from pylab import *
from sys import argv
from glob import glob
from time import time
from fys4460 import loadAtoms, boltzmannConstant
from os.path import expanduser, join, split, isdir, islink
from os import readlink
#from pylibconfig import Config
import os
import subprocess
#import os
#import h5py
#combinations = [["8","8"],
#                ["8","14"],
#                ["14","14"]]
combinations = [["1","1"]]
fileNames = argv[1]
fileNames = glob(fileNames)

fileNames.sort()

binary = "/home/svenni/Dropbox/projects/programming/emdee/build-emdee-hyperon-Desktop_Qt_5_2_1_GCC_64bit-Release/tools/radial-distribution/radial-distribution"
lib_path = "/home/svenni/Dropbox/projects/programming/emdee/build-emdee-hyperon-Desktop_Qt_5_2_1_GCC_64bit-Release/src/libs"
env = dict(os.environ)
env['LD_LIBRARY_PATH'] = lib_path

for combination in combinations:
    distBins = [0,500]
    totalBins = zeros(distBins[1])
    iFiles = 0
    temperature = 0
    
    for fileName in fileNames:
        header, lammps, atoms = loadAtoms(fileName)
        temperature += header["temperature"]
        outFileName = fileName + ".out"
        arguments = [binary, fileName, outFileName] + combination + [str(distBins[1])]
        print arguments
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
    binEdges /= 0.458
    
    plot(binEdges[:-1], totalBins, label=temperatureLabel)
    xlabel(u"r")
    ylabel("g(r)")
    xlim(0, sideLengths[0])
    ylim(0,5.0)
    legend()
    numberDensity = len(atoms) / prod(header["upperBounds"] - header["lowerBounds"])
    totalMass = sum(atoms["type"] == 8)*2.65676264126474e-26  + sum(atoms["type"] == 14) * 4.66370658657455e-26
    grid()
    print "Mass density (g/cm^3): " + str(totalMass * 1e3 / (prod(header["upperBounds"] - header["lowerBounds"]) * 1e6))
    
#        latestFigure.savefig(join(saveDir, "radial-distribution-latest.pdf"))
#    savefig(join(saveDir, "radial-distribution-" + "-".join(combination) + ".pdf"))
show()
