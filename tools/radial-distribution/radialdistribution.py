# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 08:54:28 2013

@author: svenni
"""

import sys
sys.path.append("../../analysis")
from pylab import *
from glob import glob
from fys4460 import loadAtoms, boltzmannConstant
from os.path import expanduser, join, split, isdir, islink, exists, abspath, dirname
from os import readlink
import shutil
from argparse import ArgumentParser

import os
import subprocess

parser = ArgumentParser()
parser.add_argument("--lammps_files", nargs="+", required=True)
parser.add_argument("-n", "--bins", default=200)
parser.add_argument("-m", "--max", default=inf)
parser.add_argument("--id", default="tmp")
args = parser.parse_args()

combination = ["1","1"]
fileNames = args.lammps_files
if len(fileNames) == 1:
    fileNames = glob(fileNames[0])

fileNames.sort()

current_path = dirname(os.path.realpath(__file__))
build_path = abspath(os.path.join(current_path, "..", "..", "..", "build-radial-distribution"))
project_path = abspath(os.path.join(current_path, "..", ".."))

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

distBins = [0,args.bins]
totalBins = zeros(distBins[1])
iFiles = 0
temperature = 0
pressure = 0
total_count_bins = zeros(distBins[1])

if len(fileNames) < 1:
    raise Exception("Found no files matching input")

for fileName in fileNames:
    header, lammps, atoms = loadAtoms(fileName)
    temperature += header["temperature"]
    pressure += header["pressure"]
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

    total_count_bins += binContents
    slice_number_densities = binContents / Vdiff
    newBins = slice_number_densities / system_number_density
    
    #newBins = binContents / Vdiff
    totalBins = totalBins + newBins
        
    iFiles += 1
    
temperature /= iFiles
pressure /= iFiles
totalBins /= iFiles
total_count_bins /= iFiles

binEdges /= 0.52917721092
sideLengths /= 0.52917721092

dissociation_limit = 1.832 # Angstrom

pairs_below_dissociation_limit = sum(total_count_bins[(where(binEdges < dissociation_limit))])

print "----"
print "Mean temperature: ", temperature
print "Mean pressure: ", pressure
print "Pairs below dissociation limit: ", pairs_below_dissociation_limit
print "Atoms above dissociation limit: ", atoms_count - 2*pairs_below_dissociation_limit
print "----"

figure()
plot(binEdges[:-1], totalBins)
xlabel(r"$r\ [a_0]$")
ylabel(r"$g(r)$")
max_bin_size = float(args.max)
if max_bin_size > min(sideLengths) / 2.0:
    max_bin_size = min(sideLengths) / 2.0
else:
    ax = gca()
    ax.set_xticks(range(0,int(max_bin_size+1),1))

xlim(0,max_bin_size)
ylim(0,5.0)

if args.id != "tmp":
    from sumatra.projects import load_project
    from os import makedirs
    output_dir = join(abspath(load_project().data_store.root), args.id)
    if not exists(output_dir):
        makedirs(output_dir)
    savefig(join(output_dir, "radial-distribution.pdf"))
    savefig(join(output_dir, "radial-distribution.png"))
    savetxt(join(output_dir, "radial-distribution.dat"), array([binEdges[:-1], totalBins]).transpose())
    print "Data saved to", output_dir
else:
    show()
