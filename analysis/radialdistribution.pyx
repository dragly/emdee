# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 08:54:28 2013

@author: svenni
"""
from sys import argv
from glob import glob
from time import time
from fys4460 import loadAtoms, boltzmannConstant
from os.path import expanduser, join, split
from pylibconfig import Config
from numpy import array, zeros, dot, inf, ndarray
cimport cython
from numpy cimport *
#import os
#import h5py
@cython.boundscheck(False)
def calculateDistances(ndarray[double, ndim=2] positions, ndarray[double, ndim=1] lowerBounds, ndarray[double, ndim=1] upperBounds):
    cdef ndarray[double, ndim=2] directionalVectors = zeros((27, 3), dtype=float)
    cdef int counter = 0
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):    
                 directionalVectors[counter] = (upperBounds - lowerBounds) * array([i, j, k]) 
                 counter += 1
    cdef ndarray[double, ndim=2] distances = zeros((int(len(positions) * (len(positions) / 2.0)), 3), dtype=float)
    cdef ndarray[double, ndim=1] distance = zeros(3, dtype=float)
    cdef ndarray[double, ndim=1] shortestDistance = zeros(3, dtype=float)
    counter = 0
    cdef unsigned int nPositions = len(positions)
    cdef unsigned int nDirectionalVectors = len(directionalVectors)
    for i in range(nPositions):
        for j in range(i + 1, nPositions):
            shortestDistance = array([inf, inf, inf])
            for k in range(nDirectionalVectors):
                distance = positions[j] + directionalVectors[k] - positions[i]
#                if dot(distance, distance) < dot(shortestDistance,shortestDistance):
#                    shortestDistance = distance
            distance = positions[j] - positions[i]  
            distances[counter] = distance
            counter += 1
        print i
    return distances

def loadFiles(configFilePath):
    saveDir, configFileName = split(configFilePath)
    
    config = Config()
    config.readFile(configFilePath)
    fileNames = config.value("simulation.saveFileName")[0]
    fileNames = expanduser(fileNames)
    fileNames = glob(fileNames)
    
    fileNames.sort()
    
    for fileName in fileNames[400:600]:
        print "Loading file\n" + fileName
        header, atoms = loadAtoms(fileName)
    
        lowerBounds = header["lowerBounds"]    
        upperBounds = header["upperBounds"]
        
        
        positions = atoms["position"]
        distances = calculateDistances(positions, lowerBounds[0], upperBounds[0])
        
        break