from pylab import *
from pylibconfig import Config
from os.path import expanduser, join, split
from glob import glob
import os
import errno

headerType = dtype([("rank", 'int32'),("nProcessors", 'int32'),
                  ("step", 'int32'), ("time", float), ("timeStep", float), 
                  ("nAtoms", 'int32'), ("temperature", float), ("pressure", float),
                    ("averageDisplacement", float),
                    ("averageSquareDisplacement", float),
                    ("lowerBounds", float, (3,)), ("upperBounds", float, (3,)),
                    ("kineticEnergy", float),
                    ("potentialEnergy", float),])
lammpsHeaderType = dtype([("step", 'int32'),("nAtoms", 'int32'),
                          ("bounds", float, (9,)),
                          ("nColumns", 'int32'),("nChunks", 'int32'),
                          ("chunkLength", 'int32')])
dataType = dtype([("id", float),("type", float),
                  ("position", float, (3,)),
                  ("velocity", float, (3,)),
                    ("force", float, (3,)),
                    ("displacement", float, (3,)),
                    ("potential", float), 
                    ("isPositionFixed", float),
                    ("cellID", float)])
boltzmannConstant = 1.3806503e-23

def listSaveFileNames(configFilePath):
#    saveDir, configFileName = split(configFilePath)

    config = Config()
    config.readFile(configFilePath)
    fileNames = config.value("simulation.saveFileName")[0]
    fileNames = expanduser(fileNames)
    fileNames = glob(fileNames)
    
    fileNames.sort()
    
    return fileNames
    
def saveAtoms(header, lammpsHeader, atoms, fileName):
    fileName = expanduser(fileName)
    header["nProcessors"] = 1
    header["nAtoms"] = len(atoms)
    lammpsHeader["nAtoms"] = len(atoms)
    lammpsHeader["chunkLength"] = len(atoms) * lammpsHeader["nColumns"]
    atoms["position"] *= 1e10 # Convert to LAMMPS units
    
    headerFile = open(fileName, "wb")
    print "Saving " + split(fileName)[1] 
    header.tofile(headerFile)
    headerFile.close()
    lammpsFileName = fileName.replace(".bin", ".lmp")
    print "Saving " +  str(len(atoms)) + " atoms to " + split(lammpsFileName)[1] 
    lammpsFile = open(lammpsFileName, "wb")
    lammpsHeader.tofile(lammpsFile)
    atoms.tofile(lammpsFile)
    lammpsFile.close()
    atoms["position"] *= 1e-10 # Convert back
    
def loadHeader(fileName):
    fileName = expanduser(fileName)
    headerFile = open(fileName, "rb")
    header = fromfile(headerFile, dtype=headerType, count=1)
    headerFile.close()
    
    return header
    

def loadAtoms(fileName):
    fileName = expanduser(fileName)
    header = loadHeader(fileName)
    lammpsFileName = fileName.replace(".bin", ".lmp")
    lammpsFile = open(lammpsFileName, "rb")
    lammpsHeader = fromfile(lammpsFile, dtype=lammpsHeaderType, count=1)
    atoms = fromfile(lammpsFile, dtype=dataType)
    lammpsFile.close()
    
    nProcessors = header['nProcessors'][0]
    print "Has", nProcessors, "processor(s)"
    
    if nProcessors > 1:
        for i in range(1, nProcessors):
            subFileName = fileName + (".%04d" % i)
            subLammpsFileName = subFileName.replace(".bin", ".lmp")
            print "Loading " + split(subFileName)[1]
            lammpsFile2 = open(subLammpsFileName, "rb")
            atoms2 = fromfile(lammpsFile2, dtype=dataType)
            lammpsFile2.close()
            
            atoms = concatenate([atoms, atoms2])
    
    atoms["position"] *= 1e-10
    return header, lammpsHeader, atoms
    
def makedirsSilent(directory):
    try:
        os.makedirs(directory)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(directory):
            pass
        else: 
            raise 
    
def createSymlink(source, destination):
    if os.path.islink(destination):
        os.unlink(destination)
    os.symlink(source, destination)
    
if __name__ == "__main__":
    from sys import argv
    header, lammpsHeader, atoms = loadAtoms(argv[1])
