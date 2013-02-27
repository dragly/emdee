from pylab import dtype, fromfile
import os
import errno
headerType = dtype([("step", 'int32'), ("time", float), ("timeStep", float), 
                  ("nAtoms", 'int32'), ("temperature", float), ("pressure", float),
                    ("averageDisplacement", float),
                    ("averageSquareDisplacement", float),
                    ("lowerBounds", float, (3,)), ("upperBounds", float, (3,))])
dataType = dtype([("type", "a3"),
                  ("position", float, (3,)),
                  ("velocity", float, (3,)),
                    ("force", float, (3,)),
                    ("displacement", float, (3,)),
                    ("potential", float), 
                    ("cellID", 'int32')])
boltzmannConstant = 1.3806503e-23

def loadAtoms(fileName):
    f = open(fileName, "rb")
    header = fromfile(f, dtype=headerType, count=1)
    atoms = fromfile(f, dtype=dataType)
    f.close()
    return header, atoms
    
    
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