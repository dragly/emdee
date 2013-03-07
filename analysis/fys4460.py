from pylab import dtype, fromfile, concatenate
import os
import errno
headerType = dtype([("rank", 'int32'),("nProcessors", 'int32'),
                  ("step", 'int32'), ("time", float), ("timeStep", float), 
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
    
    nProcessors = header['nProcessors'][0] 
    
    if nProcessors > 1:
        for i in range(1, nProcessors):
            subFileName = fileName + (".%04d" % i)
            print "Should load " + subFileName
            f2 = open(subFileName, "rb")
            header2 = fromfile(f2, dtype=headerType, count=1)
            atoms2 = fromfile(f2, dtype=dataType)
            f2.close()
            
            atoms = concatenate([atoms, atoms2])
            
            print header2["temperature"][0]
            
            header["nAtoms"][0] += header2["nAtoms"][0]
            header["temperature"][0] += header2["temperature"][0]
            header["pressure"][0] += header2["pressure"][0]
            header["averageDisplacement"][0] += header2["averageDisplacement"][0]
            header["averageSquareDisplacement"][0] += header2["averageSquareDisplacement"][0]
            
        header["temperature"][0] /= nProcessors
        header["pressure"][0] /= nProcessors
        header["averageDisplacement"][0] /= nProcessors
        header["averageSquareDisplacement"][0] /= nProcessors
    
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
    
if __name__ == "__main__":
    header, atoms = loadAtoms("/home/svenni/scratch/fys4460/tmp/parallel/testing/data000000.bin")
