from pylab import dtype
headerType = dtype([("step", int), ("time", float), ("timeStep", float), 
                  ("nAtoms", int), ("temperature", float),
                    ("averageDisplacementX", float), ("averageDisplacementY", float), ("averageDisplacementZ", float), 
                    ("averageSquareDisplacement", float)])
dataType = dtype([("type", "a3"),
                  ("positionX", float), ("positionY", float), ("positionZ", float), 
                  ("velocityX", float), ("velocityY", float), ("velocityZ", float), 
                    ("forceX", float), ("forceY", float), ("forceZ", float), 
                    ("potential", float), 
                    ("cellID", 'int32')])
boltzmannConstant = 1.3806503e-23                    