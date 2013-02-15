from pylab import dtype
headerType = dtype([("step", 'int32'), ("time", float), ("timeStep", float), 
                  ("nAtoms", 'int32'), ("temperature", float),
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
