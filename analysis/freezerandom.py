from pylab import *
from fys4460 import loadAtoms, listSaveFileNames, saveAtoms
from sys import argv
fileName = argv[1]
header, atoms = loadAtoms(fileName)

class Ball:
    radius = 0
    center = array([0,0,0])
    
directionalVectors = []

systemLengths = header["upperBounds"] - header["lowerBounds"]

for i in range(-1,2):
    for j in range(-1,2):
        for k in range(-1,2):
            directionalVector = systemLengths * array([i,j,k])
            directionalVectors.append(directionalVector)

balls = []

for i in range(0,10):
    ball = Ball()
    ball.center = random(3) * (header["upperBounds"] - header["lowerBounds"]) + header["lowerBounds"]
    ball.radius = (header["upperBounds"] - header["lowerBounds"]).min() * 0.2
    
    balls.append(ball)
 
nFixed = 0
for atom in atoms:    
    for ball in balls:
        for directionalVector in directionalVectors:
            if(norm(atom["position"] - ball.center + directionalVector) < ball.radius):
                atom["isPositionFixed"] = True
                nFixed += 1
    
print "Fixed", nFixed, "of", len(atoms), "atoms"
saveAtoms(header, atoms, argv[2])