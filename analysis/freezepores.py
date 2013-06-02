from pylab import *
from fys4460 import loadAtoms, listSaveFileNames, saveAtoms
from sys import argv
from os.path import expanduser
fileName = expanduser(argv[1])
header, lammpsHeader, atoms = loadAtoms(fileName)

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
nPoresSet = [8, 16, 20, 24, 28, 32, 64]
for nPores in nPoresSet:
    while len(balls) <= nPores:
        print "Creating pore", len(balls)
        ball = Ball()
        overlapping = (len(balls) > 0)
        ball.radius = random() * 0.5e-9 + 1e-9
        attempts = 0
    #    while overlapping:
        overlapping = False
        ball.center = random(3) * (header["upperBounds"] - header["lowerBounds"]) + header["lowerBounds"]
        for otherBall in balls:
            for directionalVector in directionalVectors:
                vector = (otherBall.center - ball.center + directionalVector)
                if norm(vector) < (ball.radius +  otherBall.radius):
                    overlapping = True
    
        if attempts > 1e3:
            print "Tried too many times, restarting pore list..."
            balls = []
            overlapping = False
        attempts += 1
        
        balls.append(ball)
     
    iBall = 0
    for ball in balls:
        print "Calculating pores for ball nr.", iBall
        iDirection = 0
        nFixedBall = 0
        for directionalVector in directionalVectors:
            vector = atoms["position"] - ball.center + directionalVector
            dotProducts = sum(vector * vector, axis=1)
            fixedPositions = (dotProducts < ball.radius**2)
            atoms["isPositionFixed"] = ((atoms["isPositionFixed"] + fixedPositions) > 0)
        iBall+=1
    
    nFixed = int(sum(atoms["isPositionFixed"]))
    print "Fixed", nFixed, "of", len(atoms), "atoms"
    outFileName = argv[2]
    outFileName = outFileName.replace(".bin", str(nPores) + ".bin")
    saveAtoms(header, lammpsHeader, atoms, outFileName)