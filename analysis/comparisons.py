from fys4460 import loadAtoms
from pylab import *
header1, atoms1 = loadAtoms("/home/svenni/scratch/fys4460/tmp/parallel/normal/data000080.bin")
header2, atoms2 = loadAtoms("/home/svenni/scratch/fys4460/tmp/parallel/testing/data000080.bin")

positions1 = sort(atoms1['position'])
positions2 = sort(atoms2['position'])

print "Biggest difference: ", (positions1 - positions2).max(), (positions1 - positions2).min()
print "Mean absolute difference: ", mean(sqrt((positions1 - positions2)**2))