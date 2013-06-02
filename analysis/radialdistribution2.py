# -*- coding: utf-8 -*-
"""
Created on Sat May 25 16:45:24 2013

@author: svenni
"""

from fys4460 import loadAtoms
from itertools import combinations
from pylab import *

header, lammps, atoms = loadAtoms("/home/svenni/scratch/fys4460/tmp/radial/data000800.bin")

combos = asarray(list(combinations(atoms, 2)))

atoms1 = combos[:,0]
atoms2 = combos[:,1]

#distances = atoms2["position"] - atoms1["position"]
    