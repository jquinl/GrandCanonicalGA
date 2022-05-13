#Import these for the GA code to run

from GCGA.CoreUtils.SubunitAnalysis import SubunitFinder
from ase import Atoms
import numpy as np
###########################
from ase.io import read,write
#Calculator
from ase.calculators.emt import EMT
from ase.build import connected_indices as con_ind

#Define the fitness fucntion for your atoms it must take in an atoms object as parameter and return a float for the code to work"
aaaa = Atoms("CO")

for a in aaaa:
    a.number = 200
print(aaaa.numbers)
#print(aaaa.symbols)
a = read("co_test.traj")
b = read("co.traj")
at = SubunitFinder.find_subunits(b,a)
write("reuslt.traj",at)