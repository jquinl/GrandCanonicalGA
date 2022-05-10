#Import these for the GA code to run

from GCGA.FitnessFunction.BaseFitness import GibbsFreeEnergy as BF
from ase import Atoms
import numpy as np
###########################
from ase.io import read,write
#Calculator
from ase.calculators.emt import EMT
from ase.build import connected_indices as con_ind

#Define the fitness fucntion for your atoms it must take in an atoms object as parameter and return a float for the code to work"
a = read("co_test.traj")
b = read("co.traj")
bf = BF({"pO2":5.0},{"Pt":10.0,"Au":5.0})
bf.submolecules_count(b,a)
