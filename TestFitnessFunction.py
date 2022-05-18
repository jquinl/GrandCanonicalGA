#Import these for the GA code to run

from GCGA.CoreUtils.SubunitAnalysis import SubunitFinder
from ase import Atoms
import numpy as np
###########################
from ase.io import read,write
#Calculator
from ase.calculators.emt import EMT
from GCGA.FitnessFunction.BaseFitness import GibbsFreeEnergy


def fitness_function(atoms)-> float:
    env = 0.5
    ref=read('pt.traj@:')[0].get_potential_energy()    

    au_en =read('gold_bulk.traj').get_potential_energy() / 500.0
    at_num= atoms.get_atomic_numbers()
    pt_num = np.count_nonzero(at_num == 78)
    au_num= np.count_nonzero(at_num == 79)

    fre = atoms.get_potential_energy() - ref- au_num * (au_en) - au_num*env

    return -fre
#Define the fitness fucntion for your atoms it must take in an atoms object as parameter and return a float for the code to work"


base_reference = read('pt.traj@:')[0].get_potential_energy() 
gold_reference = read('gold_bulk.traj').get_potential_energy() / 500.0
references = {"Au": gold_reference}
environment = {"mu_Au": 0.5}

gfe = GibbsFreeEnergy(variable_reference_energies=references,environmental_variables=environment, base_reference_energy= base_reference)

slab = Atoms()

atoms = read("sorted_structures.traj@:")

for a in atoms:

   print("GFE: {}".format(gfe.evaluate(slab,a)))
   print("FF : {}".format(fitness_function(a)))
   if(abs(gfe.evaluate(slab,a) - fitness_function(a)) > 0.000001 ): print("Not Equal")



