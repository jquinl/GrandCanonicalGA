#Import these for the GA code to run
import pytest
from GCGA.GCGA import GCGA
from GCGA.FitnessFunction.BaseFitness import GibbsFreeEnergy

from ase import Atoms
import numpy as np
###########################
from ase.io import read
#Calculator
from ase.optimize import BFGS
from ase.calculators.emt import EMT

#Define the fitness fucntion for your atoms it must take in an atoms object as parameter and return a float for the code to work
def fitness_function(atoms)-> float:
    env = 1.0
    ref=read('pt.traj@:')[0].get_potential_energy()    

    au_en =read('gold_bulk.traj').get_potential_energy() / 500.0
    at_num= atoms.get_atomic_numbers()
    pt_num = np.count_nonzero(at_num == 78)
    au_num= np.count_nonzero(at_num == 79)

    fre = atoms.get_potential_energy() - ref- au_num * (au_en) - au_num*env

    return -fre

def simple_fitness_function(atoms)-> float:
    return -atoms.get_potential_energy()

#Or specify an already implemented fitness function such as the gibbs free energy of the system


def test_length():

    base_reference = read('pt.traj@:')[0].get_potential_energy() 
    gold_reference = read('gold_bulk.traj').get_potential_energy() / 500.0
    references = {"Au": gold_reference}
    environment = {"mu_Au": 0.0}

    gfe = GibbsFreeEnergy(variable_reference_energies=references,environmental_variables=environment, base_reference_energy= base_reference)

    #Indicate static part of the system"
    a = 24.0
    slab = Atoms(cell=[a,a,a],
                pbc=True)

    #If slab atoms arent fixed they will be generated  for all structures, but will be relaxed!"

    #---------Generate variable part of the system----------------------"
    #Part of the system that can be relaxed and that varyies in number over the duration of the search"
    variable_types = [Atoms('Pt'),Atoms('Au')]
    variable_range = [[1],[10]]


    #Instantiating of the GCGA object with the selected parameters
    gcga = GCGA(slab,variable_types,variable_range,
    gfe,
    calculator = EMT(),
    starting_candidates_per_stc = 2,population_size=20,steps=50)

    #Calling the run function will initialize the run
    gcga.run()

    final = read("structures.traj@:")
    assert(len(final) == 52 )
