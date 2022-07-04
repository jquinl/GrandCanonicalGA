#Import these for the GA code to run

from GCGA.GCGA import GCGA
from GCGA.Operations.CrossOperation import CrossOperation as CO
from GCGA.Operations.RandomCandidateGenerator import RandomCandidateGenerator as RCG
from GCGA.Operations.PermutationOperation import PermutationOperation as PM
from GCGA.Operations.RattleOperation import RattleOperation as RT

from GCGA.Operations.RemoveOperation import RemoveOperation as RM
from GCGA.Operations.AddOperation import AddOperation as AD

from GCGA.FitnessFunction.BaseFitness import GibbsFreeEnergy

from ase import Atoms
import numpy as np
###########################
from ase.io import read
#Calculator
from ase.optimize import BFGS
from ase.calculators.emt import EMT

#Define the fitness fucntion for your atoms it must take in an atoms object as parameter and return a float for the code to work"

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


base_reference = read('pt.traj@:')[0].get_potential_energy() 
gold_reference = read('gold_bulk.traj').get_potential_energy() / 500.0
references = {"Au": gold_reference}
environment = {"mu_Au": 1.0}

gfe = GibbsFreeEnergy(variable_reference_energies=references,environmental_variables=environment, base_reference_energy= base_reference)

#Indicate static part of the system"
a = 24.0
slab = Atoms(cell=[a,a,a],
             pbc=True)

#If slab atoms arent fixed they will be generated  for all structures, but will be relaxed!"

#---------Generate variable part of the system----------------------"
#Part of the system that can be relaxed and that varyies in number over the duration of the search"
variable_types = [Atoms('Pt'),Atoms('Au')]
variable_range = [[1],list(range(1,21))]

#Indicate mutations to be performed during the run. They can be passes as an instantiated class, as a type or as a string,
# If passed as type or string, they will be instantiated with default values. A more precisse fine tuning of the run requeres the classes to be preinstantiated when
#  passed into rhe GCGA object"
candidateGenerator = RCG(slab,variable_types,variable_range,max_bond_lenght_multi=2.0)
crossing = CO(slab,variable_types,variable_range,minfrac = 0.2)
permutating = PM(slab,variable_types,variable_range)
removing = RM(slab,variable_types,variable_range)
adding = AD(slab,variable_types,variable_range)
rattling = RT(slab,variable_types,variable_range,n_to_move= 2,rattle_strength=0.1)
# Once instantiated, they have to be added to a tuple. A second tuple containing the chances of occurrence of each mutation if the crossing fails
# (Ideally, the sum of mutations must be equal to 1.0)
# The crossing can also be passed as a complimentary mutation, although its called independently during the execution.
#Supported mutations are  "random", "cross", "add", "remove", "permute" and "rattle"
# The types can be passed as:

#mutations = [candidateGenerator,permutating,rattling,removing,adding]
mutations = [RCG,PM,RT,RM,AD]
#mutations = ["random","permute","rattle","remove","add"]
chances = 0.3

#Instantiating of the GCGA object with the selected parameters
gcga = GCGA(slab,variable_types,variable_range,mutations,
gfe,
calculator = EMT(),
mutation_chance=chances,
crossing_operator=crossing,
initial_structure_generator = candidateGenerator,
starting_population = 20,population_size=20,steps=1000)

#Calling the run function will initialize the run
gcga.run()
