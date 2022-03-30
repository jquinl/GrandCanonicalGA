#Import these for the GA code to run
from GCGA.CoreUtils.StartingCandidateGenerator import StartingCandidateGenerator as SCG
from GCGA.CoreUtils.DataBaseInterface import DataBaseInterface as DBI
from GCGA.Operations.CrossOperation import CrossOperation as CO
from GCGA.Operations.AddOperation import AddOperation as AD
from GCGA.Operations.RemoveOperation import RemoveOperation as RM
from GCGA.Operations.PermutationOperation import PermutationOperation as PM

from ase import Atoms
import numpy as np
###########################
from ase.io import read,write
#Calculator
from ase.optimize import BFGS
from ase.calculators.emt import EMT

#--------------------------Define the fitness fucntion for your atoms-----------------------"

def fitness_function(atoms,env,reference = 0.0,au_energy = 0.0)-> float:

    at_num= atoms.get_atomic_numbers()
    pt_num = np.count_nonzero(at_num == 78)
    au_num= np.count_nonzero(at_num == 79)

    fre = atoms.get_potential_energy() - reference - au_num * (au_energy) - au_num*env
    return fre

#---------------------------Generate static part of the system------------------------------"
a = 12.0
slab = Atoms(cell=[a,a,a],
             pbc=True)

#If slab atoms arent fixed they will be generated  for all structures, but will be relaxed!"

#---------------------------Generate constant part of the system------------------------------"
#Part of the system that can be relaxed but that does not vary in number over the duration of the search"
constant = Atoms('Pt1')

#---------Generate variable part of the system-(Single atom types only)------------------"
#Part of the system that can be relaxed and that varyies in number over the duration of the search"
#
variable_types = [Atoms('PtAu2'),Atoms('Pt5Au'),Atoms('CO')]
for i in variable_types:
    print(i.symbols.indices())



#Considered number of variable type atoms n the search"
variable_range = [[1,2],[2,3],[4,5,6]]
lengths = 1
lengths_array = []
for i in range(len(variable_types)):
   # print(len(variable_range[i]))
    lengths_array.append(len(variable_range[i]))
    lengths = lengths * len(variable_range[i])

combiantion_matrix  = np.zeros((lengths,len(variable_range)),dtype = int)

x=0
while(x<lengths):
    for i in range(len(variable_range)):
      #  combiantion_matrix[x,i] = variable_range[i][j] 
        advancement = int(np.prod(lengths_array[i+1:len(variable_range)]))
        pos = int(x / advancement) % len(variable_range[i])
        combiantion_matrix[x,i]=variable_range[i][pos]
    x+=1
print(combiantion_matrix)
#---------------------------Define starting population--------------------------------"

candidateGenerator = SCG(slab,constant,variable_types,variable_range)
"""crossing = CO(slab,constant,variable_types,variable_range,stc_change_chance = 0.5)
adding = AD(slab,constant,variable_types,variable_range)
removing = RM(slab,constant,variable_types,variable_range)
permutating = PM(slab,constant,variable_types,variable_range)
"""
db = DBI('databaseGA.db')
population = 25

starting_pop = candidateGenerator.get_starting_population(population_size=population)

#-----------------------------------Define Calculator----------------------------------------"

calc = EMT()

#Get Reference energies##########

ref = 0.0

au_en = 0.0

#####################################



for i in starting_pop:
    db.add_unrelaxed_candidate(i )

#---------------------------------Relax initial structures-----------------------------------"
env = -0.5 #Environment chem pot 
while db.get_number_of_unrelaxed_candidates() > 0:

    atoms = db.get_first_unrelaxed()
    
    atoms.calc = calc
    dyn = BFGS(atoms)
    dyn.run(steps=100, fmax=0.05)
    atoms.get_potential_energy()
    
    atoms.info['key_value_pairs']['raw_score'] = -fitness_function(atoms,env,reference = ref, au_energy = au_en)
   
    db.update_to_relaxed(atoms.info['key_value_pairs']['dbid'],atoms)

#--------------------------------------Find better stoich to srtart eval----------------------------
#Overall fitness value
steps = 10
sub_steps = 20
for i in range(steps):
    for j in range(sub_steps):
        atomslist = db.get_better_candidates(n=10)
        #Choose two of the most stable structures to pair
        cand1 = np.random.randint(0,10)
        cand2 = np.random.randint(0,10)
        while cand1 == cand2:
            cand2 = np.random.randint(0,10)
            
        #Mate the particles
        res  = crossing.cross(atomslist[cand1],atomslist[cand2])
        if(res is not None):
            db.update_penalization(atomslist[cand1])
            db.update_penalization(atomslist[cand2])
            db.add_unrelaxed_candidate(res)
            
    
    while db.get_number_of_unrelaxed_candidates() > 0:

        atoms = db.get_first_unrelaxed()

        atoms.calc = calc
        dyn = BFGS(atoms)
        dyn.run(steps=100, fmax=0.05)
        atoms.get_potential_energy()

        atoms.info['key_value_pairs']['raw_score'] = -fitness_function(atoms,env,reference = ref, au_energy = au_en)

        db.update_to_relaxed(atoms.info['key_value_pairs']['dbid'],atoms)
        
atomslist = db.get_better_candidates_weighted_penalized(n=2)

