#Import these for the GA code to run
from GCGA.CoreUtils.DataBaseInterface import DataBaseInterface as DBI
from GCGA.Operations.CrossOperation import CrossOperation as CO
from GCGA.Operations.RandomCandidateGenerator import RandomCandidateGenerator as RCG
from GCGA.Operations.AddOperation import AddOperation as AD
from GCGA.Operations.RemoveOperation import RemoveOperation as RM
from GCGA.Operations.PermutationOperation import PermutationOperation as PM
from GCGA.Operations.RattleOperation import RattleOperation as RT

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
a = 24.0
slab = Atoms(cell=[a,a,a],
             pbc=True)

#If slab atoms arent fixed they will be generated  for all structures, but will be relaxed!"

#---------------------------Generate constant part of the system------------------------------"
#Part of the system that can be relaxed but that does not vary in number over the duration of the search"


#---------Generate variable part of the system-(Single atom types only)------------------"
#Part of the system that can be relaxed and that varyies in number over the duration of the search"
#
variable_types = [Atoms('Pt'),Atoms('Au')]
variable_range = [[1],[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20]]

#---------------------------Define starting population--------------------------------"
db = DBI('databaseGA.db')

candidateGenerator = RCG(slab,variable_types,variable_range)
crossing = CO(slab,variable_types,variable_range,minfrac = 0.2)
adding = AD(slab,variable_types,variable_range)
removing = RM(slab,variable_types,variable_range)
permutating = PM(slab,variable_types,variable_range)
rattling = RT(slab,variable_types,variable_range,n_to_move= 2,rattle_strength=0.1)

mutations = [crossing,candidateGenerator,adding,removing,permutating,rattling]
chances = [0.3,0.2,0.1,0.2,0.2,0.1]

population = 21

starting_pop = candidateGenerator.get_starting_population(population_size=population)

#-----------------------------------Define Calculator----------------------------------------"

calc = EMT()

#Get Reference energies##########


ref=read('pt.traj@:')[0].get_potential_energy()    

au_en =read('gold_bulk.traj').get_potential_energy() / 500.0

#####################################

for i in starting_pop:
    db.add_unrelaxed_candidate(i)

#---------------------------------Relax initial structures-----------------------------------"
env = 1.0 #Environment chem pot

while db.get_number_of_unrelaxed_candidates() > 0:

    atoms = db.get_first_unrelaxed()
    
    atoms.calc = calc
    dyn = BFGS(atoms)
    dyn.run(steps=100, fmax=0.05)
    atoms.get_potential_energy()
    
    atoms.info['key_value_pairs']['raw_score'] = -fitness_function(atoms,env,reference = ref, au_energy = au_en)
   
    db.update_to_relaxed(atoms)

#--------------------------------------Find better stoich to srtart eval----------------------------

#Overall fitness value
steps = 1000
maxtries = 10000

counter = 0
maxcounter = 0

while counter < steps and maxcounter < maxtries:
    maxcounter += 1

    atomslist = db.get_better_candidates_weighted(n=50)
    ranges = len(atomslist)
    #Choose two of the most stable structures to pair
    cand1 = np.random.randint(0,ranges - 1)
    cand2 = np.random.randint(0,ranges - 1)
    if(ranges >1):
        while cand1 == cand2:
            cand2 = np.random.randint(0,ranges -1)

    #Mate the particles
        res,mut  = crossing.mutate(atomslist[cand1],atomslist[cand2])
        if(res is not None):
            print("Crossed--------------------------------------------------------------------------------------------------------------------------------")
            db.update_penalization(atomslist[cand1])
            db.update_penalization(atomslist[cand2])
            db.add_unrelaxed_candidate(res)
            counter+=1
        else:
        #If it doesn't succesfully mate particles it performs the selected mutations
            rnd = np.random.rand()
            chance = 0.0
            permut = 0
            if len(mutations) != len(chances): raise Exception("Lenght of mutations array is different from the chances array")
            child = None
            for k in range(len(mutations)):
                chance+= chances[k]
                placeholder,mut = mutations[k].mutate(atomslist[cand1],atomslist[cand2])
                if(placeholder is not None):
                    child = placeholder.copy()
                    permut = mut
                if (rnd < chance and placeholder is not None):
                    break
            if child is not None:
                db.add_unrelaxed_candidate(child)
                if permut == 0:
                    pass
                elif permut == 1:
                    db.update_penalization(atomslist[cand1])
                elif permut == 2:
                    db.update_penalization(atomslist[cand1])
                    db.update_penalization(atomslist[cand2])
                counter+=1


        while db.get_number_of_unrelaxed_candidates() > 0:

            atoms = db.get_first_unrelaxed()

            atoms.calc = calc
            dyn = BFGS(atoms)
            dyn.run(steps=100, fmax=0.05)
            atoms.get_potential_energy()

            atoms.info['key_value_pairs']['raw_score'] = -fitness_function(atoms,env,reference = ref, au_energy = au_en)

            db.update_to_relaxed(atoms)

        
        atomslist = db.get_better_candidates(n=2000)
        write('structures.traj',atomslist)



