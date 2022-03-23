#Import these for the GA code to run
from StartingCandidateGenerator import StartingCandidateGenerator as SCG
from CoreUtils.DataBaseInterface import DataBaseInterface as DBI
from Operations.CrossOperation import CrossOperation as CO
from ase import Atoms
import numpy as np
###########################

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
variable_type = Atoms('Au')
#Considered number of variable type atoms n the search"
variable_range = [1,2,3,4,5,7,8,9,10]

#---------------------------Define starting population--------------------------------"

candidateGenerator = SCG(slab,constant,variable_type,variable_range)
crossing = CO(slab,constant,variable_type,variable_range,stc_change_chance = 0.5)
db = DBI('databaseGA.db')
population = 25

starting_pop = candidateGenerator.get_starting_population(population_size=population)

#-----------------------------------Define Calculator----------------------------------------"

calc = EMT()

#Get Reference energies##########
constant.cell = slab.get_cell()
variable_type.cell = slab.get_cell()

constant.calc = calc
dyn = BFGS(constant)
dyn.run(steps=100, fmax=0.05)
ref = constant.get_potential_energy()

variable_type.calc = calc
dyn = BFGS(variable_type)
dyn.run(steps=100, fmax=0.05)
au_en = variable_type.get_potential_energy()

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
            db.add_unrelaxed_candidate(res)
            
    
    while db.get_number_of_unrelaxed_candidates() > 0:

        atoms = db.get_first_unrelaxed()

        atoms.calc = calc
        dyn = BFGS(atoms)
        dyn.run(steps=100, fmax=0.05)
        atoms.get_potential_energy()

        atoms.info['key_value_pairs']['raw_score'] = -fitness_function(atoms,env,reference = ref, au_energy = au_en)

        db.update_to_relaxed(atoms.info['key_value_pairs']['dbid'],atoms)
        
atomslist = db.get_better_candidates(n=2)