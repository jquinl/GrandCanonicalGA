from ase import Atom, Atoms
from GCGA.CoreUtils.DataBaseInterface import DataBaseInterface as DBI
from GCGA.Operations.OperationsBase import OperationsBase
from GCGA.Operations.RandomCandidateGenerator import RandomCandidateGenerator as RCG
from GCGA.Operations.CrossOperation import CrossOperation as CO
from GCGA.Operations.AddOperation import AddOperation as AD
from GCGA.Operations.RemoveOperation import RemoveOperation as RM
from GCGA.Operations.PermutationOperation import PermutationOperation as PM
from GCGA.Operations.RattleOperation import RattleOperation as RT

import numpy as np
from ase.io import write

from ase.optimize import BFGS
from ase.calculators.emt import EMT

class GCGA:

    def __init__(self,calculator, slab,atomic_types,atomic_ranges,mutation_operations
                ,mutation_chances,fitness_function,
                structures_filename = 'structures.traj',db_name = 'databaseGA.db',
                starting_poupulation = 20,population_size = 50,
                stoichiometry_weight = 1.0,penalty_strength = 0.0,
                initial_structure_generator = RCG, crossing_operator = CO, 
                steps = 1000,maxtries = 10000,
                write_while_running = True,write_size = 10000):
        
        self.calc = calculator
        self.slab = slab
        self.atomic_types = atomic_types
        self.atomic_ranges = atomic_ranges
        self.mutation_operations, self.mutation_chances = self.__initialize_mutations(mutation_operations,mutation_chances)
        self.fitness_function = fitness_function

        self.db_name = db_name
        self.filename = structures_filename

        self.starting_pupulation = starting_poupulation
        self.population = population_size
        self.wt = stoichiometry_weight
        self.pts = penalty_strength

        self.initial_structure_generator = self.__initialize_generator(initial_structure_generator,RCG)
        self.crossing_operator =self.__initialize_crossing(crossing_operator,CO)

        self.steps = steps
        self.maxtries = maxtries

        self.write_running = write_while_running
        self.write_size = write_size

    def __initialize_crossing(self,crossing,ctype)-> object:

        if issubclass(crossing,CO):
            return crossing
        elif type(crossing) is type(ctype):
            return self.__default_instancing(crossing)
        else:
            raise TypeError ("Unssupported Crossing operator")

    def __initialize_generator(self,gen,gtype)-> object:
        print(type(gtype),type(gen))
        if issubclass(gen,RCG):
            return gen
        elif type(gen) is type(gtype):
            return self.__default_instancing(gen)
        else:
            raise TypeError ("Unssupported Random structure generator")

    def __initialize_mutations(self,mutations,mutation_chance):

        if(not isinstance(mutations,list)): raise TypeError("The mutations variable is not a list of mutation operator")
        if(not isinstance(mutation_chance,list)):raise TypeError("The mutations variable is not a list of floats")
        if(len(mutations)== 0): raise ValueError("The mutations list is empty")
        if(len(mutation_chance)== 0): raise ValueError("The mutations chance list is empty")
        if(not isinstance(mutation_chance[0],float)):raise TypeError("The mutations chance is not a list of floats")
        if(len(mutations)!= len(mutation_chance)): raise ValueError("The mutations list and the chances list dont have the same size")
        print(mutations)
        mut = []
        for i in mutations:
            if isinstance(i,str):
                mut.append( self.__default_instancing_string(i))
            else:
                try:
                    i.mutate()
                    mut.append(i)
                except:
                    try:
                        new_cls = self.__default_instancing(i)
                        new_cls.mutate()
                        mut.append(new_cls)
                    except:
                        raise ValueError("Provided class parameter could not be parsed to a mutation operator ")
        return list(mut),mutation_chance
    def __default_instancing(self,isClass):

            if type(isClass) is RCG:
                return RCG(self.slab,self.atomic_types,self.atomic_ranges)
            if type(isClass) is  CO:
                return CO(self.slab,self.atomic_types,self.atomic_ranges,minfrac=0.1)
            if type(isClass) is  AD:
                return AD(self.slab,self.atomic_types,self.atomic_ranges)
            if type(isClass) is  RM:
                return RM(self.slab,self.atomic_types,self.atomic_ranges)
            if type(isClass) is  PM:
                return PM(self.slab,self.atomic_types,self.atomic_ranges)
            if type(isClass) is  RT:
                return RT(self.slab,self.atomic_types,self.atomic_ranges,n_to_move = 1,rattle_strength = 0.1)

            raise TypeError("Provided mutation type is not supported",type(isClass))
    
    def __default_instancing_string(self,isClass):   
            if isClass == "random":
                return RCG(self.slab,self.atomic_types,self.atomic_ranges)
            if isClass == "cross":
                return CO(self.slab,self.atomic_types,self.atomic_ranges,minfrac=0.1)
            if isClass == "add":
                return AD(self.slab,self.atomic_types,self.atomic_ranges)
            if isClass == "remove":
                return RM(self.slab,self.atomic_types,self.atomic_ranges)
            if isClass == "permute":
                return PM(self.slab,self.atomic_types,self.atomic_ranges)
            if isClass == "rattle":
                return RT(self.slab,self.atomic_types,self.atomic_ranges,n_to_move = 1,rattle_strength = 0.1)

            raise TypeError("Provided mutation string is not supported:",isClass)
    

    def run(self):
#---------------------------Define starting population--------------------------------"
        db = DBI(self.db_name)

        calc = self.calc
        mutations = self.mutation_operations
        chances = self.mutation_chances
        population = self.starting_poupulation

        #--------------------------------Generate initial population---------------------------------"
        starting_pop = self.initial_structure_generator.get_starting_population(population_size=population)


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
            
            atoms.info['key_value_pairs']['raw_score'] = self.fitness_function(atoms)
        
            db.update_to_relaxed(atoms)

        #--------------------------------------Find better stoich to srtart eval----------------------------

        #Overall fitness value
        steps = self.steps
        maxtries = self.maxtries

        counter = 0
        maxcounter = 0

        while counter < steps and maxcounter < maxtries:
            maxcounter += 1

            atomslist = db.get_better_candidates_weighted_penalized(n=self.population,wt_strength =self.wt,penalty_strength=self.pts )
            ranges = len(atomslist)
            #Choose two of the most stable structures to pair
            cand1 = np.random.randint(0,ranges - 1)
            cand2 = np.random.randint(0,ranges - 1)
            if(ranges >1):
                while cand1 == cand2:
                    cand2 = np.random.randint(0,ranges -1)

            #Mate the particles
                res,mut  = self.crossing_operator.mutate(atomslist[cand1],atomslist[cand2])
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

                    atoms.info['key_value_pairs']['raw_score'] = -fitness_function(atoms)

                    db.update_to_relaxed(atoms)

                if(self.write_running):
                    atomslist = db.get_better_candidates(n=self.write_size)
                    write(self.structures_filename,atomslist)

        atomslist = db.get_better_candidates(max=True)
        write(self.structures_filename,atomslist)