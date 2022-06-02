import random
from GCGA.CoreUtils.DataBaseInterface import DataBaseInterface as DBI
from GCGA.FitnessFunction.BaseFitness import BaseFitness
from GCGA.Operations.RandomCandidateGenerator import RandomCandidateGenerator as RCG
from GCGA.Operations.CrossOperation import CrossOperation as CO
"""from GCGA.Operations.AddOperation import AddOperation as AD
from GCGA.Operations.RemoveOperation import RemoveOperation as RM
from GCGA.Operations.PermutationOperation import PermutationOperation as PM
from GCGA.Operations.RattleOperation import RattleOperation as RT"""

import numpy as np
from ase.io import write, Trajectory
from os import path

"Supported calculators"
from ase.calculators.singlepoint import SinglePointCalculator
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.calculators.lammpslib import LAMMPSlib


class GCGA:

    #similarity_penalty is a Prototype, use with caution, for runs with high number of steps impacts performance significantly"
    def __init__(self, slab,atomic_types,atomic_ranges,mutation_operations,
                fitness_function,mutation_chance=0.3,
                structures_filename = 'structures.traj',db_name = 'databaseGA.db',
                starting_population = 10,population_size = 2,
                stoichiometry_weight = True,similarity_penalty = False,calculator = EMT(),
                initial_structure_generator = RCG, crossing_operator = CO, 
                steps = 1000,maxtries = 10000,
                ):

        self.calc = calculator
        self.slab = slab
        self.atomic_types = atomic_types
        self.atomic_ranges = atomic_ranges
        self.mutation_operations, self.mutation_chance = self.__initialize_mutations(mutation_operations,mutation_chance)

        self.is_fitness_an_object = False
        self.fitness_function = self.__initialize_fitness_function(fitness_function)

        self.db_name = db_name

        if isinstance(structures_filename, str):
            self.filename = structures_filename
            "This performs a hard overwrite, change it later in order to include restarts"
            if path.exists(structures_filename):
                self.trajfile = Trajectory(filename=structures_filename, mode='w')
                self.trajfile.close()
                self.trajfile = Trajectory(filename=structures_filename, mode='a')
            else:
                self.trajfile = Trajectory(filename=structures_filename, mode='a')
        
        self.starting_population = starting_population
        self.population = population_size
        self.wt = stoichiometry_weight
        self.pts = similarity_penalty
        if(not stoichiometry_weight and similarity_penalty):
            print("Similarity penalty is marked true but stoichiometry weight is not. Currently only similarity penalty is not supported so sotichiometry weight will be marked true")
            self.wt = True
            self.pts = True

        self.initial_structure_generator = self.__initialize_generator(initial_structure_generator)
        self.crossing_operator =self.__cross_to_instance(crossing_operator)

        self.steps = steps
        self.maxtries = max(steps*10,maxtries)

    def __initialize_generator(self,gen)-> object:
        try:
            gen.rand_generator()
            return gen
        except:
            raise TypeError ("Unssupported Random structure generator")

    def __initialize_mutations(self,mutations,mutation_chance):

        if(not isinstance(mutation_chance,float)):      raise TypeError("Mutation chance is not a float")
        if(not isinstance(mutations,list)):             raise TypeError("The mutations variable is not a list of mutation operator")
        if(len(mutations)== 0):                         raise ValueError("The mutations list is empty")

        mut = []
        for i in mutations:
            if isinstance(i,str):
                mut.append( self.__default_instancing_string(i))
            else:
                mut.append(self.__mutation_to_instance(i))

        return list(mut),mutation_chance

    def __mutation_to_instance(self,isClass):
            try:
                isClass.mutation_instance()
                return isClass
            except:
                try:
                    isClass.mutation_class()
                    new_cls = self.__default_instancing(isClass)
                    new_cls.mutation_instance()
                    return new_cls
                except:
                    raise ValueError("Provided class parameter could not be parsed to a mutation operator ")
    def __cross_to_instance(self,isClass):
            try:
                isClass.cross_instance()
                return isClass
            except:
                try:
                    isClass.cross_class()
                    new_cls = self.__default_instancing(isClass)
                    new_cls.cross_instance()
                    return new_cls
                except:
                    raise ValueError("Provided class parameter could not be parsed to a mutation operator ")

    def __default_instancing(self,isClass):

            if issubclass(isClass,RCG):
                return RCG(self.slab,self.atomic_types,self.atomic_ranges)
            if issubclass(isClass,CO):
                return CO(self.slab,self.atomic_types,self.atomic_ranges,minfrac=0.2)
            """  if issubclass(isClass,AD):
                return AD(self.slab,self.atomic_types,self.atomic_ranges)
            if issubclass(isClass,RM):
                return RM(self.slab,self.atomic_types,self.atomic_ranges)
            if issubclass(isClass,PM):
                return PM(self.slab,self.atomic_types,self.atomic_ranges)
            if issubclass(isClass,RT):
                return RT(self.slab,self.atomic_types,self.atomic_ranges,n_to_move = 1,rattle_strength = 0.1)
            """
            raise TypeError("Provided mutation type is not supported",type(isClass))
    
    def __default_instancing_string(self,isClass):   
            if isClass == "random":
                return RCG(self.slab,self.atomic_types,self.atomic_ranges)
            if isClass == "cross":
                return CO(self.slab,self.atomic_types,self.atomic_ranges,minfrac=0.2)
            """if isClass == "add":
                return AD(self.slab,self.atomic_types,self.atomic_ranges)
            if isClass == "remove":
                return RM(self.slab,self.atomic_types,self.atomic_ranges)
            if isClass == "permute":
                return PM(self.slab,self.atomic_types,self.atomic_ranges)
            if isClass == "rattle":
                return RT(self.slab,self.atomic_types,self.atomic_ranges,n_to_move = 1,rattle_strength = 0.1)
            """
            raise TypeError("Provided mutation string is not supported:",isClass)
    
    def __initialize_fitness_function(self,function):
        if(callable(function)):
            self.is_fitness_an_object = False
            return function
        if(isinstance(function,BaseFitness)):
            self.is_fitness_an_object = True
            return function
        
        raise Exception("Provided function parameter is not a function or a Object derived from the GCGA.FitnessFunction.BaseFitness class")


    def run(self):
#---------------------------Define starting population--------------------------------"
        db = DBI(self.db_name)

        mutations = self.mutation_operations
        chance = self.mutation_chance
        population = self.starting_population

        #--------------------------------Generate initial population---------------------------------"
        starting_pop = self.initial_structure_generator.get_starting_population(population_size=population)

        for i in starting_pop:
            db.add_unrelaxed_candidate(i)

        #---------------------------------Relax initial structures-----------------------------------"

        while db.get_number_of_unrelaxed_candidates() > 0:

            atoms = db.get_first_unrelaxed()
            
            atoms = self.relax(atoms)
            if(self.is_fitness_an_object):
                atoms.info['key_value_pairs']['raw_score'] = self.fitness_function.evaluate(self.slab,atoms)
            else:    
                atoms.info['key_value_pairs']['raw_score'] = self.fitness_function(atoms)
            db.update_to_relaxed(atoms)


        steps = self.steps
        maxtries = self.maxtries

        counter = 0
        maxcounter = 0

        while counter < steps and maxcounter < maxtries:
            maxcounter += 1

            atomslist = db.get_better_candidates(n=self.population,weighted=self.wt,structure_similarity=self.pts)
            ranges = len(atomslist)
            #Choose two of the most stable structures to pair
            cand1 = np.random.randint(ranges)
            cand2 = np.random.randint(ranges)
            if(ranges >1):
                while cand1 == cand2:
                    cand2 = np.random.randint(ranges)

            #Mate the particles
                res = self.crossing_operator.cross(atomslist[cand1],atomslist[cand2])
                if(res is not None):
                    db.update_penalization(atomslist[cand1])
                    db.update_penalization(atomslist[cand2])
                    counter+=1
                    child = res.copy()

                    rnd = np.random.rand()
                    mutated = None
                    if(rnd < self.mutation_chance):
                        ran = range(len(mutations))
                        random.shuffle(ran)
                        for k in ran:
                            if(mutated == None):
                                mutated = mutations[k].mutate(res)
                    
                    if(mutated is not None):
                        child = mutated.copy()
                    

                    if child is not None:
                        db.add_unrelaxed_candidate(child)
                        counter+=1


                while db.get_number_of_unrelaxed_candidates() > 0:

                    atoms = db.get_first_unrelaxed()
                    atoms = self.relax(atoms)

                    if(self.is_fitness_an_object):
                        atoms.info['key_value_pairs']['raw_score'] = self.fitness_function.evaluate(self.slab,atoms)
                    else:    
                        atoms.info['key_value_pairs']['raw_score'] = self.fitness_function(atoms)
                    
                    db.update_to_relaxed(atoms)

        atomslist = db.get_better_candidates_raw(max_num=True)

        write("sorted_" + self.filename,atomslist)

    def append_to_file(self,atoms):
        if(self.trajfile is not None):
            self.trajfile.write(atoms)

    def relax(self,atoms):
        print("structure evaluation")
        results = None
        if(isinstance(self.calc,EMT)):
            atoms.set_calculator( self.calc)
            dyn = BFGS(atoms, trajectory=None, logfile=None)
            dyn.run(fmax=0.05, steps=100)
            E = atoms.get_potential_energy()
            F = atoms.get_forces()
            results = {'energy': E,'forces': F}
        if(isinstance(self.calc,LAMMPSlib)):
            try:
                atoms.set_calculator(self.calc)
                E = atoms.get_potential_energy()
                F = atoms.get_forces()
                results = {'energy': E,'forces': F}
            except:
                raise Exception("LAMMPS not installed")
        else:
            atoms.set_calculator(self.calc)
            atoms.get_potential_energy()
            E = atoms.get_potential_energy()
            F = atoms.get_forces()
            results = {'energy': E,'forces': F}
        if(results is not None):
            calc_sp = SinglePointCalculator(atoms, **results)
            atoms.set_calculator(calc_sp)
            self.append_to_file(atoms)
            return atoms
        else:
            return None