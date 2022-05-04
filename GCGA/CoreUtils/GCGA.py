from struct import calcsize
from GCGA.CoreUtils.DataBaseInterface import DataBaseInterface as DBI
from GCGA.Operations.RandomCandidateGenerator import RandomCandidateGenerator as RCG
from GCGA.Operations.CrossOperation import CrossOperation as CO
from GCGA.Operations.AddOperation import AddOperation as AD
from GCGA.Operations.RemoveOperation import RemoveOperation as RM
from GCGA.Operations.PermutationOperation import PermutationOperation as PM
from GCGA.Operations.RattleOperation import RattleOperation as RT

import numpy as np
from ase.io import write, Trajectory
from os import path

"Supported calculators"
from ase.calculators.singlepoint import SinglePointCalculator
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.calculators.lammpslib import LAMMPSlib


class GCGA:

    def __init__(self, slab,atomic_types,atomic_ranges,mutation_operations,
                mutation_chances,fitness_function,
                structures_filename = 'structures.traj',db_name = 'databaseGA.db',
                starting_population = 20,population_size = 50,
                stoichiometry_weight = 1.0,penalty_strength = 0.0,calculator = EMT(),
                initial_structure_generator = RCG, crossing_operator = CO, 
                steps = 1000,maxtries = 10000,
                ):
        
        self.calc = calculator
        self.slab = slab
        self.atomic_types = atomic_types
        self.atomic_ranges = atomic_ranges
        self.mutation_operations, self.mutation_chances = self.__initialize_mutations(mutation_operations,mutation_chances)
        self.fitness_function = fitness_function

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
        self.pts = penalty_strength

        self.initial_structure_generator = self.__initialize_generator(initial_structure_generator,RCG)
        self.crossing_operator =self.__initialize_crossing(crossing_operator,CO)

        self.steps = steps
        self.maxtries = max(steps*10,maxtries)

    def __initialize_crossing(self,crossing,ctype)-> object:

        if issubclass(crossing,ctype):
            return self.__mutation_to_instance(crossing)
        else:
            raise TypeError ("Unssupported Crossing operator")

    def __initialize_generator(self,gen,gtype)-> object:
        if issubclass(gen,gtype):
            return self.__mutation_to_instance(gen)
        else:
            raise TypeError ("Unssupported Random structure generator")

    def __initialize_mutations(self,mutations,mutation_chance):

        if(not isinstance(mutations,list)):             raise TypeError("The mutations variable is not a list of mutation operator")
        if(not isinstance(mutation_chance,list)):       raise TypeError("The mutations variable is not a list of floats")
        if(len(mutations)== 0):                         raise ValueError("The mutations list is empty")
        if(len(mutation_chance)== 0):                   raise ValueError("The mutations chance list is empty")
        if(not isinstance(mutation_chance[0],float)):   raise TypeError("The mutations chance is not a list of floats")
        if(len(mutations)!= len(mutation_chance)):      raise ValueError("The mutations list and the chances list dont have the same size")

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

    def __default_instancing(self,isClass):

            if issubclass(isClass,RCG):
                return RCG(self.slab,self.atomic_types,self.atomic_ranges)
            if issubclass(isClass,CO):
                return CO(self.slab,self.atomic_types,self.atomic_ranges,minfrac=0.1)
            if issubclass(isClass,AD):
                return AD(self.slab,self.atomic_types,self.atomic_ranges)
            if issubclass(isClass,RM):
                return RM(self.slab,self.atomic_types,self.atomic_ranges)
            if issubclass(isClass,PM):
                return PM(self.slab,self.atomic_types,self.atomic_ranges)
            if issubclass(isClass,RT):
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

        mutations = self.mutation_operations
        chances = self.mutation_chances
        population = self.starting_population

        #--------------------------------Generate initial population---------------------------------"
        starting_pop = self.initial_structure_generator.get_starting_population(population_size=population)

        for i in starting_pop:
            db.add_unrelaxed_candidate(i)

        #---------------------------------Relax initial structures-----------------------------------"

        while db.get_number_of_unrelaxed_candidates() > 0:

            atoms = db.get_first_unrelaxed()
            
            atoms = self.relax(atoms)

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
            cand1 = np.random.randint(ranges)
            cand2 = np.random.randint(ranges)
            if(ranges >1):
                while cand1 == cand2:
                    cand2 = np.random.randint(ranges)

            #Mate the particles
                res,mut  = self.crossing_operator.mutate(atomslist[cand1],atomslist[cand2])
                if(res is not None):
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
                    atoms = self.relax(atoms)

                    atoms.info['key_value_pairs']['raw_score'] = self.fitness_function(atoms)
                    
                    db.update_to_relaxed(atoms)

        atomslist = db.get_better_candidates(max=True)

        write("sorted_" + self.filename,atomslist)

    def append_to_file(self,atoms):
        if(self.trajfile is not None):
            self.trajfile.write(atoms)

    def relax(self,atoms):

        results = None
        if(isinstance(self.calc,EMT)):
            atoms.set_calculator( self.calc)
            dyn = BFGS(atoms)
            dyn.run(steps=100, fmax=0.05)
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