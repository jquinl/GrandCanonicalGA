import random
import time
from typing import List
from ase import Atoms

from Core.SubunitAnalysis import NonEnergyInteratomicDistanceComparator
from ase.ga.utilities import closest_distances_generator,get_all_atom_types

from GCGA.Core.DataBaseInterface import DataBaseInterface as DBI
from GCGA.Core.Population import Population
from GCGA.FitnessFunction.BaseFitness import BaseFitness
from GCGA.Operations.RandomCandidateGenerator import RandomCandidateGenerator as RCG
from GCGA.Operations.CrossOperation import CrossOperation as CO

from GCGA.Operations.RemoveOperation import RemoveOperation as RM
from GCGA.Operations.AddOperation import AddOperation as AD
from GCGA.Operations.PrepareForDB import PrepareForDB as PDB

from ase.ga.standardmutations import MirrorMutation
from ase.ga.standardmutations import RattleMutation
from ase.ga.standardmutations import PermutationMutation



import numpy as np
from ase.io import read,write, Trajectory
from os import path

"Supported calculators"
from ase.calculators.singlepoint import SinglePointCalculator
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.calculators.lammpslib import LAMMPSlib


class GCGA:

    def __init__(self, slab,atomic_types,atomic_ranges,mutation_operations,
                fitness_function,mutation_chance=0.3,
                structures_filename = 'structures.traj',db_name = 'databaseGA.db',
                starting_population = 10,population_size = 20,population_size_even = False,

                calculator = EMT(),
                initial_structure_generator = RCG, crossing_operator = CO, 
                steps = 1000,maxtries = 10000,
                restart=False,restart_filename=None
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
        
        if(population_size_even):
            while(len(self.__set_combination_matrix()) % population_size != 0):
                population_size += 1

        self.pop_size = population_size
        self.starting_population = starting_population

        self.initial_structure_generator = self.__initialize_generator(initial_structure_generator)
        self.crossing_operator =self.__cross_to_instance(crossing_operator)

        self.steps = steps
        self.maxtries = max(steps*10,maxtries)

        self.restart = restart
        self.restart_filename = restart_filename
        
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
            if issubclass(isClass,PM):
                return PM(self.slab,self.atomic_types,self.atomic_ranges)
            if issubclass(isClass,RT):
                return RT(self.slab,self.atomic_types,self.atomic_ranges,n_to_move = 1,rattle_strength = 0.1)
            if issubclass(isClass,AD):
                return AD(self.slab,self.atomic_types,self.atomic_ranges)
            if issubclass(isClass,RM):
                return RM(self.slab,self.atomic_types,self.atomic_ranges)
            
            raise TypeError("Provided mutation type is not supported",type(isClass))
    
    def __default_instancing_string(self,isClass):   
            if isClass == "random":
                return RCG(self.slab,self.atomic_types,self.atomic_ranges)
            if isClass == "cross":
                return CO(self.slab,self.atomic_types,self.atomic_ranges,minfrac=0.2)
            if isClass == "permute":
                return PM(self.slab,self.atomic_types,self.atomic_ranges)
            if isClass == "rattle":
                return RT(self.slab,self.atomic_types,self.atomic_ranges,n_to_move = 1,rattle_strength = 0.1)
            if isClass == "add":
                return AD(self.slab,self.atomic_types,self.atomic_ranges)
            
            if isClass == "remove":
                return RM(self.slab,self.atomic_types,self.atomic_ranges)
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
        population = self.starting_population
        
        
        #Here start population if size is smaller than pop size fill 
        # it w random structures

        pop = Population(self.pop_size,db)

        if(self.restart and self.restart_filename is not None):

            restart_pop = self.restart_run()
            if(restart_pop is not None):

                for atoms in restart_pop:
                    db.add_unrelaxed_candidate(atoms)
                    atoms = db.get_first_unrelaxed()
                    
                    if(self.is_fitness_an_object):
                        atoms.info['key_value_pairs']['raw_score'] = self.fitness_function.evaluate(self.slab,atoms)
                    else:    
                        atoms.info['key_value_pairs']['raw_score'] = self.fitness_function(atoms)
                    
                    
                    
                    self.append_to_file(atoms)
                    pop.update_population(atoms)
                    
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

            self.append_to_file(atoms)

            pop.update_population(atoms)


        steps = self.steps
        maxtries = self.maxtries

        counter = 0
        maxcounter = 0
        poprange = [2,5,10]

        while counter < steps and maxcounter < maxtries:
            maxcounter += 1
            succes = False
            subtries = 0
            while(not succes and subtries<3):
                atomslist = pop.get_better_candidates_singles(n=poprange[subtries])
                
                subtries +=1

                ranges = len(atomslist)
                #Choose two of the most stable structures to pair
                cand1 = 0#np.random.randint(ranges)
                cand2 = 1
                if(ranges > 1):
                    cand2 = np.random.randint(1,ranges)
                #Mate the particles
                    res = self.crossing_operator.cross(atomslist[cand1],atomslist[cand2])
                    
                    if(res is not None):
                        succes = True
                        db.update_penalization(atomslist[cand1])
                        db.update_penalization(atomslist[cand2])
                        child = res.copy()


                        rnd = np.random.rand()
                        mutated = None
                        if(rnd < self.mutation_chance):
                            ran = list(range(len(mutations)))
                            random.shuffle(ran)
                            for k in ran:
                                if(mutated == None):
                                    mutated = mutations[k].mutate(res)
                        if(mutated is not None):
                            child = mutated.copy()
                                                    
                        if child is not None:
                            db.add_unrelaxed_candidate(child)
                            print("Structure evaluation {}".format(counter))
                            
                            counter+=1
                else:
                    starting_pop = self.initial_structure_generator.get_random_candidate()
            
            while db.get_number_of_unrelaxed_candidates() > 0:
                atoms = db.get_first_unrelaxed()
                atoms = self.relax(atoms)

                if(self.is_fitness_an_object):
                    atoms.info['key_value_pairs']['raw_score'] = self.fitness_function.evaluate(self.slab,atoms)
                else:    
                    atoms.info['key_value_pairs']['raw_score'] = self.fitness_function(atoms)
           
                self.append_to_file(atoms)

                pop.update_population(atoms)

        atomslist = db.get_better_candidates_raw(max_num=True)

        write("sorted_" + self.filename,atomslist)

    def append_to_file(self,atoms):
        if(self.trajfile is not None):
            self.trajfile.write(atoms)

    def restart_run(self):
        

        atomslist = list(read(self.restart_filename + "@:"))
        if(not isinstance(atomslist,list)): return None
        pdb = PDB(self.slab,self.atomic_types,self.atomic_ranges,self.crossing_operator.ratio_of_covalent_radii,self.crossing_operator.rng)
        returnatoms = []
        for atoms in atomslist:
            if(isinstance(atoms.calc,SinglePointCalculator)):
                returnatoms.append(pdb.prepare(atoms))
        if(len(returnatoms) == 0): return None
        if(len(returnatoms) != len(atoms)): print("Not all atoms in {} were included in the run".format(self.restart_filename))
        return list(returnatoms)

    def relax(self,atoms):
        results = None
            
        if(isinstance(self.calc,EMT)):
            atoms.set_calculator( self.calc)
            dyn = BFGS(atoms, trajectory=None, logfile=None)
            dyn.run(fmax=0.05, steps=100)
            E = atoms.get_potential_energy()
            F = atoms.get_forces()
            results = {'energy': E,'forces': F}
        elif(isinstance(self.calc,LAMMPSlib)):
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
            #self.append_to_file(atoms)
            return atoms
        else:
            return None

    def __set_combination_matrix(self):
            if(len(self.atomic_ranges) != len(self.atomic_types)): raise ValueError("Variable type list length and variable range list length dont match")
            try:
                lengths = 1
                lengths_array = []
                for i in range(len(self.atomic_types)):
                    lengths_array.append(len(self.atomic_ranges[i]))
                    lengths = lengths * len(self.atomic_ranges[i])

                combiantion_matrix  = np.zeros((lengths,len(self.atomic_ranges)),dtype = int)

                for x in range(lengths):
                    for i in range(len(self.atomic_ranges)):
                        advancement = int(np.prod(lengths_array[i+1:len(self.atomic_ranges)]))
                        pos = int(x / advancement) % len(self.atomic_ranges[i])
                        combiantion_matrix[x,i]=self.atomic_ranges[i][pos]
                return combiantion_matrix
            except:
                raise Exception("Could not generate variable dictionary: Make sure variable type length and variable range length match")

    def _mantains_ordering(self,atoms):
        if(len(atoms) < len(self.slab)):
            return False
        if(self.slab.symbols.indices() != atoms[:len(self.slab)].symbols.indices()):
            return False

        return True
    def _get_var_id(self,atoms) -> int:
        if(not self.mantains_ordering(atoms)): raise Exception("Does not mantain atomic ordering")
        if(len(atoms) == len(self.slab)):
            return self.variable_dict[0]
        dict = self.atoms_to_hash(atoms[len(self.slab):])
        if(dict in  self.variable_dict):
            return self.variable_dict[dict]
        else:
            return None
  
    def _get_range(self,variable_range) -> List[List[int]]:
        if isinstance(variable_range,List) and isinstance(variable_range[0],List):
            for i in range(len(variable_range)):
                for j in range(len(variable_range[i])):
                    if(not isinstance(variable_range[i][j],int)):
                        raise Exception("variable_ range is not al List of List of integers")
            return variable_range
        else:
            raise Exception("variable_range is not al List of List of integers")
    def _get_variable_types(self,types) -> List[Atoms]:
        "Gets an atoms object based on user input"
        if isinstance(types, List):
            for i in types:
                if( not isinstance(self.__get_atoms_object(i),Atoms)):
                    raise ValueError('Cannot parse this element to Atoms object:', i)
            return types
        else:
            raise ValueError('variable_types not a list of atoms objects:', types)
    
    def __get_atoms_object(self,atoms) -> Atoms:
        "Gets an atoms object based on user input"
        if isinstance(atoms, Atoms):
            return atoms
        elif isinstance(atoms, str):
            return Atoms(atoms)
        elif isinstance(atoms,List):
            for i in atoms:
                if(i not in atomic_numbers.values()):
                    raise ValueError('Cannot parse this element {} in :'.format(i),atoms )
            return Atoms(numbers=atoms)
        else:
            raise ValueError('Cannot parse this element:', atoms)
    def _get_cell_params(self,slab,random_generation_box_size):
        "Gets cell parameters from inputed slab"
        if(random_generation_box_size < 0.0): raise ValueError("random_generation_box_size negative value")
        if(random_generation_box_size > 1.0): raise ValueError("random_generation_box_size too big")

        pos = slab.get_positions()
        cell = slab.get_cell()
        if(len(pos) == 0):
            v1 = cell[0, :] * random_generation_box_size
            v2 = cell[1, :] * random_generation_box_size
            v3 = cell[2, :] * random_generation_box_size
            p0 = np.array([0,0,0])
        else:
            p0 = np.array([0., 0., max(pos[:, 2]) + 2.])
            v1 = cell[0, :] * random_generation_box_size
            v2 = cell[1, :] * random_generation_box_size
            v3 = cell[2, :] * random_generation_box_size
            v3 = v3-p0

        return p0,v1,v2,v3

    def sort_atoms_by_type(self,atoms):
        at = Atoms()
        for i in self.variable_types:
            for k in i:
                for j in atoms:
                    if(k.symbol == j.symbol):
                        at.extend(j)
        return at

    def __set_blmin(self,slab, variable_types):
        uniques = Atoms()
        for i in variable_types:
            uniques.extend(i)
            
        unique_atom_types = get_all_atom_types(slab, uniques.numbers)

        return closest_distances_generator(atom_numbers=unique_atom_types,
                                    ratio_of_covalent_radii=self.ratio_of_covalent_radii)

    def is_structure_equal(self,atoms1,atoms2):
        if(len(atoms1) != len(atoms2)): return False
        
        comp = NonEnergyInteratomicDistanceComparator(n_top=len(atoms1), pair_cor_cum_diff=0.015,
                pair_cor_max=0.7, mic=True)
        return comp.looks_like(atoms1,atoms2)
    
    def __get_blmin(self,slab, atoms):
        uniques = Atoms()
        for i in atoms:
            if(i.symbol not in [ j.symbol for j in uniques]):
                uniques.extend(i)
            
        unique_atom_types = get_all_atom_types(slab, uniques.numbers)

        return closest_distances_generator(atom_numbers=unique_atom_types)