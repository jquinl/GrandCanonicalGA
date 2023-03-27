import time
from typing import List, Dict, Any
import hashlib
import json
import numpy as np
import copy
from ase import Atoms
from ase.data import atomic_numbers
from ase.ga.utilities import closest_distances_generator,get_all_atom_types
from ase.io import read,write, Trajectory
import os

from GCGA.Core.NumpyArrayEncoder import NumpyArrayEncoder
from GCGA.Core.SubunitAnalysis import NonEnergyInteratomicDistanceComparator

from GCGA.Core.DataBaseInterface import DataBaseInterface as DBI
from GCGA.Core.Population import Population

from GCGA.FitnessFunction.BaseFitness import BaseFitness

from GCGA.Operations.RandomCandidateGenerator import RandomCandidateGenerator as RCG
from GCGA.Operations.CrossOperation import CrossOperation as CO
from GCGA.Operations.AddRemoveOperation import AddRemoveOperation as CHG

#---------Supported calculators--------------
from ase.calculators.singlepoint import SinglePointCalculator
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.calculators.lammpsrun import LAMMPS


class GCGA:

    def __init__(self, slab,atomic_types,atomic_ranges,
                fitness_function,
                structures_filename = 'structures.traj',db_name = 'databaseGA.db',
                starting_candidates_per_stc = 4,population_size = 20,
                calculator = EMT(),
                ratio_of_covalent_radii = 0.7,
                initial_structure_generator = RCG, crossing_operator = CO, stc_change_operator = CHG,
                mutations = None,mutation_chance=0.3,
                steps = 1000,maxtries = 10000,
                restart_filename=None,debug=False
                ):

        #--------Population settings----------------
        self.calc = calculator
        self.slab = slab
        self.atomic_types = self._get_variable_types(atomic_types)
        self.atomic_ranges = atomic_ranges
        self.combination_matrix = self.__set_combination_matrix()
        self.stc_dict = self.__set_stc_dict()
        self.is_fitness_an_object = False
        self.fitness_function = self.__initialize_fitness_function(fitness_function)

        #--------Population settings----------------
        self.pop_size = population_size 
        self.starting_candidates = starting_candidates_per_stc

        #--------Operator settings----------------
        self.initial_structure_generator = self.__initialize_generator(initial_structure_generator)
        self.crossing_operator =self.__initialize_crossing(crossing_operator)
        self.change_operator =self.__initialize_change_operator(stc_change_operator)
        self.ratio_of_covalent_radii= ratio_of_covalent_radii
        self.blmin = self.__set_blmin(self.slab, self.atomic_types)

        #--------File Management----------------
        self.db_name = db_name

        if isinstance(structures_filename, str):
            self.filename = structures_filename
            "This performs a hard overwrite, change it later in order to include restarts"
            if os.path.exists(structures_filename):
                self.trajfile = Trajectory(filename=structures_filename, mode='w')
                self.trajfile.close()
                self.trajfile = Trajectory(filename=structures_filename, mode='a')
            else:
                self.trajfile = Trajectory(filename=structures_filename, mode='a')

        if(mutations is not None):
            self.mutation_operations, self.mutation_chance = self.__initialize_mutations(None,mutation_chance)

        #--------Run length parameters----------------
        self.steps = steps
        self.maxtries = max(steps*10,maxtries)
        self.restart_filename = restart_filename
        self.evalnum = 0

        self.debug = debug
        if os.path.exists('debug/'):
            pass
        else:
            os.mkdir('debug')

#--------Functions only called during initialization---------------
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

    def __set_stc_dict(self):
        symbol_dictionary = {}
        if(len(self.atomic_ranges) != len(self.atomic_types)): raise ValueError("Variable type list length and variable range list length dont match")
        try:
            for x in range(len(self.combination_matrix)):
                ats = Atoms()
                for i in range(len(self.atomic_ranges)):
                    for j in range(self.combination_matrix[x,i]):
                        ats.extend(self.atomic_types[i].copy())
                variable_id = x
                ats = self.sort_atoms_by_type(ats)
                symbol_dictionary[self.__atoms_to_hash(ats)] = variable_id
            return symbol_dictionary

        except:
            raise Exception("Could not generate variable dictionary: Make sure variable type length and variable range length match")

    def __dict_hash(self,dictionary: Dict[str, Any]) -> str:
        dhash = hashlib.sha1()
        encoded = json.dumps(dictionary, cls=NumpyArrayEncoder,sort_keys=True).encode()
        dhash.update(encoded)
        return dhash.hexdigest()

    def __atoms_to_hash(self,atoms):
        if(not isinstance(atoms,Atoms)):
            raise ValueError("Tried to hash object of type different than atoms Object")
        hashed_dict = self.__dict_hash(atoms.symbols.indices())
        return hashed_dict

    def __initialize_generator(self,gen)-> object:
        try:
            gen.rand_generator()
            try:
                gen.rand_instance()
                return gen
            except:
                return RCG(atom_spread = 4.0)
        except:
            raise TypeError ("Unssupported Random structure generator")

    def __initialize_crossing(self,crs)-> object:
        try:
            crs.cross_class()
            try:
                crs.cross_instance()
                return crs
            except:
                return CO()
        except:
            raise TypeError ("Unssupported Crossing operator")

    def __initialize_change_operator(self,chg)-> object:
        try:
            chg.chg_class()
            try:
                chg.chg_instance()
                return chg
            except:
                return CHG()
        except:
            raise TypeError ("Unssupported Addition operator")

    def __initialize_fitness_function(self,function):
        if(callable(function)):
            self.is_fitness_an_object = False
            return function
        if(isinstance(function,BaseFitness)):
            self.is_fitness_an_object = True
            return function
        raise Exception("Provided function parameter is not a function or a Object derived from the GCGA.FitnessFunction.BaseFitness class")

    def __initialize_mutations(mutation_operations,mutation_chance):
        return mutation_operations,mutation_chance

    def __set_blmin(slab, variable_types):
        uniques = Atoms()
        for i in variable_types:
            uniques.extend(i)
        unique_atom_types = get_all_atom_types(slab, uniques.numbers)
        return closest_distances_generator(atom_numbers=unique_atom_types,
                                    ratio_of_covalent_radii=0.7)

#------------Actual run of the EA-------------------------
    def run(self):
        self.evalnum = 0
        db = DBI(self.db_name)
        pop = Population(self.pop_size,db,self.steps)
        starting_pop = []
        for i in self.combination_matrix:
            for j in range(self.starting_candidates):
                new_cand = self.initial_structure_generator.new_candidate(self.slab,i,self.atomic_types,self.blmin)
                if(self.debug): print("Random structure {0} {1} generated ".format(new_cand.symbols, j))
                starting_pop.append(new_cand)

        #--------------------------------Generate initial population---------------------------------"

        count = 0
        for i in starting_pop:
            at = self.prepare_candidate(i)
            db.add_unrelaxed_candidate(at)

        #---------------------------------Relax initial structures-----------------------------------"
            while db.get_number_of_unrelaxed_candidates() > 0:
                atoms = db.get_first_unrelaxed()
                atoms = self.relax(atoms)
                if(self.is_fitness_an_object):
                    atoms.info['key_value_pairs']['raw_score'] = self.fitness_function.evaluate(self.slab,atoms)
                else:
                    atoms.info['key_value_pairs']['raw_score'] = self.fitness_function(atoms)

                pop.update_population(atoms)
                self.append_to_file(atoms)
                if(self.debug): write("debug/random_eval_{}_{}.traj".format(count,atoms.symbols),atoms)
                count+=1

        #Restart
        if self.restart_filename is not None and os.path.isfile(self.restart_filename):
            def restart():
                atomslist = list(read(self.restart_filename + "@:"))
                if(not isinstance(atomslist,list)): 
                    print("No restart_filename is not a list of atoms objects")
                    return
                if(len(atomslist) == 0 ):
                    print("No atoms in restart_filename")
                    return
                for i,atoms in enumerate(atomslist):
                    if atoms.cell == self.slab.cell:
                        try:
                            E = atoms.get_potential_energy()
                            F = atoms.get_forces()
                            results = {'energy': E,'forces': F}
                            atoms = self.prepare_candidate(atoms)
                            db.add_unrelaxed_candidate(at)
                            calc_sp = SinglePointCalculator(atoms, **results)
                            atoms.calc = calc_sp
                            if(self.is_fitness_an_object):
                                atoms.info['key_value_pairs']['raw_score'] = self.fitness_function.evaluate(self.slab,atoms)
                            else:
                                atoms.info['key_value_pairs']['raw_score'] = self.fitness_function(atoms)
                            pop.update_population(atoms)
                            self.append_to_file(atoms)
                        except:
                            print("Could not add supplied atoms #{} for restart".format(i))
            restart()

        #------------Time-----------------
        cross_time_total = 0
        cross_time = 0
        chg_time_total = 0
        chg_time = 0
        #------------Run cycle---------------

        steps = self.steps
        maxtries = self.maxtries

        counter = 0
        maxcounter = 0
        while counter < steps and maxcounter < maxtries:
            maxcounter += 1
            a1=None
            a2=None
            if(pop.should_change_stc()):
                tm = time.time()
                debug_optput = ""
                if(self.debug):print("--Changing sotichiometry at step {} try {}--".format(counter, maxcounter))
                res = None
                if(np.random.random() < 0.5):
                    if(self.debug):print("----With cross operation")
                    a2,disc = pop.get_better_candidates()
                    a1,disc = pop.get_better_candidates_stc()
                    subtries = 0
                    while(res is None and subtries < 10):
                        res = self.crossing_operator.cross(self.slab,a1,a2,self.blmin)
                        res = self.prepare_candidate(res)
                        debug_optput += " {}".format(subtries)
                        subtries+=1
                else:
                    if(self.debug):print("----By addition/removal")
                    #target_stc = pop.target_stc()
                    #a1 = pop.get_better_candidate()
                    a1,disc = pop.get_better_candidates()
                    a2,disc = pop.get_better_candidates_stc()
                    target_stc = a2.info['key_value_pairs']['var_stc']
                    subtries = 0
                    while(res is None and subtries < 10):
                        current_stc = a1.info['key_value_pairs']['var_stc']
                        res = self.change_operator.change(self.slab,a1,self.combination_matrix[current_stc],
                            self.combination_matrix[target_stc],self.atomic_types,self.blmin)
                        res = self.prepare_candidate(res)
                        debug_optput += " {}".format(subtries)
                        subtries+=1

                chg_time = time.time() -tm
                chg_time_total += chg_time
                if(res is not None):
                    if(self.debug):
                        stri = "--Stc change SUCCES\n"
                        stri += "   at attempt {0}\n   New structure symbol {1}\n".format(debug_optput,res.symbols)
                        stri += "   writen to stc_chg_op_{0}_{1}_{2}.traj\n".format(res.symbols,counter,maxcounter)
                        stri += "   took {} s".format(chg_time)
                        print(stri)
                        if(a2== None):
                            write("debug/stc_chg_op_{0}_{1}_{2}.traj".format(res.symbols,counter,maxcounter),[a1,res])
                        else:
                            write("debug/stc_chg_op_{0}_{1}_{2}.traj".format(res.symbols,counter,maxcounter),[a1,a2,res])

                    #pop.change_current_stc(res.info['stc'])
                    db.add_unrelaxed_candidate(res)
                else:
                    if(self.debug):
                        print("""Crossing at step FAILURE at {0} try {1}  took {2} s""".format(counter,maxcounter,chg_time))
            else:
                debug_optput = ""
                a1,a2 = pop.get_better_candidates()
                res = None
                subtries = 0
                tm = time.time()
                if(self.debug):print("--Crossing at step {} try {}".format(counter, maxcounter))
                while(res is None and subtries < 10):
                    debug_optput += " {}".format(subtries)
                    subtries+=1
                    res = self.crossing_operator.cross(self.slab,a1,a2,self.blmin)
                    res = self.prepare_candidate(res)
                cross_time = time.time() -tm
                cross_time_total += cross_time
                if(res is not None):
                    if(self.debug):
                        stri = "--Crossing operation SUCCES\n"
                        stri += "   at attempt {0}\n   new structure symbol {1}\n".format(debug_optput,res.symbols)
                        stri += "   writen to cross_op_{0}_{1}_{2}.traj\n".format(res.symbols,counter,maxcounter)
                        stri += "   took {} s".format(cross_time)
                        print(stri)
                        if(a2== None):
                            write("debug/cross_op_{0}_{1}_{2}.traj".format(res.symbols,counter,maxcounter),[a1,res])
                        else:
                            write("debug/cross_op_{0}_{1}_{2}.traj".format(res.symbols,counter,maxcounter),[a1,a2,res])

                    succes = True
                    pop.update_penalization(a1,a2)
                    db.add_unrelaxed_candidate(res)
                else:
                    if(self.debug):
                        print("""Crossing at step FAILURE at {0} try {1}  took {2} s""".format(counter,maxcounter,cross_time))

            while db.get_number_of_unrelaxed_candidates() > 0:
                atoms = db.get_first_unrelaxed()
                atoms = self.relax(atoms)
                counter+=1
                if(self.is_fitness_an_object):
                    atoms.info['key_value_pairs']['raw_score'] = self.fitness_function.evaluate(self.slab,atoms)
                else:
                    atoms.info['key_value_pairs']['raw_score'] = self.fitness_function(atoms)
                if(self.debug):print("candidate evaluated")
                self.append_to_file(atoms)
                if(self.debug):print("candidate appended to file")
                pop.update_population(atoms)
                if(self.debug):print("candidate added to population")
                if(self.debug):
                    print("New candidate evaluated with a score of {}  at step {} try {}".format(atoms.info['key_value_pairs']['raw_score'],counter,maxcounter))
                    write("debug/run_eval_{}_{}.traj".format(counter,maxcounter),atoms)

        print("Exited run loop after {0} steps and {1} tries".format(counter,maxcounter))
        if(self.debug):print("Crossing operations took a total of {} s:".format(cross_time_total))
        atomslist = db.get_better_candidates_raw(max_num=True)
        write("sorted_" + self.filename,atomslist)

    #---Methos employed during the run---
    def prepare_candidate(self,atoms):
        if(atoms is None ):return None
        cand = self.slab.copy()
        ats = self.sort_atoms_by_type(atoms[len(self.slab):])
        cand.extend(ats)

        if(not self._mantains_ordering(cand)): return None
        if(self._check_overlap_all_atoms(cand,self.blmin)): return None

        var_id = self._get_var_id(cand)
        if(var_id is None): return None
        cand.info['stc']= var_id
        return cand

    def append_to_file(self,atoms):
        if(self.trajfile is not None):
            self.trajfile.write(atoms)

    def relax(self,atoms):
        results = None
        stri =  "------------------------------\n"
        stri += "----Candidate Evaluation----\n"
        stri += "------------------------------\n"
        stri += "- Atoms symbol: {}\n".format(atoms.symbols)
        stri += "- Calculator : {}\n".format(self.calc)
        c_tm = time.time()
        if(isinstance(self.calc,EMT)):
            stri += "- Using BFGS optimizer\n"
            atoms.calc = self.calc
            dyn = BFGS(atoms, trajectory=None, logfile=None)
            dyn.run(fmax=0.05, steps=100)
            E = atoms.get_potential_energy()
            F = atoms.get_forces()
            results = {'energy': E,'forces': F}
        elif(isinstance(self.calc,LAMMPS)):
            try:
                stri += "- Using BFGS optimizer\n"
                atoms.calc = self.calc
                dyn = BFGS(atoms, trajectory=None, logfile=None)
                dyn.run(fmax=0.05, steps=100)
                E = atoms.get_potential_energy()
                F = atoms.get_forces()
                results = {'energy': E,'forces': F}
            except:
                raise Exception("LAMMPS not installed")
        else:
            atoms.calc = self.calc
            atoms.get_potential_energy()
            E = atoms.get_potential_energy()
            F = atoms.get_forces()
            results = {'energy': E,'forces': F}

        stri += "- Elapsed time: {}\n".format(time.time() - c_tm)
        if(results is not None):
            stri += "- SUCCESS!\n"
            stri += "- Evaluated Energy: {}".format(results['energy'])
            calc_sp = SinglePointCalculator(atoms, **results)
            atoms.calc = calc_sp
            if(self.debug): print(stri)
            print("Structure number {} evaluated\n".format(self.evalnum))
            self.evalnum+=1
            return atoms
        else:
            stri += "- FAILURE\n"
            if(self.debug): print(stri)
            return None

    def _mantains_ordering(self,atoms):
        if(len(atoms) < len(self.slab)):
            return False
        if(self.slab.symbols.indices() != atoms[:len(self.slab)].symbols.indices()):
            return False
        return True

    def _get_var_id(self,atoms) -> int:
        if(not self._mantains_ordering(atoms)): raise Exception("Does not mantain atomic ordering")
        dict = self.__atoms_to_hash(atoms[len(self.slab):])
        if(dict in  self.stc_dict):
            return self.stc_dict[dict]
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
        for i in self.atomic_types:
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

    def prepare(self, atoms):
            if( not isinstance(atoms.calc,SinglePointCalculator)):
                return None

            E = atoms.get_potential_energy()
            F = atoms.get_forces()
            results = {'energy': E,'forces': F}
            if(results is None):
                return None

            return_atoms = self.slab.copy()
            return_atoms.extend(self.sort_atoms_by_type(atoms[len(self.slab):]))
            if(self._check_overlap_all_atoms(return_atoms,self.blmin)):
                return None

            var_id = self.get_var_id(return_atoms)
            if(var_id is None):
                return None

            calc_sp = SinglePointCalculator(atoms, **results)
            return_atoms.calc =  calc_sp
            return_atoms.info['stc']= var_id
            return return_atoms

    def _check_overlap_all_atoms(self,atoms,blmin):
        indices = np.array([ a for a in np.arange(len(atoms))])
        for i in indices:
            for j in indices:
                if(i != j):
                    if(not self._check_overlap(atoms[i],atoms[j],blmin[(atoms[i].number,atoms[j].number)])):
                        return True
        return False

    def _check_overlap(self,atom1,atom2,dist):
        return dist*dist < ((atom1.position[0]-atom2.position[0]) * (atom1.position[0]-atom2.position[0]) +
                            (atom1.position[1]-atom2.position[1]) * (atom1.position[1]-atom2.position[1]) +
                            (atom1.position[2]-atom2.position[2]) * (atom1.position[2]-atom2.position[2]))