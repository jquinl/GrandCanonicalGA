from asyncio import constants
from typing import List
import numpy as np
from ase import Atom, Atoms
from ase.data import atomic_numbers
from ase.ga.cutandsplicepairing import Positions
from ase.geometry import find_mic
from ase.ga.utilities import (closest_distances_generator, get_all_atom_types,atoms_too_close, atoms_too_close_two_sets)
class CrossOperation:
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    """
    def __init__(self, slab,constant,variable,variable_range,ratio_of_covalent_radii=0.7,
                minfrac = None,rng=np.random,stc_change_chance = 0.1):
        self.slab = slab
        self.constant = self.__get_atoms_object(constant)
        self.variable = self.__get_atoms_object(variable)
        self.variable_number = self.variable.numbers[0]
        self.variable_range = self.__get_range(variable_range)
        self.ratio_of_covalent_radii = ratio_of_covalent_radii
        self.stc_change_chance = stc_change_chance
        self.rng = rng
        self.minfrac = minfrac

        uniques = self.constant.copy()
        uniques.extend(self.variable)
            
        unique_atom_types = get_all_atom_types(self.slab, uniques.numbers)

        self.blmin = closest_distances_generator(atom_numbers=unique_atom_types,
                                    ratio_of_covalent_radii=self.ratio_of_covalent_radii)


    def cross(self, a1, a2):
        """Crosses the two atoms objects and returns one"""

        allowed_stc1 = a1.info['key_value_pairs']['var_stc'] 
        allowed_stc2 = a2.info['key_value_pairs']['var_stc'] 
        if (len(a1)-len(self.slab)-len(self.constant) not in self.variable_range):
            raise ValueError('Wrong size of structure a1 to optimize')
        if (len(a2)-len(self.slab)-len(self.constant) not in self.variable_range):
            raise ValueError('Wrong size of structure a2 to optimize')

        #check that a1 and a2 share a cell with initialized slef.slab
        if(self.slab.get_cell().all() != a1.get_cell().all() or self.slab.get_cell().all() != a2.get_cell().all() ):
            raise ValueError('Different cell sizes found for slab and inputed structures')

        # Only consider the atoms to optimize
        a1 = a1[len(self.slab) :len(a1)]
        a2 = a2[len(self.slab) :len(a2)]
        
        invalid = True
        counter = 0
        maxcount = 1000
        a1_copy = a1.copy()
        a2_copy = a2.copy()
        cell = self.slab.get_cell()

        results_array = [None,None]
        # Run until a valid pairing is made or maxcount pairings are tested.
        while invalid and results_array[0] is None and counter < maxcount:
            counter += 1
            
            # Choose direction of cutting plane normal
            # Will be generated entirely at random
            theta = np.pi * self.rng.rand()
            phi = 2. * np.pi * self.rng.rand()
            cut_n = np.array([np.cos(phi) * np.sin(theta),
                                np.sin(phi) * np.sin(theta), np.cos(theta)])
           
            # Randomly translate parent structures
            for a_copy, a in zip([a1_copy, a2_copy], [a1, a2]):
                a_copy.set_positions(a.get_positions())
                for i in range(len(cell)):
                    a_copy.positions += self.rng.rand() * cell[i]
                a_copy.wrap()

            # Generate the cutting point in scaled coordinates
            cosp1 = np.average(a1_copy.get_scaled_positions(), axis=0)
            cosp2 = np.average(a2_copy.get_scaled_positions(), axis=0)
            cut_p = np.zeros((1, 3))
            for i in range(3):
                cut_p[0, i] = 0.5 * (cosp1[i] + cosp2[i])
            child = self.get_pairing(a1_copy, a2_copy, cut_p, cut_n)
           
            if child is None:
                continue

            atoms  = self.slab.copy()

            atoms.extend(child)

            if atoms_too_close(atoms, self.blmin):
                continue
            if(not self.mantains_ordering(atoms)):
                continue
            if(self.get_var_stc(atoms) not in self.variable_range):
                continue
            # Passed all the tests
            atoms.wrap()
            var_stc = self.get_var_stc(atoms)
            if(self.get_var_stc(atoms) not in self.variable_range):
                continue
            if(var_stc != allowed_stc1 and var_stc != allowed_stc2):
                    atoms.info['stc']= self.get_var_stc(atoms)
                    results_array[1] = atoms
                    continue
        
            atoms.info['stc']= self.get_var_stc(atoms)
            results_array[0] = atoms

        #Cross loop ends, check results
        if(results_array[0] is not None and results_array[1] is not None):
            if(self.rng.rand() < self.stc_change_chance):
                return results_array[1]
            else:
                return results_array[0]
        elif(results_array[0] is not None and results_array[1] is None):
            return results_array[0]
        elif(results_array[0] is None and results_array[1] is not None):
            if(self.rng.rand() < self.stc_change_chance):
                return results_array[1]
            else:
                return None
        else:
            return None

    def mantains_ordering(self,atoms):
        for i in range(len(self.slab)):
            if(atoms[i].symbol != self.slab[i].symbol):
                print("Eror in ordering")
                return False
        for i in range(len(self.constant)):
            if(atoms[len(self.slab)+i].symbol != self.constant[i].symbol):
                return False
        return True

    def get_var_stc(self,atoms) -> int:
        var_stc = len(atoms)-len(self.slab)-len(self.constant)
        if(var_stc >= 0 ):
            for i in atoms[(len(self.slab)+len(self.constant)):len(atoms)]:
                if(i.symbol != self.variable[0].symbol):
                    raise Exception("Variable type of atoms does not match stored type")
        else:
            raise Exception("Negative numer of variable atoms")
        return var_stc
    def __get_range(self,variable_range) -> List[int]:
        if isinstance(variable_range,List) and isinstance(variable_range[0],int):
            return variable_range
        else:
            raise Exception("variable_ range is not al ist of integers")

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

    def get_pairing(self,a1,a2,cutting_point, cutting_normal):

        """Creates a child from two parents using the given cut.

        Returns None if the generated structure does not contain
        a large enough fraction of each parent (see self.minfrac).

        Does not check whether atoms are too close.

        Assumes the 'slab' parts have been removed from the parent
        structures. Stoichiometry agnostic"""
        atoms_result = Atoms()
        atoms_result.set_cell(self.slab.get_cell())
        a1_copy = a1.copy()
        a2_copy = a2.copy()
        "Generate constant part first"
        for num in self.constant.numbers:
            has_been_added = False
            for atoms,value in zip([a1_copy,a2_copy,a1_copy,a2_copy],[1,-1,-1,1]):
                for atom in atoms:
                    if(atom.number == num and not has_been_added):
                        at_vector =  atom.position - cutting_point
                        if(np.dot(at_vector,cutting_normal)[0] * value < 0 ):
                            atoms_result.append(atom)
                            atom.number = 200
                            has_been_added = True
        
        atoms_result.wrap()

        #Check wether the constant part has been correctly paired
        for x,y in zip(self.constant.numbers, atoms_result.numbers):
            if(x != y):
                return None
        for atoms,value in zip([a1_copy[len(atoms_result):],a2_copy[len(atoms_result):]],[1,-1]):
            for atom in atoms:
                if(atom.number == self.variable_number):
                    at_vector =  atom.position - cutting_point
                    if(np.dot(at_vector,cutting_normal)[0] * value < 0 ):
                        atoms_result.append(atom)
                        atom.number = 200
                        has_been_added = True
        atoms_result.wrap()
        return atoms_result
