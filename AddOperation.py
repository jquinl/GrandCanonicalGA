from asyncio import constants
from typing import List
import numpy as np
from ase import Atom, Atoms
from ase.data import atomic_numbers
from ase.ga.cutandsplicepairing import Positions
from ase.geometry import find_mic
from ase.ga.utilities import (closest_distances_generator, get_all_atom_types,atoms_too_close, atoms_too_close_two_sets)
from scipy import rand
class CrossOperation:
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    """
    def __init__(self, slab,constant,variable,variable_range,ratio_of_covalent_radii=0.7,
                rng=np.random):
        self.slab = slab
        self.constant = self.__get_atoms_object(constant)
        self.variable = self.__get_atoms_object(variable)
        self.variable_number = self.variable.numbers[0]
        self.variable_range = self.__get_range(variable_range)
        self.ratio_of_covalent_radii = ratio_of_covalent_radii
        self.rng = rng

        uniques = self.constant.copy()
        uniques.extend(self.variable)
            
        unique_atom_types = get_all_atom_types(self.slab, uniques.numbers)

        self.blmin = closest_distances_generator(atom_numbers=unique_atom_types,
                                    ratio_of_covalent_radii=self.ratio_of_covalent_radii)

    
    def add(self, a1, a2):
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

        if((allowed_stc1 + 1) not in self.variable_range):
            return None

        # Only consider the atoms to optimize
        a1 = a1[len(self.slab) :len(a1)]
        a2 = a2[len(self.slab) :len(a2)]
        
        invalid = True
        counter = 0
        maxcount = 1000
        a1_copy = a1.copy()
        a2_copy = a2.copy()

        # Run until a valid pairing is made or maxcount pairings are tested.
        while invalid and counter < maxcount:
            counter += 1
        
            rand_displacement = int(self.rng.rand(0,allowed_stc2))
            for i in range(allowed_stc2):
                child = self.get_addition_by_pairing(a1_copy, a2_copy,place = (i+rand_displacement)%allowed_stc2)
                if atoms_too_close(child, self.blmin):
                    child = None

            if child is None:
                child = self.get_addition_by_random(a1_copy)

            if child in None:
                continue

            atoms  = self.slab.copy()

            atoms.extend(child)

            if atoms_too_close(atoms, self.blmin):
                continue
            if(not self.mantains_ordering(atoms)):
                continue
            if(self.__get_var_stc(atoms) not in self.variable_range):
                continue
            # Passed all the tests
            atoms.wrap()
            var_stc = self.__get_var_stc(atoms)
            if(self.__get_var_stc(atoms) not in self.variable_range):
                continue
            if(var_stc != allowed_stc1 and var_stc != allowed_stc2):
                if(self.rng.rand() > self.stc_change_chance):
                    continue
        
            atoms.info['stc']= self.__get_var_stc(atoms)
            return atoms
        return None

    def mantains_ordering(self,atoms):
        try:
            for i in range(len(self.slab)):
                if(atoms[i].symbol != self.slab[i].symbol):
                    print("Eror in ordering")
                    return False
            for i in range(len(self.constant)):
                if(atoms[len(self.slab)+i].symbol != self.constant[i].symbol):
                    return False
        except:
            return False
        return True


    def get_addition_by_pairing(self,a1,a2,place = 0):

        """Adds a variable atom object to a structure. It gets its position from a second.
        

        Does not check whether atoms are too close.

        Assumes the 'slab' parts have been removed from the parent
        structures. Stoichiometry agnostic"""
        atoms_result= a1.copy()
        atoms_result.set_cell(self.slab.get_cell())

        a2_copy = a2.copy()
        variable_atoms = Atom()

        #Check wether the constant part has been correctly created
        for x,y in zip(self.constant.numbers, atoms_result.numbers):
            if(x != y):
                return None

        for atom in a2_copy:
            if(atom.number == self.variable_number):
                variable_atoms.append(atom)
                atom.number = 200
        
        if(len(variable_atoms) == 0):
            return None
        atoms_result.extend(variable_atoms[place])
        atoms_result.wrap()
        return atoms_result

    def get_addition_by_random(self,a1):

        """Adds a variable atom object to a structure. It gets its position from a second.
        

        Does not check whether atoms are too close.

        Assumes the 'slab' parts have been removed from the parent
        structures. Stoichiometry agnostic"""
        atoms_result= a1.copy()
        atoms_result.set_cell(self.slab.get_cell())

        #Check wether the constant part has been correctly created
        for x,y in zip(self.constant.numbers, atoms_result.numbers):
            if(x != y):
                return None
        x,y,z = self.rng.rand(),self.rng.rand(),self.rng.rand()
        atom = Atoms(self.variable_number)
        atom.set_cell(self.slab.get_cell())
        atom.set_scaled_positions(np.array([[x,y,z]]))
        atoms_result.append(atom)

        return atoms_result

    def __get_var_stc(self,atoms) -> int:
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
