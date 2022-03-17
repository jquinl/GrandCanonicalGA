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


    def remove(self, a1):
        """Crosses the two atoms objects and returns one"""

        allowed_stc1 = a1.info['key_value_pairs']['var_stc'] 

        
        if (len(a1)-len(self.slab)-len(self.constant) not in self.variable_range):
            raise ValueError('Wrong size of structure a1 to optimize')

        #check that a1 and a2 share a cell with initialized slef.slab
        if(self.slab.get_cell().all() != a1.get_cell().all()):
            raise ValueError('Different cell sizes found for slab and inputed structures')

        if((allowed_stc1 - 1) not in self.variable_range):
            return None

        # Only consider the atoms to optimize
        a1 = a1[len(self.slab) :len(a1)]
        
        invalid = True
        counter = 0
        maxcount = 1000
        a1_copy = a1.copy()

        poppable_indices = np.array([ a for a in np.arange(len(a1)) if a1.numbers[a] == self.variable_number])
       
        # Run until a valid pairing is made or maxcount pairings are tested.
        while invalid and counter < maxcount:
            counter += 1
            
            child = a1_copy().copy()

            rand_atm = int(self.rng.rand(0,len(poppable_indices)))
            
            child.pop(poppable_indices[rand_atm])

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


