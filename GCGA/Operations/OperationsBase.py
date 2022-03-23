
from typing import List
import numpy as np
from ase import Atoms
from ase.data import atomic_numbers
from ase.ga.utilities import (closest_distances_generator, get_all_atom_types)

class OperationsBase:
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    """
    def __init__(self, slab,constant,variable,variable_range,ratio_of_covalent_radii=0.7,
                rng=np.random):
        self.slab = slab
        self.constant = self.__get_atoms_object(constant)
        self.variable = self.__get_atoms_object(variable)
        self.variable_number = int(self.variable.numbers[0])
        self.variable_range = self.__get_range(variable_range)
        self.ratio_of_covalent_radii = ratio_of_covalent_radii
        self.rng = rng
      
        uniques = self.constant.copy()
        uniques.extend(self.variable)
            
        unique_atom_types = get_all_atom_types(self.slab, uniques.numbers)

        self.blmin = closest_distances_generator(atom_numbers=unique_atom_types,
                                    ratio_of_covalent_radii=self.ratio_of_covalent_radii)

    def mantains_ordering(self,atoms):
        try:
            for i in range(len(self.slab)):
                if(atoms[i].symbol != self.slab[i].symbol):
                    return False
            for i in range(len(self.constant)):
                if(atoms[len(self.slab)+i].symbol != self.constant[i].symbol):
                    return False
        except:
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