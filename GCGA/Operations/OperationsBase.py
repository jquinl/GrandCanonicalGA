
from tkinter import Variable
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
    def __init__(self, slab,constant,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                rng=np.random):
        self.slab = slab
        self.constant = self.__get_atoms_object(constant)
        self.variable_types = self.__get_variable_types(variable_types)
        self.variable_range = self.__get_range(variable_range)
        self.variable_dict = self.__set_variable_dict()
        self.ratio_of_covalent_radii = ratio_of_covalent_radii
        self.rng = rng
      
        uniques = self.constant.copy()
        uniques.extend(self.variable)
            
        unique_atom_types = get_all_atom_types(self.slab, uniques.numbers)

        self.blmin = closest_distances_generator(atom_numbers=unique_atom_types,
                                    ratio_of_covalent_radii=self.ratio_of_covalent_radii)

    def mantains_ordering(self,atoms):
        if(len(atoms) < len(self.slab)+ len(self.constant)):
            return False
        if(self.slab.symbols.indices() != atoms[:len(self.slab)].symbols.indices()):
            return False
        if(self.constant.symbols.indices() != atoms[len(self.slab):len(self.slab)+len(self.constant)].symbols.indices()):
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

    def get_var_id(self,atoms) -> int:
        if(not self.mantains_ordering(atoms)): raise Exception("Does not mantain ordering of constant part")
        if(len(atoms) == len(self.slab)+ len(self.constant)):
            return self.variable_dict[0]
        try:
            var_id = self.variable_dict[atoms[len(self.slab)+ len(self.constant):].symbols.indices()]
            return var_id
        except:
            raise Exception("var_id not found in dictionary")
       

    def __set_variable_dict(self):
        symbol_dictionary = {0:0}
        if(len(self.variable_range) != len(self.variable_types)): raise ValueError("Variable type list length and variable range list length dont match")
        try:
            lengths = 1
            lengths_array = []
            for i in range(len(self.variable_types)):
                lengths_array.append(len(self.variable_range[i]))
                lengths = lengths * len(self.variable_range[i])

            combiantion_matrix  = np.zeros((lengths,len(self.variable_range)),dtype = int)

            for x in range(lengths):
                for i in range(len(self.variable_range)):
                    advancement = int(np.prod(lengths_array[i+1:len(self.variable_range)]))
                    pos = int(x / advancement) % len(self.variable_range[i])
                    combiantion_matrix[x,i]=self.variable_range[i][pos]
                
            for x in range(lengths):
                ats = Atoms()
                for i in range(len(self.variable_range)):
                    ats.extend(self.variable_types[i].copy() * combiantion_matrix[x,i])
                variable_id = x+1
                symbol_dictionary[ats.symbols.indices()] = variable_id
        except:
            raise Exception("Could not generate variable dictionary: Make sure variable type length and variable range length match")
        print(symbol_dictionary)



    def __get_range(self,variable_range) -> List[List[int]]:
        if isinstance(variable_range,List) and isinstance(variable_range[0],List):
            for i in range(len(variable_range)):
                for j in range(len(variable_range[i])):
                    if(not isinstance(variable_range[i][j],int)):
                        raise Exception("variable_ range is not al List of List of integers")
            return variable_range
        else:
            raise Exception("variable_ range is not al List of List of integers")

    def __get_variable_types(self,types) -> List[Atoms]:
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