

from abc import ABC
from turtle import Turtle
from typing import List, Dict, Any
import hashlib
import json
from ..CoreUtils.NumpyArrayEncoder import NumpyArrayEncoder
from ..CoreUtils.SubunitAnalysis import NonEnergyInteratomicDistanceComparator
import numpy as np
from ase import Atoms
from ase.data import atomic_numbers
from ase.ga.utilities import closest_distances_generator,get_all_atom_types


class OperationsBase(ABC):
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    """
    def __init__(self, slab,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                rng=np.random,gen_box_size = 0.8):
        self.slab = slab
        self.variable_types = self.__get_variable_types(variable_types)
        self.variable_range = self.__get_range(variable_range)
        self.combination_matrix = self.__set_combination_matrix()
        self.variable_dict = self.__set_variable_dict()
        self.ratio_of_covalent_radii = ratio_of_covalent_radii
        self.rng = rng

        self.blmin = self.__set_blmin(slab, variable_types)
        self.p0,self.v1,self.v2,self.v3 = self._get_cell_params(slab,gen_box_size)

    def mantains_ordering(self,atoms):
        if(len(atoms) < len(self.slab)):
            return False
        if(self.slab.symbols.indices() != atoms[:len(self.slab)].symbols.indices()):
            return False

        return True

    def get_var_id(self,atoms) -> int:
        if(not self.mantains_ordering(atoms)): raise Exception("Does not mantain atomic ordering")
        if(len(atoms) == len(self.slab)):
            return self.variable_dict[0]
        dict = self.atoms_to_hash(atoms[len(self.slab):])
        if(dict in  self.variable_dict):
            return self.variable_dict[dict]
        else:
            return None

    def atoms_to_hash(self,atoms):
        if(not isinstance(atoms,Atoms)):
            raise ValueError("Tried to hash object of type different than atoms Object")
        hashed_dict = self.dict_hash(atoms.symbols.indices())
        return hashed_dict

    def dict_hash(self,dictionary: Dict[str, Any]) -> str:
        
        dhash = hashlib.sha1()
        encoded = json.dumps(dictionary, cls=NumpyArrayEncoder,sort_keys=True).encode()
        dhash.update(encoded)
        return dhash.hexdigest()

    def __set_variable_dict(self):
        symbol_dictionary = {}
        if(len(self.variable_range) != len(self.variable_types)): raise ValueError("Variable type list length and variable range list length dont match")
        try:
            for x in range(len(self.combination_matrix)):
                ats = Atoms()
                for i in range(len(self.variable_range)):
                    for j in range(self.combination_matrix[x,i]):
                        ats.extend(self.variable_types[i].copy())
                variable_id = x
                symbol_dictionary[self.atoms_to_hash(ats)] = variable_id
            return symbol_dictionary

        except:
            raise Exception("Could not generate variable dictionary: Make sure variable type length and variable range length match")

    def __set_combination_matrix(self):
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
            return combiantion_matrix
        except:
            raise Exception("Could not generate variable dictionary: Make sure variable type length and variable range length match")

    def __get_range(self,variable_range) -> List[List[int]]:
        if isinstance(variable_range,List) and isinstance(variable_range[0],List):
            for i in range(len(variable_range)):
                for j in range(len(variable_range[i])):
                    if(not isinstance(variable_range[i][j],int)):
                        raise Exception("variable_ range is not al List of List of integers")
            return variable_range
        else:
            raise Exception("variable_range is not al List of List of integers")

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

    def __set_blmin(self,slab, variable_types):
        uniques = Atoms()
        for i in variable_types:
            uniques.extend(i)
            
        unique_atom_types = get_all_atom_types(slab, uniques.numbers)

        return closest_distances_generator(atom_numbers=unique_atom_types,
                                    ratio_of_covalent_radii=self.ratio_of_covalent_radii)

    def sort_atoms_by_type(self,atoms):
        at = Atoms()
        for i in self.variable_types:
            for k in i:
                for j in atoms:
                    if(k.symbol == j.symbol):
                        at.extend(j)
        return at

    def is_structure_equal(self,atoms1,atoms2):
        if(len(atoms1) != len(atoms2)): return False
        
        comp = NonEnergyInteratomicDistanceComparator(n_top=len(atoms1), pair_cor_cum_diff=0.015,
                pair_cor_max=0.7, mic=True)
        return comp.looks_like(atoms1,atoms2)
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

    def _overlaps(self,atoms,atomstoadd,blmin):
        for at in atoms:
            for nat in atomstoadd:
                if(not self._check_overlap(at,nat,blmin[(at.number,nat.number)])):
                    return True
        return False
    
    def _normalize(self,vector):
        return vector / np.linalg.norm(vector)
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
    def _get_lens(self):
        lens = []

        for i in range(len(self.combination_matrix)):
            sums  = 0
            for k in range(len(self.combination_matrix[i])):
                sums+=self.combination_matrix[i][k] * len(self.variable_types[k])
            lens.append(sums)
        return list(lens)