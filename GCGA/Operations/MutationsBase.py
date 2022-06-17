

from abc import ABC,abstractmethod
from turtle import Turtle
from typing import List, Dict, Any
import hashlib
import json
from ..CoreUtils.NumpyArrayEncoder import NumpyArrayEncoder
from ..CoreUtils.SubunitAnalysis import NonEnergyInteratomicDistanceComparator
import numpy as np
from ase import Atoms
from ase.data import atomic_numbers


from .OperationsBase import OperationsBase
class MutationsBase(OperationsBase):
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    """
    def __init__(self, slab,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                rng=np.random,mut_box_size = 0.8):
         super().__init__(slab,variable_types,variable_range,ratio_of_covalent_radii,rng,gen_box_size=mut_box_size)

    @classmethod
    def mutation_class(cls):
        return True
    def mutation_instance(self):
        return True

    @abstractmethod
    def mutate(self, a1):
        pass

    def _check_index_overlaps(self,index,atoms,blmin):
        index_overlaps = []
        indices = np.array([ a for a in np.arange(len(atoms))])
        for i in indices:
            if(index != i):
                if(not self._check_overlap(atoms[i],atoms[index],blmin[(atoms[i].number,atoms[index].number)])):
                    index_overlaps.append(i)
        return list(index_overlaps)

    def _displace(self,index,overlap_indices,atoms):
        if(len(overlap_indices) == 0): return
        displacement_vector  = np.zeros(3)
        for i in overlap_indices:
            displacement_vector += atoms[index].position - atoms[i].position
        
        atoms[index].position += -displacement_vector
    """
    def _check_overlap_all_atoms(self,atoms,blmin):
        indices = np.array([ a for a in np.arange(len(atoms))])
        for i in indices:
            for j in indices:
                if(i != j):
                    if(not self._check_overlap(atoms[i],atoms[j],blmin[(atoms[i].number,atoms[j].number)])):
                        return True
        return False"""
    