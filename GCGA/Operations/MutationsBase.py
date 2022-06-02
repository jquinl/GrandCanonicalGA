

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
from ase.ga.utilities import (closest_distances_generator, get_all_atom_types)


from .OperationsBase import OperationsBase
class MutationsBase(OperationsBase):
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    """
    def __init__(self, slab,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                rng=np.random):
         super().__init__(slab,variable_types,variable_range,ratio_of_covalent_radii,rng)

    @classmethod
    def mutation_class(cls):
        return True
    def mutation_instance(self):
        return True

    @abstractmethod
    def mutate(self, a1):
        pass