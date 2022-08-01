

from abc import ABC,abstractmethod
from turtle import Turtle
from typing import List, Dict, Any
import hashlib
import json
import numpy as np
from ase import Atoms
from ase.data import atomic_numbers


from .OperationsBase import OperationsBase

class CrossBase(OperationsBase):
    """
    Abstract class for the CrossOver operations, if a different implementation is needed inherit from this and implement the cross method
    """
    def __init__(self, slab,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                rng=np.random,box_size = 0.8):
        super().__init__(slab,variable_types,variable_range,ratio_of_covalent_radii,rng,gen_box_size=box_size)

    @classmethod
    def cross_class(cls):
        pass
    def cross_instance(self):
        pass

    @abstractmethod
    def cross(self, a1,a2):
        pass

