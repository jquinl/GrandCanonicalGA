from abc import ABC
import numpy as np
from ase.data import atomic_numbers


class OperationsBase(ABC):
    """
    Base class for random generation, crossover, addition, and removal operations, contains some functionalities common to all 
    """
    def __init__(self, ratio_of_covalent_radii=0.7,
                rng=np.random):
        self.ratio_of_covalent_radii = ratio_of_covalent_radii
        self.rng = rng
    
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
    
    def _get_lens(self):
        lens = []

        for i in range(len(self.combination_matrix)):
            sums  = 0
            for k in range(len(self.combination_matrix[i])):
                sums+=self.combination_matrix[i][k] * len(self.variable_types[k])
            lens.append(sums)
        return list(lens)