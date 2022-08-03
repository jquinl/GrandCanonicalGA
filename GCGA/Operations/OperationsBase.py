from abc import ABC
import numpy as np
from ase import Atoms
from ase.data import atomic_numbers
from ase.ga.utilities import closest_distances_generator,get_all_atom_types


class OperationsBase(ABC):
    """
    Base class for random generation, crossover, addition, and removal operations, contains some functionalities common to all 
    """
    def __init__(self, rng=np.random):
        self.rng = rng
    
    @classmethod
    def operation_class(cls):
        pass
    def operation_instance(self):
        pass

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

    def _get_blmin(self,slab, atoms):
        uniques = Atoms()
        for i in atoms:
            if(i.symbol not in [ j.symbol for j in uniques]):
                uniques.extend(i)
            
        unique_atom_types = get_all_atom_types(slab, uniques.numbers)

        return closest_distances_generator(atom_numbers=unique_atom_types)
    
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

    def _get_cell_params(self,slab,generation_box_size = 0.8):
        "Gets cell parameters from inputed slab"
        if(generation_box_size < 0.0): raise ValueError("generation_box_size negative value")
        if(generation_box_size > 1.0): raise ValueError("generation_box_size too big")

        pos = slab.get_positions()
        cell = slab.get_cell()
        if(len(pos) == 0):
            v1 = cell[0, :] * generation_box_size
            v2 = cell[1, :] * generation_box_size
            v3 = cell[2, :] * generation_box_size
            p0 = np.array([0,0,0])
        else:
            p0 = np.array([0., 0., max(pos[:, 2]) + 2.])
            v1 = cell[0, :] * generation_box_size
            v2 = cell[1, :] * generation_box_size
            v3 = cell[2, :] * generation_box_size
            v3 = v3-p0

        return p0,v1,v2,v3

    def _random_position_in_box(self,p0,v1,v2,v3):
        return p0 + (self.rng.random() * v1 + self.rng.random() * v2 + self.rng.random() * v3)

    def _random_position_from_atom(self,atom,distance,p0,v1,v2,v3):
        
        new_vec = self.rng.normal(size=3)
        while(new_vec[0] == 0.0 and new_vec[1] == 0.0 and new_vec[2] == 0.0):
            new_vec = self.rng.normal(size=3)

        norm = self._normalize(new_vec)

        unif = self.rng.uniform(distance,distance*self.max_blen)
        norm *= unif
        pos = atom.position + norm
        pos[0] = max(min(p0[0]+v1[0]+v2[0]+v3[0], pos[0]),p0[0])
        pos[1] = max(min(p0[1]+v1[1]+v2[1]+v3[1], pos[1]),p0[1])
        pos[2] = max(min(p0[2]+v1[2]+v2[2]+v3[2], pos[2]),p0[2])
        
        return pos
    def _compare_stc(self,slab,a1,a2):
        s1 = a1[len(slab):].get_atomic_numbers()
        s2 = a2[len(slab):].get_atomic_numbers()
        if(len(s1)!=len(s2)):return False
        s1.sort()
        s2.sort()
        for a,b in zip(s1,s2):
            if(a != b): return False
        return True