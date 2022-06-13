import numpy as np
from ase.ga.utilities import (closest_distances_generator)
import random


from .MutationsBase import MutationsBase
class PermutationOperation(MutationsBase):

    def __init__(self,slab,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                   rng=np.random):
        super().__init__(slab,variable_types,variable_range,ratio_of_covalent_radii,rng)

        

    def mutate(self, a1):
        super().mutate(a1)
        #check that a1 and a2 share a cell with initialized slef.slab
        if(self.slab.get_cell().all() != a1.get_cell().all()):
            raise ValueError('Different cell sizes found for slab and inputed structures')

        unique_types =[]
        for a in a1.numbers:
            if(a not in unique_types):
                unique_types.append(a)
        # Only consider the atoms to optimize
        a1_copy = a1.copy()
        a1_copy = a1_copy[len(self.slab) :len(a1_copy)]

        indices = np.array([ a for a in np.arange(len(a1_copy))])

        if(len(unique_types) < 2):
            print("Provided atoms object only contains one type of atom, unable to permute,returning initial structure")
            return a1

        counter = 0
        maxcount = 1000
        
        # Run until a valid pairing is made or maxcount pairings are tested.
        while counter < maxcount:
            counter += 1
            i1 = self.rng.randint(0,len(indices)-1)
            i2 = self.rng.randint(0,len(indices)-1)
            while(a1[indices[i1]].number == a1[indices[i2]].number):
                i2 = self.rng.randint(0,len(indices)-1)
            
            ind1 = indices[i1]
            ind2 = indices[i2]
            child = self.swap(a1,ind1,ind2,self.blmin)
            if(child is None):
                continue

            return_atoms = self.slab.copy()

            return_atoms.extend(self.sort_atoms_by_type(child))
            
            if(self._check_overlap_all_atoms(return_atoms,self.blmin)):
                continue
            
            var_id = self.get_var_id(return_atoms)
            if(var_id is not None):
                return_atoms.info['stc']= var_id
                return return_atoms
            else:
                raise Exception("Provided atomic combination is not present in combination matrix")


    def swap(self,at,i1,i2,blmin):
        atoms = at.copy()
        pos = atoms.get_positions()
    
        atoms[i1].position = pos[i2]
        atoms[i2].position = pos[i1]

        overlap = True
        tries = 0
        while(overlap and tries<10):
            tries += 1
            
            a1overlap = self._check_index_overlaps(i1,atoms,blmin)
            a2overlap = self._check_index_overlaps(i2,atoms,blmin)
            overlap = not (len(a1overlap) != 0 and len(a2overlap) != 0)

            self._displace(i1,a1overlap,atoms)
            self._displace(i2,a2overlap,atoms)
        
        if(self._check_overlap_all_atoms(atoms,blmin)):
            return None

        return atoms


   