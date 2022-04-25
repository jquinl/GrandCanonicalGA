import numpy as np
from ase.ga.utilities import atoms_too_close

from .OperationsBase import OperationsBase
class PermutationOperation(OperationsBase):
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    """
    def __init__(self, slab,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                rng=np.random,delete_chance = 0.0,):
        super().__init__(slab,variable_types,variable_range,ratio_of_covalent_radii,rng)
 
        self.delete_chance = delete_chance



    def mutate(self, a1, a2):
        super().mutate(a1,a2)

        """Swaps two atoms of different type
        """

        chg_chance = self.delete_chance
        
        #check that a1 and a2 share a cell with initialized slef.slab
        if(self.slab.get_cell().all() != a1.get_cell().all()):
            raise ValueError('Different cell sizes found for slab and inputed structures')

        # Only consider the atoms to optimize
        a1 = a1[len(self.slab) :len(a1)]

        indices = np.array([ a for a in np.arange(len(a1))])

        unique_types =[]
        for a in a1.numbers:
            if(a not in unique_types):
                unique_types.append(a)
        
        if(len(unique_types) < 2):
            raise Exception("Provided atoms object only contains one type of atom, unable to permute")

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
            child = self.swap(a1,ind1,ind2) 
            
            if child is None:
                continue

            #Random chance of deleting one of the swapped atoms
            if(self.rng.rand() < chg_chance):
                rnd = self.rng.rand()
                if(rnd > 0.5):
                    child.pop(ind1)
                else:
                    child.pop(ind2)

            if child is None:
                continue
            atoms  = self.slab.copy()
            atoms.extend(child)
            atoms.wrap()

            #Does not check for minimum distance on swap
            #if atoms_too_close(atoms, self.blmin):
            #   continue
            var_id = self.get_var_id(atoms)
            if( var_id  is None):
                continue
            # Passed all the tests
        
            atoms.info['stc']= var_id 
            return atoms,1

        return None,1

    def swap(self,at,i1,i2):
        atoms = at.copy()
        pos = atoms.get_positions()
    
        atoms[i1].position = pos[i2]
        atoms[i2].position = pos[i1]
        return atoms


