import numbers
from tokenize import Number
import numpy as np
from ase import Atoms
from ase.ga.utilities import atoms_too_close

from .OperationsBase import OperationsBase
class PermutationOperation(OperationsBase):
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    """
    def __init__(self, slab,constant,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                rng=np.random,delete_chance = 0.0,):
        super().__init__(slab,constant,variable_types,variable_range,ratio_of_covalent_radii,rng)
 
        self.delete_chance = delete_chance

    def permutate(self, a1):
        """Swaps two atoms of different type
            If variable_by_constant is True, will swap variable by constant positions
            If False, might swap atomic position between 
        """

        allowed_stc1 = a1.info['key_value_pairs']['var_stc'] 
        chg_chance = self.delete_chance
        if(allowed_stc1 -1 not in self.variable_range):
            chg_chance = 0.0
        if (len(a1)-len(self.slab)-len(self.constant) not in self.variable_range):
            raise ValueError('Wrong size of structure a1 to optimize')

        #check that a1 and a2 share a cell with initialized slef.slab
        if(self.slab.get_cell().all() != a1.get_cell().all()):
            raise ValueError('Different cell sizes found for slab and inputed structures')

        # Only consider the atoms to optimize
        a1 = a1[len(self.slab) :len(a1)]

        all_indices = np.array([ a for a in np.arange(len(a1))])
        variable_indices = np.array([ a for a in np.arange(len(a1)) if a1.numbers[a] == self.variable_number])
        

        unique_types =[]
        for a in a1.numbers:
            if(a not in unique_types):
                unique_types.append(a)
        if(len(unique_types) < 2):
            raise Exception("Provided atoms object only contains one type of atom, unable to permute")
        invalid = True
        counter = 0
        maxcount = 1000
        
        # Run until a valid pairing is made or maxcount pairings are tested.
        while invalid and counter < maxcount:
            counter += 1
            
            i1 = self.rng.randint(0,len(all_indices)-1)
            i2 = self.rng.randint(0,len(all_indices)-1)
             
            while(a1[all_indices[i1]].number == a1[all_indices[i2]].number):
                i2 = self.rng.randint(0,len(all_indices)-1)

            ind1 = all_indices[i1]
            ind2 = all_indices[i2]
            child = self.swap(a1,ind1,ind2) 
            
            if child is None:
                continue
            
            if(ind1 in variable_indices or ind2 in variable_indices):
                if(self.rng.rand() < chg_chance):
                    if(ind1 in variable_indices):
                        child.pop(ind1)
                    elif(ind2 in variable_indices):
                        child.pop(ind2)
            if child is None:
                continue
            atoms  = self.slab.copy()

            atoms.extend(child)

            atoms.wrap()

            #Does not check for minimum distance on swap
            #if atoms_too_close(atoms, self.blmin):
            #   continue
            if(not self.mantains_ordering(atoms)):
                continue
            var_stc = self.get_var_stc(atoms)
            if(var_stc not in self.variable_range):
                continue
            # Passed all the tests
        
            atoms.info['stc']= self.get_var_stc(atoms)
            return atoms

        return None


    def swap(self,at,i1,i2):
        atoms = at.copy()
        pos = atoms.get_positions()
    
        atoms[i1].position = pos[i2]
        atoms[i2].position = pos[i1]
        return atoms


