
import numpy as np

from .OperationsBase import OperationsBase
class RemoveOperation(OperationsBase):
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    """
    def __init__(self, slab,constant,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                rng=np.random):
        super().__init__(slab,constant,variable_types,variable_range,ratio_of_covalent_radii,rng)

    def remove(self, a1,index=None):
        """If index is not provided removes a random atom of the variable type, If it is provided removes the atom at that index as long as its of variable type """

        if (len(a1)-len(self.slab)-len(self.constant) not in self.variable_range):
            raise ValueError('Wrong size of structure a1 to optimize')

        #check that a1 and a2 share a cell with initialized slef.slab
        if(self.slab.get_cell().all() != a1.get_cell().all()):
            raise ValueError('Different cell sizes found for slab and inputed structures')

        if((a1.info['key_value_pairs']['var_stc']  - 1) not in self.variable_range):
            return None

        # Only consider the atoms to optimize
        a1 = a1[len(self.slab) :len(a1)]
        
        counter = 0
        maxcount = 1000
        a1_copy = a1.copy()

        poppable_indices = np.array([ a for a in np.arange(len(a1)) if a1.numbers[a] == self.variable_number])
        # Run until a valid pairing is made or maxcount pairings are tested.
        while counter < maxcount:
            counter += 1
            
            child = a1_copy.copy()
            pop_int = self.rng.randint(0,len(poppable_indices)-1)
            rand_atm = poppable_indices[pop_int]
            
            if(index is not None and index not in poppable_indices):
                raise ValueError("Provided index does not belong to a removable atom")
            elif(index is not None and index  in poppable_indices):
                child.pop(poppable_indices[index])
            else:
                child.pop(poppable_indices[rand_atm])

            atoms  = self.slab.copy()

            atoms.extend(child)

            if(not self.mantains_ordering(atoms)):
                continue

            # Passed all the tests
            atoms.wrap()

            atoms.info['stc']= self.get_var_stc(atoms)
            return atoms
        return None
