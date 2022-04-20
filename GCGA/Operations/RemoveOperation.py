
import numpy as np

from .OperationsBase import OperationsBase
class RemoveOperation(OperationsBase):
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    """
    def __init__(self, slab,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                rng=np.random):
        super().__init__(slab,variable_types,variable_range,ratio_of_covalent_radii,rng)

    def mutate(self, a1,a2):
        super().mutate( a1,a2)
        """If index is not provided removes a random atom of the variable type, If it is provided removes the atom at that index as long as its of variable type """

        #check that a1 and a2 share a cell with initialized slef.slab
        if(self.slab.get_cell().all() != a1.get_cell().all()):
            raise ValueError('Different cell sizes found for slab and inputed structures')

        # Only consider the atoms to optimize
        a1 = a1[len(self.slab) :len(a1)]
        
        counter = 0
        maxcount = 1000
        a1_copy = a1.copy()

        poppable_indices = np.array([ a for a in np.arange(len(a1))])
        # Run until a valid pairing is made or maxcount pairings are tested.
        while counter < maxcount:
            counter += 1
            
            child = a1_copy.copy()
            pop_int = self.rng.randint(0,len(poppable_indices)-1)
            rand_atm = poppable_indices[pop_int]
            
            child.pop(poppable_indices[rand_atm])

            atoms  = self.slab.copy()

            atoms.extend(child)

            # Passed all the tests
            atoms.wrap()
            var_id = self.get_var_id(atoms)
            if(var_id is None):
                continue
            atoms.info['stc']= var_id
            return atoms,1
        return None,1
