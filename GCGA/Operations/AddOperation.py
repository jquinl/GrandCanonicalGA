import numpy as np
from ase import  Atoms
from ase.ga.utilities import atoms_too_close

from .OperationsBase import OperationsBase
class AddOperation(OperationsBase):
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    """
    def __init__(self, slab,constant,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                rng=np.random):
        super().__init__(slab,constant,variable_types,variable_range,ratio_of_covalent_radii,rng)
    
    def add(self, a1, a2):
        """Crosses the two atoms objects and returns one"""

        allowed_stc1 = a1.info['key_value_pairs']['var_stc'] 
        allowed_stc2 = a2.info['key_value_pairs']['var_stc']

        #check that a1 and a2 share a cell with initialized slef.slab
        if(self.slab.get_cell().all() != a1.get_cell().all() or self.slab.get_cell().all() != a2.get_cell().all() ):
            raise ValueError('Different cell sizes found for slab and inputed structures')

        # Only consider the atoms to optimize
        a1 = a1[len(self.slab) :len(a1)]
        a2 = a2[len(self.slab)+len(self.constant) :len(a2)]
        
        invalid = True
        counter = 0
        maxcount = 1000
        a1_copy = a1.copy()
        a2_copy = a2.copy()

        # Run until a valid pairing is made or maxcount pairings are tested.
        while invalid and counter < maxcount:
            counter += 1
        
            rand_displacement = self.rng.randint(0,allowed_stc2)
            for i in range(len(a2_copy)):
                child = self.get_addition_by_pairing(a1_copy, a2_copy,place = (i+rand_displacement)%len(a2_copy))
                if atoms_too_close(child, self.blmin):
                    child = None

            if child is None:
                child = self.get_addition_random(a1_copy)

            if child is None:
                continue

            atoms  = self.slab.copy()

            atoms.extend(child)

            if atoms_too_close(atoms, self.blmin):
                continue
            if(not self.mantains_ordering(atoms)):
                continue

            # Passed all the tests
            atoms.wrap()

            if(self.get_var_id(atoms) is None):
                continue
            atoms.info['stc']= self.get_var_id(atoms)
            
            return atoms
        return None

    def get_addition_by_pairing(self,a1,a2,place = 0):

        """Adds a variable atom object to a structure. It gets its position from a second.
        

        Does not check whether atoms are too close.

        Assumes the 'slab' parts have been removed from the parent
        structures. Stoichiometry agnostic"""
        atoms_result= a1.copy()
        atoms_result.set_cell(self.slab.get_cell())

        a2_copy = a2.copy()

        if(place > len(a2_copy)):
            return None
        #Check wether the constant part has been correctly created
        if(self.constant.symbols.indices() != atoms_result.self.symbols.indices()):
            return None

      
        atoms_result.append(a2_copy[place])
        
        atoms_result.wrap()
        return atoms_result

    def get_addition_random(self,a1):

        """Adds a variable atom object to a structure. It gets its position from a second

        Does not check whether atoms are too close.

        Assumes the 'slab' parts have been removed from the parent
        structures. Stoichiometry agnostic"""
        atoms_result= a1.copy()
        atoms_result.set_cell(self.slab.get_cell())

        #Check wether the constant part has been correctly created
        if(self.constant.symbols.indices() != atoms_result.self.symbols.indices()):
            return None

        x,y,z = self.rng.rand(),self.rng.rand(),self.rng.rand()

        at = self.rng.randint(len(self.variable_types)-1)
        atom = self.variable_types[at].copy()
        
        atom.set_cell(self.slab.get_cell())
        atom.set_scaled_positions(np.array([[x,y,z]]))
        atoms_result.extend(atom)

        return atoms_result
