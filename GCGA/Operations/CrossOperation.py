
import numpy as np
from ase import Atoms
from ase.ga.utilities import atoms_too_close

from .OperationsBase import OperationsBase
class CrossOperation(OperationsBase):
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    """
    def __init__(self, slab,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                rng=np.random,stc_change_chance = 0.1,minfrac = None):
        super().__init__(slab,variable_types,variable_range,ratio_of_covalent_radii,rng)
 
        self.minfrac = self.__get_minfrac(minfrac)
        self.stc_change_chance = stc_change_chance
    
    def mutate(self, a1, a2):
        super().mutate( a1,a2)

        """Crosses the two atoms objects and returns one"""

        allowed_stc1 = a1.info['key_value_pairs']['var_stc'] 
        allowed_stc2 = a2.info['key_value_pairs']['var_stc'] 
        
        #check that a1 and a2 share a cell with initialized slef.slab
        if(self.slab.get_cell().all() != a1.get_cell().all() or self.slab.get_cell().all() != a2.get_cell().all() ):
            raise ValueError('Different cell sizes found for slab and inputed structures')

        # Only consider the atoms to optimize
        a1 = a1[len(self.slab) :len(a1)]
        a2 = a2[len(self.slab) :len(a2)]
        
        counter = 0
        maxcount = 1000
        a1_copy = a1.copy()
        a2_copy = a2.copy()
        cell = self.slab.get_cell()

        while counter < maxcount:
            counter += 1
            
            # Choose direction of cutting plane normal
            # Will be generated entirely at random
            theta = np.pi * self.rng.rand()
            phi = 2. * np.pi * self.rng.rand()
            cut_n = np.array([np.cos(phi) * np.sin(theta),
                                np.sin(phi) * np.sin(theta), np.cos(theta)])
           
            # Randomly translate parent structures
            for a_copy, a in zip([a1_copy, a2_copy], [a1, a2]):
                a_copy.set_positions(a.get_positions())
                for i in range(len(cell)):
                    a_copy.positions += self.rng.rand() * cell[i]
                a_copy.wrap()
            
            # Generate the cutting point in scaled coordinates
            cosp1 = np.average(a1_copy.get_scaled_positions(), axis=0)
            cosp2 = np.average(a2_copy.get_scaled_positions(), axis=0)
            cut_p = np.zeros((1, 3))
            for i in range(3):
                cut_p[0, i] = 0.5 * (cosp1[i] + cosp2[i])
            child = self.get_pairing(a1_copy, a2_copy, cut_p, cut_n)
           
            if child is None:
                continue

            atoms  = self.slab.copy()

            atoms.extend(child)

            if atoms_too_close(atoms, self.blmin):
                continue

            atoms.wrap()
            var_id = self.get_var_id(atoms)
            if(var_id is None):
                continue
            # Passed all the tests if it generates a valid var_id its in the possible combination list
            if(var_id != allowed_stc1 and var_id != allowed_stc2):
                if(self.rng.rand() > self.stc_change_chance):
                    continue
        
            atoms.info['stc']= var_id
            return atoms,2

        return None,2


    def get_pairing(self,a1,a2,cutting_point, cutting_normal):

        """Creates a child from two parents using the given cut.

        Returns None if the generated structure does not contain
        a large enough fraction of each parent (see self.minfrac).

        Does not check whether atoms are too close.

        Assumes the 'slab' parts have been removed from the parent
        structures. Stoichiometry agnostic"""
        atoms_result = Atoms()
        atoms_result.set_cell(self.slab.get_cell())
        a1_copy = a1.copy()
        a2_copy = a2.copy()

        #Minfrac checkers 
        len_sys1 = 0
        len_sys2 = 0
        
        var_numbers = []
        for  i in self.variable_types:
            var_numbers.extend(i.numbers)

        for atoms,value,sys in zip([a1_copy,a2_copy],[1,-1],[1,2]):
            for atom in atoms:
                if(atom.number in var_numbers):
                    at_vector =  atom.position - cutting_point
                    if(np.dot(at_vector,cutting_normal)[0] * value < 0 ):
                        atoms_result.append(atom)
                        atom.number = 200
                        has_been_added = True
                        if(sys == 1): len_sys1 += 1
                        if(sys == 2): len_sys2 += 1
                        
        if(len(atoms_result) == 0): return None
        if(self.minfrac is not None):
            if(self.minfrac > float(float(len_sys1)/len(atoms_result))): return None
            if(self.minfrac > float(float(len_sys2)/len(atoms_result))): return None
        atoms_result.wrap()
        return atoms_result 

    def __get_minfrac(self,minfrac):
        if minfrac is not None:
            if(isinstance(minfrac,float)):
                return minfrac
            else:
                raise ValueError("Specified minfrac value not a float")
        else:
            return None
