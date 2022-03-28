
import numpy as np
from ase import Atoms
from ase.ga.utilities import atoms_too_close

from .OperationsBase import OperationsBase
class CrossOperation(OperationsBase):
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    """
    def __init__(self, slab,constant,variable,variable_range,ratio_of_covalent_radii=0.7,
                rng=np.random,stc_change_chance = 0.1,minfrac = None):
        super().__init__(slab,constant,variable,variable_range,ratio_of_covalent_radii,rng)
 
        self.minfrac = self.__get_minfrac(minfrac)
        self.stc_change_chance = stc_change_chance


    def cross(self, a1, a2):
        """Crosses the two atoms objects and returns one"""

        allowed_stc1 = a1.info['key_value_pairs']['var_stc'] 
        allowed_stc2 = a2.info['key_value_pairs']['var_stc'] 
        if (len(a1)-len(self.slab)-len(self.constant) not in self.variable_range):
            raise ValueError('Wrong size of structure a1 to optimize')
        if (len(a2)-len(self.slab)-len(self.constant) not in self.variable_range):
            raise ValueError('Wrong size of structure a2 to optimize')

        #check that a1 and a2 share a cell with initialized slef.slab
        if(self.slab.get_cell().all() != a1.get_cell().all() or self.slab.get_cell().all() != a2.get_cell().all() ):
            raise ValueError('Different cell sizes found for slab and inputed structures')

        # Only consider the atoms to optimize
        a1 = a1[len(self.slab) :len(a1)]
        a2 = a2[len(self.slab) :len(a2)]
        
        invalid = True
        counter = 0
        maxcount = 1000
        a1_copy = a1.copy()
        a2_copy = a2.copy()
        cell = self.slab.get_cell()

        # Run until a valid pairing is made or maxcount pairings are tested.
        while invalid and counter < maxcount:
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
            if(not self.mantains_ordering(atoms)):
                continue
            # Passed all the tests
            atoms.wrap()
            var_stc = self.get_var_stc(atoms)
            if(var_stc not in self.variable_range):
                continue
            if(var_stc != allowed_stc1 and var_stc != allowed_stc2):
                if(self.rng.rand() > self.stc_change_chance):
                    continue
        
            atoms.info['stc']= self.get_var_stc(atoms)
            print(a1.info['key_value_pairs']['parent_penalty'])
            return atoms

        return None


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
        "Generate constant part first"

        for num in self.constant.numbers:
            has_been_added = False
            #for atoms,value in zip([a1_copy,a2_copy,a1_copy,a2_copy],[1,-1,-1,1]):
            len_sys1 = 0
            len_sys2 = 0
            for atoms,value,sys in zip([a1_copy,a2_copy],[1,-1],[1,2]):
                for atom in atoms:
                    if(atom.number == num and not has_been_added):
                        at_vector =  atom.position - cutting_point
                        if(np.dot(at_vector,cutting_normal)[0] * value < 0 ):
                            atoms_result.append(atom)
                            atom.number = 200
                            has_been_added = True
                            if(sys == 1): len_sys1 += 1
                            if(sys == 2): len_sys2 += 1
        
        if(self.minfrac is not None):
            if(self.minfrac > float(float(len_sys1)/len(self.constant))): return None
            if(self.minfrac > float(float(len_sys2)/len(self.constant))): return None
        
        atoms_result.wrap()

        #Check wether the constant part has been correctly paired
        for x,y in zip(self.constant.numbers, atoms_result.numbers):
            if(x != y):
                return None
        for atoms,value in zip([a1_copy[len(atoms_result):],a2_copy[len(atoms_result):]],[1,-1]):
            for atom in atoms:
                if(atom.number == self.variable_number):
                    at_vector =  atom.position - cutting_point
                    if(np.dot(at_vector,cutting_normal)[0] * value < 0 ):
                        atoms_result.append(atom)
                        atom.number = 200
                        has_been_added = True
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
