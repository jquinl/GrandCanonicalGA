
import numpy as np
import random
from ase import Atoms
from ase.ga.utilities import atoms_too_close

from .CrossBase import CrossBase
class CrossOperation(CrossBase):
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    random_translation: Applies random translation to both parent structures before mating
    """
    def __init__(self, slab,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                rng=np.random,stc_change_chance = 0.1,minfrac = 0.2,random_translation =False):
        super().__init__(slab,variable_types,variable_range,ratio_of_covalent_radii,rng)
 
        self.minfrac = self.__get_minfrac(minfrac)
        self.stc_change_chance = stc_change_chance
        self.random_translation = random_translation
        self.combination_lens =  self.__get_lens()
    
    def cross(self, a1, a2):
        super().cross( a1,a2)
        """Crosses the two atoms objects and returns one"""
        
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
           
            if(self.random_translation):
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
            atoms.extend(self.sort_atoms_by_type(child))
            if atoms_too_close(atoms, self.blmin):

                continue
            atoms.wrap()
            var_id = self.get_var_id(atoms)
            if(var_id is None):
                new_ats = self.reassign_atoms(atoms[len(self.slab):])
                if(new_ats is None):
                    continue
                atoms  = self.slab.copy()
                atoms.extend(self.sort_atoms_by_type(new_ats))
                
            atoms.info['stc']= var_id
            return atoms

        return None
    def reassign_atoms(self,atoms):
        
        candidates = []
        for i in range(len(self.combination_lens)):
            if len(atoms) == self.combination_lens[i]:
                candidates.append(i)
        if(len(candidates)== 0): return None
        cand = random.choice(candidates)
        new_atoms = Atoms()
        for i in range(len(self.combination_matrix[cand])):
            for k in range(self.combination_matrix[cand][i]):
                for i in self.variable_types[i]:
                    new_atoms.append(i)
        
        if(len(atoms)!= len(new_atoms)): return None

        nums = new_atoms.get_atomic_numbers()
        random.shuffle(nums)

        ratoms = atoms.copy()
        ratoms.set_atomic_numbers(nums)

        return ratoms



        


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

        half1 = Atoms()
        half2 = Atoms()
        half1.set_cell(self.slab.get_cell())
        half2.set_cell(self.slab.get_cell())

        "Create two halves of the system"
        for atom in a1_copy:
            at_vector =  atom.position - cutting_point
            if(np.dot(at_vector,cutting_normal)[0] > 0 ):
                half1.append(atom)
        for atom in a2_copy:
            at_vector =  atom.position - cutting_point
            if(np.dot(at_vector,cutting_normal)[0] * -1 > 0 ):
                half2.append(atom)

        if(len(half1)+len(half2) == 0):return None
        if(self.minfrac is not None):
            if(self.minfrac > float(float(len(half1))/float(len(half1)+len(half2)))): return None
            if(self.minfrac > float(float(len(half2))/float(len(half1)+len(half2)))): return None

        half1.wrap()
        half2.wrap()
        comb = half1.copy()
        comb.extend(half2.copy())
        tries = 0
        while atoms_too_close(comb,self.blmin) and tries < 10:
            comb = half1.copy()
            comb.extend(half2.copy())
            tries += 1
            half1.positions += self.rng.rand() * cutting_normal
            half1.wrap()

        if(atoms_too_close(comb,self.blmin)):
            return None

        atoms_result.extend(half1)
        atoms_result.extend(half2)
        return atoms_result 

    def __get_minfrac(self,minfrac):
        if minfrac is not None:
            if(isinstance(minfrac,float)):
                return minfrac
            else:
                raise ValueError("Specified minfrac value not a float")
        else:
            return None

    def __get_lens(self):
        lens = []
        for i in self.combination_matrix:
            sums = 0
            for k in i:
                sums+k
            lens.append(sums)
        return lens
