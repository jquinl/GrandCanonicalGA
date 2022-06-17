
from itertools import count,chain
import numpy as np
import random
from ase import Atoms

from .CrossBase import CrossBase
class CrossOperation(CrossBase):
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    random_translation: Applies random translation to both parent structures before mating
    """
    def __init__(self, slab,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                rng=np.random,stc_change_chance = 0.3,minfrac = 0.2,random_translation =False,box_size=0.8):
        super().__init__(slab,variable_types,variable_range,ratio_of_covalent_radii,rng)
 
        self.minfrac = self.__get_minfrac(minfrac)
        self.stc_change_chance = stc_change_chance
        self.random_translation = random_translation
        self.combination_lens =  self._get_lens()
        self.max_blen = 2.0
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
        stcchg = False
        if(self.rng.random()<self.stc_change_chance): stcchg = True

        while counter < maxcount:
            counter += 1
            # Choose direction of cutting plane normal
            # Will be generated entirely at random
            
            child = self.get_new_candidate(a1_copy, a2_copy)
           
            if child is None:
                continue

            atoms  = self.slab.copy()
            atoms.extend(self.sort_atoms_by_type(child))
            if self._check_overlap_all_atoms(atoms, self.blmin):
                continue
            atoms.wrap()
            var_id = self.get_var_id(atoms)
            

            if(var_id is None):
                new_ats = self.reassign_atoms(atoms[len(self.slab):])
                if(new_ats is None):
                    continue
                atoms  = self.slab.copy()
                atoms.extend(self.sort_atoms_by_type(new_ats))
            if(var_id != a1.info['key_value_pairs']['var_stc'] and stcchg == False):
                
                cand = self.add(atoms[len(self.slab):],a1.info['key_value_pairs']['var_stc'])
                atoms  = self.slab.copy()
                atoms.extend(self.sort_atoms_by_type(cand))


            if self._check_overlap_all_atoms(atoms, self.blmin):
                continue
            var_id = self.get_var_id(atoms)
            
            if(var_id is None):
                continue
            if(var_id != a1.info['key_value_pairs']['var_stc'] and stcchg == False):
                continue
            print(var_id)
            print(a1.info['key_value_pairs']['var_stc'])
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
                for z in self.variable_types[i]:
                    new_atoms.append(z)
        
        
        if(len(atoms)!= len(new_atoms)): return None

        nums = new_atoms.get_atomic_numbers()
        random.shuffle(nums)

        ratoms = atoms.copy()
        ratoms.set_atomic_numbers(nums)
        return ratoms


    def get_new_candidate(self,a1,a2):
        a1_copy = a1.copy()
        a2_copy = a2.copy()

        cm = (a1_copy.get_center_of_mass() + a2_copy.get_center_of_mass())/2.0
        theta = self.rng.rand() * 2 * np.pi 
        phi = self.rng.rand() * np.pi  
        e = np.array((np.sin(phi) * np.cos(theta),
                      np.sin(theta) * np.sin(phi),
                      np.cos(phi)))

        a1_copy.translate(-a1_copy.get_center_of_mass())
        a2_copy.translate(-a2_copy.get_center_of_mass())



        fmap = [np.dot(x, e) for x in a1_copy.get_positions()]
        mmap = [-np.dot(x, e) for x in a2_copy.get_positions()]
        ain = sorted([i for i in chain(fmap, mmap) if i > 0],
                     reverse=True)
        aout = sorted([i for i in chain(fmap, mmap) if i < 0],
                      reverse=True)

        try:
            if(len(ain)-min(self.combination_lens) < 0):
                off = len(ain)-random.choice(self.combination_lens)
                dist = (abs(aout[abs(off) - 1]) + abs(aout[abs(off)])) * .5
                a1_copy.translate(e * dist)
                a2_copy.translate(-e * dist)
            elif(len(ain)-max(self.combination_lens) >0):
                off = len(ain)-random.choice(self.combination_lens)
                dist = (abs(aout[abs(off) - 1]) + abs(aout[abs(off)])) * .5
                a1_copy.translate(e * dist)
                a2_copy.translate(-e * dist)
        except:
            pass
        

        fmap = [np.dot(x, e) for x in a1_copy.get_positions()]
        mmap = [-np.dot(x, e) for x in a2_copy.get_positions()]
        ain = sorted([i for i in chain(fmap, mmap) if i > 0],
                     reverse=True)
        aout = sorted([i for i in chain(fmap, mmap) if i < 0],
                      reverse=True)

        if(len(ain) not in self.combination_lens): return None
        tmpf, tmpm = Atoms(), Atoms()
        for atom in a1_copy:
            if np.dot(atom.position, e) > 0:
                tmpf.append(atom)
        for atom in a2_copy:
            if np.dot(atom.position, e) < 0:
                tmpm.append(atom)

        ratoms = Atoms()
        ratoms.set_cell(self.slab.get_cell())
        ratoms.extend(tmpf)
        ratoms.extend(tmpm)
        if(len(ratoms) == 0): return None
        ratoms.translate(cm)

        return ratoms


    def __get_minfrac(self,minfrac):
        if minfrac is not None:
            if(isinstance(minfrac,float)):
                return minfrac
            else:
                raise ValueError("Specified minfrac value not a float")
        else:
            return None


    def _get_lens(self):
        lens = []

        for i in range(len(self.combination_matrix)):
            sums  = 0
            for k in range(len(self.combination_matrix[i])):
                sums+=self.combination_matrix[i][k] * len(self.variable_types[k])
            lens.append(sums)
        return list(lens)

    def add(self,atoms,target_comb):
        combination = None
        try:
            combination = self.combination_matrix[target_comb]
        except:
            return None

        if(combination is None): return None

        numbers = []
        
        for k in range(len(combination)):
            for i in range(combination[k]):
                numbers.extend(self.variable_types[k].numbers)

        at_nums = atoms.get_atomic_numbers()
        for i in range(len(numbers)):
            for k in range(len(at_nums)):
                if(numbers[i] == at_nums[k] and numbers[i] != 200):
                    at_nums[k] = 200
                    numbers[i] = 200
        final_numbers = [i for i in numbers if i!= 200]
        
        additional_atoms = []
        newAtoms = Atoms(numbers = final_numbers)

        tries = 0
        added = False
        cand = None
        while tries < 10 and added == False:
            cand = atoms.copy()
            random.shuffle(additional_atoms)
            for i in newAtoms:
                candidate = random.choice(range(len(cand)))
                i.position  = self._random_position_from_atom(cand[candidate],self.blmin[(cand[candidate].number,i.number)])

            if(not self._overlaps(cand,newAtoms,self.blmin)):
                    cand.extend(newAtoms)
                    added = True
        if not added:
            return None
        return cand

    def _random_position_from_atom(self,atom,distance):
        new_vec = self.rng.normal(size=3)
        while(new_vec[0] == 0.0 and new_vec[1] == 0.0 and new_vec[2] == 0.0):
            new_vec = self.rng.normal(size=3)

        norm = self._normalize(new_vec)

        unif = self.rng.uniform(distance,distance*self.max_blen)
        norm *= unif
        pos = atom.position + norm
        pos[0] = max(min(self.p0[0]+self.v1[0]+self.v2[0]+self.v3[0], pos[0]),self.p0[0])
        pos[1] = max(min(self.p0[1]+self.v1[1]+self.v2[1]+self.v3[1], pos[1]),self.p0[1])
        pos[2] = max(min(self.p0[2]+self.v1[2]+self.v2[2]+self.v3[2], pos[2]),self.p0[2])
        
        return pos

