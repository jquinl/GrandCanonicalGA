from audioop import add
import random
from ase import Atoms,Atom
import numpy as np
from ase.io import write

from .MutationsBase import MutationsBase
class AddOperation(MutationsBase):

    def __init__(self, slab,variable_types,variable_range,ratio_of_covalent_radii=0.7,
            rng=np.random,max_bond_multiplier = 2.0,jump_to_any = True,addition_box_size = 0.8):
        super().__init__(slab,variable_types,variable_range,ratio_of_covalent_radii,rng,mut_box_size=addition_box_size)
        
        self.combination_lens =  self._get_lens()
        self.jmp = jump_to_any
        self.max_blen = max_bond_multiplier

    def mutate(self, a1):
        super().mutate(a1)

        resulting_lens = []
        smallest_diff = 100
        for i in range(len(self.combination_lens)):
            if(len(a1[len(self.slab):]) < self.combination_lens[i]):
                resulting_lens.append(i)
                smallest_diff = min( smallest_diff, self.combination_lens[i]-len(a1[len(self.slab):]))

        if(len(resulting_lens) == 0): return None

        if (not self.jmp):
            for i in resulting_lens:
                resulting_lens = [i for i in resulting_lens if smallest_diff == len(a1[len(self.slab):]) -self.combination_lens[i] ]

        a1_copy = a1.copy()
        a1_copy = a1_copy[len(self.slab):]

        counter = 0
        maxcount = 1000
        # Run until a valid pairing is made or maxcount pairings are tested.
        while counter < maxcount:
            counter += 1

            random.shuffle(resulting_lens)

            child = self.add(a1_copy,resulting_lens[0])

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

    