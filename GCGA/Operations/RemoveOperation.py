import random
import numpy as np

from .MutationsBase import MutationsBase
class RemoveOperation(MutationsBase):

    def __init__(self, slab,variable_types,variable_range,ratio_of_covalent_radii=0.7,
            rng=np.random,jump_to_any = True):
        super().__init__(slab,variable_types,variable_range,ratio_of_covalent_radii,rng)
        
        self.combination_lens =  self._get_lens()
        self.jmp = jump_to_any

    def mutate(self, a1):
        super().mutate(a1)

        resulting_lens = []
        smallest_diff = 100
        for i in range(len(self.combination_lens)):
            if(len(a1[len(self.slab):]) > self.combination_lens[i]):
                resulting_lens.append(i)
                smallest_diff = min(smallest_diff,len(a1[len(self.slab):]) -self.combination_lens[i] )

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
            child = self.remove(a1_copy,resulting_lens[0])
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
            

    def remove(self,atoms,target_comb):
        combination = None
        try:
            combination = self.combination_matrix[target_comb]
        except:
            return None

        if(combination is None): return None
        indices = np.array([ a for a in np.arange(len(atoms))])

        new_numbers = []
        
        for k in range(len(combination)):
            for i in range(combination[k]):
                new_numbers.extend(self.variable_types[k].numbers)
                    
        final_indices = []
        random.shuffle(indices)
        for i in new_numbers:
            added = False
            for j in indices:
                if j not in final_indices and not added:
                    if(i == atoms[j].number):
                        final_indices.append(j)
                        added = True
        
        if(self.__get_len(combination) != len(new_numbers)): return None
        final_indices = np.array(final_indices)
        return atoms[final_indices].copy()


    def __get_len(self,comb):
            sums  = 0
            for k in range(len(comb)):
                sums+=comb[k] * len(self.variable_types[k])
            return sums
