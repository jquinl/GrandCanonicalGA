import random
import numpy as np

from GCGA.Operations.OperationsBase import OperationsBase
class RemoveOperation(OperationsBase):

    def __init__(self,rng=np.random):
        super().__init__(rng)

    @classmethod
    def remove_class(cls):
        pass
    def remove_instance(self):
        pass

    def remove(self,slab, a1,combination,atom_symbols,blmin):

        a1_copy = a1.copy()
        a1_copy = a1_copy[len(slab):]

        counter = 0
        maxcount = 1000
        # Run until a valid pairing is made or maxcount pairings are tested.
        while counter < maxcount:
            counter += 1

            child = self.remove_atoms(a1_copy,combination,atom_symbols)
            if(child is None):
                continue
            
            return_atoms = slab.copy()

            return_atoms.extend(child)
            
            if(self._check_overlap_all_atoms(return_atoms,blmin)):
                continue

            return return_atoms

    def remove_atoms(self,atoms,combination,atom_symbols):

        indices = np.array([ a for a in np.arange(len(atoms))])

        new_numbers = []
        for k in range(len(combination)):
            for i in range(combination[k]):
                new_numbers.extend(atom_symbols[k].numbers)

        final_indices = []
        random.shuffle(indices)
        for i in new_numbers:
            added = False
            for j in indices:
                if j not in final_indices and not added:
                    if(i == atoms[j].number):
                        final_indices.append(j)
                        added = True

        if(self.__get_len(combination,atom_symbols) != len(new_numbers)): return None
        final_indices = np.array(final_indices)
        return atoms[final_indices].copy()


    def __get_len(self,comb,atom_symbols):
        sums  = 0
        for k in range(len(comb)):
            sums+=comb[k] * len(atom_symbols[k])
        return sums
