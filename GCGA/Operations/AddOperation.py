from audioop import add
import random
from ase import Atoms,Atom
import numpy as np
from ase.io import write

from GCGA.Operations.OperationsBase import OperationsBase
class AddOperation(OperationsBase):

    def __init__(self, ratio_of_covalent_radii=0.7,
            rng=np.random,spread = 2.0,addition_box_size = 0.8):
        super().__init__(ratio_of_covalent_radii,rng)
        
        self.max_blen = spread
        self.box_size = addition_box_size

    @classmethod
    def add_class(cls):
        pass
    def add_instance(self):
        pass


    def add(self,slab, a1,combination,atom_symbols,blmin):

        a1_copy = a1.copy()
        a1_copy = a1_copy[len(slab):]

        counter = 0
        maxcount = 1000
        # Run until a valid pairing is made or maxcount pairings are tested.
        while counter < maxcount:
            counter += 1

            child = self.add_atoms(slab,a1_copy,combination,atom_symbols,blmin)

            if(child is None):
                continue
            
            return_atoms = slab.copy()

            return_atoms.extend(child)
            
            if(self._check_overlap_all_atoms(return_atoms,blmin)):
                continue
            
            return return_atoms

    def add_atoms(self,slab,atoms,combination,atom_symbols,blmin):

        if len(atom_symbols) != len(combination): raise Exception("Different array lengths in combiantion and atom_symbols arrays for the random")
        numbers = []
        
        p0,v1,v2,v3 = self._get_cell_params(slab,generation_box_size = self.box_size)

        for k in range(len(combination)):
            for i in range(combination[k]):
                numbers.extend(atom_symbols[k].numbers)

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
                i.position  = self._random_position_from_atom(cand[candidate],blmin[(cand[candidate].number,i.number)],p0,v1,v2,v3)

            if(not self._overlaps(cand,newAtoms,blmin)):
                    cand.extend(newAtoms)
                    added = True
        if not added:
            return None
        return cand


    