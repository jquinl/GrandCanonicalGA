from audioop import add
from os import remove
import random
from ase import Atoms,Atom
import numpy as np
from ase.io import write
from GCGA.Operations.RemoveOperation import RemoveOperation as RM
from GCGA.Operations.AddOperation import AddOperation as AD

from GCGA.Operations.OperationsBase import OperationsBase
class AddRemoveOperation(OperationsBase):

    def __init__(self,rng=np.random,spread = 2.0,addition_box_size = 0.8):
        super().__init__(rng)
        
        self.max_blen = spread
        self.box_size = addition_box_size
        self.add = AD(spread=spread,addition_box_size=addition_box_size,rng=rng)
        self.remove = RM(rng=rng)

    @classmethod
    def chg_class(cls):
        pass
    def chg_instance(self):
        pass

    def change(self,slab, a1,current_comb,target_combination,atom_symbols,blmin):

        a1_copy = a1.copy()
        a1_copy = a1_copy[len(slab):]

        target_add = []
        target_remove = []
        add=False
        rem=False

        for i in range(len(current_comb)):
            num = target_combination[i]-current_comb[i]
            if(num>0):
                add = True
                target_add.append(current_comb[i]+num)
                target_remove.append(current_comb[i])
            elif(num<0):
                rem=True
                target_remove.append(current_comb[i]+num)
                target_add.append(current_comb[i])
            else:
                target_remove.append(current_comb[i])
                target_add.append(current_comb[i])

        counter = 0
        maxcount = 1000
        # Run until a valid pairing is made or maxcount pairings are tested.
        while counter < maxcount:
            counter += 1
            child = a1_copy.copy()
            if(rem):
                child = self.remove.remove(slab,child,target_remove,atom_symbols,blmin)
            if(child is None):
                continue
            if(add):
                child = self.add.add(slab,child,target_add,atom_symbols,blmin)
            if(child is None):
                continue
            
            
            return_atoms = slab.copy()

            return_atoms.extend(child)
            
            if(self._check_overlap_all_atoms(return_atoms,blmin)):
                continue
            
            return return_atoms

   