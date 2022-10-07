import copy
import random
import numpy as np
from ase import Atoms

from GCGA.Operations.OperationsBase import OperationsBase
class RandomCandidateGenerator(OperationsBase):
    """
    Class instantiated when passing:

    ratio_of_covalent_radii:
        Ratio of covalent radii passed to generate the minimum distance dictionary

    random_generation_box_size:
        The volume available for the random positioning of atoms, percentage of the unoccupied space in the slab cell atoms object
        - 1.0 to use the entirety of the unoccupied cell in the slab atoms object
        - default 0.8
    rng:
        random_generator
    atom_spread:
        maximum distance from an existing atom chosen at random of a newly added atom as a multiple of the minumom bond distance. Direct correlation with the compactness of
        the generated structures
    """

    def __init__(self,random_generation_box_size = 0.8,atom_spread=2.0,rng=np.random):
        super().__init__(rng)
        self.max_blen = atom_spread
        self.random_generation_box_size = random_generation_box_size

    @classmethod
    def rand_generator(cls):
        pass
    def rand_instance(self):
        pass


    def new_candidate(self,slab,combination,atom_symbols,blmin) -> Atoms:

        if len(atom_symbols) != len(combination): raise Exception("Different array lengths in combiantion and atom_symbols arrays for the random")

        atoms = Atoms()
        new_atoms = combination
        for i in range(len(new_atoms)):
            for j in range(new_atoms[i]):
                atoms.extend(atom_symbols[i])

        if(len(atoms) == 0):
            raise Exception("Empty atoms being generated at random, revise stoichiometies")
        if(len(atoms) == 1):
            return_atoms = atoms.copy()
            return_atoms.set_cell(slab.get_cell())

        atoms_numbers  = atoms.numbers

        return self.__generate(slab,atom_numbers=atoms_numbers,blmin=blmin)


    def __generate(self,slab, atom_numbers, blmin):

        p0,v1,v2,v3 = self._get_cell_params(slab,generation_box_size = self.random_generation_box_size)

        cand = slab.copy()
        nums = copy.deepcopy(atom_numbers)
        random.shuffle(nums)
        if(len(atom_numbers) == 0):
            raise ValueError("Empty atom_numbers arrays")

        maxtries = 100
        if(len(cand)==0):
            newAtoms = Atoms(numbers = [nums[-1]])
            nums = nums[:-1]
            tries = 0
            done = False
            while tries< maxtries and not done:
                tries +=1
                for at in newAtoms:
                    at.position = self._random_position_in_box(p0,v1,v2,v3)
                cand.extend(at)
                cand.wrap()
                if(len(cand)>0):
                    done = True
            if not done:
                return None

        for i in nums:
            tries = 0
            done = False
            newAtoms = Atoms(numbers = [i])
            while tries< maxtries and not done:
                for at in newAtoms:
                    candidate = random.choice(range(len(cand)))
                    at.position = self._random_position_from_atom(cand[candidate],blmin[(cand[candidate].number,at.number)],p0,v1,v2,v3)
                
                if(not self._overlaps(cand,newAtoms,blmin)):
                    cand.extend(newAtoms)
                    done = True
                tries +=1
            if(not done):
                return None
        final_atoms = slab.copy()
        final_atoms.extend(cand[len(slab):])
        return final_atoms