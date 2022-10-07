
from itertools import count,chain
import numpy as np
import random
from ase import Atoms
from ase.io import write

from GCGA.Operations.OperationsBase import OperationsBase
class CrossOperation(OperationsBase):
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    random_translation: Applies random translation to both parent structures before mating
    """
    def __init__(self,rng=np.random,minfrac = 0.2,generation_box_size=0.8,new_atom_spread = 2.0):
        super().__init__(rng)

        self.minfrac = self.__get_minfrac(minfrac)
        self.max_blen = new_atom_spread
        self.gen_box_size = generation_box_size

    @classmethod
    def cross_class(cls):
        pass
    def cross_instance(self):
        pass


    def cross(self, slab, a1, a2,blmin):
        """Crosses the two atoms objects and returns one"""

        #check that a1 and a2 share a cell with initialized slef.slab
        if(slab.get_cell().all() != a1.get_cell().all() or slab.get_cell().all() != a2.get_cell().all() ):
            raise ValueError('Different cell sizes found for slab and inputed structures')

        # Only consider the atoms to optimize
        a1 = a1[len(slab) :len(a1)]
        a2 = a2[len(slab) :len(a2)]

        counter = 0
        maxcount = 1000
        a1_copy = a1.copy()
        a2_copy = a2.copy()
        while counter < maxcount:
            counter += 1
            # Choose direction of cutting plane normal
            # Will be generated entirely at random
            child = self.get_new_candidate(slab,a1_copy, a2_copy)
            if child is None:
                continue

            atoms  = slab.copy()
            atoms.extend(child)
            if self._check_overlap_all_atoms(atoms, blmin):
                continue
            atoms.wrap()
            if self._check_overlap_all_atoms(atoms, blmin):
                continue
            if(not self._compare_stc(slab,a1_copy,child)):
                continue

            return atoms

        return None

    def get_new_candidate(self,slab,a1,a2):
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
            if(len(ain)-len(a1) < 0):
                off = len(ain)-random.choice([len(a1),len(a2)])
                dist = (abs(aout[abs(off) - 1]) + abs(aout[abs(off)])) * .5
                a1_copy.translate(e * dist)
                a2_copy.translate(-e * dist)
            elif(len(ain)-len(a1) >0):
                off = len(ain)-random.choice([len(a1),len(a2)])
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

        tmpf, tmpm = Atoms(), Atoms()
        for atom in a1_copy:
            if np.dot(atom.position, e) > 0:
                tmpf.append(atom)
        for atom in a2_copy:
            if np.dot(atom.position, e) < 0:
                tmpm.append(atom)
        ratoms = Atoms()
        ratoms.set_cell(slab.get_cell())
        ratoms.extend(tmpf)
        ratoms.extend(tmpm)
        if(len(ratoms) == 0): return None
        if(self.minfrac > len(tmpf)/len(ratoms)): return None
        if(self.minfrac > len(tmpm)/len(ratoms)): return None

        ratoms.translate(cm)

        tmpm.set_cell(slab.get_cell())
        tmpf.set_cell(slab.get_cell())
        tmpm.translate(cm)
        tmpf.translate(cm)
        return ratoms

    def __get_minfrac(self,minfrac):
        if minfrac is not None:
            if(isinstance(minfrac,float)):
                return minfrac
            else:
                raise ValueError("Specified minfrac value not a float")
        else:
            return None