from asyncio import constants
from typing import List
import numpy as np
from ase import Atom, Atoms
from ase.data import atomic_numbers
from ase.ga.cutandsplicepairing import Positions
from ase.geometry import find_mic
from ase.ga.utilities import (closest_distances_generator, get_all_atom_types,atoms_too_close, atoms_too_close_two_sets)
class CrossOperation:
    """
    Modified cross operation found in the Atomic Simulation Environment (ASE) GA package. ga.cutandspliceparing.py
    Modified in order to allow the cut and splice pairing to happen between neighboring stoichiometries
    """
    def __init__(self, slab,constant,variable,variable_range,ratio_of_covalent_radii=0.7,
                minfrac = None,rng=np.random):
        self.slab = slab
        self.constant = self.__get_atoms_object(constant)
        self.variable = self.__get_atoms_object(variable)
        self.variable_range = self.__get_range(variable_range)
        self.ratio_of_covalent_radii = ratio_of_covalent_radii
        self.rng = rng
        self.minfrac = minfrac

        uniques = self.constant.copy()
        uniques.extend(self.variable)
            
        unique_atom_types = get_all_atom_types(self.slab, uniques.numbers)

        self.blmin = closest_distances_generator(atom_numbers=unique_atom_types,
                                    ratio_of_covalent_radii=self.ratio_of_covalent_radii)


    def cross(self, a1, a2):
        """Crosses the two atoms objects and returns one"""

        if (len(a1)-len(self.slab)-len(self.constant) not in self.variable_range):
            raise ValueError('Wrong size of structure a1 to optimize')
        if (len(a2)-len(self.slab)-len(self.constant) not in self.variable_range):
            raise ValueError('Wrong size of structure a2 to optimize')

       
        # Only consider the atoms to optimize
        a1 = a1[len(self.slab) :len(a1)]
        a2 = a2[len(self.slab) :len(a2)]

        invalid = True
        counter = 0
        maxcount = 1000
        a1_copy = a1.copy()
        a2_copy = a2.copy()

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
                cell = a_copy.get_cell()
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
            
            # Perform the pairing:
            #child = self._get_pairing(a1_copy, a2_copy, cut_p, cut_n, self.slab.get_cell())
            if child is None:
                continue

            # Verify whether the atoms are too close or not:
            

            #if len(self.slab) > 0:
            #    if atoms_too_close_two_sets(self.slab, child, self.blmin):
            #        continue

            atoms  = self.slab.copy()

            atoms.extend(child)

            if atoms_too_close(atoms, self.blmin):
                continue
            if(not self.mantains_ordering(atoms)):
                continue
            if(self.get_var_stc(atoms) not in self.variable_range):
                continue
            # Passed all the tests
            atoms.wrap()
            atoms.info['stc']= self.get_var_stc(atoms)
            return atoms

        return None
    
    def mantains_ordering(self,atoms):
        for i in range(len(self.slab)):
            if(atoms[i].symbol != self.slab[i].symbol):
                print("Eror in ordering")
                return False
        for i in range(len(self.constant)):
            if(atoms[len(self.slab)+i].symbol != self.constant[i].symbol):
                return False
        return True

        
                
    def get_var_stc(self,atoms) -> int:
        var_stc = len(atoms)-len(self.slab)-len(self.constant)
        if(var_stc >= 0 ):
            for i in atoms[(len(self.slab)+len(self.constant)):len(atoms)]:
                if(i.symbol != self.variable[0].symbol):
                    raise Exception("Variable type of atoms does not match stored type")
        else:
            raise Exception("Negative numer of variable atoms")
        return var_stc
    def __get_range(self,variable_range) -> List[int]:
        if isinstance(variable_range,List) and isinstance(variable_range[0],int):
            return variable_range
        else:
            raise Exception("variable_ range is not al ist of integers")


    def __get_atoms_object(self,atoms) -> Atoms:
        "Gets an atoms object based on user input"
        if isinstance(atoms, Atoms):
            return atoms
        elif isinstance(atoms, str):
            return Atoms(atoms)
        elif isinstance(atoms,List):
            for i in atoms:
                if(i not in atomic_numbers.values()):
                    raise ValueError('Cannot parse this element {} in :'.format(i),atoms )
            return Atoms(numbers=atoms)
        else:
            raise ValueError('Cannot parse this element:', atoms)

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
        cell = self.slab.get_cell()
        "Generate constant part first"
        for num in self.constant.numbers:
            has_been_added = False
            for atom in a1_copy:
                if(atom.symbol == num and not has_been_added):
                    at_vector =  np.linalg.norm(atom.position - cutting_point)
                    if(np.dot(at_vector,cutting_normal)<0):
                        atoms_result.append(atom)
                        atom.symbol = 200
                        has_been_added = True      
            for atom in a2_copy:
                if(atom.symbol == num and not has_been_added):
                    at_vector =  np.linalg.norm(atom.position - cutting_point)
                    if(np.dot(at_vector,cutting_normal)>0):
                        atoms_result.append(atom)
                        atom.symbol = 200
                        has_been_added = True
            for atom in a1_copy:
                if(atom.symbol == num and not has_been_added):
                    at_vector =  np.linalg.norm(atom.position - cutting_point)
                    if(np.dot(at_vector,cutting_normal)>0):
                        atoms_result.append(atom)
                        atom.symbol = 200
                        has_been_added = True
            for atom in a2_copy:
                if(atom.symbol == num and not has_been_added):
                    at_vector =  np.linalg.norm(atom.position - cutting_point)
                    if(np.dot(at_vector,cutting_normal)<0):
                        atoms_result.append(atom)
                        atom.symbol = 200
                        has_been_added = True
            if(not has_been_added):
                atoms_result.append(Atom(symbol=num,position = (self.rng.rand()* cell[0],self.rng.rand()*cell[1],self.rng.rand()*cell[2])))
        
        atoms_result.wrap()
        return atoms_result


    def _get_pairing(self, a1, a2, cutting_point, cutting_normal, cell):
        """Creates a child from two parents using the given cut.

        Returns None if the generated structure does not contain
        a large enough fraction of each parent (see self.minfrac).

        Does not check whether atoms are too close.

        Assumes the 'slab' parts have been removed from the parent
        structures and that these have been checked for equal
        lengths, stoichiometries, and tags (if self.use_tags).

        Parameters:

        cutting_normal: int or (1x3) array

        cutting_point: (1x3) array
            In fractional coordinates

        cell: (3x3) array
            The unit cell for the child structure
        """
        symbols = a1.get_chemical_symbols()
        tags = np.arange(len(a1))

        # Generate list of all atoms / atom groups:
        p1, p2, sym = [], [], []
        for i in np.unique(tags):
            indices = np.where(tags == i)[0]
            s = ''.join([symbols[j] for j in indices])
            sym.append(s)

            for i, (a, p) in enumerate(zip([a1, a2], [p1, p2])):
                print(i, a, p)
                c = a.get_cell()
                cop = np.mean(a.positions[indices], axis=0)
                cut_p = np.dot(cutting_point, c)
                if isinstance(cutting_normal, int):
                    vecs = [c[j] for j in range(3) if j != cutting_normal]
                    cut_n = np.cross(vecs[0], vecs[1])
                else:
                    cut_n = np.dot(cutting_normal, c)
                d = np.dot(cop - cut_p, cut_n)
                spos = a.get_scaled_positions()[indices]
                scop = np.mean(spos, axis=0)
                p.append(Positions(spos, scop, s, d, i))

        all_points = p1 + p2
        unique_sym = np.unique(sym)
        types = {s: sym.count(s) for s in unique_sym}

        # Sort these by chemical symbols:
        all_points.sort(key=lambda x: x.symbols, reverse=True)

        # For each atom type make the pairing
        unique_sym.sort()
        use_total = dict()
        for s in unique_sym:
            used = []
            not_used = []
            # The list is looked trough in
            # reverse order so atoms can be removed
            # from the list along the way.
            for i in reversed(range(len(all_points))):
                # If there are no more atoms of this type
                if all_points[i].symbols != s:
                    break
                # Check if the atom should be included
                if all_points[i].to_use():
                    used.append(all_points.pop(i))
                else:
                    not_used.append(all_points.pop(i))

            assert len(used) + len(not_used) == types[s] * 2

            # While we have too few of the given atom type
            while len(used) < types[s]:
                index = self.rng.randint(len(not_used))
                used.append(not_used.pop(index))

            # While we have too many of the given atom type
            while len(used) > types[s]:
                # remove randomly:
                index = self.rng.randint(len(used))
                not_used.append(used.pop(index))

            use_total[s] = used

        n_tot = sum([len(ll) for ll in use_total.values()])
        assert n_tot == len(sym)

        # check if the generated structure contains
        # atoms from both parents:
        count1, count2, N = 0, 0, len(a1)
        for x in use_total.values():
            count1 += sum([y.origin == 0 for y in x])
            count2 += sum([y.origin == 1 for y in x])

        nmin = 1 if self.minfrac is None else int(round(self.minfrac * N))
        if count1 < nmin or count2 < nmin:
            return None

        # Construct the cartesian positions and reorder the atoms
        # to follow the original order
        newpos = []
        pbc = a1.get_pbc()
        for s in sym:
            p = use_total[s].pop()
            c = a1.get_cell() if p.origin == 0 else a2.get_cell()
            pos = np.dot(p.scaled_positions, c)
            cop = np.dot(p.cop, c)
            vectors, lengths = find_mic(pos - cop, c, pbc)
            newcop = np.dot(p.cop, cell)
            pos = newcop + vectors
            for row in pos:
                newpos.append(row)

        newpos = np.reshape(newpos, (N, 3))
        num = a1.get_atomic_numbers()
        child = Atoms(numbers=num, positions=newpos, pbc=pbc, cell=cell,
                      tags=tags)
        child.wrap()
        return child
