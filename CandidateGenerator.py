from typing import List
import numpy as np
from ase import Atom, Atoms
from ase.data import atomic_numbers
from ase.build import molecule
from ase.ga.utilities import (closest_distances_generator, get_all_atom_types)
from ase.ga.startgenerator import StartGenerator

class CandidateGenerator:
    """
    Class instantiated when passing:

    slab: 
        Atoms Object which will be used as a fixed cell size template on top of which the atoms are generated

    atomsList: 
        Array containing stoichiometry of atoms to be generated accepts:
            atoms_array = ['Pt6O5','Pt6O6','Pt6O7']
            atoms_array = [[78]*6+[8]*5,[[78]*6+[8]*6, [78]*6+[8]*7]
            atoms_array = [Atoms('CO', positions=[[0, 0, 0], [1.4, 0, 0]])]

    ratio_of_covalent_radii:
        Ratio of covalent radii passed to generate the minimum distance dictionary

    sort_atoms_by_quantity:
        Whether to sort the stored atoms list by quantity of atoms or not. Useful for accesing each stoichiometry by system size.
        Leave False if the inputed order is deliberate
    random_generation_box_size:
        The volume available for the random positioning of atoms, percentage of the unoccupied space in the slab cell atoms object
        - 1.0 to use the entirety of the unoccupied cell in the slab atoms object
        - default 0.8
    """

    def __init__(self,slab,atoms_array,ratio_of_covalent_radii=0.7,
                    sort_atoms_by_quantity = False, random_generation_box_size = 0.8):

        self.slab = slab
        self.atoms_array = self.__get_atoms_array(atoms_array,sort_atoms_by_quantity)
        self.ratio_of_covalent_radii = ratio_of_covalent_radii
        self.p0,self.v1,self.v2,self.v3 = self.__get_cell_params(slab,random_generation_box_size)

    def get_candidate_by_position(self,position,maxiter=None) -> Atoms:
        
        unique_atom_types = get_all_atom_types(self.slab, self.atoms_array[position].numbers)

        blmin = closest_distances_generator(atom_numbers=unique_atom_types,
                                    ratio_of_covalent_radii=self.ratio_of_covalent_radii)

        sg = StartGenerator(self.slab, self.atoms_array[position].numbers, blmin,
                    box_to_place_in=[self.p0, [self.v1, self.v2, self.v3]])
        
        return sg.get_new_candidate()

    def get_candidate_by_atoms(self,atoms) -> Atoms:
        "Returns a random structure matching the stoichiometry of the inputed atoms object"
        matches = [self.atoms_array.index(x) for x in self.atoms_array 
                    if self.__get_atoms_object(x).get_atomic_numbers().sort()  
                    == self.__get_atoms_object(atoms).get_atomic_numbers().sort()]
        if(len(matches)==0):
            raise Exception("No match for the inputed atoms object")
        return self.get_candidate_by_position(matches[0])

    def get_random_candidate(self,maxiter=None) -> Atoms:
        "Returns a random structure from all the possible stoichiometries"
        return  self.get_candidate(position = np.random(low = 0,high = len(self.atoms_array)),maxiter = maxiter)
        
    "Private Methods do not touch"

    def __get_cell_params(self,slab,random_generation_box_size):
        "Gets cell parameters from inputed slab"
        if(random_generation_box_size < 0.0): raise ValueError("random_generation_box_size negative value")
        if(random_generation_box_size > 1.0): raise ValueError("random_generation_box_size too big")
        pos = slab.get_positions()
        cell = slab.get_cell()
        if(len(pos) == 0):
            v1 = cell[0, :] * random_generation_box_size
            v2 = cell[1, :] * random_generation_box_size
            v3 = cell[2, :]
            v3[2] = 3.
            p0 = np.array([0,0,0])

        else:
            p0 = np.array([0., 0., max(pos[:, 2]) + 2.])
            v1 = cell[0, :] * random_generation_box_size
            v2 = cell[1, :] * random_generation_box_size
            v3 = cell[2, :]
            v3[2] = 3.
        return p0,v1,v2,v3

    def __get_atoms_array(self,atoms_array,sort_atoms_by_quantity) -> List[Atoms]:
        "Gets an array of atoms object based on user input"
        atoms_result_array = []
        if len(atoms_array) == 0:
            raise Exception("Empty atoms_array")
        if isinstance(atoms_array[0],Atom):
            raise Exception("Inputed atoms_array list not deep enough")
        for i in atoms_array:
            atoms_result_array.append(self.__get_atoms_object(i))
        if(sort_atoms_by_quantity):
            atoms_result_array.sort(key = lambda x: len(x))
        return atoms_result_array

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
    
    
        


    
    