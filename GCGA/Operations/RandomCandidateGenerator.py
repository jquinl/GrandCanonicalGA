from typing import List
import numpy as np
from ase import Atoms
from ase.data import atomic_numbers
from ase.build import molecule
from ase.ga.utilities import (closest_distances_generator, get_all_atom_types)
from ase.ga.startgenerator import StartGenerator

from .OperationsBase import OperationsBase
class RandomCandidateGenerator(OperationsBase):
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

    def __init__(self,slab,constant,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                    random_generation_box_size = 0.8,rng=np.random):
        super().__init__(slab,constant,variable_types,variable_range,ratio_of_covalent_radii,rng)
        self.p0,self.v1,self.v2,self.v3 = self.__get_cell_params(slab,random_generation_box_size)

    def get_candidate_by_number(self,number,maxiter=None) -> Atoms:
        if(number > len(self.combination_matrix)):
            raise Exception("Provided number higher than possible combinations")

        atoms = self.constant.copy()
        var = Atoms()
        new_atoms = self.combination_matrix[number]

        for i in range(len(new_atoms)):
            for j in range(new_atoms[i]):
                print(new_atoms[i])
                var.extend(self.variable_types[i])
        
        atoms.extend(var)
        unique_atom_types = get_all_atom_types(self.slab, atoms.numbers)
        blmin = closest_distances_generator(atom_numbers=unique_atom_types,
                                    ratio_of_covalent_radii=self.ratio_of_covalent_radii)
        
        atoms_numbers  = atoms.numbers
        sg = StartGenerator(self.slab, atoms_numbers, blmin,
                    box_to_place_in=[self.p0, [self.v1, self.v2, self.v3]])
        return_atoms = sg.get_new_candidate()
        return_atoms.info['stc']= self.get_var_id(return_atoms)
        return return_atoms

    def get_random_candidate(self,maxiter=None) -> Atoms:
        "Returns a random structure from all the possible stoichiometries"
        return  self.get_candidate_by_number(number = np.randnint(len(self.combination_matrix),maxiter=maxiter))

    def get_starting_population(self,population_size=20,maxiter=None):
        starting_population = []
        single_population_size = int(population_size/len(self.combination_matrix))
        for i in range(len(self.combination_matrix)):
            for j in range(single_population_size):
                atoms = self.get_candidate_by_number(i,maxiter=maxiter)
                #atoms.info['stc'] = i
                starting_population.append(atoms)
        return starting_population

    def get_starting_population_single(self,variable_number,population_size=20,maxiter=None,):
        starting_population = []
        if(variable_number < len(self.combination_matrix)):
            for j in range(population_size):
                atoms = self.get_candidate_by_number(variable_number,maxiter=maxiter)
                starting_population.append(atoms)
            return starting_population
        else:
            raise Exception("Provided variable number not in range")
    #"Private Methods do not touch"
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
    
    
        


    
    