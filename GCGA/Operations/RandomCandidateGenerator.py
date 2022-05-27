from configparser import RawConfigParser
import random
import numpy as np
from ase import Atoms
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

    def __init__(self,slab,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                    random_generation_box_size = 0.8,rng=np.random,test_too_far=True):
        super().__init__(slab,variable_types,variable_range,ratio_of_covalent_radii,rng)
        self.p0,self.v1,self.v2,self.v3 = self.__get_cell_params(slab,random_generation_box_size)
        self.test_too_far = test_too_far
    def get_candidate_by_number(self,number,maxiter=None) -> Atoms:
        if(number > len(self.combination_matrix)):
            raise Exception("Provided number higher than possible combinations")

        atoms = Atoms()
        new_atoms = self.combination_matrix[number]
        for i in range(len(new_atoms)):
            for j in range(new_atoms[i]):
                atoms.extend(self.variable_types[i])

        #Consider using self.blmin here see if it works

        if(len(atoms) == 0):
            raise Exception("Empty atoms being generated at random, revise stoichiometies")
        if(len(atoms) == 1):
            return_atoms = atoms.copy()
            return_atoms.set_cell(self.slab.get_cell())
            var_id = self.get_var_id(return_atoms)
            if(var_id is not None):
                return_atoms.info['stc']= var_id
                return return_atoms
            else:
                raise Exception("Provided atomic combination is not present in combination matrix")
        unique_atom_types = get_all_atom_types(self.slab, atoms.numbers)
        blmin = closest_distances_generator(atom_numbers=unique_atom_types,
                                    ratio_of_covalent_radii=self.ratio_of_covalent_radii)
        
        atoms_numbers  = atoms.numbers
        sg = StartGenerator(self.slab, atoms_numbers, blmin,
                    box_to_place_in=[self.p0, [self.v1, self.v2, self.v3]],test_too_far=self.test_too_far)
        return_atoms = sg.get_new_candidate()
        var_id = self.get_var_id(return_atoms)
        if(var_id is not None):
            return_atoms.info['stc']= var_id
            return return_atoms
        else:
            raise Exception("Provided atomic combination is not present in combination matrix")

    def get_random_candidate(self,maxiter=None) -> Atoms:
        "Returns a random structure from all the possible stoichiometries"
        return  self.get_candidate_by_number(number = np.random.randint(len(self.combination_matrix)),maxiter=maxiter)
    

    "The mutate override to be used in runtime is placed here, it calls the get_random_candidate method"
    def mutate(self, a1, a2):
        super().mutate(a1,a2)
        return self.get_random_candidate(), 0

    def get_starting_population(self,population_size=20,maxiter=None):
        starting_population = []
        single_population_size = max(1,int(population_size/len(self.combination_matrix)))
        for i in range(len(self.combination_matrix)):
            for j in range(single_population_size):
                atoms = self.get_candidate_by_number(i,maxiter=maxiter)
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
            v3 = cell[2, :] * random_generation_box_size
            p0 = np.array([0,0,0])
        else:
            p0 = np.array([0., 0., max(pos[:, 2]) + 2.])
            v1 = cell[0, :] * random_generation_box_size
            v2 = cell[1, :] * random_generation_box_size
            v3 = cell[2, :] * random_generation_box_size

        return p0,v1,v2,v3
    ############Work In Progress##################
    def __generate(self,slab, atom_numbers, blmin):
        cand = slab.copy()
        
        random.shuffle(atom_numbers)
        for i in atom_numbers:
            maxtries = 100
            tries = 0
            done = False
            while tries< maxtries and not done:
                self.__check_overlap()
                done = True

        return
    
    def __check_overlap(self,atom1,atom2,dist):
        return dist*dist < ((atom1.position[0]-atom2.postion[0]) * (atom1.position[0]-atom2.postion[0]) +
                            (atom1.position[1]-atom2.postion[1]) * (atom1.position[1]-atom2.postion[1]) +
                            (atom1.position[2]-atom2.postion[2]) * (atom1.position[2]-atom2.postion[2]))
    def __random_position_in_box(self):
        
        return  
    def __random_position_from_atom(self,atom,distance):
        
        new_vec = self.rng.normal(size=3)
        while(new_vec == np.array([0.0,0.0,0.0])):
            new_vec = self.rng.normal(size=3)

        norm = np.linalg.norm(new_vec)
        norm *= self.rng.uniform(distance,distance*2.0)

        pos = atom.position + norm

        pos[0] = max(min(self.p0[0]+self.v1[0]+self.v2[0]+self.v3[0], pos[0]),self.p0[0])
        pos[1] = max(min(self.p0[1]+self.v1[1]+self.v2[1]+self.v3[1], pos[1]),self.p0[1])
        pos[2] = max(min(self.p0[2]+self.v1[2]+self.v2[2]+self.v3[2], pos[2]),self.p0[2])

        return pos



    
    