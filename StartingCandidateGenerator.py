from typing import List
import numpy as np
from ase import Atom, Atoms
from ase.data import atomic_numbers
from ase.build import molecule
from ase.ga.utilities import (closest_distances_generator, get_all_atom_types)
from ase.ga.startgenerator import StartGenerator

class StartingCandidateGenerator:
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

    def __init__(self,slab,constant,variable,variable_range,ratio_of_covalent_radii=0.7,
                    random_generation_box_size = 0.8):

        self.slab = slab
        self.constant = self.__get_atoms_object(constant)
        self.variable = self.__get_atoms_object(variable)
        self.variable_range = self.__get_range(variable_range)
        self.ratio_of_covalent_radii = ratio_of_covalent_radii
        self.p0,self.v1,self.v2,self.v3 = self.__get_cell_params(slab,random_generation_box_size)

    def get_candidate_by_number(self,number,maxiter=None) -> Atoms:
        if(number not in self.variable_range):
            raise Exception("Provided amount of variable systems not present in range")
        atoms = self.constant.copy()
        var = Atoms()
        for j in range(number):
            var.extend(self.variable)
        atoms.extend(var)
        
        unique_atom_types = get_all_atom_types(self.slab, atoms.numbers)

        blmin = closest_distances_generator(atom_numbers=unique_atom_types,
                                    ratio_of_covalent_radii=self.ratio_of_covalent_radii)

        atoms_numbers  = self.constant.numbers

        for i in range(number):
            atoms_numbers = np.concatenate((atoms_numbers,self.variable.numbers),axis=None)
        
        sg = StartGenerator(self.slab, atoms_numbers, blmin,
                    box_to_place_in=[self.p0, [self.v1, self.v2, self.v3]])
        atoms = sg.get_new_candidate()
        atoms.info['stc']= self.get_var_stc(atoms)

        return atoms
    def get_var_stc(self,atoms) -> int:
        var_stc = len(atoms)-len(self.slab)-len(self.constant)
        if(var_stc >= 0 ):
            for i in atoms[(len(self.slab)+len(self.constant)):len(atoms)]:
                if(i.symbol != self.variable[0].symbol):
                    raise Exception("Variable type of atoms does not match stored type")
        else:
            raise Exception("Negative numer of variable atoms")
        return var_stc
    def get_random_candidate(self,maxiter=None) -> Atoms:
        "Returns a random structure from all the possible stoichiometries"
        return  self.get_candidate_by_number(number = np.random.choice(self.variable_range,size=1)[0],maxiter=maxiter)

    def get_starting_population(self,population_size=20,maxiter=None):
        starting_population = []
        starting_population_numbers = []
        single_population_size = int(population_size/len(self.variable_range))
        for i in self.variable_range:
            for j in range(single_population_size):
                atoms = self.get_candidate_by_number(i,maxiter=maxiter)
                #atoms.info['stc'] = i
                starting_population.append(atoms)
        return starting_population
    def get_starting_population_single(self,variable_number,population_size=20,maxiter=None,):
        starting_population = []
        starting_population_numbers = []
        if(variable_number in self.variable_range):
            for j in range(population_size):
                atoms = self.get_candidate_by_number(variable_number,maxiter=maxiter)
                starting_population.append(atoms)
            return starting_population
        else:
            raise Exception("Provided variable number not in range")
    #"Private Methods do not touch"
    def __get_range(self,variable_range) -> List[int]:
        if isinstance(variable_range,List) and isinstance(variable_range[0],int):
            return variable_range
        else:
            raise Exception("variable_ range is not al ist of integers")

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
    
    
        


    
    