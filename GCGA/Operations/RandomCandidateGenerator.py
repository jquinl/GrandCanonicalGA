import copy
import random
import numpy as np
from ase import Atoms



from .MutationsBase import MutationsBase
class RandomCandidateGenerator(MutationsBase):
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
                    random_generation_box_size = 0.8,rng=np.random,max_bond_lenght_multi=4.0):
        super().__init__(slab,variable_types,variable_range,ratio_of_covalent_radii,rng,mut_box_size=random_generation_box_size)
        self.max_blen = max_bond_lenght_multi

    @classmethod
    def rand_generator(cls):
        pass


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
        
        atoms_numbers  = atoms.numbers
       
        child = self.__generate(self.slab,atom_numbers=atoms_numbers,blmin=self.blmin)
        return_atoms = self.slab.copy()
        return_atoms.extend(self.sort_atoms_by_type(child[len(self.slab):]))
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
    def mutate(self, a1):
        super().mutate(a1)
        return self.get_random_candidate()

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


    ############Work In Progress##################
    def __generate(self,slab, atom_numbers, blmin):
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
                    at.position = self.__random_position_in_box()
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
                    at.position = self.__random_position_from_atom(cand[candidate],blmin[(cand[candidate].number,at.number)])
                
                if(not self._overlaps(cand,newAtoms,blmin)):
                    cand.extend(newAtoms)
                    done = True
                tries +=1
            if(not done):
                return None
        final_atoms = slab.copy()
        final_atoms.extend(self.sort_atoms_by_type(cand[len(slab):]))
        return final_atoms


    def __random_position_in_box(self):
        return self.p0 + (self.rng.random() * self.v1 + self.rng.random() * self.v2 + self.rng.random() * self.v3)

    def __random_position_from_atom(self,atom,distance):
        
        new_vec = self.rng.normal(size=3)
        while(new_vec[0] == 0.0 and new_vec[1] == 0.0 and new_vec[2] == 0.0):
            new_vec = self.rng.normal(size=3)

        norm = self._normalize(new_vec)

        unif = self.rng.uniform(distance,distance*self.max_blen)
        norm *= unif
        pos = atom.position + norm
        pos[0] = max(min(self.p0[0]+self.v1[0]+self.v2[0]+self.v3[0], pos[0]),self.p0[0])
        pos[1] = max(min(self.p0[1]+self.v1[1]+self.v2[1]+self.v3[1], pos[1]),self.p0[1])
        pos[2] = max(min(self.p0[2]+self.v1[2]+self.v2[2]+self.v3[2], pos[2]),self.p0[2])
        
        return pos



    
    