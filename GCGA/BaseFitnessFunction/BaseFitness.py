from enum import unique
from turtle import distance
from ase import Atoms
from ase.build import connected_indices as con_ind

from abc import ABC,abstractmethod

import numpy as np
class BaseFitness(ABC):

    def __init__(self,environmental_variables,references):
        self.env = self.__get_env(environmental_variables)
        self.references = self.__get_refs(references)
        pass

    @classmethod
    def fitness_class(self):
        return True
    
    def fitness_instance(self):
        return True

    @abstractmethod
    def evaluate(self,atoms)-> float:
        if(not issubclass(atoms,Atoms)):
            raise TypeError("Fitness function tried to evaluate something other than an ASE atoms object")
        
        pass
    
    def __get_refs(self,**references):
        if(type(references) is not dict):raise TypeError("references is not a string:float dictionary")
        for a in references:
            if(not (type(a.value()) == float and type(a.key()) == float)):
                raise TypeError("references is not a string:float dictionary")
            try:
                at = Atoms(a.key())
            except:
                raise ValueError("Cannot create atoms with key: " + a.key())

        return references

    def __get_env(self,env):
        if(type(env) is not dict):raise TypeError("environment is not a string:float dictionary")
        for e in env:
                if(not (type(e.value()) == float and type(e.key()) == float)):
                    raise TypeError("environment is not a string:float dictionary")
        return env

    def __submolecules_count(self,str,atoms,dmax=None,scale=1.5):
        if(type(str) is not str): raise TypeError("Passed argument to __submolecule is not a string")
        try:
            at = Atoms(str)
        except:
            raise ValueError("Not a valid molecule: " + str)
        symbols = [a.symbol for a in at]

        count = 0
        symbols.sort()
        comb_dict = {}

      
        for a in atoms:
            if(a.symbol in symbols):
                comb_dict[a.index] = con_ind(atoms,a.index,dmax=dmax,scale=scale).sort()

        clusters = []
        for a in comb_dict.keys:
            clusters.append(self.__traverse_tree(atoms,comb_dict,a,symbols).sort())
        
        clusters = unique(clusters)

        for i in clusters:
            if len(i) == len(at):
                count += 1
            if(len(i)% len(at) == 0):
                self.__subdivide_clusters(i,atoms,at,dmax)
        return count

    def __subdivide_clusters(self,cluster,atoms,subunit,dmax):
        from itertools import combinations
        distance_matrix = np.array(len(cluster),len(cluster))
        for i in cluster:
            for j in cluster:
                distance_matrix[i,j] = abs(np.linalg.norm(atoms[i].position -atoms[j].position))
                if(i==j): distance_matrix[i,j] = 1000.0

        nums = subunit.numbers.sort()
        comb = combinations(cluster,len(subunit))
        remove_from_list = []
        for i in range(len(comb)):
            at = Atoms()
            for j in comb[i]:
                at.extend(atoms[i])
            if(not at.numbers.sort() == nums):
                remove_from_list.append(i)
        
        full_combination = [comb[i] for i in range(len(comb)) if i!= remove_from_list]
        




    def __traverse_tree(self,atoms,dictionary,start_index,symbols):
        if(not atoms[start_index].symbol in symbols): raise ValueError("Initial index not in symbols")

        visited_nodes = []
        node_queue = []
        current_node = dictionary[start_index]
        node_queue.append(current_node)
        visited_nodes.append(start_index)

        
        while(len(node_queue)>0):
            for i in node_queue[0]:
                if(atoms[i].symbol in symbols and i not in visited_nodes):
                    node_queue.append(dictionary[i])
                    visited_nodes.append(i)
            node_queue.remove(0)

        return visited_nodes



class GibbsFreeEnergy(BaseFitness):
    def __init__(self, environmental_variables, references):
        super().__init__(environmental_variables, references)

    def evaluate(self, atoms) -> float:
        super().evaluate(atoms)
        
        atoms_energy = atoms.get_potential_energy()

        at_num= atoms.get_atomic_numbers()


        self.references[]
        pt_num = np.count_nonzero(at_num == 78)
        au_num= np.count_nonzero(at_num == 79)
        return self.references



def fitness_function(atoms)-> float:
    env = 1.0
    ref=read('pt.traj@:')[0].get_potential_energy()    

    au_en =read('gold_bulk.traj').get_potential_energy() / 500.0
    at_num= atoms.get_atomic_numbers()
    pt_num = np.count_nonzero(at_num == 78)
    au_num= np.count_nonzero(at_num == 79)

    fre = a - ref- au_num * (au_en) - au_num*env

    return -fre
