from enum import unique
from itertools import combinations
from operator import index

from ase import Atoms
from ase.build import connected_indices as con_ind

from abc import ABC,abstractmethod
from matplotlib.pyplot import sca

import numpy as np
class BaseFitness(ABC):

    def __init__(self,environmental_variables,references):
        self.env = self.__get_env(environmental_variables)
        self.references = self.__get_refs(**references)
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
        for a,k in references.items():
            if(not (type(a) == str and type(k) == float)):
                raise TypeError("references is not a string:float dictionary")
            try:
                at = Atoms(a)
            except:
                raise ValueError("Cannot create atoms with key: " + a)

        return references

    def __get_env(self,env):
        if(type(env) is not dict):raise TypeError("environment is not a string:float dictionary")
        for e,k in env.items():
                if(not (type(e) == str and type(k) == float)):
                    raise TypeError("environment is not a string:float dictionary")
        return env

    def submolecules_count(self,submolecule,atoms,dmax=None,scale=1.5):
       # if(not issubclass(submolecule, Atoms)): raise TypeError("Passed argument to __submolecule is not an atoms objects")
        
        symbols = [a.symbol for a in submolecule]

        count = 0
        symbols.sort()
        comb_dict = {}
        for a in atoms:
            if(a.symbol in symbols):
                new = [i for i in con_ind(atoms,a.index,dmax=dmax,scale=scale) if i != a.index]
                comb_dict[a.index] = new
        clusters = []
        for a in comb_dict.keys():
            atree = self.__traverse_tree(atoms,comb_dict,a,symbols)
            atree.sort()
            clusters.append(atree)
        
        clusters = clusters

        
        clusters = list(set(frozenset(item) for item in clusters))
        atypes = ["H","He","Ni","Pt","Au","Xe"]
        counter = 0
        for i in clusters:
            
            for j in self.__subdivide_clusters(i,atoms,submolecule,dmax=dmax,scale=scale):
                for k in j:
                    atoms[k].symbol = atypes[counter]
            counter+=1
        return count

    def __subdivide_clusters(self,cluster,atoms,subunit,dmax,scale):
 
        #Gets  unique combinations of the atoms
        nums = subunit.numbers
        nums.sort()
        comb = []
        for i in  combinations(cluster,len(subunit)):
            comb.append(list(i))
      
        #Gets combinaitons that contain only the atoms in the subunit
        remove_from_list = []
        for i in range(len(comb)):
            at = Atoms()
            for j in comb[i]:
                at.extend(atoms[j])
            newnums = at.numbers
            newnums.sort()
            for a,b in zip(newnums,nums):
                if(a != b):
                    remove_from_list.append(i)

        
        full_combinations = [comb[i] for i in range(len(comb)) if i not in remove_from_list]

        #Creates a sorted matrix containing unique atom pairings and bond distances in the subunit
        subunit_bonds = []
        pairs = []
        for i in range(len(subunit)):
            for a in con_ind(subunit,i,dmax=dmax,scale=scale):
                if([i,a] not in pairs):
                    pairs.append([i,a])
                    pairs.append([a,i])
                    subunit_bonds.append([subunit[i].symbol,subunit[a].symbol,subunit.get_distance(i,a)])

        subunit_bonds.sort(key=lambda x:  x[2])

        #Checks wether any of the combinations share the same atom
        accepted_combinations = []
        for k in full_combinations:
            pairs = []
            test_bonds = []
            test_bonds = self.__bond_type_array(k,atoms,dmax=dmax,scale=scale)
            
            test_bonds.sort(key=lambda x: x[2])
            if len(test_bonds) != len(subunit_bonds):
                continue
            for i in range(len(test_bonds)):
                sub = [subunit_bonds[i][0],subunit_bonds[i][1]]
                tst = [test_bonds[i][0],test_bonds[i][1]]
                sub.sort()
                tst.sort()
                if(tst != sub): continue
            
            for i in range(len(test_bonds)):
                if(subunit_bonds[i][2] * scale < test_bonds[i][2]):
                    continue
           
            #Check that the accepted combinations are unique (dont share atoms) if they do they prioritize the one with the shortest sum of bond lengths
            for i in range(len(accepted_combinations)):
                print(accepted_combinations[i]) 
                print(k)
                if any(x in accepted_combinations[i] for x in k):
                    compare = self.__bond_type_array(accepted_combinations[i],atoms,dmax=dmax,scale=scale)
                    compare.sort(key=lambda x : x[2])
                    sum_new = 0.0
                    sum_old = 0.0

                    for z,al in zip(test_bonds,compare):
                        sum_new += z[2]
                        sum_old += al[2]
                   
                    if(sum_new < sum_old):
                        accepted_combinations[i] = k
                  
                        continue
                    else:

                        continue
            accepted_combinations.append(k)

        return accepted_combinations
    
    def __bond_type_array(self,comb,atoms,dmax,scale):
        pairs = []
        bond_array = []
        for i in comb:
            con = [a for a in con_ind(atoms,i,dmax=dmax,scale=scale) if a != i]
            for a in con:
                if([i,a] not in pairs):
                    pairs.append([i,a])
                    pairs.append([a,i])
                    bond_array.append([atoms[i].symbol,atoms[a].symbol,atoms.get_distance(i,a)])
        return bond_array
    #Bond Type array but return only the disntance between the given indices
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
            node_queue.pop(0)

        return visited_nodes



class GibbsFreeEnergy(BaseFitness):
    def __init__(self, environmental_variables, references):
        super().__init__(environmental_variables, references)

    def evaluate(self, atoms) -> float:
        super().evaluate(atoms)
    
        return self.references
