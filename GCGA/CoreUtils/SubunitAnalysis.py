from ase.build import connected_indices as con_ind
from itertools import combinations

from ase import Atoms

import numpy as np

class SubunitFinder():

    def __main__(self):
        pass

    @classmethod
    def find_subunits(self,subunit,atoms,dmax=None,scale=1.5):
       # if(not issubclass(submolecule, Atoms)): raise TypeError("Passed argument to __submolecule is not an atoms objects")
        
        symbols = [a.symbol for a in subunit]

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
        print(clusters)
        atypes = ["H","He","Ni","Pt","Au","Xe","C"]
        counter = 0
        for i in clusters:
            print(i)
            for j in self.__subdivide_clusters(i,atoms,subunit,dmax=dmax,scale=scale):
                print(j)
                for k in j:
                    atoms[k].symbol = atypes[counter]
                counter+=1
        print(atoms)
        return atoms.copy()
    @classmethod
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

    @classmethod
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
        subunit_bonds = self.__bond_type_length(range(len(subunit)),subunit)

        subunit_bonds.sort(key=lambda x:  x[2])

        #Checks wether any of the combinations share the same atom
        accepted_combinations = []
        for k in full_combinations:
            test_bonds = []
            test_bonds = self.__bond_type_length(k,atoms)
            print(test_bonds)
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
           
            add_to_accepted = True
            #Check that the accepted combinations are unique (dont share atoms) if they do they prioritize the one with the shortest sum of bond lengths
            for i in range(len(accepted_combinations)):
              
                if any(x in accepted_combinations[i] for x in k):
                    sum_new = self.__total_bond_lenght(k,atoms)
                    sum_old = self.__total_bond_lenght(accepted_combinations[i],atoms)

                    if(sum_new < sum_old):
                        accepted_combinations[i] = k
                        add_to_accepted = False
                        continue
                    else:
                        add_to_accepted = False
                        continue
            if(add_to_accepted):
                accepted_combinations.append(k)
        return accepted_combinations
    
    @classmethod
    def __bond_type_length(self,indices,atoms):
        pairs = []
        total_bond = []
        for i in indices:
            for j in indices:
                if(i != j and [j,i] not in pairs):
                    pairs.append([j,i])
                    pairs.append([i,j])
                    total_bond.append([atoms[i].symbol,atoms[j].symbol,atoms.get_distance(i,j)])
        return total_bond
    
    @classmethod
    def __total_bond_lenght(self,indices,atoms):
        pairs = []
        total_bond_length = 0.0
        for i in indices:
            for j in indices:
                if(i != j and [j,i] not in pairs):
                    pairs.append([j,i])
                    pairs.append([i,j])
                    total_bond_length += atoms.get_distance(i,j)
        return total_bond_length