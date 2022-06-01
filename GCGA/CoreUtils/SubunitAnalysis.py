from ase.build import connected_indices as con_ind
from itertools import combinations

from ase import Atoms

import numpy as np

class SubunitFinder():

    def __main__(self):
        pass

    @classmethod
    def find_subunits(cls,subunit,atoms,dmax=None,scale=1.5):
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
            atree = cls.__traverse_tree(atoms,comb_dict,a,symbols)
            atree.sort()
            clusters.append(atree)
        
        clusters = clusters

        
        clusters = list(set(frozenset(item) for item in clusters))
        atypes = ["H","He","Ni","Pt","Au","Xe","C"]
        counter = 0
        for i in clusters:
            for j in cls.__subdivide_clusters(i,atoms,subunit,dmax=dmax,scale=scale):
                for k in j:
                    atoms[k].symbol = atypes[counter]
                counter+=1
        return atoms.copy()
    @classmethod
    def __traverse_tree(cls,atoms,dictionary,start_index,symbols):
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
    def __subdivide_clusters(cls,cluster,atoms,subunit,dmax,scale):
 
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
        subunit_bonds = cls.__bond_type_length(range(len(subunit)),subunit)

        subunit_bonds.sort(key=lambda x:  x[2])

        #Checks wether any of the combinations share the same atom
        accepted_combinations = []
        for k in full_combinations:
            test_bonds = []
            test_bonds = cls.__bond_type_length(k,atoms)
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
                    sum_new = cls.__total_bond_lenght(k,atoms)
                    sum_old = cls.__total_bond_lenght(accepted_combinations[i],atoms)

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
    def __bond_type_length(cls,indices,atoms):
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
    def __total_bond_lenght(cls,indices,atoms):
        pairs = []
        total_bond_length = 0.0
        for i in indices:
            for j in indices:
                if(i != j and [j,i] not in pairs):
                    pairs.append([j,i])
                    pairs.append([i,j])
                    total_bond_length += atoms.get_distance(i,j)
        return total_bond_length
class NonEnergyInteratomicDistanceComparator:

    """ An implementation of the comparison criteria described in
          L.B. Vilhelmsen and B. Hammer, PRL, 108, 126101 (2012)

        Parameters:

        n_top: The number of atoms being optimized by the GA.
            Default 0 - meaning all atoms.

        pair_cor_cum_diff: The limit in eq. 2 of the letter.
        pair_cor_max: The limit in eq. 3 of the letter
        dE: The limit of eq. 1 of the letter
        mic: Determines if distances are calculated
        using the minimum image convention

        Modified to not check for the energy of structures
    """
    def __init__(self, n_top=None, pair_cor_cum_diff=0.015,
                 pair_cor_max=0.7, mic=False):
        self.pair_cor_cum_diff = pair_cor_cum_diff
        self.pair_cor_max = pair_cor_max
        self.n_top = n_top or 0
        self.mic = mic

    def looks_like(self, a1, a2):
        """ Return if structure a1 or a2 are similar or not. """
        if len(a1) != len(a2):
            raise Exception('The two configurations are not the same size')

        # then we check the structure
        a1top = a1[-self.n_top:]
        a2top = a2[-self.n_top:]
        cum_diff, max_diff = self.__compare_structure__(a1top, a2top)

        return (cum_diff < self.pair_cor_cum_diff
                and max_diff < self.pair_cor_max)

    def __compare_structure__(self, a1, a2):
        """ Private method for calculating the structural difference. """
        p1 = get_sorted_dist_list(a1, mic=self.mic)
        p2 = get_sorted_dist_list(a2, mic=self.mic)
        numbers = a1.numbers
        total_cum_diff = 0.
        max_diff = 0
        for n in p1.keys():
            cum_diff = 0.
            c1 = p1[n]
            c2 = p2[n]

            if( len(c1) != len(c2)): return (self.pair_cor_cum_diff + 1.0,self.pair_cor_max +1)
            if len(c1) == 0:
                continue
            t_size = np.sum(c1)
            d = np.abs(c1 - c2)
            cum_diff = np.sum(d)
            max_diff = np.max(d)
            ntype = float(sum([i == n for i in numbers]))
            total_cum_diff += cum_diff / t_size * ntype / float(len(numbers))
        return (total_cum_diff, max_diff)


def get_sorted_dist_list(atoms, mic=False):
    """ Utility method used to calculate the sorted distance list
        describing the cluster in atoms. """
    numbers = atoms.numbers
    unique_types = set(numbers)
    pair_cor = dict()
    for n in unique_types:
        i_un = [i for i in range(len(atoms)) if atoms[i].number == n]
        d = []
        for i, n1 in enumerate(i_un):
            for n2 in i_un[i + 1:]:
                d.append(atoms.get_distance(n1, n2, mic))
        d.sort()
        pair_cor[n] = np.array(d)
    return pair_cor

     
