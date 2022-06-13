from ase.ga.standard_comparators import InteratomicDistanceComparator
from math import tanh, sqrt, exp
class Population:
    def __init__(self, population_size, database_interface):
        self.pop_size = population_size
        self.dbi = database_interface
        self.pop = []
        self.confid = 0

    
    def update_population(self,atoms):
        isSimilar = False
        
        if(len(self.pop) == 0):
            atoms.info['key_value_pairs']['confid'] =  self.confid
            self.confid += 1
            self.pop.append(atoms)
        else:
            for i in range(len(self.pop)):
                if self.get_similarity(atoms,self.pop[i]):
                    isSimilar = True
                    atoms.info['key_value_pairs']['confid'] = self.pop[i].info['key_value_pairs']['confid']
                    if(atoms.info['key_value_pairs']['raw_score'] > self.pop[i].info['key_value_pairs']['confid']):
                        self.pop[i] = atoms
        
            if(not isSimilar and len(self.pop) <= self.pop_size):

                atoms.info['key_value_pairs']['confid'] =  self.confid
                self.pop.append(atoms)
            elif(not isSimilar and len(self.pop)>self.pop_size):

                confids = [a.info['key_value_pairs']['confid'] for a in self.pop ]
                atoms.info['key_value_pairs']['confid'] =  self.__check_other_confids(atoms,confids)
                if(self.pop[-1].info['key_value_pairs']['raw_score'] < atoms.info['key_value_pairs']['raw_score']):
                    self.pop[-1] = atoms
            
        
        self.pop.sort(key=lambda x:(x.info['key_value_pairs']['raw_score']),reverse = True)
        self.dbi.update_to_relaxed(atoms)
    
    def refresh_populations(self):
        ids = [x.info['key_value_pairs']['dbid'] for x in self.pop]
        self.pop = []
        for i in ids:
            self.pop.append(self.dbi.get_atoms_from_id(i))

        self.pop.sort(key=lambda x:(x.info['key_value_pairs']['raw_score']),reverse = True)


    def __check_other_confids(self,atoms,popconfids):
        compatoms = self.dbi.get_other_confids_atoms(popconfids)
        if(len(compatoms) == 0): 
            self.confid += 1
            return self.confid
        for at in compatoms:
            if self.get_similarity(atoms, at):
                return at.info['key_value_pairs']['confid']
        self.confid += 1
        return self.confid

    def get_similarity(self,a1,a2):
        if(len(a1)!= len(a2)): return False

        comp = InteratomicDistanceComparator(n_top=len(a1), pair_cor_cum_diff=0.030,
                pair_cor_max=1.0, dE=0.5, mic=True)
        if comp.looks_like(a1,a2): return True
        return False


    def get_better_candidates(self,n=2):

        """Calculates the fitness using the formula from
            L.B. Vilhelmsen et al., JACS, 2012, 134 (30), pp 12807-12816
            Applied to the fitness function instead of the energy of the particles

            If weighted is selected it adds a variation to the fitness function that accounts for how many times the strucute has been chose for pairing
            and for how many times structures with the same stoichiometry exist in the population 

             
        """
        self.refresh_populations()

        "How many times the same stoichiometry appears in the run"
        
        raw_scores = [ x.info['key_value_pairs']['raw_score'] for x in self.pop]
        max_score = max(raw_scores)
        min_score = min(raw_scores)

        atoms = []
        for at in self.pop:
            atoms.append(at.copy())

        T = min_score - max_score

        
        atoms.sort(key=lambda x: (0.5 * (1. - tanh(2. * (x.info['key_value_pairs']['raw_score']-max_score)/ T - 1.))) *
            1.0/sqrt(1.0 + self.dbi.stc_count(x.info['key_value_pairs']['var_stc'])) *
            1.0/sqrt(1.0 + x.info['key_value_pairs']['parent_penalty']) * 
            1.0/sqrt(1.0 + self.dbi.confid_count(x.info['key_value_pairs']['confid'])),reverse = True)

        """elif(weighted and not structure_similarity):
            atoms.sort(key=lambda x: (0.5 * (1. - tanh(2. * (x.info['key_value_pairs']['raw_score']-max_score)/ T - 1.))) *
             1.0/sqrt(1.0+self.__times_paired(x.info['key_value_pairs']['var_stc'])) *
             1.0/sqrt(1.0+x.info['key_value_pairs']['parent_penalty']),reverse = True)
        elif(not weighted and not structure_similarity):
            atoms.sort(key=lambda x: (0.5 * (1. - tanh(2. * (x.info['key_value_pairs']['raw_score']-max_score)/ T - 1.))),reverse = True)
        """

        
        if(len(atoms)>n):
            return list(atoms[:n])
        else:
            return list(atoms[:len(atoms)-1])
        