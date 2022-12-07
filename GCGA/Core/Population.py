from ase.ga.standard_comparators import InteratomicDistanceComparator
from math import tanh, sqrt, exp
import random
import numpy as np
class Population:
    def __init__(self, population_size, database_interface, stc_hopping=0.3):
        self.pop_size = population_size
        self.dbi = database_interface
        self.run_stcs = []
        self.pop = []
        self.pop_stc = []
        self.confid = 0
        self.current_stc = 0
        self.max_stc = 0
        self.min_stc = 100
        self.pairs = []

        self.stc_hop = stc_hopping
        self.stc_atempts = 0
        self.oldstc = None
    def initialize_population(self):
        self.change_current_stc(self.pop_stc[0].info['key_value_pairs']['var_stc'])

    def update_population(self,atoms):
        isSimilar = False

        if(len(self.pop) == 0):
            atoms.info['key_value_pairs']['confid'] =  self.confid
            self.confid += 1
            self.pop.append(atoms)
        elif(atoms.info['key_value_pairs']['var_stc'] == self.current_stc):
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
        else:
            confids = [a.info['key_value_pairs']['confid'] for a in self.pop ]
            atoms.info['key_value_pairs']['confid'] =  self.__check_other_confids(atoms,confids)

        self.pop.sort(key=lambda x:(x.info['key_value_pairs']['raw_score']),reverse = True)
        if(len(self.pop_stc) == 0):
            self.pop_stc.append(atoms)
        else:
            stcs = [i.info['key_value_pairs']['var_stc'] for i in self.pop_stc]
            if(atoms.info['key_value_pairs']['var_stc'] not in stcs):
                self.pop_stc.append(atoms)
            else:
                for i in range(len(self.pop_stc)):
                    if(self.pop_stc[i].info['key_value_pairs']['var_stc'] == atoms.info['key_value_pairs']['var_stc'] and
                    atoms.info['key_value_pairs']['raw_score']>self.pop_stc[i].info['key_value_pairs']['raw_score'] ):
                        self.pop_stc[i] = atoms

        self.pop_stc.sort(key=lambda x:x.info['key_value_pairs']['raw_score'],reverse=True)
        stcs = [i.info['key_value_pairs']['var_stc'] for i in self.pop_stc]

        self.dbi.update_to_relaxed(atoms)
        if(self.current_stc not in stcs):
            self.current_stc = self.pop_stc[0].info['key_value_pairs']['var_stc']

    def should_change_stc(self):
        if(self.oldstc is None):
            self.oldstc  =  self.current_stc
        if(self.oldstc  ==  self.current_stc):
            self.stc_atempts += 1
        else:
            self.stc_atempts = 0

        if(len(self.pop_stc)<2):return False

        if (np.random.random() < (self.stc_hop/ 0.001 + abs(self.pop[0].info['key_value_pairs']['raw_score'] - self.pop[-1].info['key_value_pairs']['raw_score']))):
            return True

        if(abs(self.pop[0].info['key_value_pairs']['raw_score'] - self.pop[-1].info['key_value_pairs']['raw_score']) < 0.01 and self.stc_atempts > 10):
            return True

        return False

    def target_stc(self):
        stcs = [i.info['key_value_pairs']['var_stc'] for i in self.pop_stc]
        if(np.random.random() < 0.1):
                random.shuffle(stcs)
                if(len(stcs) == 1):return stcs[0]
                if(stcs[0]== self.current_stc): return stcs[1]
                return stcs[0]
        a1,a2  = self.get_better_candidates_stc()
        if(a1.info['key_value_pairs']['var_stc'] == self.current_stc  ): return a2.info['key_value_pairs']['var_stc']
        return a1.info['key_value_pairs']['var_stc']

    def change_current_stc(self,stc):
        self.current_stc = stc
        new_pop =  self.dbi.get_atoms_with_stc(stc)
        new_pop.sort(key=lambda x:(x.info['key_value_pairs']['raw_score']),reverse = True)
        if(len(new_pop) < self.pop_size):
            self.pop= new_pop
        else:
            self.pop = new_pop[:self.pop_size]

    def get_similarity(self,a1,a2):
        if(a1.info['key_value_pairs']['var_stc'] != a2.info['key_value_pairs']['var_stc']) :return False

        comp = InteratomicDistanceComparator(n_top=len(a1), pair_cor_cum_diff=0.030,
                pair_cor_max=1.0, dE=0.5, mic=True)
        if comp.looks_like(a1,a2): return True
        return False

    def refresh_populations(self):
        ids = [x.info['key_value_pairs']['dbid'] for x in self.pop]
        self.pop = []
        for i in ids:
            self.pop.append(self.dbi.get_atoms_from_id(i))

        self.pop.sort(key=lambda x:(x.info['key_value_pairs']['raw_score']),reverse = True)

        ids = [x.info['key_value_pairs']['dbid'] for x in self.pop_stc]
        self.pop_stc = []
        for i in ids:
            self.pop_stc.append(self.dbi.get_atoms_from_id(i))

    def get_better_candidate(self):
        self.refresh_populations()
        atoms = [at.copy() for at in self.pop]
        return atoms[0]

    def _get_fitness(self):
        """Calculates the fitness using the formula from
            L.B. Vilhelmsen et al., JACS, 2012, 134 (30), pp 12807-12816
            Applied to the fitness function instead of the energy of the particles
        """
        self.refresh_populations()

        "How many times the same stoichiometry appears in the run"

        raw_scores = [ x.info['key_value_pairs']['raw_score'] for x in self.pop]
        max_score = max(raw_scores)
        min_score = min(raw_scores)

        atoms = [at.copy() for at in self.pop]

        T = min_score - max_score

        fit = []
        for at in atoms:
            fit.append((0.5 * (1. - tanh(2. * (at.info['key_value_pairs']['raw_score']-max_score)/ T - 1.))) *
            1.0/sqrt(1.0 + at.info['key_value_pairs']['parent_penalty']) * 
            1.0/sqrt(1.0 + self.dbi.confid_count(at.info['key_value_pairs']['confid'])))
        return fit

    def get_better_candidates(self):
        fit = self._get_fitness()
        fmax = max(fit)

        atoms = [at.copy() for at in self.pop]

        c1 = atoms[0]
        c2 = atoms[0]
        pairs = (0,0)
        used_before = False
        while c1.info['key_value_pairs']['confid'] == c2.info['key_value_pairs']['confid'] and not used_before:
            nnf = True
            while nnf:
                t = np.random.randint(len(atoms))
                if fit[t] > np.random.rand() * fmax:
                    c1 = atoms[t]
                    nnf = False
            nnf = True
            while nnf:
                t = np.random.randint(len(atoms))
                if fit[t] > np.random.rand() * fmax:
                    c2 = atoms[t]
                    nnf = False

            c1id = c1.info['key_value_pairs']['confid']
            c2id = c2.info['key_value_pairs']['confid']
            pairs = (min([c1id, c2id]), max([c1id, c2id]))
            used_before = pairs in self.pairs

        self.pairs.append(pairs)
        return c1,c2


    def get_better_candidates_stc(self):
        """Calculates the fitness using the formula from
            L.B. Vilhelmsen et al., JACS, 2012, 134 (30), pp 12807-12816
            Applied to the fitness function instead of the energy of the particles
            This selects from a population containing the best candidate for each stoichiometry
        """
        self.refresh_populations()

        en_dic = {}
        max_stc = 0
        max_score = 0.0
        current_stc_pos = 0
        for j,i in enumerate(self.pop_stc):
            en_dic[i.info['key_value_pairs']['var_stc']] = i.info['key_value_pairs']['raw_score']
            if(i.info['key_value_pairs']['raw_score']> max_score):
                max_score = i.info['key_value_pairs']['raw_score']
                max_stc = i.info['key_value_pairs']['var_stc']

            if(i.info['key_value_pairs']['var_stc'] == self.current_stc):
                current_stc_pos = j

        raw_scores = [ x.info['key_value_pairs']['raw_score'] for x in self.pop_stc]
        max_score = max(raw_scores)
        min_score = min(raw_scores)

        atoms = [at.copy() for at in self.pop_stc]

        T = min_score - max_score

        atoms.sort(key=lambda x: (0.5 * (1. - tanh(2. * (x.info['key_value_pairs']['raw_score']-max_score)/ T - 1.))) *
            0.5/sqrt(1.0 + self.stc_distance(x,self.pop_stc[max_stc])) * 
            sqrt(1.0 + self.stc_distance(x,self.pop_stc[current_stc_pos])),reverse = True)
        if(atoms[0].info['key_value_pairs']['var_stc'] == self.current_stc and len(atoms)>2):
            return atoms[1],atoms[2]
        elif(atoms[0].info['key_value_pairs']['var_stc'] == self.current_stc):
            return atoms[1],atoms[0]
        else:
            return atoms[0],atoms[1]

        #return atoms[0],atoms[1]

    def stc_distance(self,at1,at2):
        import math
        d1 = at1.symbols.indices()
        d2 = at2.symbols.indices()
        d1_dist = {}
        d2_dist = {}
        for i,j in d1.items():
            d1_dist[i] = len(j)
            d2_dist[i] = 0
        for i,j in d2.items():
            d2_dist[i] = len(j)
            if(i not in d1_dist.keys()):
                d1_dist[i] = 0
        suma = 0
        for i in d1_dist.keys():
            dd = (d1_dist[i] - d2_dist[i])**2
            suma += dd
        return math.sqrt(suma)

    def update_penalization(self,atoms1,atoms2):
        self.dbi.update_penalization(atoms1)
        self.dbi.update_penalization(atoms2)

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

