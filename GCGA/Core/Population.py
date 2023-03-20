from math import tanh, sqrt, exp
import random
import numpy as np
class Population:
    def __init__(self, population_size, database_interface,steps):
        self.pop_size = population_size
        self.dbi = database_interface
        self.run_stcs = []
        self.pop = []
        self.pop_stc = []

        self.current_stc = 0
        self.pairs = []

        self.stc_atempts = 0
        self.ctr = 0
        self.steps = steps

    def update_population(self,atoms):
        stc = atoms.info['key_value_pairs']['var_stc']
        if(self.current_stc != stc):
            self.change_current_stc(stc)

        atoms.info['key_value_pairs']['confid'] = self.dbi.get_atom_confid(atoms)

        if(len(self.pop)<self.pop_size):
            self.pop.append(atoms)
        else:
            if(self.pop[-1].info['key_value_pairs']['raw_score'] < atoms.info['key_value_pairs']['raw_score']):
                self.pop[-1] = atoms

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
        self.ctr+=1

    def should_change_stc(self):
        if(len(self.pop)<2): return False

        self.refresh_populations()
        
        max_stc = 0
        max_score = 0.0
        for j,i in enumerate(self.pop_stc):
            if(i.info['key_value_pairs']['raw_score']> max_score):
                max_score = i.info['key_value_pairs']['raw_score']
                max_stc = i.info['key_value_pairs']['var_stc']

        if(len(self.pop_stc) < 2): return False

        max_stc_dist = float(max([self.stc_distance(x,self.pop_stc[max_stc]) for x in self.pop_stc]))

        dist = self.stc_distance(self.pop[0],self.pop_stc[max_stc]) / max_stc_dist
        proportion = float(len(self.dbi.get_ids_from_stc(self.current_stc))) / float(len(self.dbi.get_relaxed_candidates())+1)

        should_change = np.random.random() * (2.0-dist) * (1.0 + self.ctr/self.steps) < (proportion + (self.stc_atempts * 0.01))
        if(self.stc_atempts <5):
            should_change =  False
        if(should_change):
            self.stc_atempts = 0
        else:
            self.stc_atempts += 1

        return should_change

    def change_current_stc(self,stc):
        self.current_stc = stc
        all_atoms =  self.dbi.get_atoms_with_stc(stc)
        all_atoms.sort(key=lambda x:(x.info['key_value_pairs']['raw_score']),reverse = True)

        if(len(all_atoms) < self.pop_size):
            self.pop= all_atoms
        else:
            self.pop = all_atoms[:self.pop_size]

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

    def _get_fitness(self):
        """Calculates the fitness using the formula from
            L.B. Vilhelmsen et al., JACS, 2012, 134 (30), pp 12807-12816
            Applied to the fitness function instead of the energy of the particles
        """
        self.refresh_populations()

        indices = list(range(len(self.pop)))

        raw_scores = [ x.info['key_value_pairs']['raw_score'] for x in self.pop]
        max_score = max(raw_scores)
        min_score = min(raw_scores)

        T = min_score - max_score

        if(T == 0.0000): T = 1.0

        f = [float(0.5 * (1. - tanh(2. * (raw_scores[i] - max_score) / T - 1.))) for i in indices]
        P = [float(self.pop[i].info['key_value_pairs']['parent_penalty']) for i in indices]
        S = [float(len(self.dbi.get_ids_from_confid(self.pop[i].info['key_value_pairs']['confid']))) for i in indices]

        fit = [(f[i] * 1. / sqrt(1. + P[i]) * 1. / sqrt(1. + S[i])) for i in range(len(f))]
        return fit

    def get_better_candidates(self):
        fit = self._get_fitness()
        fmax = max(fit)

        atoms = [at.copy() for at in self.pop]

        c1 = atoms[0]
        c2 = atoms[0]
        pairs = (0,0)
        used_before = False

        tries = 0
        while c1.info['key_value_pairs']['dbid'] == c2.info['key_value_pairs']['dbid'] or used_before:
            tries +=1
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

            c1id = c1.info['key_value_pairs']['dbid']
            c2id = c2.info['key_value_pairs']['dbid']
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

        all_cands = float(len(self.dbi.get_relaxed_candidates())+1)
        proportions ={}
        for at in self.pop_stc:
            proportions[at.info['key_value_pairs']['var_stc']] = float(len(self.dbi.get_ids_from_stc(at.info['key_value_pairs']['var_stc'])))/all_cands 

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

        max_stc_dist = max([self.stc_distance(x,self.pop_stc[max_stc]) for x in self.pop_stc])

        raw_scores = [ x.info['key_value_pairs']['raw_score'] for x in self.pop_stc]
        max_len = max([len(x) for x in self.pop_stc])
        min_len = min([len(x) for x in self.pop_stc])
        M = max_len - min_len
        if(M == 0): M = 1.0
        max_score = max(raw_scores)
        min_score = min(raw_scores)

        atoms = list([at.copy() for at in self.pop_stc])

        T = max_score - min_score
        if(T == 0.0): T= 1.0

        mean = self.pop_stc[max_stc].info['key_value_pairs']['raw_score']
        atoms.sort(key=lambda x:( ((x.info['key_value_pairs']['raw_score']-min_score )/ T) +
                (1.0 - (self.stc_distance(x,self.pop_stc[max_stc])/max_stc_dist)) +
               (len(x) - min_len)/M),reverse=True)

        scl = (1.4 - (self.ctr/self.steps))*len(self.pop_stc)
        next_stc = self.current_stc
        while next_stc == self.current_stc:
            next_stc = abs(int(np.random.normal(0,scale=scl)))%len(self.pop_stc)

       
        return atoms[next_stc],atoms[(next_stc+1) % len(self.pop_stc)]

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
