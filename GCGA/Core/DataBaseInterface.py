from pathlib import Path
from math import tanh, sqrt, exp
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.db import connect
class DataBaseInterface:

    "Object to interface between the databese and the GCGA code in order to decouple database handling from the GA algortihm"
    def __init__(self,db_name = "database.db"):

        self.db_name =self.__get_db_name(db_name)
        self.db = connect(self.db_name)
        self.confid = 0

    def add_unrelaxed_candidate(self,atoms):
        try:
            atoms.info['stc']
        except:
            raise Exception("stc parameter not set in unrelaxed candidate")
        dbid = self.db.write(atoms,relaxed = False, var_stc = atoms.info['stc'])
        atoms.info['dbid'] = dbid
        self.db.update(dbid,atoms = atoms,dbid=dbid)

        try:
            atoms =  self.db.get_atoms(dbid=dbid,add_additional_information=True)
        except:
            raise Exception("No atoms with that dbid")
        try:
                atoms.info['key_value_pairs']['dbid']
        except KeyError:
            raise Exception("No dbid Parameter found in fetched atoms")


    def update_to_relaxed(self, atoms):
        try:
            atoms.info['key_value_pairs']['dbid']
        except:
            raise Exception("No dbid Parameter found in fetched atoms")
        try:
            atoms.info['key_value_pairs']['var_stc']
        except KeyError:
            raise Exception("No var_stc Parameter found in fetched atoms")
        try:
            atoms.info['key_value_pairs']['raw_score']
        except KeyError:
           raise Exception("Relaxed candidate does not have raw_score")
        try:
            atoms.info['key_value_pairs']['confid']
        except KeyError:
           raise Exception("Relaxed candidate does not have confid")

        self.db.update(atoms.info['key_value_pairs']['dbid'],atoms=atoms,relaxed = True,raw_score = atoms.info['key_value_pairs']['raw_score'],
            parent_penalty = 0,confid = atoms.info['key_value_pairs']['confid'])


    def update_penalization(self,atoms):
        try:
            atoms.info['key_value_pairs']['dbid']
        except:
            raise Exception("No atoms object was not included in database before")
        dbid = atoms.info['key_value_pairs']['dbid']
        self.db.update(dbid,parent_penalty = atoms.info['key_value_pairs']['parent_penalty'] + 1)

    """def update_similarity(self,atoms):
        try:
            atoms.info['key_value_pairs']['dbid']
        except:
            raise Exception("No atoms object was not included in database before")
        dbid = atoms.info['key_value_pairs']['dbid']
        self.db.update(dbid,similarity = atoms.info['key_value_pairs']['similarity'] + 1)"""


    def get_atoms_from_id(self,dbid):
        try:
            atoms =  self.db.get_atoms(dbid=dbid,add_additional_information=True)
        except:
            raise Exception("No atoms with that dbid")
        try:
            atoms.info['key_value_pairs']['dbid']
        except KeyError:
            raise Exception("No dbid Parameter found in fetched atoms")
        try:
            atoms.info['key_value_pairs']['var_stc']
        except KeyError:
            raise Exception("No var_stc Parameter found in fetched atoms")
        return atoms

    def get_atoms_with_stc(self,stc):
        ids = self.get_ids_from_stc(stc)
        if(len(ids) == 0 ): return list([])
        atoms = []
        for i in ids:
            atoms.append(self.get_atoms_from_id(i))
        return list(atoms)

    def get_atom_confid(self,atom):
        compare_ids = []
        for i in range(self.confid + 1):
            ids = self.get_ids_from_confid(i)
            if(len(ids)>0):
                compare_ids.append(ids[0])
        for i in compare_ids:
            cmp = self.get_atoms_from_id(i)
            if(self.get_similarity(atom,cmp)):
                return cmp.info['key_value_pairs']['confid']

        return self.new_confid()

    def new_confid(self):
        self.confid += 1
        return self.confid


    def get_all_unrelaxed_ids(self):
        return list(set([t.dbid for t in self.db.select(relaxed = False)]))
    def get_all_relaxed_ids(self):
        return list(set([t.dbid for t in self.db.select(relaxed = True)]))
    def get_ids_from_stc(self,stc):
        return list(set([t.dbid for t in self.db.select(var_stc = stc,relaxed = True)]))

    def get_ids_from_confid(self,confid):
        return list(set([t.dbid for t in self.db.select(confid = confid,relaxed = True)]))
    def get_number_of_unrelaxed_candidates(self):
        return len(set([t.dbid for t in self.db.select(relaxed = False)]))
    def get_number_of_relaxed_candidates(self):
        return len(set([t.dbid for t in self.db.select(relaxed = True)]))


    def get_first_unrelaxed(self):
        id_list = self.get_all_unrelaxed_ids()

        if(len(id_list) == 0): raise Exception("No unrelaxed candidates")

        atom_id = id_list[0]

        atoms = self.get_atoms_from_id(atom_id)

        return atoms

    def get_relaxed_candidates(self):
        atoms = []
        for i in self.get_all_relaxed_ids():
            at = self.get_atoms_from_id(i)
            try:
                at.info['key_value_pairs']['dbid']
            except KeyError:
                raise Exception("No dbid parameter found in fetched atoms")
            try:
                at.info['key_value_pairs']['var_stc']
            except KeyError:
                raise Exception("No var_stc parameter found in fetched atoms")
            try:
                at.info['key_value_pairs']['raw_score']
            except KeyError:
                raise Exception("No var_stc parameter found in fetched atoms")
            try:
                at.info['key_value_pairs']['parent_penalty']
            except KeyError:
                raise Exception("No var_stc parameter found in fetched atoms")
            atoms.append(at)

        return list(atoms)

    def confid_count(self,confid):
        atoms = []
        for i in self.get_ids_from_confid(confid):
            at = self.get_atoms_from_id(i)
            try:
                at.info['key_value_pairs']['dbid']
            except KeyError:
                raise Exception("No dbid parameter found in fetched atoms")
            try:
                at.info['key_value_pairs']['var_stc']
            except KeyError:
                raise Exception("No var_stc parameter found in fetched atoms")
            try:
                at.info['key_value_pairs']['raw_score']
            except KeyError:
                raise Exception("No var_stc parameter found in fetched atoms")
            atoms.append(at)

        return len(list(atoms))

    def stc_count(self,stc):
        atoms = []
        for i in self.get_ids_from_stc(stc):
            at = self.get_atoms_from_id(i)
            try:
                at.info['key_value_pairs']['dbid']
            except KeyError:
                raise Exception("No dbid parameter found in fetched atoms")
            try:
                at.info['key_value_pairs']['var_stc']
            except KeyError:
                raise Exception("No var_stc parameter found in fetched atoms")
            try:
                at.info['key_value_pairs']['raw_score']
            except KeyError:
                raise Exception("No var_stc parameter found in fetched atoms")
            atoms.append(at)

        return len(list(atoms))

    def get_better_candidates_raw(self,n=1,max_num=False):

        atoms = self.get_relaxed_candidates()
        atoms.sort(key=lambda x: x.info['key_value_pairs']['raw_score'],reverse = True)

        if max_num == False:
            if(len(atoms)>n):
                return list(atoms[:n])
            else:
                return list(atoms[:len(atoms)-1])
        else:
            return list(atoms[:len(atoms)-1])

    def __get_db_name(self,db_name):
        if Path(db_name ).is_file():
            Path(db_name ).unlink()
        return db_name

    def get_similarity(self,a1,a2):
        if(a1.info['key_value_pairs']['var_stc'] != a2.info['key_value_pairs']['var_stc']) :return False

        comp = InteratomicDistanceComparator(n_top=len(a1), pair_cor_cum_diff=0.015,
                pair_cor_max=0.5, dE=0.15, mic=True)
        if comp.looks_like(a1,a2): return True
        return False