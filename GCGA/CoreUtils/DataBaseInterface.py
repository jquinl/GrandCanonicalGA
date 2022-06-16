from distutils.command.config import config
from pathlib import Path
from math import tanh, sqrt, exp

from ase.db import connect
class DataBaseInterface:

    "Object to interface between the databese and the GCGA code in order to decouple database handling from the GA algortihm"
    def __init__(self,db_name = "database.db"):

        self.db_name =self.__get_db_name(db_name)
        self.db = connect(self.db_name)

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


    def get_other_confids_atoms(self,confids):
        popids = []
        for i in confids:
            popids.extend(self.get_ids_from_confid(i))

        allids = self.get_all_relaxed_ids()
        finalids = [i for i in allids if i not in popids]
        atoms = []
        for i in finalids:
            atoms.append(self.get_atoms_from_id(i))
        return list(atoms)


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
    

    """def get_better_candidates(self,n=1,max_num=False,weighted = False,structure_similarity = False):

        Calculates the fitness using the formula from
            L.B. Vilhelmsen et al., JACS, 2012, 134 (30), pp 12807-12816
            Applied to the fitness function instead of the energy of the particles

            If weighted is selected it adds a variation to the fitness function that accounts for how many times the strucute has been chose for pairing
            and for how many times structures with the same stoichiometry exist in the population 

             
        
        atoms = []
        if(structure_similarity and not weighted):
            weighted = True
        
        
        atoms = self.get_relaxed_candidates()


        "How many times the same stoichiometry appears in the run"
        
        raw_scores = [ x.info['key_value_pairs']['raw_score'] for x in atoms]
        max_score = max(raw_scores)
        min_score = min(raw_scores)

        T = min_score - max_score

        if len(atoms) > 100:
            atoms = self.get_relaxed_candidates()[:100]


        if( weighted and structure_similarity ):
            atoms.sort(key=lambda x: (0.5 * (1. - tanh(2. * (x.info['key_value_pairs']['raw_score']-max_score)/ T - 1.))) *
             1.0/sqrt(1.0 + self.__times_paired(x.info['key_value_pairs']['var_stc'])) *
             1.0/sqrt(1.0 + x.info['key_value_pairs']['parent_penalty']) * 
             1.0/sqrt(1.0 + x.info['key_value_pairs']['similarity']),reverse = True)
        elif(weighted and not structure_similarity):
            atoms.sort(key=lambda x: (0.5 * (1. - tanh(2. * (x.info['key_value_pairs']['raw_score']-max_score)/ T - 1.))) *
             1.0/sqrt(1.0+self.__times_paired(x.info['key_value_pairs']['var_stc'])) *
             1.0/sqrt(1.0+x.info['key_value_pairs']['parent_penalty']),reverse = True)
        elif(not weighted and not structure_similarity):
            atoms.sort(key=lambda x: (0.5 * (1. - tanh(2. * (x.info['key_value_pairs']['raw_score']-max_score)/ T - 1.))),reverse = True)


        if max_num == False:
            if(len(atoms)>n):
                return list(atoms[:n])
            else:
                return list(atoms[:len(atoms)-1])
        else:
            return list(atoms[:len(atoms)-1]) """

    def __get_db_name(self,db_name):
        if Path(db_name ).is_file():
            Path(db_name ).unlink()
        return db_name