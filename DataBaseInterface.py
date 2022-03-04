from pathlib import Path
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

    

    def update_to_relaxed(self, dbid, atoms):
        try:
            atoms.info['key_value_pairs']['dbid']
        except:
            raise Exception("No dbid Parameter found in fetched atoms")
        try:
            atoms.info['key_value_pairs']['var_stc']
        except KeyError:
            raise Exception("No var_stc Parameter found in fetched atoms")
        if(dbid != atoms.info['key_value_pairs']['dbid']):
            raise Exception("Atoms object dbid and database dbid do not coincide")
        try:
            atoms.info['key_value_pairs']['raw_score']
        except KeyError:
           raise Exception("Relaxed candidate does not have raw_score")

        self.db.update(dbid,atoms=atoms,relaxed = True,raw_score = atoms.info['key_value_pairs']['raw_score'])


    "THIS NEEDS FIXING"
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

    def get_all_unrelaxed_ids(self):
        return list(set([t.dbid for t in self.db.select(relaxed = False)]))
    def get_all_relaxed_ids(self):
        return list(set([t.dbid for t in self.db.select(relaxed = True)]))

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
                raise Exception("No dbid Parameter found in fetched atoms")
            try:
                at.info['key_value_pairs']['var_stc']
            except KeyError:
                raise Exception("No var_stc Parameter found in fetched atoms")
            try:
                at.info['key_value_pairs']['raw_score']
            except KeyError:
                raise Exception("No var_stc Parameter found in fetched atoms")

            atoms.append(at)
        return list(atoms)
    
    def get_better_candidates(self,n=1):

        atoms = self.get_relaxed_candidates()
        atoms.sort(key=lambda x: x.info['key_value_pairs']['raw_score'],reverse = True)

        return list(atoms[:n])

    def __get_db_name(self,db_name):
        if Path(db_name ).is_file():
            Path(db_name ).unlink()
        return db_name