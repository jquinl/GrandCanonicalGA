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
        print(dbid)
        self.db.update(dbid,atoms,dbid=dbid)
    

    def update_to_relaxed(self, dbid, atoms):
        if(dbid != atoms.info['dbid']):
            raise Exception("Atoms objec dbid and database dbid do not coincide")
        try:
            atoms.info['raw_score']
            self.db.update(dbid,atoms=atoms,relaxed = True)

        except KeyError:
           raise Exception("Relaxed candidate doe snot have raw_score")


    "THIS NEEDS FIXING"
    def get_atoms_from_id(self,dbid):
        rows = self.db.select(dbid=dbid)
        for row in rows:
            return row.toatoms()

    def get_all_unrelaxed_ids(self):
        return list(set([t.dbid for t in self.db.select(relaxed = False)]))
    def get_all_relaxed_ids(self):
        return list(set([t.dbid for t in self.db.select(relaxed = True)]))

    def get_number_of_unrelaxed_candidates(self):
        return len(set([t.dbid for t in self.db.select(relaxed = False)]))
    def get_number_of_relaxed_candidates(self):
        return len(set([t.dbid for t in self.db.select(relaxed = True)]))

    def get_relaxed_candidates(self):
        atoms = []
        for i in self.get_all_relaxed_ids():
            atoms.append(self.get_atoms_from_id(i))
        return atoms
    
    def get_better_candidates(self,n=1):

        atoms = self.get_relaxed_candidates()
        atoms.sort(key=lambda x: x.info['raw_score'],reverse = True)

        return list(atoms[:n])

    def __get_db_name(self,db_name):
        if Path(db_name ).is_file():
            Path(db_name ).unlink()
        return db_name