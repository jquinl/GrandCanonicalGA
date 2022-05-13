from itertools import count
from tkinter.tix import Tree
from ase import Atoms
from abc import ABC,abstractmethod
from CoreUtils.SubunitAnalysis import SubunitFinder
import numpy as np

class BaseFitness(ABC):

    def __init__(self):
        pass

    @classmethod
    def fitness_class(cls):
        return True
    
    def fitness_instance(self):
        return True

    @abstractmethod
    def evaluate(self,atoms)-> float:
        if(not issubclass(atoms,Atoms)):
            raise TypeError("Fitness function tried to evaluate something other than an ASE atoms object")
        pass
    
    def __mantains_ordering(self,slab, atoms):
        if(len(atoms) < len(slab)):
            return False
        if(slab.symbols.indices() != atoms[:len(slab)].symbols.indices()):
            return False

        return True


class GibbsFreeEnergy(BaseFitness):

    supported_environment_variable_types = ["mu"]
    
    def __init__(self, variable_reference_energies,environmental_variables,base_reference_energy= 0.0):
        super().__init__()
        self.multi_ref = False
        self.base_ref = float(base_reference_energy)
        self.references = self.__get_refs(variable_reference_energies)
        self.references.sort()
        self.env = self.__get_env(environmental_variables)

        if len(self.references) == 0: raise ValueError("No reference energies found")
        if len(self.env) == 0: raise ValueError("No environmental variables found")
    
    def evaluate(self,slab, atoms) -> float:
        super().evaluate(atoms)
        energy = atoms.get_potential_energy()
        fitness = 0.0
        if(slab is not None):
            if(not self.__mantains_ordering(slab,atoms)):
                raise ValueError("Slab atoms")
        
        at = atoms[len(slab):]
        indices = at.symbols.indices()
        self.__compare_indices_to_ref(indices)

        if(not self.multi_ref):
            fitness = energy - self.base_ref
            for i,j in indices.items():
                fitness -=  (self.references[i] * len(j))
            for i,j in self.env[0].items():
                fitness -= (len(indices[i]) * j)
        else:
            fitness = energy - self.base_ref
            counts = self.__get_molecule_counts(indices)
            for i,j in indices.items():
                fitness -=  (self.references[i] * counts[i])
            for i,j in self.env[0].items():
                fitness -= (len(indices[i]) * counts[j])

        return fitness
        
    def __get_molecule_counts(self,atoms,indices):
        atoms_obj = atoms.copy()
        ats = [Atoms(i) for i in self.references.keys()]
        ats.sort(reverse=True,key= lambda x: len(x))
        print(ats)
        counts = {}
        for i in range(len(ats)):
            mat = []
            for j in  ats[i].numbers:
                sym_array = []
                for k in range(len(atoms_obj.numbers)):
                    if(j == atoms_obj.numbers[k]):
                        sym_array.append(k)
                sym_array.sort(reverse=True)
                mat.append(sym_array)
            shortest = min([len(j) for j in mat])
            counts[self.references.keys()[i]] = shortest
            counts.append(shortest)
            for k in mat:
                for j in range(shortest):
                    atoms_obj[k[j]].number = 200
        return counts

    def fitness_function(atoms)-> float:
        env = 1.0
        ref=read('pt.traj@:')[0].get_potential_energy()    

        au_en =read('gold_bulk.traj').get_potential_energy() / 500.0
        at_num= atoms.get_atomic_numbers()
        pt_num = np.count_nonzero(at_num == 78)
        au_num= np.count_nonzero(at_num == 79)

        fre = atoms.get_potential_energy() - ref- au_num * (au_en) - au_num*env

        return -fre


    def __compare_indices_to_ref(self,indices):
        indices.sort()
        for i in indices.keys():
            if (i not in self.references.keys()):
                self.references[i] = 0.0
            if(i not in self.env[0].keys()):
                self.env[0][i] = 0.0
        

                

    def __get_refs(self,references):
        if(type(references) is not dict):raise TypeError("references is not a string:float dictionary")
        for a,k in references.items():
            if(not (type(a) == str and type(k) == float)):
                raise TypeError("references is not a string: float dictionary")
            try:
                if(a is not None):
                    at = Atoms(a)
                    if(len(at)>1): self.multi_ref = True
            except:
                raise ValueError("Cannot create atoms with key: " + a)
        return references

    def __get_env(self,env):
        if(type(env) is not dict):raise TypeError("environment is not a string:float dictionary")

        for e,k in env.items():
                if(not (type(e) == str and type(k) == float)):
                    raise TypeError("environment is not a string:float dictionary")
        env_var = []
    
        for i in env.keys():
            x = i.split("_")
            if(len(x) != 2): raise ValueError("Environment variable: " + x + " Does not follow the rule [var]_[atoms] p.e = mu_CO")
            if(x[0] not in self.supported_environment_variable_types ): raise ValueError(" unsupported environment variable in : " + x + " Supported types are: ".join(str(i) for i in self.supported_environment_variable_types))
            try:
                Atoms(x[1])
            except:
                raise ValueError(" unsupported atom types in : " + x + " Supported types are: ".join(str(i) for i in self.supported_environment_variable_types))
            if(x[1] not in self.references.keys()):
                raise ValueError("Environmental variable " + x[1] + " not found in references")
        for i in self.supported_environment_variable_types:
            env_var.append({})
        for i,j in env.items():
            x = i.split("_")
            if(x[0] == "mu"):
                env_var[0][x[0]] = x[1]
        

        return list(env_var)