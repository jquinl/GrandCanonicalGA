from ase import Atoms
from abc import ABC,abstractmethod

class BaseFitness(ABC):

    def __init__(self,environmental_variables,references):
        self.env = self.__get_env(environmental_variables)
        self.references = self.__get_refs(**references)
        pass

    @classmethod
    def fitness_class(self):
        return True
    
    def fitness_instance(self):
        return True

    @abstractmethod
    def evaluate(self,atoms)-> float:
        if(not issubclass(atoms,Atoms)):
            raise TypeError("Fitness function tried to evaluate something other than an ASE atoms object")
        
        pass
    
    def __get_refs(self,**references):
        if(type(references) is not dict):raise TypeError("references is not a string:float dictionary")
        for a,k in references.items():
            if(not (type(a) == str and type(k) == float)):
                raise TypeError("references is not a string:float dictionary")
            try:
                at = Atoms(a)
            except:
                raise ValueError("Cannot create atoms with key: " + a)

        return references

    def __get_env(self,env):
        if(type(env) is not dict):raise TypeError("environment is not a string:float dictionary")
        for e,k in env.items():
                if(not (type(e) == str and type(k) == float)):
                    raise TypeError("environment is not a string:float dictionary")
        return env

 
class GibbsFreeEnergy(BaseFitness):
    def __init__(self, environmental_variables, references):
        super().__init__(environmental_variables, references)

    def evaluate(self, atoms) -> float:
        super().evaluate(atoms)
    
        return self.references
