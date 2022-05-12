from ase import Atoms
from abc import ABC,abstractmethod
from CoreUtils.SubunitAnalysis import SubunitFinder

class BaseFitness(ABC):



    def __init__(self):
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


 
class GibbsFreeEnergy(BaseFitness):

    supported_environment_variable_types = ["mu"]

    def __init__(self,environmental_variables, references):
        super().__init__()
        self.env = self.__get_env(environmental_variables)
        self.references = self.__get_refs(references)
        
    
    def evaluate(self, atoms) -> float:
        super().evaluate(atoms)

        return self.references
        
    def __get_refs(self,references):
        if(type(references) is not dict):raise TypeError("references is not a string:float dictionary")
        for a,k in references.items():
            if(not (type(a) == str and type(k) == float)):
                raise TypeError("references is not a string: float dictionary")
            try:
                if(a is not None):
                    at = Atoms(a)
            except:
                raise ValueError("Cannot create atoms with key: " + a)

        return references

    def __get_env(self,env):
        if(type(env) is not dict):raise TypeError("environment is not a string:float dictionary")

        for e,k in env.items():
                if(not (type(e) == str and type(k) == float)):
                    raise TypeError("environment is not a string:float dictionary")

        for i in env.keys():
            x = i.split("_")
            if(len(x) != 2): raise ValueError("Environment variable: " + x + " Does not follow the rule [var]_[atoms] p.e = mu_CO")
            if(x[0] not in self.supported_environment_variable_types ): raise ValueError(" unsupported environment variable in : " + x + " Supported types are: ".join(str(i) for i in self.supported_environment_variable_types))
            try:
                Atoms(x[1])
            except:
                raise ValueError(" unsupported atom types in : " + x + " Supported types are: ".join(str(i) for i in self.supported_environment_variable_types))



        return env