import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.ga.utilities import (closest_distances_generator)
import random


from .MutationsBase import OperationsBase
class PrepareForDB(OperationsBase):

    def __init__(self,slab,variable_types,variable_range,ratio_of_covalent_radii=0.7,
                   rng=np.random):
        super().__init__(slab,variable_types,variable_range,ratio_of_covalent_radii,rng)

        

    def prepare(self, atoms):

        if( not isinstance(atoms.calc,SinglePointCalculator)):
            return None        
                
        E = atoms.get_potential_energy()
        F = atoms.get_forces()
        results = {'energy': E,'forces': F}

        if(results is None):
            return None

        return_atoms = self.slab.copy()

        return_atoms.extend(self.sort_atoms_by_type(atoms[len(self.slab):]))

        if(self._check_overlap_all_atoms(return_atoms,self.blmin)):
            return None

        var_id = self.get_var_id(return_atoms)

        if(var_id is None):
            return None

    
        calc_sp = SinglePointCalculator(atoms, **results)
        return_atoms.set_calculator(calc_sp)
       
        return_atoms.info['stc']= var_id
        return return_atoms