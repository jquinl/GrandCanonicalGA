import numpy as np
from ase.ga.utilities import atoms_too_close

from .OperationsBase import OperationsBase
class RattleOperation(OperationsBase):
    """Rattle operation
    
        If Strict is true, will only return an atoms object if the number of atoms indicated in n_to_move have been succesfully rattled
        Else will return an atoms object if at least one of the atoms has been rattled
    """
    def __init__(self, slab,variable_types,variable_range,n_to_move =1,ratio_of_covalent_radii=0.7,
            rng=np.random,strict = False,rattle_strength = 1.0):
        super().__init__(slab,variable_types,variable_range,ratio_of_covalent_radii,rng)

        self.n_to_move = n_to_move
        self.strict = strict
        self.rattle_strength = rattle_strength
    

    def mutate(self):
        return super().mutate()
    def mutate(self, a1, a2):
        super().mutate(a1,a2)

        if(self.slab.get_cell().all() != a1.get_cell().all()):
            raise ValueError('Different cell sizes found for slab and inputed structures')

        # Only consider the atoms to optimize
        a1 = a1[len(self.slab) :len(a1)]

        return_atoms = a1
        
        maxtries = 100
        n_moved = self.n_to_move
        if(self.n_to_move > len(a1)):
            print("Rattle operation being performed on all atoms of structure")
            n_moved = len(a1)

        displaced_atoms = []

        while len(displaced_atoms) != n_moved:
            pos = self.rng.randint(0,len(a1)-1)
            if(pos not in displaced_atoms):
                displaced_atoms.append(pos)
        success = 0 

        for i in displaced_atoms:
            atoms = None
            count = 0
            while atoms == None and count < maxtries: 
                count += 1
                atoms = self.rattle_operation(return_atoms, i)

                ats = self.slab.copy()

                atoms.extend(self.slab.copy())

                if atoms_too_close(ats, self.blmin):
                    atoms = None

            if atoms != None:
                success += 1
                return_atoms = atoms.copy()

        if(self.strict):
            if(success == n_moved):
                return return_atoms,1
            return None,1
        else:
            if(success == 0):
                return None,1
            return return_atoms,1
       
    def rattle_operation(self,atoms, position):

        at = atoms.copy()
        
        if position > len(atoms): raise ValueError("Rattle position higher than number of atoms")

        cm = at.get_center_of_mass(scaled = True)
        at_pos = at[position].scaled_position
    
        x = (self.rattle_strength * self.rng.rand()+at_pos[0]) % 1.0
        y = (self.rattle_strength * self.rng.rand()+at_pos[1]) % 1.0
        z = (self.rattle_strength * self.rng.rand()+at_pos[2]) % 1.0

        new_pos = np.array([x,y,z])
        if( -np.dot(np.linalg.norm(at_pos - cm ),np.linalg.norm(new_pos-at_pos)) > self.rng.rand() ):
            print("None2")
            return None
        
        if (len(self.slab )> 0 ):
            slab_cm = self.slab.get_center_of_mass(scaled = True)
            if not (self.rng.rand() > np.dot(np.linalg.norm(slab_cm - cm),np.linalg.norm(new_pos -cm))):
                print("None")
                return None
            at[position].scaled_position = new_pos
            return at

        at[position].scaled_position = new_pos
        return at

        