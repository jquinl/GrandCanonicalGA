import random
import numpy as np
from ase.io import write


from .MutationsBase import MutationsBase
class RattleOperation(MutationsBase):


    def __init__(self, slab,variable_types,variable_range,n_to_move =1,ratio_of_covalent_radii=0.7,
            rng=np.random,strict = False,rattle_strength = 1.0):
        super().__init__(slab,variable_types,variable_range,ratio_of_covalent_radii,rng)

        self.n_to_move = n_to_move
        self.strict = strict
        self.rattle_strength = rattle_strength

    def mutate(self, a1):
        super().mutate(a1)
       

        unique_types =[]
        for a in a1.numbers:
            if(a not in unique_types):
                unique_types.append(a)

        a1_copy = a1.copy()
        a1_copy = a1_copy[len(self.slab) :len(a1_copy)]


        n_moved = self.n_to_move
        if(self.n_to_move > len(a1)):
            print("Rattle operation being performed on all atoms of structure")
            n_moved = len(a1)

        counter = 0
        maxcount = 1000
        
        # Run until a valid pairing is made or maxcount pairings are tested.
        while counter < maxcount:
            counter += 1

            displaced_atoms = []
            while len(displaced_atoms) != n_moved:
                pos = self.rng.randint(0,len(a1))
                if(pos not in displaced_atoms):
                    displaced_atoms.append(pos)

            random.shuffle(displaced_atoms)

            displaced = False
            child = a1_copy.copy()
            for i in displaced_atoms:
                rattled = False
                tries = 0
                while tries < 10 and rattled == False:
                    tries+= 1
                    atoms = self.rattle_operation(child,i)
                    if(atoms is not None):
                        atoverlaps = self._check_index_overlaps(i,atoms,self.blmin)
                        self._displace(i,atoverlaps,atoms)
                        if(len(atoverlaps) == 0): 
                            rattled = True
                            displaced = True
                            child = atoms.copy()
            
            if(displaced == False): continue
            return_atoms = self.slab.copy()

            return_atoms.extend(self.sort_atoms_by_type(child))
            
            if(self._check_overlap_all_atoms(return_atoms,self.blmin)):
                continue
            
            var_id = self.get_var_id(return_atoms)
            if(var_id is not None):
                return_atoms.info['stc']= var_id
                return return_atoms
            else:
                raise Exception("Provided atomic combination is not present in combination matrix")
                

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
            return None
        
        if (len(self.slab )> 0 ):
            slab_cm = self.slab.get_center_of_mass(scaled = True)
            if not (self.rng.rand() > np.dot(np.linalg.norm(slab_cm - cm),np.linalg.norm(new_pos -cm))):
                return None
            at[position].scaled_position = new_pos
            return at

        at[position].scaled_position = new_pos
        return at

    