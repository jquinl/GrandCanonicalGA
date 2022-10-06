#Import these for the GA code to run

from GCGA.GCGA import GCGA
from GCGA.FitnessFunction.BaseFitness import GibbsFreeEnergy
from ase.calculators.lammpsrun import LAMMPS

from ase import Atoms
import numpy as np
###########################
from ase.io import read
#Calculator
from ase.calculators.lammpslib import LAMMPSlib
#Define the fitness fucntion for your atoms it must take in an atoms object as parameter and return a float for the code to work"
base_reference = 1.0
references = {"Zn": 1.0,"O": 2.0,"H": 3.0}
environment = {"mu_Zn": 1.0,"mu_O":5.0,"mu_H":2.0}
gfe = GibbsFreeEnergy(variable_reference_energies=references,environmental_variables=environment, base_reference_energy= base_reference)
#Indicate static part of the system"
a = 24.0
slab = Atoms(cell=[a,a,a],
             pbc=True)

#---------Generate variable part of the system----------------------"
variable_types = [Atoms('Zn'),Atoms('O'),Atoms('H')]
variable_range = [[1],[1,2,3,4],[1,2]]

#Single point evaluation of the sistem, followed by an optimization with BFGS (Provisional)
params={"clear":"",
    "units":"real",
    "atom_style":"charge",
    "pair_style": "reaxff NULL",
    "pair_coeff": ["* * ffield.reax.PdO O Pd"],
    "fix":["1 all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff"]
}
files = ["ffield.reax.PdO"]
lammps= LAMMPS(parameters=params, files=files,keep_alive=False)

#Instantiating of the GCGA object with the selected parameters
gcga = GCGA(slab,variable_types,variable_range,
gfe,
calculator = lammps,
starting_candidates_per_stc = 2,population_size=20,steps=100)

#Calling the run function will initialize the run
gcga.run()

