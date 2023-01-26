#Import these for the GA code to run
from GCGA.GCGA import GCGA
from ase.build import bulk
from ase.build import molecule
from GCGA.FitnessFunction.BaseFitness import GibbsFreeEnergy

from ase import Atoms
###########################
from ase.io import read
#Calculator
from ase.calculators.lammpsrun import LAMMPS
#Define the fitness fucntion for your atoms it must take in an atoms object as parameter and return a float for the code to work"
base_reference = 1.0
a = 24.0


#Indicate static part of the system"
slab = Atoms(cell=[a,a,a],
             pbc=True)

#---------Generate variable part of the system----------------------"
variable_types = [Atoms('Pd'),Atoms('O')]
variable_range = [[3],list(range(7,11))]

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

references = {
    "O": read("references/Oref.traj").get_potential_energy()/2.0,
    }

environment = {"mu_O":ENVIR} #Change envir vaariable for desired mu_O

gfe = GibbsFreeEnergy(variable_reference_energies=references,environmental_variables=environment, base_reference_energy= base_reference)


#Instantiating of the GCGA object with the selected parameters
gcga = GCGA(slab,variable_types,variable_range,
gfe,
calculator = lammps,
starting_candidates_per_stc = 2,population_size=20,steps=1000)

#Calling the run function will initialize the run
gcga.run()

