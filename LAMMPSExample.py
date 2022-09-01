#Import these for the GA code to run

from GCGA.GCGA import GCGA
from GCGA.FitnessFunction.BaseFitness import GibbsFreeEnergy

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


#Define calculator for LAMMPS
header = [
    "units		real",
    "atom_style	charge"]

masses = {'Zn':65.3900,'O':19.9990,'H':1.0080}
cmds = [
    "mass 1 65.3900",
    "mass 2 15.9990",
    "mass 3 1.0080",
    "pair_style	reax/c lmp_control",
    "pair_coeff	* * ffield.reax.ZnOH H O Zn",
    "neighbor	2 bin",
    "neigh_modify	every 10 delay 0 check no",
    "fix		1 all nve",
    "fix             2 all qeq/reax 1 0.0 10.0 1e-6 param.qeq",
    "fix             3 all temp/berendsen 500.0 500.0 100.0",
    "timestep	0.25",
    "run		3000"]
        

lammps = LAMMPSlib(lammps_header = header,lmpcmds=cmds, log_file='test.log')

#Instantiating of the GCGA object with the selected parameters
gcga = GCGA(slab,variable_types,variable_range,
gfe,
calculator = lammps,
starting_candidates_per_stc = 2,population_size=20,steps=100)

#Calling the run function will initialize the run
gcga.run()

