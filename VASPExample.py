#Import these for the GA code to run
from GCGA.GCGA import GCGA
from ase.build import bulk
from ase.build import molecule
from GCGA.FitnessFunction.BaseFitness import GibbsFreeEnergy

from ase import Atoms
###########################
from ase.io import read
#Calculator
from ase.calculators.vasp import Vasp
#Define the fitness fucntion for your atoms it must take in an atoms object as parameter and return a float for the code to work"

calc = Vasp(
	prec='Low',
	lreal = True,
        xc='PW91',
        nwrite = 2,
        istart= 1,
        icharg = 1,
        ispin = 1,
        voskown = 1,
        lmaxmix = 6,

        nelm = 80,
        nelmin = 4,
        ediff = 0.1E-03,

        ediffg = -0.1,
        ibrion = 2,
        potim = 0.3,
        iwavpr = 1,
        nsw = 200,

        ialgo = 38,

        ismear= 0,
        sigma= 0.1,
        lorbit= 11,
        )

a = 24.0

#Indicate static part of the system"
slab = Atoms(cell=[a,a,a],
             pbc=True)

#---------Generate variable part of the system----------------------"
variable_types = [Atoms('Pt'),Atoms('O')]
variable_range = [[6],list(range(0,12))]



o_ref = read("references/Oref.traj")
o_ref.calc = calc


references = {
    "O": o_ref.get_potential_energy()/2.0,
    }

environment = {"mu_O":ENVIR} #Change envir vaariable for desired mu_O

gfe = GibbsFreeEnergy(variable_reference_energies=references,environmental_variables=environment)


#Instantiating of the GCGA object with the selected parameters
gcga = GCGA(slab,variable_types,variable_range,
gfe,
calculator = calc,
starting_candidates_per_stc = 2,population_size=20,steps=1000)

#Calling the run function will initialize the run
gcga.run()

