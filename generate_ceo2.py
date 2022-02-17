def create_ceo2_111(number_of_layers = 3,
                    lattice_parameter = 5.399,
                    vacuum = 10., 
                    repetitions = (1,1,1),
                    ):
    from ase import Atoms, Atom
    from ase.build import surface
    from ase import geometry



    number_of_layers = number_of_layers # this corresponds to the number of trilayers (Ce layers)

    a = lattice_parameter # lattice parameter

    #Coordinates of the unit-cell atoms in fractional 
    #coordinates (the fractional coordinates can be 
    #found in the link given in external link section)
    CeO2 = Atoms([
                  Atom('Ce', (0., 0., 0.)),
                  Atom('Ce', (0., 0.5, 0.5)),
                  Atom('Ce', (0.5, 0., 0.5)),
                  Atom('Ce', (0.5, 0.5, 0.)),
                  Atom('O', (0.75, 0.25, 0.25)),
                  Atom('O', (0.25, 0.75, 0.75)),
                  Atom('O', (0.75, 0.75, 0.75)),
                  Atom('O', (0.25, 0.25, 0.25)),
                  Atom('O', (0.25, 0.25, 0.75)),
                  Atom('O', (0.75, 0.75, 0.25)),
                  Atom('O', (0.25, 0.75, 0.25)),
                  Atom('O', (0.75, 0.25, 0.75))
                                               ])

    #Defining lattice size in a, b and c directions. 
    #In this case, they are all equal.
    cell = [(a, 0., 0.),
            (0., a, 0.),
            (0., 0., a)]

    #Scales the atomic positions with the unit cell
    CeO2.set_cell(cell, scale_atoms=True)
    #(1,1,1) is the slab type. There are 2 unit cells along 
    #z direction
    slab = surface(CeO2, (1, 1, 1), number_of_layers+1)
    cell = slab.get_cell()
    cell[0]*= 1/2
    cell[1]*= 1/2
    slab.set_cell(cell)
    slab.wrap()
    geometry.get_duplicate_atoms(slab, cutoff=0.1, delete=True)
    #view(slab)
    del slab[[0]]
    del slab[[0]]
    del slab[[-1]]
    #view(slab)
    slab = slab.repeat(repetitions)
    slab.center(vacuum=vacuum/2., axis=2)
    return slab
