#ifndef _CRYSTAL_H
#define _CRYSTAL_H

typedef struct
{
    float lattice_vectors[3][3];
    float *Xcord;
    float *Ycord;
    float *Zcord;
    char *atoms;
    int spg;
    int wyckoff_position;
    int num_atoms_in_molecule;
    int Z;
    int Zp;

}crystal;

/*
Crystal structure

 -  lattice_vectors is a 3x3 array. Lattice vectors are stored in a
    row major way.
    Eg. a = (2, 4.4, 5), b = (4, 5, 6.7) c = ( 8.1, 9, 10) are stored as
    lattice_vectors[3][3] = {{2.0, 4.4, 5.0}, {4.0, 5.0, 6.7}, {8.1, 9.0, 10.0}}

 -  Xcord, Ycord, Zcord store the positions of the atoms within the
    crystal. Generally, co-rdinates are stored using in units of Angstroms
    and using a cartesian system.
    length of Xcord/Ycord/Zcord = Z * num_atom_in_molecule

    IMPORTANT ASSUMPTION:
    1) Atoms are stored in blocks. Each block represents a molecule.
    2) The order of the atoms within each molecule is same across all
    molecules.
    Eg:

    X   Y   Z    atom
    -----------------

    0.0 2.0 0.0  H  |
    2.0 0.0 0.0  H  |  Molecule 1
    0.0 0.0 0.0  O  |

    7.0 2.0 0.0  H  |
    9.0 0.0 0.0  H  |  Molecule 2
    7.0 0.0 0.0  O  |

    Then,
    Xcord = {0.0, 2.0, 0.0, 7.0, 9.0, 7.0};

 -  atoms stores the information of species at the corresponding site.
    it a char array of 2 x length of coordinate array. If the atom symbol
    occupies only a single character, the second character is a space or ' '.

    Eg:
    Oxygen    - 'O '   // Notice the space after O
    Bromine   - 'Br'

    For the above example, the atoms array should be:
    atoms = {'H', ' ', 'H', ' ', 'O', ' ', 'H', ' ', 'H', ' ', 'O', ' '};

 -  spg stores the space group number according to ITC. If unknown, use 0.

 -  wyckoff position stores the position within the spacegroup according
    to ITC. If unknown, use 0.

 -  num_atoms_in_molecule stores the number of atoms within each molecule.
    The assumption is that the crystal is homomolecular. i.e., not a
    co-crystal.

 -  Z is the number of molecule in the unit cell.

 -  Zp is Z'. Not used. Set it to 0.
*/


#endif
