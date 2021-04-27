#ifndef _COCRYSTAL_H_
#define _COCRYSTAL_H_

typedef struct
{
    float lattice_vectors[3][3];
    float *Xcord;
    float *Ycord;
    float *Zcord;
    char *atoms;
    int *num_atoms_in_molecule;
    int spg;
    int *wyckoff_position;
    int Z;
    int Zp;
}cocrystal;


#endif  // Cocrystal.h
