#ifndef _COCRYSTAL_UTILS_H_
#define _COCRYSTAL_UTILS_H_

#include "cocrystal.h"
#include "input_settings.h"

void cxtal_init(cocrystal *cxtal, int *stoic, int *n_atoms_in_mol, int n_mol_types, int Z);
void cxtal_allocate(cocrystal *cxtal, int total_atoms);
void cxtal_print(cocrystal *cxtal, FILE* out, int fractional);
int cxtal_check_structure(cocrystal *cxtal, Settings *set);
float cxtal_get_cell_volume(cocrystal *cxtal);
int cxtal_check_structure(cocrystal *cxtal, Settings *set);
void cxtal_bring_molecules_first_cell(cocrystal *cxtal);

#endif