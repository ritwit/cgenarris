#ifndef READ_INPUT_H
#define READ_INPUT_H

#include <stdio.h>
#include <stdlib.h>

#include "molecule.h"
#include "crystal.h"

void read_control(int* num_structures, int* Z, float* Zp_max,
                 float* volume_mean,float* volume_std,
                 float *sr,long *max_attempts, char *spg_dist_type,
                 int *vol_attempts,int *random_seed,
                 int *crystal_generation,
                 float* interface_area_mean,float* interface_area_std,
                 int* volume_multiplier,
                 float lattice_vector_2d[2][3], float* norm_dev,
                 float* angle_std);

void read_geometry(molecule* mol);

void read_crystal_geometry(crystal* xtal, char* geo_name);

void print_input_geometry(molecule* mol);

void print_molecule(molecule *mol);

void print_input_settings(int* num_structures,
                          int* Z,
                          float* Zp_max,
                          float* volume_mean,
                          float* volume_std,
                          float *sr,
                          long *max_attempts,
                          char * spg_dist_type,
                          int *vol_attempt,
                          int *random_seed,
                          float *norm_dev,
                          float *angle_std);

void print_input_settings_layer(int* num_structures,
                          int* Z,
                          float* Zp_max,
                          float* volume_mean,
                          float* volume_std,
                          float* interface_area_mean,
                          float* interface_area_std,
                          int* volume_multiplier,
                          float *sr,
                          float global_lattice_vector_2d[2][3],
                          long *max_attempts,
                          char * spg_dist_type,
                          int *vol_attempt,
                          int *random_seed);

#endif
