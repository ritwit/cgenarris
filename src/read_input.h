#ifndef READ_INPUT_H
#define READ_INPUT_H

#include <stdio.h>
#include <stdlib.h>

#include "molecule.h"

void read_control(int* num_structures,
                  int* Z,
                  float* Zp_max,
                  float* volume_mean,
                  float* volume_std,
                  float *sr,
                  long *max_attempts,
                  char *spg_dist_type,
                  int *vol_attempt,
                  int *random_seed);

void read_control_layer(int* num_structures, int* Z, float* Zp_max, float* volume_mean,
				 float* volume_std,float* interface_area_mean,
				 float* interface_area_std,int* volume_multiplier, 
                 float *sr, long *max_attempts, char *spg_dist_type, 
                 float lattice_vector_2d[2][3])

void read_geometry(molecule* mol);

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
                          int *random_seed);

#endif
