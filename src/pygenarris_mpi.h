#ifndef PYGENARRIS_MPI_H
#define PYGENARRIS_MPI_H

#include "mpi.h"
#include "crystal.h"

void print_time(void);

void mpi_generate_molecular_crystals_with_vdw_cutoff_matrix(
    float *vdw_matrix,
    int dim1,
    int dim2,
    int num_structures,
    int Z,
    double volume_mean1,
    double volume_std1,
    double tol1,
    long max_attempts,
    char *spg_dist_type,
    int vol_attempt,
    int random_seed,
    float norm_dev,
    float angle_std,
    MPI_Comm world_comm);

void mpi_generate_layer_with_vdw_cutoff_matrix(
	float *vdw_matrix,
	int dim1,
	int dim2,
	int num_structures,
	int Z,
	double volume_mean1,
	double volume_std1,
	double interface_area_mean1,
	double interface_area_std1,
	int volume_multiplier,
	double tol1,
	long max_attempts,
	char *spg_dist_type,
	float lattice_vector_2d_from_geo[2][3],
        int vol_attempt,
        int random_seed,
	MPI_Comm world_comm);

void send_xtal(MPI_Comm comm, int destination, crystal* xtal, int total_atoms);
void receive_xtal(MPI_Comm comm, int source, crystal* xtal, int total_atoms);

int num_compatible_spacegroups(int Z, double tolerance);
int num_compatible_layergroups(int Z, double tolerance,float volume,float lattice_vector_2d_from_geo[2][3]);

void generate_molecular_crystals(char *filename,int num_structures, int Z,
    double volume_mean1, double volume_std1, double sr1, double tol1,
    int max_attempts);

int c_check_structure(crystal xtal, double sr);

void create_crystal_from_array(crystal *xtal, double lattice_vector[3][3],
    double *Xc,int total_atoms1, double *Yc,int total_atoms2,
    double *Zc, int total_atoms3, char* atoms, int total_atoms,
    int Z, int spg);

#endif
