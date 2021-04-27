#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "crystal.h"
#include "crystal_utils.h"
#include "molecule.h"
#include "read_input.h"
#include "pygenarris_mpi.h"

enum status{FAIL, SUCCESS};

// Declarations
int regression_test(MPI_Comm world_comm);
void create_vdw_matrix_from_sr(molecule *mol,
                                float *vdw_matrix,
                                float sr,
                                int Z);
int compare_geometry_out();
char* get_metadata_line(char* line);
int read_geometry_out(FILE *fptr, crystal *xtal);

int regression_test(MPI_Comm world_comm)
{
    printf("Running regression test...\n");
    // settings
    molecule mol;
    float sr = 0.85;
    double volume_mean = 600;
    double volume_std = 40;
    double tol = 0.1;
    int num_structures = 2;
    int Z = 2;
    int seed = 223344;
    int volume_attempts = 10000;
    long max_attempts = 100000;
    char spg_dist_type[10] = "uniform";
    float norm_dev = 0.4;
    float angle_std = 8;

    // Read molecule
    read_geometry(&mol);

    // Create vdw distance cutoff matrix
    int num_atoms_in_molecule = mol.num_of_atoms;
    int dim_vdw_matrix = num_atoms_in_molecule * Z ;
    float *vdw_cutoff_matrix = (float *) malloc( dim_vdw_matrix *
                                dim_vdw_matrix *
                                sizeof(float) ); //square matrix
    create_vdw_matrix_from_sr(&mol, vdw_cutoff_matrix, sr, Z);

    //call the generator from pygenarris_mpi
    mpi_generate_molecular_crystals_with_vdw_cutoff_matrix(
        vdw_cutoff_matrix,
        dim_vdw_matrix,
        dim_vdw_matrix,
        num_structures,
        Z,
        volume_mean,
        volume_std,
        tol,
        max_attempts,
        spg_dist_type,
        volume_attempts,
        seed,
        norm_dev,
        angle_std,
        world_comm);

    free(vdw_cutoff_matrix);
    int status = compare_geometry_out();
    return status;
}

int compare_geometry_out()
{
    crystal xtal, xtal_ref;
    FILE *fptr, *fptr_ref;
    fptr = fopen("geometry.out", "r");
    fptr_ref = fopen("geometry.out.test", "r");

    while(read_geometry_out(fptr, &xtal))
    {
        read_geometry_out(fptr_ref, &xtal_ref);
        int result = is_equal_xtal(&xtal, &xtal_ref, 0.001);
        if(result != 1)
        {
            printf("Test failed!\n");
            return FAIL;
        }
    }

    fclose(fptr);
    fclose(fptr_ref);
    return SUCCESS;
}

int main(int argc, char **argv)
{
    //Initialise MPI
    MPI_Init(&argc, &argv);
    MPI_Comm world_comm = MPI_COMM_WORLD;

    int status;
    status = regression_test(world_comm);
    if(status == SUCCESS)
    {
        printf("Regression tests passed\n");
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
    else if(status == FAIL)
    {
        printf("Regression test failed\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    MPI_Finalize();
    exit(EXIT_FAILURE);
}


int read_geometry_out(FILE* fptr, crystal *xtal)
{
  char *line = NULL;
  char *sub_line = NULL;
  size_t len = 0;
  int num_atom = 0;

  int spg = 0;
  int Z = 0;
  int wyckoff_position = 0;
  int Zp = 0;
  int i = 0;
  int counter = 0;
  int num_atoms_in_molecule = 0;

  if(!fptr)
  {
       printf("***ERROR: Cannot open the file \n");
       exit(EXIT_FAILURE);
  }

    while(getline(&line,&len,fptr) != -1)
    {
        // Starts with # then
        if(strstr(line,"#") != NULL || *line== '\n')
        {
            if(strstr(line, "attempted_spacegroup") != NULL)
            {
                 sub_line = get_metadata_line(line);
                 spg = atoi(sub_line);
                 // Set temporary counting index for new xtal structure
                 i = 0;
                 counter = 0;
            }

            else if (strstr(line, "Z") != NULL)
            {
                sub_line = get_metadata_line(line);
                Z = atoi(sub_line);
            }

            else if (strstr(line,"number_of_atoms_in_molecule") != NULL)
            {
                sub_line = get_metadata_line(line);
                num_atoms_in_molecule = atof(sub_line);
            }

            else if (strstr(line,"END") != NULL)
            {
                 xtal->spg = spg;
                 xtal->wyckoff_position = wyckoff_position;
                 xtal->num_atoms_in_molecule = num_atoms_in_molecule;
                 xtal->Z = Z;
                 xtal->Zp = Zp;
                 return 1;
            }

            continue;
        }

        sub_line = strtok(line, "    ");
        if(strcmp(sub_line, "lattice_vector") == 0)
        {
            sub_line = strtok(NULL, "    ");
            xtal->lattice_vectors[i][0] = atof(sub_line);
            sub_line = strtok(NULL, "    ");
            xtal->lattice_vectors[i][1] = atof(sub_line);
            sub_line = strtok(NULL, "    ");
            xtal->lattice_vectors[i][2] = atof(sub_line);
            if(!i)
            {
                num_atom = num_atoms_in_molecule *Z;
                xtal->Xcord = (float *) malloc (num_atom*sizeof(float));
                xtal->Ycord = (float *) malloc (num_atom*sizeof(float));
                xtal->Zcord = (float *) malloc (num_atom*sizeof(float));
                xtal->atoms = (char * ) malloc (2*num_atom*sizeof(char));
             }
            i++;
        }

        if(strcmp(sub_line, "atom") == 0)
        {
            xtal->Xcord[counter] = atof(strtok(NULL, "    "));
            xtal->Ycord[counter] = atof(strtok(NULL, "    "));
            xtal->Zcord[counter] = atof(strtok(NULL, "    "));

            sub_line = strtok(NULL, " ");
            xtal->atoms[2*counter] = *sub_line;
            if(*(sub_line+1) == '\n' || *(sub_line+1) == ' ' || *(sub_line+1) == '\0' )
                xtal->atoms[2*counter+1] = ' ';
            else
                xtal->atoms[2*counter+1] = *(sub_line+1);
            counter++;
        }
    }

    return 0;
}

// extract the meta data line like Z, spg etc
char* get_metadata_line(char* line)
{
    char* sub_line = NULL;
    sub_line = strtok(line," ");
    sub_line = strtok(NULL," ");
    sub_line = strtok(NULL," ");
    return sub_line;
}