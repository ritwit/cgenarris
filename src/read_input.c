#include "read_input.h"
#include "molecule_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


float TOL;
void read_control(int* num_structures, int* Z, float* Zp_max,
                 float* volume_mean,float* volume_std,
                 float *sr,long *max_attempts, char *spg_dist_type,
                 int *vol_attempts,int *random_seed,
                 int *crystal_generation,
                 float* interface_area_mean,float* interface_area_std,
                 int* volume_multiplier,
                 float lattice_vector_2d[2][3], float* norm_dev,
                 float* angle_std)
{
	FILE *fileptr;
	size_t len = 0;
	char *line = NULL;
	char *sub_line = NULL;
	int read;
	fileptr = fopen("control.in","r");

	if(!fileptr)
	{
		printf("***ERROR: no control.in file \n");
		exit(EXIT_FAILURE);
	}

	//defaults
	*sr = 0.85;
    *norm_dev = 0.4;
    *angle_std = 8;
	*vol_attempts = 100000;
    *random_seed = 0;
	*volume_multiplier = 3;
    *interface_area_mean = 0;
    *interface_area_std = 0;
    lattice_vector_2d[0][0] = 0;
    lattice_vector_2d[0][1] = 0;
    lattice_vector_2d[0][2] = 0;
    lattice_vector_2d[1][0] = 0;
    lattice_vector_2d[1][2] = 0;
    lattice_vector_2d[1][3] = 0;
    *crystal_generation = 1;

	//read from control
	while ((read = getline(&line, &len, fileptr)) != -1)
	{

	    //if comment
        if (strstr(line, "#") != NULL)
                continue;

	    sub_line=strtok(line," ");
            if(strcmp(sub_line, "Z") == 0)
		{
		    sub_line = strtok(NULL," ");
		    *Z = atoi(sub_line);
		    continue;
		}
	    if (strcmp(sub_line,"lattice_vector_a")==0)
		{

		    sub_line = strtok(NULL,"        ");
		    lattice_vector_2d[0][0] = atof(sub_line);
		    sub_line = strtok(NULL,"        ");
		    lattice_vector_2d[0][1] = atof(sub_line);
		    sub_line = strtok(NULL,"        ");

            lattice_vector_2d[0][2] = atof(sub_line);
                continue;

		}
	    if (strcmp(sub_line,"lattice_vector_b")==0)
		{
		    sub_line = strtok(NULL,"        ");
            lattice_vector_2d[1][0] = atof(sub_line);
            sub_line = strtok(NULL,"        ");
            lattice_vector_2d[1][1] = atof(sub_line);
            sub_line = strtok(NULL,"        ");
            lattice_vector_2d[1][2] = atof(sub_line);
            continue;
		}
            if (strcmp(sub_line,"crystal_generation")==0)
		{
		    sub_line = strtok(NULL," ");
	            *crystal_generation = atof(sub_line);
		    continue;
		}
	    if(strcmp(sub_line, "sr") == 0)
		{
		    sub_line = strtok(NULL," ");
		    *sr = atof(sub_line);
		    continue;
		}

	    if(strcmp(sub_line, "volume_mean") == 0)
		{
		    sub_line = strtok(NULL," ");
		    *volume_mean = atof(sub_line);
		    continue;
		}

	    if(strcmp(sub_line, "volume_std") == 0)
		{
		    sub_line = strtok(NULL," ");
		    *volume_std = atof(sub_line);
		    continue;
		}

	    if(strcmp(sub_line, "interface_area_mean") == 0)
		{
		    sub_line = strtok(NULL," ");
		    *interface_area_mean = atof(sub_line);
		    continue;
		}

	    if(strcmp(sub_line, "interface_area_std") == 0)
		{
		    sub_line = strtok(NULL," ");
		    *interface_area_std = atof(sub_line);
		    continue;
		}

	    if(strcmp(sub_line, "volume_multiplier") == 0)
        {
            sub_line = strtok(NULL," ");
            *volume_multiplier = atof(sub_line);
            continue;
        }
	    if(strcmp(sub_line, "number_of_structures") == 0)
		{
		    sub_line = strtok(NULL," ");
		    *num_structures = atoi(sub_line);
		    continue;
		}

	    if(strcmp(sub_line, "tolerance") == 0)
		{
		    sub_line = strtok(NULL," ");
		    TOL = atof(sub_line);
		    continue;
		}

	    if(strcmp(sub_line, "max_attempts") == 0)
		{
		    sub_line = strtok(NULL," ");
		    *max_attempts = atol(sub_line);
		    continue;
		}

        if(strcmp(sub_line, "lattice_angle_std") == 0)
        {
            sub_line = strtok(NULL," ");
            *angle_std = atof(sub_line);
            continue;
        }

        if(strcmp(sub_line, "lattice_norm_dev") == 0)
        {
            sub_line = strtok(NULL," ");
            *norm_dev = atof(sub_line);
            continue;
        }

	    if(strcmp(sub_line, "spg_distribution_type") == 0)
		{
		    sub_line = strtok(NULL," ");
		    sub_line = strtok(sub_line, "\n");
		    strcpy(spg_dist_type, sub_line);
		    if(!( strcmp(spg_dist_type, "standard") ||
			      strcmp(spg_dist_type, "uniform")  ||
			      strcmp(spg_dist_type, "chiral")   ||
			      strcmp(spg_dist_type, "csd")       )
               )
			{
				printf("***ERROR: read_input: bad value of spg_distribution_type %s", spg_dist_type);
				exit(EXIT_FAILURE);
			}
			continue;
		}

        if(strcmp(sub_line, "random_seed") == 0)
        {
        	sub_line = strtok(NULL," ");
        	*random_seed = atoi(sub_line);
        	continue;
       	}

        if(strcmp(sub_line, "volume_attempts") == 0)
        {
           	sub_line = strtok(NULL," ");
            int inp = atoi(sub_line);
            if(inp != 0)
                *vol_attempts = inp;
            continue;
        }

   	}

	/*
    *num_structures = 1;
    *volume_mean = 1200;
    *volume_std  = 100;
	*Z = 6;
	*sr = 0.75;
	*/
	fclose(fileptr);
    *Zp_max = 192;
}


void read_geometry(molecule* mol)
{
    FILE *fileptr;
    size_t len = 0;
    int read;
    char *line = NULL;
    char *sub_line = NULL;
    int i = 0;
    int atom_count = 0;

    //find_number_of atoms
    fileptr = fopen("geometry.in","r");
    //check if file exits
    if(!fileptr)
    {
        printf("***ERROR: no geometry.in file \n");
        exit(EXIT_FAILURE);
    }

     while ((read = getline(&line, &len, fileptr)) != -1)
    {
        if (strstr(line, "#") != NULL)
            continue;

        sub_line=strtok(line," ");
        //printf("%s \n" , sub_line);
        if(strcmp(sub_line, "atom") == 0)
            atom_count++;
        else
            continue;
    }
    fclose(fileptr);

    //printf("Total number of atoms in molecule = %d\n", atom_count);
    int N = atom_count;
    //memory allocation
    (*mol).atoms = (char *)malloc(2*N*sizeof(char));
    (*mol).X = (float *)malloc(N*sizeof(float));
    (*mol).Y = (float *)malloc(N*sizeof(float));
    (*mol).Z = (float *)malloc(N*sizeof(float));

    //creates a file pointer and opens geometry.in
    fileptr = fopen("geometry.in","r");
    while ((read = getline(&line, &len, fileptr)) != -1)
    {
        //printf("%s", line);
        if (strstr(line, "#") != NULL)
            continue;

        sub_line=strtok(line," ");
        if(strcmp(sub_line, "atom") == 0)
            atom_count++;
        else
            continue;
        sub_line=strtok(NULL," ");
        (*mol).X[i]=atof(sub_line);
        sub_line=strtok(NULL," ");
        (*mol).Y[i]=atof(sub_line);
        sub_line=strtok(NULL," ");
        (*mol).Z[i]=atof(sub_line);
        sub_line=strtok(NULL," ");
        (*mol).atoms[2*i]=*sub_line;
        if(*(sub_line+1) == '\n' || *(sub_line+1) == ' ' || *(sub_line+1) == '\0' )
            (*mol).atoms[2*i+1]=' ';
        else
            (*mol).atoms[2*i+1]=*(sub_line+1);
        i++;

    }
    fclose(fileptr);
   // printf("atoms = %d \n", atom_count);
    mol->num_of_atoms = N;


}

//added here, read a geometty.in for crystal, convert to xtal obj
void read_crystal_geometry(crystal* xtal, char* geo_name)
{
    FILE *fileptr;
    size_t len = 0;
    int read;
    char *line = NULL;
    char *sub_line = NULL;
    int i = 0;
    int lv_count = 0;

    //find_number_of atoms
    fileptr = fopen(geo_name,"r");
    //check if file exits
    if(!fileptr)
    {
        printf("***ERROR: no geometry.in file \n");
        exit(EXIT_FAILURE);
    }

     while ((read = getline(&line, &len, fileptr)) != -1)
    {
        if (strstr(line, "#") != NULL)
            continue;

        sub_line=strtok(line," ");
		if (strcmp(sub_line,"lattice_vector")==0)
		{
			sub_line = strtok(NULL," ");
		    xtal->lattice_vectors[lv_count][0] = atof(sub_line);
			printf("sub_line: %s",sub_line);
		    fflush(stdout);
		    sub_line = strtok(NULL," ");
		    printf("sub_line: %s",sub_line);
		    fflush(stdout);
		    xtal->lattice_vectors[lv_count][1] = atof(sub_line);
		    sub_line = strtok(NULL," ");
		    printf("sub_line: %s",sub_line);
		    fflush(stdout);
            xtal->lattice_vectors[lv_count][2] = atof(sub_line);
			lv_count ++;
            continue;
		}
        //printf("%s \n" , sub_line);
        if(strcmp(sub_line, "atom") == 0)
		{
			sub_line=strtok(NULL," ");
        	xtal->Xcord[i]=atof(sub_line);
        	printf("%f \t",xtal->Xcord[i]);
        	sub_line=strtok(NULL," ");
        	xtal->Ycord[i]=atof(sub_line);
        	printf("%f \t",xtal->Ycord[i]);
        	sub_line=strtok(NULL," ");
        	xtal->Zcord[i]=atof(sub_line);
        	printf("%f \t",xtal->Zcord[i]);
        	sub_line=strtok(NULL," ");
        	xtal->atoms[2*i]=*sub_line;
        	if(*(sub_line+1) == '\n' || *(sub_line+1) == ' ' || *(sub_line+1) == '\0' )
            	xtal->atoms[2*i+1]=' ';
        	else
            	xtal->atoms[2*i+1]=*(sub_line+1);
        	i++;
		}
       
    }
    fclose(fileptr);


}


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
                          float *angle_std)
{
    *Zp_max = 192; //useless argument
    printf("INPUT SETTINGS:\n");
    printf("-----------------------------\n");
    printf("Number of structures per space group:         %d \n", *num_structures);
    printf("Number of molecules in the cell:              %d\n", *Z);
    printf("Mean volume of unit cell:                     %f\n", *volume_mean);
    printf("Standard deviation of unit cell volume:       %f\n", *volume_std);
    printf("Spacegroup distribution type:                 %s\n", spg_dist_type);
    //printf("Specific radius proportion:                   %f\n", *sr );
    printf("Maximum attempts per space group:             %ld\n", *max_attempts);
    printf("Volume attempts:                              %d\n", *vol_attempt);
    printf("Random seed:                                  %d\n", *random_seed);
    printf("Lattice angle standard deviation:             %f\n", *angle_std);
    printf("Lattice principal component deviation:        %f\n", *norm_dev);
    printf("Tolerance:                                    %f\n", TOL);
    printf("-----------------------------\n\n");

}

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
                          	int *random_seed)
{
    *Zp_max = 192; //useless argument
    printf("INPUT SETTINGS:\n");
    printf("-----------------------------\n");
    printf("Number of structures per layer group:         %d \n", *num_structures);
    printf("Number of molecules in the cell:              %d\n", *Z);
    printf("Mean volume of unit cell:                     %f\n", *volume_mean);
    printf("Standard deviation of unit cell volume:       %f\n", *volume_std);
    printf("Mean interface area:                          %f\n", *interface_area_mean);
    printf("Standard deviation of interface area:         %f\n", *interface_area_std);
    printf("Volume Multiplier:                            %d\n", *volume_multiplier);
    printf("Layergroup distribution type:                 %s\n", spg_dist_type);
    printf("Lattice vector:                               [[%f,%f,%f],[%f,%f,%f]]\n",global_lattice_vector_2d[0][0],
                                                global_lattice_vector_2d[0][1],global_lattice_vector_2d[0][2],
                                                global_lattice_vector_2d[1][0],global_lattice_vector_2d[1][1],
                                                global_lattice_vector_2d[1][2]  );
    printf("Maximum attempts per layer group:             %ld\n", *max_attempts);
    printf("Volume attempts:                              %d\n", *vol_attempt);
    printf("Random seed:                                  %d\n", *random_seed);
    printf("Tolerance:                                    %f\n", TOL);
    printf("-----------------------------\n\n");

}


void print_input_geometry(molecule* mol)
{
    printf("MOLECULAR GEOMETRY:\n");
    printf("-----------------------------\n");
    print_molecule(mol);
    printf("-----------------------------\n");
    printf("\n");
}


