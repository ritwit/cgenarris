#include <stdio.h>
#include <math.h>
#include <string.h>
#include "read_input.h"
#include "spg_generation.h"
#include "crystal_utils.h"
#include "crystal.h"
#include "molecule.h"
#include "combinatorics.h"
#include "molecule_placement.h"
#include "molecule_utils.h"
#include "layer_group_position_database.h"
#include "algebra.h"
#include "spglib.h"
#include "check_structure_layer_group.h"

int check_layer_group(crystal* xtal)
{

    int attempted_lg = xtal->spg;
    printf("The attempted layer group is %d\n",attempted_lg);
    int spglib_spg = detect_spg_using_spglib(xtal);
    printf( "#SPGLIB_detected_spacegroup = %d\n", spglib_spg);
    int spglib_lg[6]; //max 5 from detected and 1 of its original
    int match = 0;
    for (int i =0;i<5;i++)
    {
        spglib_lg[i] = spg_to_lg[spglib_spg-1][i];
        if (spglib_lg[i] == attempted_lg)
        {
            match = 1;
	    //printf("match = 1\n");
        }
        printf("the lg is %d\n",spglib_lg[i]);
    }

    printf("after check all lg,match = %d\n",match);

    //if result in a higher symmetry
    if (match == 0)
    {
	  //add the attemtped lg to test its operations
	  spglib_lg[5] = attempted_lg;
	  // if only one lg accosiated with the spg,just return the lg
	  //if (spglib_lg[1] == 0)  
	  //{
		  //return spglib_lg[0];
	  //}
	  for (int j = 0; j<6;j++)
	    {
		  int lg = spglib_lg[j];
		  
		  if (lg != 0)
		    {
				printf("lg is %d\n",lg);
			    float Z = xtal->Z;
    			float N = xtal->num_atoms_in_molecule;
		        crystal *tmp_xtal = (crystal *)malloc(sizeof(crystal) );
    			allocate_xtal(tmp_xtal,Z,N);
			    copy_xtal(tmp_xtal, xtal);
		        check_symmet_equiv_atoms(lg,tmp_xtal,Z,N);	
		


				free_xtal(tmp_xtal);		
		    }

	    }
		  return -1;

    }
	else
	{
		return attempted_lg;
	}

}

void check_symmet_equiv_atoms(int lg,crystal* tmp_xtal,int Z, int N)
{    
    
    float lattice_vec_a[3];
    float lattice_vec_b[3];
    float lattice_vec_c[3];
    float norm_a = 0;
    float norm_b = 0;
    float norm_c = 0;

	//convert coords to frac
	convert_xtal_to_fractional(tmp_xtal);


    int hall_number = lg;
    double translations[192][3];
    int rotations[192][3][3];
    int num_of_operations = get_lg_symmetry(hall_number,translations,rotations);
    //get lattice_vector inverse and transpose
    float inverse_lattice_vectors[3][3];
    inverse_mat3b3(inverse_lattice_vectors, tmp_xtal->lattice_vectors);
    float lattice_vectors_transpose[3][3];
    copy_mat3b3_mat3b3(lattice_vectors_transpose, tmp_xtal->lattice_vectors);
    mat3b3_transpose(lattice_vectors_transpose, lattice_vectors_transpose);
    for(int i = 0; i < num_of_operations; i++)
	{
		printf("this is lg %d, operation %d\n",lg,i+1);
		FILE *outfile;
		char fname[80];
		for (int j =0 ; j<80;j++)
		{
			fname[j] = '';
		}
		char lg_char[5];
		char i_char[5];
		sprintf(lg_char, "%d", lg);
		sprintf(i_char, "%d", i);
		strcat(fname,"lg_");
		strcat(fname,lg_char);
		strcat(fname,"_operation_");
		strcat(fname,i_char);
		strcat(fname,".in");
		outfile = fopen(fname,"w");
	    //rot and trans are ith symmetry operation
	    int rot[3][3];
	    copy_intmat3b3_intmat3b3bN(rot, rotations, i);
	    float trans[3];
	    trans[0] = translations[i][0];
	    trans[1] = translations[i][1];
	    trans[2] = translations[i][2];

		//create crystal for each operation
		crystal *tmp_xtal_each_op = (crystal *)malloc(sizeof(crystal) );
    	allocate_xtal(tmp_xtal_each_op,Z,N);
		copy_xtal(tmp_xtal_each_op, tmp_xtal);
		tmp_xtal_each_op->wyckoff_position = 3;


	    // loop over all atoms in crystal
	    for (int j = 0; j < N * Z; j++)
		{
		   float atomj_array[3];
		   atomj_array[0] = tmp_xtal->Xcord[j];
		   atomj_array[1] = tmp_xtal->Ycord[j];
		   atomj_array[2] = tmp_xtal->Zcord[j];
		   //apply symmetry operations
		   vector3_intmat3b3_multiply(rot, atomj_array, atomj_array);
		   vector3_add(trans, atomj_array, atomj_array);
		   //printf("after operat, x_frac is%f\n",atomj_array[0]);
		   //printf("after operat, y_frac is%f\n",atomj_array[1]);
		   //printf("after operat, z_frac is%f\n",atomj_array[2]);

		   //convert back to cartesian
		   vector3_mat3b3_multiply(lattice_vectors_transpose,
				atomj_array, atomj_array);
		   tmp_xtal_each_op->Xcord[j] = atomj_array[0];
		   tmp_xtal_each_op->Ycord[j] = atomj_array[1];
		   tmp_xtal_each_op->Zcord[j] = atomj_array[2];
		   

		   //printf("new_Xcord[j] is :%f\n",new_Xcord[j]);
		   //printf("new_Ycord[j] is :%f\n",new_Ycord[j]);
		   //printf("new_Zcord[j] is :%f\n",new_Zcord[j]);

		}
		
		bring_all_molecules_to_first_cell(tmp_xtal_each_op);
		// check after bring to first cell
		/*
		for (int j = 0;j< N*Z;j++)
		{
			printf("after bring first Xcord[j] %f\n",tmp_xtal_each_op->Xcord[j]);
			printf("after bring first Ycord[j] %f\n",tmp_xtal_each_op->Ycord[j]);
			printf("after bring first Zcord[j] %f\n",tmp_xtal_each_op->Zcord[j]);
		}
		*/
		print_layer2file(tmp_xtal_each_op,outfile);
		
	    //compare_all_atoms_distance(tmp_xtal->Xcord,tmp_xtal->Ycord,tmp_xtal->Zcord,
													//new_Xcord,new_Ycord,new_Zcord);
		free_xtal(tmp_xtal_each_op);
		fclose(outfile);


   



	}
}

void compare_all_atoms_distance(float* Xcord,float* Ycord,float* Zcord,
								float* new_Xcord,float* new_Ycord,float* new_Zcord)
{


}

