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

#define frac_epsilon 0.000001

int check_layer_group(crystal* xtal)
{

    int attempted_lg = xtal->spg;
    int spglib_spg = detect_spg_using_spglib(xtal);
    int spglib_lg[6]; //max 5 from detected and 1 of its original
    int match = 0;
    for (int i =0;i<5;i++)
    {
        spglib_lg[i] = spg_to_lg[spglib_spg-1][i];
        if (spglib_lg[i] == attempted_lg)
        {
            match = 1;
        }
    }



    //if result in a higher symmetry
    if (match == 0)
    {
	  printf("The attempted layer group is %d\n",attempted_lg);
	  printf( "#SPGLIB_detected_spacegroup = %d\n", spglib_spg);
	  printf("after check all lg,match = 0 \n");
	  fflush(stdout);

	  //add the attemtped lg to test its operations
	  spglib_lg[5] = attempted_lg;
	
	  /*
	  // for check all lg
	  int test = 1;
	  //printf("I am here!!!!!!!!!!!!!!!!");
	  for (int i=5; i<85;i++)
	  {
		  spglib_lg[i] = test;
		  test ++;
	  }
	  */

	  int array_compatible_lg[6];
	  int counter = 0;

	  for (int j = 0; j<6;j++)
	    {
		  int lg = spglib_lg[j];
		  int compatible = 0;
		  
		  if (lg != 0)
		    {
			    float Z = xtal->Z;
    			float N = xtal->num_atoms_in_molecule;
		        crystal *tmp_xtal = (crystal *)malloc(sizeof(crystal) );
    			allocate_xtal(tmp_xtal,Z,N);
			    copy_xtal(tmp_xtal, xtal);
		        compatible = check_symmet_equiv_atoms(lg,tmp_xtal,Z,N);	
				free_xtal(tmp_xtal);
				if (compatible == 1)  //layer group is compatible
				{
					array_compatible_lg[counter] = lg;
					counter ++;
				}
		
						
		    }

	    }
		if (counter == 0)   // no compatible lg, error
		{
			printf("**error: no layer group found compatiblen return is -1**\n");
			return -1;
		}
		else if (counter == 1)
		{
			return array_compatible_lg[0];
		}
		else
		{
			int max_num_ops = 0;
			int matched_lg = 0;
			for (int i =0;i<counter; i++) //go through each compatible lg, return lg of max num operations
			{
				double translations[192][3];
    			int rotations[192][3][3];
				int num_of_operations = get_lg_symmetry(array_compatible_lg[i],translations,rotations);
				if (num_of_operations > max_num_ops)
				{
					max_num_ops = num_of_operations;
					matched_lg = array_compatible_lg[i];
				}
			}

			return matched_lg;


		}


    }
	else
	{
		return attempted_lg;
	}

}

int check_symmet_equiv_atoms(int lg,crystal* tmp_xtal,int Z, int N)
{    
    
    float lattice_vec_a[3];
    float lattice_vec_b[3];
    float lattice_vec_c[3];
    float norm_a = 0;
    float norm_b = 0;
    float norm_c = 0;

	//convert original coords to frac and bring them to positive
	convert_xtal_to_fractional(tmp_xtal);
	bring_all_coords_to_pos(tmp_xtal);
	//FILE *out_file;
	//char f_orig_xtal_name[20] = "original_frac.in";
	//out_file = fopen(f_orig_xtal_name,"w");
	//print_layer2file(1,tmp_xtal,out_file);
	/*
	for (int i =0;i<20;i++)
	{
		printf("%f,%f,%f\n",tmp_xtal->Xcord[i],tmp_xtal->Ycord[i],tmp_xtal->Zcord[i]);
	}
	fclose(out_file);
    */

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
		/*
		// for print operations to file
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
		
		*/

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
		int symme_equivalent = 0;

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

		   //convert back to cartesian
		   vector3_mat3b3_multiply(lattice_vectors_transpose,
				atomj_array, atomj_array);
		   tmp_xtal_each_op->Xcord[j] = atomj_array[0];
		   tmp_xtal_each_op->Ycord[j] = atomj_array[1];
		   tmp_xtal_each_op->Zcord[j] = atomj_array[2];
		}
		
		bring_all_molecules_to_first_cell(tmp_xtal_each_op);
		convert_xtal_to_fractional(tmp_xtal_each_op);
		bring_all_coords_to_pos(tmp_xtal_each_op);
		//print_layer2file(1,tmp_xtal_each_op,outfile);
	    symme_equivalent = compare_all_atoms_distance(tmp_xtal,tmp_xtal_each_op);
		//printf("lg %d operation %d symme_equivalent=%d\n",lg,i+1,symme_equivalent);

		if (symme_equivalent == 0)  //one operation failed to be equivalent
		{
			//printf("lg %d is incompatible\n",lg);
			return 0;
		}

		free_xtal(tmp_xtal_each_op);
		//fclose(outfile);


   



	}

	return 1;
}

int compare_all_atoms_distance(crystal* tmp_xtal,crystal* tmp_xtal_each_op)
{
	int Z = tmp_xtal->Z;
	int N = tmp_xtal->num_atoms_in_molecule;
	int total_atoms = N * Z;
	int num_matched_atoms = 0;
	for (int i =0;i< total_atoms ;i++)
	{
		float atomi_array[3];
		atomi_array[0] = tmp_xtal->Xcord[i];
		atomi_array[1] = tmp_xtal->Ycord[i];
		atomi_array[2] = tmp_xtal->Zcord[i];
		
		//printf("atom checking is %f	%f %f\n",atomi_array[0],atomi_array[1],atomi_array[2]);
		for (int j =0;j< total_atoms  ; j++)
		{
			float atomj_array[3];
			atomj_array[0] = tmp_xtal_each_op->Xcord[j];
			atomj_array[1] = tmp_xtal_each_op->Ycord[j];
			atomj_array[2] = tmp_xtal_each_op->Zcord[j];

			float dist = sqrt(pow((atomi_array[0]-atomj_array[0]),2) + pow((atomi_array[1]-atomj_array[1]),2)
				   + pow((atomi_array[2]-atomj_array[2]),2));

			//printf("dist is %f\n",dist);
			if (dist < 0.000015)  //each coords diff by less than 0.00001, distance diff less than 0.000015
			{
				num_matched_atoms ++;
				//printf("atom %f	%f %f is found\n ",atomi_array[0],atomi_array[1],atomi_array[2]);
				break;
			}

		}

	}

	if (num_matched_atoms == total_atoms)
	{
		return 1;
	}
	else
	{
		return 0;
	}
	


}

//input crystal should be in frac coords
void bring_all_coords_to_pos(crystal* xtal)
{
	int Z = xtal->Z;
	int N = xtal->num_atoms_in_molecule;
	for (int k = 0; k < N * Z;k++)
		{
			float atomj_array[3];
			atomj_array[0] = xtal->Xcord[k];
		    atomj_array[1] = xtal->Ycord[k];
		    atomj_array[2] = xtal->Zcord[k];
			//check negative
			if (atomj_array[0] < 0 && atomj_array[0] >= -1 )
			{
				//atomj_array[0] += 1;
				xtal->Xcord[k] += 1;
			}
			if (atomj_array[1] < 0 && atomj_array[1] >= -1)
			{
				//atomj_array[1] += 1;
				xtal->Ycord[k] += 1;
			}
			if (atomj_array[2] < 0 && atomj_array[2] >= -1)
			{
				//atomj_array[2] += 1;
				xtal->Zcord[k] += 1;
			}


			//check 1
			
			if (atomj_array[0] < -1)
			{
				int int_part = (int)(atomj_array[0]-1);
				xtal->Xcord[k] = atomj_array[0] - int_part;
			}
			if (atomj_array[1] < -1)
			{
				//atomj_array[1] -= 1.0;
				int int_part = (int)(atomj_array[1]-1);
				xtal->Ycord[k] = atomj_array[1]- int_part;
			}
			if (atomj_array[2] < -1)
			{
				//atomj_array[2] -= 1.0;
				int int_part = (int)(atomj_array[2] -1);
				xtal->Zcord[k] = atomj_array[2]- int_part;
			}
			

			//check > 1
			if (atomj_array[0] >=1)
			{
				int int_part = (int)atomj_array[0];
				xtal->Xcord[k] = atomj_array[0] - int_part;
			}
			if (atomj_array[1] >= 1 )
			{
				int int_part = (int)atomj_array[1];
				xtal->Ycord[k] = atomj_array[1] - int_part;
			}
			if (atomj_array[2] >= 1 )
			{
				int int_part = (int)atomj_array[2];
				xtal->Zcord[k] = atomj_array[2] - int_part;
			}


		}
}