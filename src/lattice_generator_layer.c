#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "lattice_generator_layer.h"
#include "randomgen.h"
#include "algebra.h"
//#include "spg_generation.h"
//#include "read_input.h"
#define LOWB 2.5  //lower bound for length of lattice vector
#define PI 3.141592653
#define epsilon_length 0.001
#define epsilon 0.1


//min and max angles are in radians

void generate_oblique_triclinic(float lattice_vector[3][3], 
	float target_volume, float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,
	float gamma)


{	

	//printf("I am in generate_oblique_triclinic");
	int found =0;
	float cz = target_volume / (tmp_lattice_vec_a_norm * tmp_lattice_vec_b_norm * sin(gamma));
	float c;
	float coeff_cx;
	float coeff_cy;

	float random02 = uniform_dist_01();
	float beta = random02*(max_angle - min_angle) + min_angle;

	while (found ==0)
	{
		float random01 = uniform_dist_01();
		float alpha = random01*(max_angle - min_angle) + min_angle;
		float cz_square = pow(cz,2);
		coeff_cx=cos(beta);
		coeff_cy=(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);    //cosbeta2
		float coeff_c=1-pow(coeff_cx,2)-pow(coeff_cy,2);

		if (coeff_c > 0.2 && fabs(coeff_cy) < 1)   //set coeff_c > 0.2 to avoid lattice too oblique
		{
			c = sqrt(cz_square/coeff_c);
			found = 1;
		}

	}

	float cx = coeff_cx * c;
	float cy = coeff_cy * c;

	lattice_vector[0][0] = tmp_lattice_vec_a[0];
	lattice_vector[1][1] = tmp_lattice_vec_b[1];
	lattice_vector[2][2] = cz;

	lattice_vector[0][1] = tmp_lattice_vec_a[1];
	lattice_vector[0][2] = tmp_lattice_vec_a[2];
	lattice_vector[1][2] = tmp_lattice_vec_b[2];

	lattice_vector[1][0] = tmp_lattice_vec_b[0];
	lattice_vector[2][0] = cx;
	lattice_vector[2][1] = cy;


	

}


void generate_oblique_monoclinic(float lattice_vector[3][3], 
	float target_volume, float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,float gamma)


{
	float c = target_volume / (tmp_lattice_vec_a_norm * tmp_lattice_vec_b_norm * sin(gamma));
	//printf("I am in generate_oblique_monoclinic and c is %f",c);
	
	lattice_vector[0][0] = tmp_lattice_vec_a[0];
	lattice_vector[1][1] = tmp_lattice_vec_b[1];
	lattice_vector[2][2] = c;

	lattice_vector[0][1] = tmp_lattice_vec_a[1];
	lattice_vector[0][2] = tmp_lattice_vec_a[2];
	lattice_vector[1][2] = tmp_lattice_vec_b[2];

	lattice_vector[1][0] = tmp_lattice_vec_b[0];
	lattice_vector[2][0] = 0;
	lattice_vector[2][1] = 0;
	//printf("I am at end of generate_oblique_monoclinic");


}

void generate_rectangle_monoclinic(float lattice_vector[3][3], 
	float target_volume,float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,float gamma)


{

	float random01 = uniform_dist_01();
	float alpha = random01*(max_angle - min_angle) + min_angle;
	float c = target_volume / (tmp_lattice_vec_a_norm * tmp_lattice_vec_b_norm * sin(alpha));
	float cy = c * cos(alpha);
	float cz = c * sin(alpha);

	lattice_vector[0][0] = tmp_lattice_vec_a [0];
	lattice_vector[1][1] = tmp_lattice_vec_b [1];
	lattice_vector[2][2] = cz;

	lattice_vector[0][1] = tmp_lattice_vec_a [1];
	lattice_vector[0][2] = tmp_lattice_vec_a [2];
	lattice_vector[1][2] = tmp_lattice_vec_b [2];

	lattice_vector[1][0] = tmp_lattice_vec_b [0];
	lattice_vector[2][0] = 0;
	lattice_vector[2][1] = cy;
}


void generate_rectangle_orthorhombic(float lattice_vector[3][3], 
	float target_volume,float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,	      float gamma)
{


	float c = target_volume / (tmp_lattice_vec_a_norm * tmp_lattice_vec_b_norm);

	lattice_vector[0][0] = tmp_lattice_vec_a [0];
	lattice_vector[1][1] = tmp_lattice_vec_b [1];
	lattice_vector[2][2] = c;

	lattice_vector[0][1] = tmp_lattice_vec_a [1];
	lattice_vector[0][2] = tmp_lattice_vec_a [2];
	lattice_vector[1][2] = tmp_lattice_vec_b [2];

	lattice_vector[1][0] = tmp_lattice_vec_b [0];
	lattice_vector[2][0] = 0;
	lattice_vector[2][1] = 0;
}

void generate_hexagonal_hexagonal(float lattice_vector[3][3], 
	float target_volume,float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,        float gamma)
{
	float c = target_volume / (sqrt(3)/2 * pow(tmp_lattice_vec_a_norm,2));
	//printf("c generated is %f\n",c);
	
	lattice_vector[0][0] = tmp_lattice_vec_a [0];
	lattice_vector[1][1] = tmp_lattice_vec_b [1];
	lattice_vector[2][2] = c;

	lattice_vector[0][1] = tmp_lattice_vec_a [1];
	lattice_vector[0][2] = tmp_lattice_vec_a [2];
	lattice_vector[1][2] = tmp_lattice_vec_b [2];

	lattice_vector[1][0] = tmp_lattice_vec_b [0];
	lattice_vector[2][0] = 0;
	lattice_vector[2][1] = 0;
}

void generate_square_tetragonal(float lattice_vector[3][3], 
	float target_volume,float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,        float gamma)
{
	//printf("I am in square tetragonal");
	//fflush(stdout);
	float c = target_volume / (pow(tmp_lattice_vec_a_norm,2));

	lattice_vector[0][0] = tmp_lattice_vec_a [0];
	lattice_vector[1][1] = tmp_lattice_vec_b [1];
	lattice_vector[2][2] = c;

	lattice_vector[0][1] = tmp_lattice_vec_a [1];
	lattice_vector[0][2] = tmp_lattice_vec_a [2];
	lattice_vector[1][2] = tmp_lattice_vec_b [2];

	lattice_vector[1][0] = tmp_lattice_vec_b [0];
	lattice_vector[2][0] = 0;
	lattice_vector[2][1] = 0;
}


int generate_substrate_lattice_combs(int *all_substrate_combo, float lattice_vector_2d[2][3],
	float target_volume, float max_angle,float min_angle )
{

	//printf("I am in generate_substrate_lattice_combs \n");
	//fflush(stdout);

	float prim_lv[3][3] = {{lattice_vector_2d[0][0],lattice_vector_2d[0][1],lattice_vector_2d[0][2] },
						   {lattice_vector_2d[1][0],lattice_vector_2d[1][1],lattice_vector_2d[1][2] },
						   {0.0,0.0,1.0}};
	float max_allowed_area = 0.75 * target_volume;
	//set the max coefficient to the square root of max area
	int max_coeff_1_4 = (int) sqrt(max_allowed_area ); 
	int max_coeff_2_3 = (int) sqrt(max_allowed_area );
	int num_combo = 0;
	//printf("max_coeff_1_4 is %d\n", max_coeff_1_4);
	//printf("max_coeff_2_3 is %d\n", max_coeff_2_3);

	//int max_num_combo = 100000;
	//int *all_substrate_combo= (int *) malloc(max_num_combo * 4 * sizeof (int));   //max 10000 combinations

	//float transform_martix[3][3];
	/*
	int num_combo = 0;

	float tmp_lattice_vec_a[3];
	float tmp_lattice_vec_b[3];
	float norm_tmp_lattice_vec_a;
	float norm_tmp_lattice_vec_b;
	float gamma;
	float gamma_in_deg;
	int counter = 0;
	float lower_lattice_vector[3][3];
	*/

	// get all substrate combinations
	for (int factor_1 = -max_coeff_1_4; factor_1 <max_coeff_1_4; factor_1++)
	{
		for (int factor_2 = -max_coeff_2_3 + 1; factor_2 <max_coeff_2_3; factor_2++)
		{
			for (int factor_3 = -max_coeff_2_3 + 1; factor_3 <max_coeff_2_3; factor_3++)
			{
				for (int factor_4 = -max_coeff_1_4; factor_4 <max_coeff_1_4; factor_4++)
				{
					// if determinant is zero, skip, enforce right hand chirality now
					if (factor_1 * factor_4 == factor_2 * factor_3 || factor_1 * factor_4 < factor_2 * factor_3)
					//if (factor_1 * factor_4 == factor_2 * factor_3)
					{ 
						continue;
					}

					//else
					float transform_martix[3][3] = {{(float)factor_1,(float)factor_2,0.0},{(float)factor_3,(float)factor_4,0.0},{0.0,0.0,1.0}};
					float new_vec[3][3];
					mat3b3_mat3b3_multiply(transform_martix, prim_lv, new_vec);
					float new_vec_a[3];
					float new_vec_b[3];
					new_vec_a[0] = new_vec[0][0];
					new_vec_a[1] = new_vec[0][1];
					new_vec_a[2] = new_vec[0][2];
					new_vec_b[0] = new_vec[1][0];
					new_vec_b[1] = new_vec[1][1];
					new_vec_b[2] = new_vec[1][2];
					
					float norm_new_vec_a = sqrt(pow(new_vec_a[0],2)+pow(new_vec_a[1],2)+pow(new_vec_a[2],2));
					float norm_new_vec_b = sqrt(pow(new_vec_b[0],2)+pow(new_vec_b[1],2)+pow(new_vec_b[2],2));

					// if the ratio netween the length of two lattice vector is more than 3, skip
					if (norm_new_vec_a / norm_new_vec_b >3 || norm_new_vec_b / norm_new_vec_a >3)
					{
						continue;
					}
					float cross[3];
					cross_vector3_vector3(cross, new_vec_a,new_vec_b);
					float new_area = sqrt(pow(cross[0],2)+pow(cross[1],2)+pow(cross[2],2));
					float dot_product = new_vec_a[0] * new_vec_b[0] + new_vec_a[1] * new_vec_b[1] 
								+new_vec_a[2] * new_vec_b[2] ;
					float angle = acos(dot_product/(norm_new_vec_a * norm_new_vec_b)) * 180/PI;					

					//if film area is greater than max_allowed, break the loop
					if (new_area > max_allowed_area )
					{
						continue;

					}

					//if angle is not between 30 and 50, skip
					if (angle < min_angle * 180/PI || angle > max_angle * 180/PI)
					{
						continue;
					}

					*(all_substrate_combo + num_combo*4 + 0) = factor_1;
					*(all_substrate_combo + num_combo*4 + 1) = factor_2;
					*(all_substrate_combo + num_combo*4 + 2) = factor_3;
					*(all_substrate_combo + num_combo*4 + 3) = factor_4;

					num_combo ++;
					/*
					if (fabs(norm_new_vec_a-norm_new_vec_b)<1)
					{
						printf("fac_1, fac_2, fac_3, fac_4: %d %d %d %d\n",factor_1,factor_2,factor_3,factor_4);
						printf("norm_new_vec_a,norm_new_vec_b: %f %f\n",norm_new_vec_a,norm_new_vec_b);
						printf("new_area: %f\n",new_area);
						printf("angle: %f\n",angle);
						fflush(stdout);
					}
					*/

				}


			}
		}

	}

	return num_combo;


}


void generate_layer_lattice(int *all_substrate_combo,float lattice_vector[3][3], int spg,
	float max_angle, float min_angle, float target_volume, float lattice_vector_2d[2][3],int num_combo,float interface_area_mean,float interface_area_std,
	int volume_multiplier,int SET_INTERFACE_AREA)
{	
	
	//printf("num combo is %d\n",num_combo);
	//fflush(stdout);
	float tmp_lattice_vec_a[3];
	float tmp_lattice_vec_b[3];
	float norm_tmp_lattice_vec_a = 0;
	float norm_tmp_lattice_vec_b = 0;
	float gamma = 0;
	float gamma_in_deg = 0;
	int counter = 0;
	float new_area = 0;
	float lower_lattice_vector[3][3];
	float prim_lv[3][3] = {{lattice_vector_2d[0][0],lattice_vector_2d[0][1],lattice_vector_2d[0][2] },
						   {lattice_vector_2d[1][0],lattice_vector_2d[1][1],lattice_vector_2d[1][2] },
						   {0.0,0.0,1.0}};
	//printf("volume_multiplier is: %d\n",volume_multiplier );
	//printf("new volume is %f\n",target_volume * volume_multiplier);
	//fflush(stdout);
	if (SET_INTERFACE_AREA == 0)
		{
			interface_area_mean = target_volume;
			interface_area_std = target_volume;

		}
	
	//choose a combo based on lattice type of spg
	if (spg <=2)  //oblique_triclinic
	{
	  //printf("the layer group is oblique_triclinic");
	  //fflush(stdout);
	  do 
	    {
			
			int rand_index = genrand_int32()% num_combo;
			int chosen_factor_1= all_substrate_combo[rand_index * 4 +0 ];
			int chosen_factor_2= all_substrate_combo[rand_index * 4 +1 ];
			int chosen_factor_3= all_substrate_combo[rand_index * 4 +2 ];
			int chosen_factor_4= all_substrate_combo[rand_index * 4 +3 ];

			float transform_martix[3][3] = {{chosen_factor_1,chosen_factor_2,0},{chosen_factor_3,chosen_factor_4,0},{0,0,1}};
			float new_vec[3][3];
			mat3b3_mat3b3_multiply(transform_martix, prim_lv, new_vec);
			get_lower_triangle(new_vec,lower_lattice_vector);

			tmp_lattice_vec_a[0] = lower_lattice_vector[0][0];
			tmp_lattice_vec_a[1] = lower_lattice_vector[0][1];
			tmp_lattice_vec_a[2] = lower_lattice_vector[0][2];
			tmp_lattice_vec_b[0] = lower_lattice_vector[1][0];
			tmp_lattice_vec_b[1] = lower_lattice_vector[1][1];
			tmp_lattice_vec_b[2] = lower_lattice_vector[1][2];


			norm_tmp_lattice_vec_a = sqrt(pow(tmp_lattice_vec_a[0],2)+pow(tmp_lattice_vec_a[1],2)+pow(tmp_lattice_vec_a[2],2));
			norm_tmp_lattice_vec_b = sqrt(pow(tmp_lattice_vec_b[0],2)+pow(tmp_lattice_vec_b[1],2)+pow(tmp_lattice_vec_b[2],2));
			float dot_product = tmp_lattice_vec_a[0] * tmp_lattice_vec_b[0] + tmp_lattice_vec_a[1] * tmp_lattice_vec_b[1] 
								+tmp_lattice_vec_a[2] * tmp_lattice_vec_b[2] ;
			gamma = acos(dot_product/(norm_tmp_lattice_vec_a * norm_tmp_lattice_vec_b));
			gamma_in_deg = gamma * 180/PI;
			counter++;
			if (counter > 10000000)
			{
				for (int i = 0; i < 3; i++)
					{
						for (int j =0; j<3; j++)
							{
								lattice_vector[i][j] = 0;
							}
					}

				return;
			}

			//get area and check it if the interface mean and std is set
			float cross[3];
			cross_vector3_vector3(cross, tmp_lattice_vec_a,tmp_lattice_vec_b);
			new_area = sqrt(pow(cross[0],2)+pow(cross[1],2)+pow(cross[2],2));

	    }	
		//remove the constraint that a,b should be different
	  while (fabs(gamma_in_deg - 90) < epsilon || fabs(gamma_in_deg - 120) < epsilon ||
	        new_area > (interface_area_mean + interface_area_std ) || 
		new_area < (interface_area_mean - interface_area_std));
	  
	  
	  generate_oblique_triclinic( lattice_vector,target_volume * volume_multiplier,max_angle, min_angle, tmp_lattice_vec_a,tmp_lattice_vec_b,norm_tmp_lattice_vec_a,norm_tmp_lattice_vec_b,gamma);
	}


	else if (spg <= 7)  //oblique_monoclinic
	{
	  //printf("the layer group is oblique_monoclinic");
	  //fflush(stdout);
      do
        {
			int rand_index = genrand_int32()% num_combo;
			int chosen_factor_1= all_substrate_combo[rand_index * 4 +0 ];
			int chosen_factor_2= all_substrate_combo[rand_index * 4 +1 ];
			int chosen_factor_3= all_substrate_combo[rand_index * 4 +2 ];
			int chosen_factor_4= all_substrate_combo[rand_index * 4 +3 ];

			float transform_martix[3][3] = {{chosen_factor_1,chosen_factor_2,0},{chosen_factor_3,chosen_factor_4,0},{0,0,1}};
			float new_vec[3][3];
			mat3b3_mat3b3_multiply(transform_martix, prim_lv, new_vec);
			get_lower_triangle(new_vec,lower_lattice_vector);

			tmp_lattice_vec_a[0] = lower_lattice_vector[0][0];
			tmp_lattice_vec_a[1] = lower_lattice_vector[0][1];
			tmp_lattice_vec_a[2] = lower_lattice_vector[0][2];
			tmp_lattice_vec_b[0] = lower_lattice_vector[1][0];
			tmp_lattice_vec_b[1] = lower_lattice_vector[1][1];
			tmp_lattice_vec_b[2] = lower_lattice_vector[1][2];

			norm_tmp_lattice_vec_a = sqrt(pow(tmp_lattice_vec_a[0],2)+pow(tmp_lattice_vec_a[1],2)+pow(tmp_lattice_vec_a[2],2));
			norm_tmp_lattice_vec_b = sqrt(pow(tmp_lattice_vec_b[0],2)+pow(tmp_lattice_vec_b[1],2)+pow(tmp_lattice_vec_b[2],2));
			float dot_product = tmp_lattice_vec_a[0] * tmp_lattice_vec_b[0] + tmp_lattice_vec_a[1] * tmp_lattice_vec_b[1] 
								+tmp_lattice_vec_a[2] * tmp_lattice_vec_b[2] ;
			gamma = acos(dot_product/(norm_tmp_lattice_vec_a * norm_tmp_lattice_vec_b));
			gamma_in_deg = gamma * 180/PI;
			counter++;
			if (counter > 10000000)
			{
				for (int i = 0; i < 3; i++)
					{
						for (int j =0; j<3; j++)
							{
								lattice_vector[i][j] = 0;
							}
					}

				return;
			}

			//get area and check it if the interface mean and std is set
			float cross[3];
			cross_vector3_vector3(cross, tmp_lattice_vec_a,tmp_lattice_vec_b);
			new_area = sqrt(pow(cross[0],2)+pow(cross[1],2)+pow(cross[2],2));

		
			
        }
		 while (fabs(gamma_in_deg - 90) < epsilon || fabs(gamma_in_deg - 120) < epsilon || 
	     new_area > (interface_area_mean + interface_area_std ) || 
	     new_area < (interface_area_mean - interface_area_std));

        //while (fabs(norm_tmp_lattice_vec_a - norm_tmp_lattice_vec_b) < epsilon_length ||
	    // fabs(gamma_in_deg - 90) < epsilon || fabs(gamma_in_deg - 120) < epsilon || 
	     //new_area > (interface_area_mean + interface_area_std ) || 
	    // new_area < (interface_area_mean - interface_area_std));

	  generate_oblique_monoclinic( lattice_vector,target_volume * volume_multiplier,max_angle, min_angle, tmp_lattice_vec_a,tmp_lattice_vec_b,norm_tmp_lattice_vec_a,norm_tmp_lattice_vec_b,gamma);
        }   
	

	else if (spg <= 18)   //rectangle_monoclinic
    {
		  //printf("the layer group is rectangle_monoclinic");
	  	  //fflush(stdout);
     do
        {
			int rand_index = genrand_int32()% num_combo;
			int chosen_factor_1= all_substrate_combo[rand_index * 4 +0 ];
			int chosen_factor_2= all_substrate_combo[rand_index * 4 +1 ];
			int chosen_factor_3= all_substrate_combo[rand_index * 4 +2 ];
			int chosen_factor_4= all_substrate_combo[rand_index * 4 +3 ];
			float transform_martix[3][3] = {{chosen_factor_1,chosen_factor_2,0},{chosen_factor_3,chosen_factor_4,0},{0,0,1}};
			//float transform_martix[3][3] = {{6,3,0},{0,5,0},{0,0,1}};
			float new_vec[3][3];
			mat3b3_mat3b3_multiply(transform_martix, prim_lv, new_vec);
			get_lower_triangle(new_vec,lower_lattice_vector);

			tmp_lattice_vec_a[0] = lower_lattice_vector[0][0];
			tmp_lattice_vec_a[1] = lower_lattice_vector[0][1];
			tmp_lattice_vec_a[2] = lower_lattice_vector[0][2];
			tmp_lattice_vec_b[0] = lower_lattice_vector[1][0];
			tmp_lattice_vec_b[1] = lower_lattice_vector[1][1];
			tmp_lattice_vec_b[2] = lower_lattice_vector[1][2];


			norm_tmp_lattice_vec_a = sqrt(pow(tmp_lattice_vec_a[0],2)+pow(tmp_lattice_vec_a[1],2)+pow(tmp_lattice_vec_a[2],2));
			norm_tmp_lattice_vec_b = sqrt(pow(tmp_lattice_vec_b[0],2)+pow(tmp_lattice_vec_b[1],2)+pow(tmp_lattice_vec_b[2],2));
			float dot_product = tmp_lattice_vec_a[0] * tmp_lattice_vec_b[0] + tmp_lattice_vec_a[1] * tmp_lattice_vec_b[1] 
								+tmp_lattice_vec_a[2] * tmp_lattice_vec_b[2] ;
			gamma = acos(dot_product/(norm_tmp_lattice_vec_a * norm_tmp_lattice_vec_b));
			gamma_in_deg = gamma * 180/PI;
			counter++;
			if (counter > 10000000)
			{
				for (int i = 0; i < 3; i++)
					{
						for (int j =0; j<3; j++)
							{
								lattice_vector[i][j] = 0;
							}
					}

				return;
			}

			//get area and check it if the interface mean and std is set
			float cross[3];
			cross_vector3_vector3(cross, tmp_lattice_vec_a,tmp_lattice_vec_b);
		    new_area = sqrt(pow(cross[0],2)+pow(cross[1],2)+pow(cross[2],2));


        }

          while (fabs(norm_tmp_lattice_vec_a - norm_tmp_lattice_vec_b) < epsilon_length || 
		 fabs(gamma_in_deg - 90) > epsilon || 
	         new_area > (interface_area_mean + interface_area_std ) || 
		 new_area < (interface_area_mean - interface_area_std));
          
		  
		  generate_rectangle_monoclinic(lattice_vector,target_volume * volume_multiplier,max_angle, min_angle, tmp_lattice_vec_a,tmp_lattice_vec_b,norm_tmp_lattice_vec_a,norm_tmp_lattice_vec_b,gamma);
    }

	 else if (spg <= 48)   //rectangle_orthorombic
    {
		  //printf("the layer group is rectangle_orthorombic");
	  	  //fflush(stdout);
		  
    do
        {
			int rand_index = genrand_int32()% num_combo;
			int chosen_factor_1= all_substrate_combo[rand_index * 4 +0 ];
			int chosen_factor_2= all_substrate_combo[rand_index * 4 +1 ];
			int chosen_factor_3= all_substrate_combo[rand_index * 4 +2 ];
			int chosen_factor_4= all_substrate_combo[rand_index * 4 +3 ];

			float transform_martix[3][3] = {{chosen_factor_1,chosen_factor_2,0},{chosen_factor_3,chosen_factor_4,0},{0,0,1}};
			float new_vec[3][3];
			mat3b3_mat3b3_multiply(transform_martix, prim_lv, new_vec);
			get_lower_triangle(new_vec,lower_lattice_vector);

			tmp_lattice_vec_a[0] = lower_lattice_vector[0][0];
			tmp_lattice_vec_a[1] = lower_lattice_vector[0][1];
			tmp_lattice_vec_a[2] = lower_lattice_vector[0][2];
			tmp_lattice_vec_b[0] = lower_lattice_vector[1][0];
			tmp_lattice_vec_b[1] = lower_lattice_vector[1][1];
			tmp_lattice_vec_b[2] = lower_lattice_vector[1][2];


			norm_tmp_lattice_vec_a = sqrt(pow(tmp_lattice_vec_a[0],2)+pow(tmp_lattice_vec_a[1],2)+pow(tmp_lattice_vec_a[2],2));
			norm_tmp_lattice_vec_b = sqrt(pow(tmp_lattice_vec_b[0],2)+pow(tmp_lattice_vec_b[1],2)+pow(tmp_lattice_vec_b[2],2));
			float dot_product = tmp_lattice_vec_a[0] * tmp_lattice_vec_b[0] + tmp_lattice_vec_a[1] * tmp_lattice_vec_b[1] 
								+tmp_lattice_vec_a[2] * tmp_lattice_vec_b[2] ;
			gamma = acos(dot_product/(norm_tmp_lattice_vec_a * norm_tmp_lattice_vec_b));
			gamma_in_deg = gamma * 180/PI;
			counter++;
			if (counter > 10000000)
			{
		
				for (int i = 0; i < 3; i++)
					{
						for (int j =0; j<3; j++)
							{
								lattice_vector[i][j] = 0;
							}
					}

				return;
			}

			//get area and check it if the interface mean and std is set
			float cross[3];
			cross_vector3_vector3(cross, tmp_lattice_vec_a,tmp_lattice_vec_b);
			new_area = sqrt(pow(cross[0],2)+pow(cross[1],2)+pow(cross[2],2));

			
        }

          while (fabs(norm_tmp_lattice_vec_a - norm_tmp_lattice_vec_b) < epsilon_length || 
		 fabs(gamma_in_deg - 90) > epsilon || 
	         new_area > (interface_area_mean + interface_area_std ) || 
	  	 new_area < (interface_area_mean - interface_area_std));
          
		  generate_rectangle_orthorhombic( lattice_vector,target_volume * volume_multiplier,max_angle, min_angle, tmp_lattice_vec_a,tmp_lattice_vec_b,norm_tmp_lattice_vec_a,norm_tmp_lattice_vec_b,gamma);
    }

	else if (spg <= 64)  //square_tetragonal
    {
		  //printf("the layer group is square_tetragonal");
	  	  //fflush(stdout);
        do
        {
			int rand_index = genrand_int32()% num_combo;
			int chosen_factor_1= all_substrate_combo[rand_index * 4 +0 ];
			int chosen_factor_2= all_substrate_combo[rand_index * 4 +1 ];
			int chosen_factor_3= all_substrate_combo[rand_index * 4 +2 ];
			int chosen_factor_4= all_substrate_combo[rand_index * 4 +3 ];

			float transform_martix[3][3] = {{chosen_factor_1,chosen_factor_2,0},{chosen_factor_3,chosen_factor_4,0},{0,0,1}};
			float new_vec[3][3];
			mat3b3_mat3b3_multiply(transform_martix, prim_lv, new_vec);
			get_lower_triangle(new_vec,lower_lattice_vector);

			tmp_lattice_vec_a[0] = lower_lattice_vector[0][0];
			tmp_lattice_vec_a[1] = lower_lattice_vector[0][1];
			tmp_lattice_vec_a[2] = lower_lattice_vector[0][2];
			tmp_lattice_vec_b[0] = lower_lattice_vector[1][0];
			tmp_lattice_vec_b[1] = lower_lattice_vector[1][1];
			tmp_lattice_vec_b[2] = lower_lattice_vector[1][2];
			norm_tmp_lattice_vec_a = sqrt(pow(tmp_lattice_vec_a[0],2)+pow(tmp_lattice_vec_a[1],2)+pow(tmp_lattice_vec_a[2],2));
			norm_tmp_lattice_vec_b = sqrt(pow(tmp_lattice_vec_b[0],2)+pow(tmp_lattice_vec_b[1],2)+pow(tmp_lattice_vec_b[2],2));
			float dot_product = tmp_lattice_vec_a[0] * tmp_lattice_vec_b[0] + tmp_lattice_vec_a[1] * tmp_lattice_vec_b[1] 
								+tmp_lattice_vec_a[2] * tmp_lattice_vec_b[2] ;
			gamma = acos(dot_product/(norm_tmp_lattice_vec_a * norm_tmp_lattice_vec_b));
			gamma_in_deg = gamma * 180/PI;
			counter++;
			if (counter > 10000000)
			{
				
				for (int i = 0; i < 3; i++)
					{
						for (int j =0; j<3; j++)
							{
								lattice_vector[i][j] = 0;
							}
					}

				return;
			}

			//get area and check it if the interface mean and std is set
			float cross[3];
			cross_vector3_vector3(cross, tmp_lattice_vec_a,tmp_lattice_vec_b);
			new_area = sqrt(pow(cross[0],2)+pow(cross[1],2)+pow(cross[2],2));


        }
          while (fabs(norm_tmp_lattice_vec_a - norm_tmp_lattice_vec_b) >= epsilon_length || 
		 fabs(gamma_in_deg - 90) > epsilon || 
	         new_area > (interface_area_mean + interface_area_std ) || 
	         new_area < (interface_area_mean - interface_area_std));
          
		  
		  generate_square_tetragonal( lattice_vector,target_volume * volume_multiplier,max_angle, min_angle, tmp_lattice_vec_a,tmp_lattice_vec_b,norm_tmp_lattice_vec_a,norm_tmp_lattice_vec_b,gamma);
    }

	else if (spg <= 80)  //hexagonal_trigonnal, hexagonal_hexagonal
    {
		  //printf("the layer group is hexagonal_trigonnal, hexagonal_hexagonal");
	  	  //fflush(stdout);
        do
        {
			int rand_index = genrand_int32()% num_combo;
			int chosen_factor_1= all_substrate_combo[rand_index * 4 +0 ];
			int chosen_factor_2= all_substrate_combo[rand_index * 4 +1 ];
			int chosen_factor_3= all_substrate_combo[rand_index * 4 +2 ];
			int chosen_factor_4= all_substrate_combo[rand_index * 4 +3 ];

			float transform_martix[3][3] = {{chosen_factor_1,chosen_factor_2,0},{chosen_factor_3,chosen_factor_4,0},{0,0,1}};
			float new_vec[3][3];
			mat3b3_mat3b3_multiply(transform_martix, prim_lv, new_vec);
			get_lower_triangle(new_vec,lower_lattice_vector);

			tmp_lattice_vec_a[0] = lower_lattice_vector[0][0];
			tmp_lattice_vec_a[1] = lower_lattice_vector[0][1];
			tmp_lattice_vec_a[2] = lower_lattice_vector[0][2];
			tmp_lattice_vec_b[0] = lower_lattice_vector[1][0];
			tmp_lattice_vec_b[1] = lower_lattice_vector[1][1];
			tmp_lattice_vec_b[2] = lower_lattice_vector[1][2];
			norm_tmp_lattice_vec_a = sqrt(pow(tmp_lattice_vec_a[0],2)+pow(tmp_lattice_vec_a[1],2)+pow(tmp_lattice_vec_a[2],2));
			norm_tmp_lattice_vec_b = sqrt(pow(tmp_lattice_vec_b[0],2)+pow(tmp_lattice_vec_b[1],2)+pow(tmp_lattice_vec_b[2],2));
			float dot_product = tmp_lattice_vec_a[0] * tmp_lattice_vec_b[0] + tmp_lattice_vec_a[1] * tmp_lattice_vec_b[1] 
								+tmp_lattice_vec_a[2] * tmp_lattice_vec_b[2] ;
			gamma = acos(dot_product/(norm_tmp_lattice_vec_a * norm_tmp_lattice_vec_b));
			gamma_in_deg = gamma * 180/PI;
			counter++;
			if (counter > 10000000)
			{
	
				for (int i = 0; i < 3; i++)
					{
						for (int j =0; j<3; j++)
							{
								lattice_vector[i][j] = 0;
							}
					}

				return;
			}


			//get area and check it if the interface mean and std is set
			float cross[3];
			cross_vector3_vector3(cross, tmp_lattice_vec_a,tmp_lattice_vec_b);
			new_area = sqrt(pow(cross[0],2)+pow(cross[1],2)+pow(cross[2],2));


        }
          while (fabs(norm_tmp_lattice_vec_a - norm_tmp_lattice_vec_b) >= epsilon_length ||
		 fabs(gamma_in_deg - 120) > epsilon || 
		 new_area > (interface_area_mean + interface_area_std ) || 
	      	 new_area < (interface_area_mean - interface_area_std));


          generate_hexagonal_hexagonal(lattice_vector,target_volume * volume_multiplier,max_angle, min_angle, tmp_lattice_vec_a,tmp_lattice_vec_b,norm_tmp_lattice_vec_a,norm_tmp_lattice_vec_b,gamma);
    }



	//free(all_substrate_combo);



}

//create a large volume lattice for testing compatiility
void generate_fake_layer_lattice(float lattice_vector[3][3], int spg)
{
	const float ax = 15;
	
	
	if(spg < 1 || spg > 80)
	printf("***ERROR: generate_lattice: spg out of bounds***");
	
	else if (spg <= 2)      //oblique_triclinic
	{	
		lattice_vector[0][0] = ax;
		lattice_vector[1][1] = 0.8*ax;
		lattice_vector[2][2] = 0.5*ax;

		lattice_vector[0][1] = 0;
		lattice_vector[0][2] = 0;
		lattice_vector[1][2] = 0;

		lattice_vector[1][0] = 0.8*ax*tan(10*PI/180);
		lattice_vector[2][0] = 0.5*ax*cos(85*PI/180);
		lattice_vector[2][1] = 0.5*ax*cos(70*PI/180);
	}
	
	else if ( spg <= 7  )    //oblique_monoclinic, alpha=beta=90
        {
                lattice_vector[0][0] = ax;
                lattice_vector[1][1] = 0.8*ax;
                lattice_vector[2][2] = 0.5*ax;

                lattice_vector[0][1] = ax * tan(70*PI/180);
                lattice_vector[0][2] = 0;
                lattice_vector[1][2] = 0;

                lattice_vector[1][0] = 0;
                lattice_vector[2][0] = 0;
                lattice_vector[2][1] = 0;
        }


	else if (  spg <= 18 )    //rectangle_monoclinic, gamma = beta = 90
	{
		lattice_vector[0][0] = ax;
		lattice_vector[1][1] = 0.8*ax;
		lattice_vector[2][2] = 0.5*ax;

		lattice_vector[0][1] = 0;
		lattice_vector[0][2] = 0;
		lattice_vector[1][2] = 0;

		lattice_vector[1][0] = 0;
		lattice_vector[2][0] = 0;
		lattice_vector[2][1] = 0.5*ax/tan(70*PI/180);
	}
	
	else if (spg <= 48)       //rectangle_orthorombic
	{
		lattice_vector[0][0] = ax;
		lattice_vector[1][1] = 0.8*ax;
		lattice_vector[2][2] = 0.5*ax;

		lattice_vector[0][1] = 0;
		lattice_vector[0][2] = 0;
		lattice_vector[1][2] = 0;

		lattice_vector[1][0] = 0;
		lattice_vector[2][0] = 0;
		lattice_vector[2][1] = 0;
	}
	
	else if (spg <= 64)     //square_tetragonal
	{
		lattice_vector[0][0] = ax;
		lattice_vector[1][1] = ax;
		lattice_vector[2][2] = 0.5*ax;

		lattice_vector[0][1] = 0;
		lattice_vector[0][2] = 0;
		lattice_vector[1][2] = 0;

		lattice_vector[1][0] = 0;
		lattice_vector[2][0] = 0;
		lattice_vector[2][1] = 0;
	}
	
	else if (spg <= 80)   //hexagonal_trigonnal, hexagonal_hexagonal
	{
		lattice_vector[0][0] = ax;
		lattice_vector[1][1] = ax*sin(120*PI/180);
		lattice_vector[2][2] = 0.5*ax;

		lattice_vector[0][1] = 0;
		lattice_vector[0][2] = 0;
		lattice_vector[1][2] = 0;

		lattice_vector[1][0] = ax*cos(120*PI/180);
		lattice_vector[2][0] = 0;
		lattice_vector[2][1] = 0;
	}
	
}

static inline float fmodulo (float n, float d)
{
	long q =n/d;
	float r = n - q*d ;
	return r;	
}
//need to change standardise lattice!!!!!
void standardise_layer_lattice( float lattice[3][3], int spg)
{
	//triclinic
	if (spg == 1 || spg == 2)
	{
		float bx = lattice[1][0];
		float ax = lattice[0][0];
		float bx_new = fmodulo (bx, ax);
		lattice[1][0] = bx_new;
		
		float by = lattice[1][1];
		float cy = lattice[2][1];
		float cy_new = fmodulo(cy, by);
		lattice[2][1] = cy_new;
		
		float cx = lattice[2][0];
		float cx_new = fmodulo(cx, ax);
		lattice[2][0] = cx_new;
	}
	
	/*
	else if (spg > 2 && spg < 16)
	{
		float bx = lattice[1][0];
		float ax = lattice[0][0];
		float bx_new = fmodulo (bx, ax);
		lattice[1][0] = bx_new;
		
		float cx = lattice[2][0];
		float cx_new = fmodulo(cx, ax);
		lattice[2][0] = cx_new;
		
	}
	*/
	
	else if (spg <= 0 || spg > 80)
		printf("***ERROR: lattice_generator: standardise_lattice invalid spg");
	
}

void get_lower_triangle(float lattice_vector[3][3],float lower_lattice_vector[3][3])
{
	
	float norm_lattice_vectors_a = sqrt(pow(lattice_vector[0][0],2)+pow(lattice_vector[0][1],2)+pow(lattice_vector[0][2],2));
	float norm_lattice_vectors_b = sqrt(pow(lattice_vector[1][0],2)+pow(lattice_vector[1][1],2)+pow(lattice_vector[1][2],2));
	//float norm_lattice_vectors_c = sqrt(pow(tmp_lattice_vectors[2][0],2)+pow(tmp_lattice_vectors[2][1],2)+
	//																		pow(tmp_lattice_vectors[2][2],2));
	float dot_product_ab = lattice_vector[0][0] * lattice_vector[1][0] + lattice_vector[0][1] * lattice_vector[1][1]+lattice_vector[0][2] * lattice_vector[1][2];
	//float dot_product_bc = tmp_lattice_vectors[1][0] * tmp_lattice_vectors[2][0] + tmp_lattice_vectors[1][1] * 
	//										tmp_lattice_vectors[2][1]+tmp_lattice_vectors[1][2] * tmp_lattice_vectors[2][2];
	//float dot_product_ac = tmp_lattice_vectors[0][0] * tmp_lattice_vectors[2][0] + tmp_lattice_vectors[0][1] * 
	//										tmp_lattice_vectors[2][1]+tmp_lattice_vectors[0][2] * tmp_lattice_vectors[2][2];

	//float alpha = acos(dot_product_bc/(norm_lattice_vectors_b * norm_lattice_vectors_c));
	//float beta = acos(dot_product_ac/(norm_lattice_vectors_a * norm_lattice_vectors_c));
	float gamma = acos(dot_product_ab/(norm_lattice_vectors_a * norm_lattice_vectors_b));

	lower_lattice_vector[0][0] = norm_lattice_vectors_a;
	lower_lattice_vector[0][1] = 0;
	lower_lattice_vector[0][2] = 0;
	lower_lattice_vector[1][0] = cos(gamma) * norm_lattice_vectors_b;
	lower_lattice_vector[1][1] = sin(gamma) * norm_lattice_vectors_b;
	lower_lattice_vector[1][2] = 0;
	//random_crystal->lattice_vectors[2][0] = cos(beta) * norm_lattice_vectors_c;
	//random_crystal->lattice_vectors[2][1] = (norm_lattice_vectors_b*norm_lattice_vectors_c*cos(alpha)-
	//	random_crystal->lattice_vectors[1][0] * random_crystal->lattice_vectors[2][0])/random_crystal->lattice_vectors[1][1];
    //random_crystal->lattice_vectors[2][2] = pow((pow(norm_lattice_vectors_c,2)-pow(random_crystal->lattice_vectors[2][0],2)- 
	//													pow(random_crystal->lattice_vectors[2][1],2)),0.5);


}











