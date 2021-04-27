#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "lattice_generator.h"
#include "randomgen.h"
//#include "niggli.h"


////////// MACROS AND DECLARATIONS //////////////////
#define LOWB        3  //lower bound for length of lattice vector
#define PI          3.141592653
#define MAX_ANGLE   150 * PI / 180
#define MIN_ANGLE   30 * PI / 180
#define EPS         0.1
enum COMP_TYPE {DISTINCT, TWO_EQUAL, ALL_EQUAL};

static float gen_angle(float angle_std);
static void set_lattice_vectors_zero(float lattice_vector[3][3]);
static void gen_principal_comps(float *ax,
                                float *by,
                                float *cz,
                                float target_volume,
                                float norm_std,
                                int type);


void generate_lattice(float lattice_vector[3][3],
                      int spg,
                      float norm_std,
                      float angle_std,
                      float target_volume)
{
    static int first_time = 1;
    static FILE *lattice_file = NULL;
    static float const_lattice[3][3];

    // If there's a lattice file, donot generate random lattices.
    if(first_time)
    {
        if(access("lattice.dat", F_OK) == 0)
        {
            lattice_file = fopen("lattice.dat", "r");
            int nmatches = 0;
            nmatches = fscanf(lattice_file, "%f %f %f", &const_lattice[0][0],
                &const_lattice[0][1], &const_lattice[0][2]);
            nmatches += fscanf(lattice_file, "%f %f %f", &const_lattice[1][0],
                &const_lattice[1][1], &const_lattice[1][2]);
            nmatches += fscanf(lattice_file, "%f %f %f", &const_lattice[2][0],
                &const_lattice[2][1], &const_lattice[2][2]);
            // Check if read correctly
            if(nmatches != 9)
            {
                printf("***ERROR: bad lattice.dat file; unable to read\n");
                exit(EXIT_FAILURE);
            }
            fclose(lattice_file);
        }
        first_time = 0;
    }

    if(lattice_file != NULL)
    {
        memcpy(lattice_vector, const_lattice, 9*sizeof(lattice_vector[0][0]));
        return;
    }

    if(spg < 1 || spg > 230)
    {
        printf("***ERROR: generate_lattice: spg out of bounds***");
        exit(EXIT_FAILURE);
    }

    else if (spg <= 2)
        gen_triclinic_lattice(lattice_vector, target_volume, norm_std, angle_std);

    else if (spg <= 15)
        gen_monoclinic_lattice(lattice_vector, target_volume, norm_std, angle_std);

    else if (spg <= 74)
        gen_orthorhombic_lattice(lattice_vector, target_volume, norm_std);

    else if (spg <= 142)
        gen_tetragonal_lattice(lattice_vector, target_volume, norm_std);

    else if (spg <= 167)
        gen_hexagonal_lattice(lattice_vector, target_volume, norm_std);
        //same as hexagonal?

    else if (spg <= 194)
        gen_hexagonal_lattice(lattice_vector, target_volume, norm_std);
    else if (spg <= 230)
        gen_cubic_lattice(lattice_vector, target_volume);

    standardise_lattice(lattice_vector, spg);

    return;
}

////////////////////// FUNCTIONS FOR GENERATING DIFFERENT LATTICES ///////////////

void gen_triclinic_lattice(float lattice_vector[3][3],
                           float target_volume,
                           float norm_std,
                           float angle_std)
{
    float ax, by, cz;

    //select principal components
    gen_principal_comps(&ax, &by, &cz, target_volume, norm_std, DISTINCT);

    //generate random angles
    float beta = gen_angle(angle_std);
    float gamma = gen_angle(angle_std);

    // Not all values of alpha are allowed.
    // Loop until a valid value is selected.
    float alpha;
    float cosbeta2;
    float cosbeta3;
    do
    {
        alpha = gen_angle(angle_std);
        cosbeta2 = (cos(alpha) - cos(gamma)*cos(beta))/sin(gamma);
        cosbeta3 = sqrt(1 - cos(beta)*cos(beta) - cosbeta2*cosbeta2);
    }
    while( fabs(cosbeta2) >= 1 || (cos(beta)*cos(beta) + cosbeta2*cosbeta2) > 1 );

    float bx = by/tan(gamma);
    float modc = cz/cosbeta3;
    float cx = modc*cos(beta);
    float cy = modc*cosbeta2;

    set_lattice_vectors_zero(lattice_vector);

    lattice_vector[0][0] = ax;
    lattice_vector[1][1] = by;
    lattice_vector[2][2] = cz;

    lattice_vector[1][0] = bx;
    lattice_vector[2][0] = cx;
    lattice_vector[2][1] = cy;
}


void gen_monoclinic_lattice(float lattice_vector[3][3],
                            float target_volume,
                            float norm_std,
                            float angle_std)
{
    float ax, by, cz;
    gen_principal_comps(&ax, &by, &cz, target_volume, norm_std, DISTINCT);

    float beta = gen_angle(angle_std);
    float cx = cz / tan(beta);

    set_lattice_vectors_zero(lattice_vector);
    lattice_vector[0][0] = ax;
    lattice_vector[1][1] = by;
    lattice_vector[2][2] = cz;
    lattice_vector[2][0] = cx;
    return;
}

void gen_orthorhombic_lattice(float lattice_vector[3][3],
                              float target_volume,
                              float norm_std)
{
    float ax, by, cz;
    gen_principal_comps(&ax, &by, &cz, target_volume, norm_std, DISTINCT);

    set_lattice_vectors_zero(lattice_vector);
    lattice_vector[0][0] = ax;
    lattice_vector[1][1] = by;
    lattice_vector[2][2] = cz;
    return;
}


void gen_tetragonal_lattice(float lattice_vector[3][3],
                            float target_volume,
                            float norm_std)
{
    float ax, by, cz;
    gen_principal_comps(&ax, &by, &cz, target_volume, norm_std, TWO_EQUAL);

    set_lattice_vectors_zero(lattice_vector);
    lattice_vector[0][0] = ax;
    lattice_vector[1][1] = by;
    lattice_vector[2][2] = cz;

    return;
}

void gen_hexagonal_lattice(float lattice_vector[3][3],
                           float target_volume,
                           float norm_std)
{
    float ax, by, cz;
    gen_principal_comps(&ax, &by, &cz, target_volume, norm_std, TWO_EQUAL);

    float gamma = 120 * PI/180;
    by = ax * sin(gamma);
    float bx = ax * cos(gamma);

    set_lattice_vectors_zero(lattice_vector);
    lattice_vector[0][0] = ax;
    lattice_vector[1][1] = by;
    lattice_vector[2][2] = cz;
    lattice_vector[1][0] = bx;

    return;
}

void gen_cubic_lattice(float lattice_vector[3][3],
                       float target_volume)
{
    float a = cbrt(target_volume);

    set_lattice_vectors_zero(lattice_vector);
    lattice_vector[0][0] = a;
    lattice_vector[1][1] = a;
    lattice_vector[2][2] = a;

    return;
}

//create a large volume lattice for testing compatiility
void generate_fake_lattice(float lattice_vector[3][3], int spg)
{
    const float ax = 15;
    set_lattice_vectors_zero(lattice_vector);

    if(spg < 1 || spg > 230)
    {
        printf("***ERROR: generate_lattice: spg out of bounds***");
        exit(EXIT_FAILURE);
    }

    else if (spg <= 2)
    {
        lattice_vector[0][0] = ax;
        lattice_vector[1][1] = 0.8*ax;
        lattice_vector[2][2] = 0.5*ax;

        lattice_vector[1][0] = 0.8*ax*tan(10*PI/180);
        lattice_vector[2][0] = 0.5*ax*cos(85*PI/180);
        lattice_vector[2][1] = 0.5*ax*cos(70*PI/180);
    }

    else if ( spg <= 15 )
    {
        lattice_vector[0][0] = ax;
        lattice_vector[1][1] = 0.8*ax;
        lattice_vector[2][2] = 0.5*ax;
        lattice_vector[2][0] = 0.5*ax/tan(70*PI/180);

    }

    else if (spg <= 74)
    {
        lattice_vector[0][0] = ax;
        lattice_vector[1][1] = 0.8*ax;
        lattice_vector[2][2] = 0.5*ax;
    }

    else if (spg <= 142)
    {
        lattice_vector[0][0] = ax;
        lattice_vector[1][1] = ax;
        lattice_vector[2][2] = 0.5*ax;
    }

    else if (spg <= 194)
    {
        lattice_vector[0][0] = ax;
        lattice_vector[1][1] = ax*sin(120*PI/180);
        lattice_vector[2][2] = 0.5*ax;
        lattice_vector[1][0] = ax*cos(120*PI/180);
    }

    else if (spg <= 230)
    {
        lattice_vector[0][0] = ax;
        lattice_vector[1][1] = ax;
        lattice_vector[2][2] = ax;
    }
}

static inline float fmodulo(float n, float d)
{
    long q =n/d;
    float r = n - q*d ;
    return r;
}

void standardise_lattice( float lattice[3][3], int spg)
{
/*
    double nlattice[9] = {lattice[0][0], lattice[0][1], lattice[0][2],
                         lattice[1][0], lattice[1][1], lattice[1][2],
                         lattice[2][0], lattice[2][1], lattice[2][2]
                        };

    int status = niggli_reduce(nlattice, 0.0001);

    if(status)
    {
        lattice[0][0] = nlattice[0];
        lattice[0][1] = nlattice[1];
        lattice[0][2] = nlattice[2];

        lattice[1][0] = nlattice[3];
        lattice[1][1] = nlattice[4];
        lattice[1][2] = nlattice[5];

        lattice[2][0] = nlattice[6];
        lattice[2][1] = nlattice[7];
        lattice[2][2] = nlattice[8];
    }
*/
    //triclinic
    if (spg == 1 || spg == 2)
    {
        float by = lattice[1][1];
        float cy = lattice[2][1];
        float cx = lattice[2][0];
        float bx = lattice[1][0];
        long q = cy/by;
        float cy_new = cy - q * by;
        float cx_new = cx - q * bx;
        lattice[2][1] = cy_new;
        lattice[2][0] = cx_new;

        bx = lattice[1][0];
        float ax = lattice[0][0];
        float bx_new = fmodulo (bx, ax);
        lattice[1][0] = bx_new;

        cx = lattice[2][0];
        cx_new = fmodulo(cx, ax);
        lattice[2][0] = cx_new;
    }

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

    else if (spg <= 0 || spg > 230)
        printf("***ERROR: lattice_generator: standardise_lattice invalid spg");

}

static float gen_angle(float angle_std)
{
    float ang;
    do{
        ang = normal_dist_ab(PI/2, angle_std*PI/180);
    }
    while(ang > MAX_ANGLE || ang < MIN_ANGLE);

    return ang;
}

static void set_lattice_vectors_zero(float lattice_vector[3][3])
{
    lattice_vector[0][0] = 0;
    lattice_vector[0][1] = 0;
    lattice_vector[0][2] = 0;

    lattice_vector[1][0] = 0;
    lattice_vector[1][1] = 0;
    lattice_vector[1][2] = 0;

    lattice_vector[2][0] = 0;
    lattice_vector[2][1] = 0;
    lattice_vector[2][2] = 0;
    return;
}

static void gen_principal_comps(float *ax,
                                float *by,
                                float *cz,
                                float target_volume,
                                float norm_std,
                                int type)
{
    float x, y, z;
    do{
        do {x = normal_dist_ab(1, norm_std); } while(x < EPS);

        if (type == DISTINCT)
            do {y = normal_dist_ab(1, norm_std); } while(y < EPS);
        else
            y = x;

        do {z = normal_dist_ab(1, norm_std); } while(z < EPS);

        float factor =  cbrt(target_volume / (x*y*z) );
        *ax = x*factor;
        *by = y*factor;
        *cz = z*factor;
    }
    while(*ax < LOWB ||*by < LOWB ||*cz < LOWB);

    return;
}











