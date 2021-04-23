#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "read_input.h"
#include "spg_generation.h"
#include "crystal_utils.h"
#include "check_structure.h"
#include "algebra.h"

#define CONSTANT_TOL 0.001
#define MAXVDW 2.1

#define NORM2(a,b) ( (a[0]-b[0])*(a[0]-b[0]) +\
                     (a[1]-b[1])*(a[1]-b[1]) +\
                     (a[2]-b[2])*(a[2]-b[2]) )

#define NORMSQR(a) ( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] )
#define SQR(x) (x*x)

static void convert_atom2atom_vdw(char *atom,float *atom_vdw, int num_atoms);
static int fast_screener_vdw(crystal xtal, float *vdw_matrix);
static int check_pairwise_if_mol_close( float *vdw_matrix, int total_atoms,
    xtal_molecule_pair *mol_pair, float max_dist);
static int check_pairwise(float *vdw_matrix , int total_atoms,
    xtal_molecule_pair *mol_pair, float l_disp[3]);
static int check_lower_triangular(float L[3][3]);
static inline int int_floor(float x);
static inline int xseq_generator(int reset);
static inline int yseq_generator(int reset);
static inline int zseq_generator(int reset);
static inline int check_lower_triangular(float L[3][3]);
static inline float maxof(float a, float b);

// Public Functions

void create_vdw_matrix_from_sr( molecule *mol,
                                float *vdw_matrix,
                                float sr,
                                int Z)
{   int num_atoms_in_molecule = mol->num_of_atoms;
    float *atom_vdw_vector = (float *) malloc( num_atoms_in_molecule*
                            sizeof(float) );
    //create vector with vdw radii of size num of atoms in mol
    convert_atom2atom_vdw( mol->atoms, atom_vdw_vector,
                            num_atoms_in_molecule);

    //now create matrix
    int dim_vdw_matrix = num_atoms_in_molecule * Z ;
    for(int i = 0; i < dim_vdw_matrix; i++)
        for(int j = 0; j < dim_vdw_matrix; j++)
            {
                *(vdw_matrix + i*dim_vdw_matrix +j)
                    =  (*(atom_vdw_vector + i%num_atoms_in_molecule) +
                        *(atom_vdw_vector + j%num_atoms_in_molecule) ) *
                        sr;
            }
    free(atom_vdw_vector);
}


int check_structure_with_vdw_matrix(crystal xtal,
    float *vdw_matrix,
    int dim1,
    int dim2)
{
    int total_atoms = xtal.Z * xtal.num_atoms_in_molecule;
    int mol_id[xtal.Z];

    if( !fast_screener_vdw(xtal, vdw_matrix) )
        return 0;

    mol_id[0] = 0;
    for(int i = 1; i < xtal.Z; i++)
    {
        mol_id[i] = mol_id[i-1] + xtal.num_atoms_in_molecule;
        //printf("mol_id = %d\n", mol_id[i]);
    }

    return structure_checker(
                             &xtal,
                             vdw_matrix,
                             total_atoms,
                             mol_id,
                             xtal.Z);

}

/*
 * Assume lattice to be lower triangular in row major form
 * More generalized function designed for co-crystals
 */
int structure_checker(crystal *xtal,
    float *vdw_matrix,
    int total_atoms,
    int *mol_id,
    int num_mols
    )
{
    //preliminary checks
    /*
    if(dim1 != dim2 || dim2 != dim3 || dim3 != dim1)
    {
        printf("***ERROR: length of coordinate arrays don't match\n");
        return -1;
    }

    if(dim4 != dim5)
    {
        printf("***ERROR: vdw cutoff matrix should be a square matrix\n");
        return -1;
    }

    if(!check_lower_triangular(L))
    {
        printf("***ERROR: lattice_vectors should be a lower triangular matrix in row major form\n");
        return -1;
    }
    */

    // Unpack
    float *X = xtal->Xcord;
    float *Y = xtal->Ycord;
    float *Z = xtal->Zcord;

    //find number of atoms in each molecule
    int num_atoms_in_molecule[num_mols];
    for(int i = 0; i < num_mols -1; i++)
    {
        num_atoms_in_molecule[i] = mol_id[i+1] - mol_id[i];
        //printf("%s = %d\n", "num_atoms_molecule = ", num_atoms_in_molecule[i]);
    }
    //for last mol:
    num_atoms_in_molecule[num_mols-1] = total_atoms - mol_id[num_mols-1];

    // estimate molecule lengths
    float mol_len[num_mols];
    for(int i = 0 ; i < num_mols; i++)
    {
        mol_len[i] = find_mol_len(X + mol_id[i],
                                  Y + mol_id[i],
                                  Z + mol_id[i],
                                  num_atoms_in_molecule[i]);
        //printf("molLen = %f\n", mol_len[i]);
    }

    // find com
    float com[num_mols*3];
    for(int i = 0; i < num_mols; i++)
    {
        //printf("num_atoms_in_molecule = %d\n", num_atoms_in_molecule[i]);
        find_mol_com(X + mol_id[i], Y + mol_id[i], Z + mol_id[i],
            num_atoms_in_molecule[i], com + 3*i);
        //printf("i = %d , com  === %f, %f, %f\n",i, com[3*i], com[3*i +1], com[3*i+2] );
    }

    // Loop over molecules in a unit cell
    for(int i = 0; i < num_mols; i++)
    {
        for(int j = i; j < num_mols; j++)
        {
            //check pair of selcted molecules.
            //printf("i, j = %d, %d \n", i, j );
            float max_dist = mol_len[i] + mol_len[j];

            xtal_molecule_pair mol_pair = {.L = xtal->lattice_vectors, .com1 = com + 3*i,
                .com2 = com + 3*j, .X = X, .Y = Y, .Z = Z, .index1 = mol_id[i],
                .num_atoms1 = num_atoms_in_molecule[i], .index2 = mol_id[j],
                .num_atoms2 = num_atoms_in_molecule[j]};

            int verdict = check_pairwise_if_mol_close(vdw_matrix, total_atoms, &mol_pair, max_dist);

            if (!verdict)
                return 0 ;
        }
    }

    return 1;
}




// Private functions to support the above functions.

/* convert the atoms char array into vdw radii information for structure checking
 * uses Bondii radii. If you need to add/change vdw radii of an atom
 * to the database, do it here.
 */
static void convert_atom2atom_vdw(char *atom,float *atom_vdw, int num_atoms)
{

    for (int i = 0; i < num_atoms; i++)
    {
        if      (atom[2*i] == 'C' && atom[2*i+1] == ' ')
        atom_vdw[i]=1.7;
        else if (atom[2*i] == 'H' && atom[2*i+1] == ' ')
        atom_vdw[i]=1.1;
        else if (atom[2*i] == 'N' && atom[2*i+1] == ' ')
        atom_vdw[i] = 1.55 ;
        else if (atom[2*i] == 'O' && atom[2*i+1] == ' ')
        atom_vdw[i] = 1.52;
        else if (atom[2*i] == 'F' && atom[2*i+1] == ' ')
        atom_vdw[i] = 1.47;
        else if (atom[2*i] == 'P' && atom[2*i+1] == ' ')
        atom_vdw[i] = 1.8;
        else if (atom[2*i] == 'S' && atom[2*i+1] == ' ')
        atom_vdw[i] = 1.8;
        else if (atom[2*i] == 'C' && atom[2*i+1] == 'l')
        atom_vdw[i] = 1.75;
        else if (atom[2*i] == 'B' && atom[2*i+1] == 'r')
        atom_vdw[i] = 1.85;
        else if (atom[2*i] == 'I' && atom[2*i+1] == ' ')
        atom_vdw[i] = 1.98;
        else if (atom[2*i] == 'B' && atom[2*i+1] == ' ')
        atom_vdw[i] = 1.92;
        else if (atom[2*i] == 'H' && atom[2*i+1] == 'e')
        atom_vdw[i] = 1.40;
        else if (atom[2*i] == 'N' && atom[2*i+1] == 'e')
        atom_vdw[i] = 1.54;
        else if (atom[2*i] == 'K' && atom[2*i+1] == 'r')
        atom_vdw[i] = 2.02;
        else if (atom[2*i] == 'S' && atom[2*i+1] == 'i')
        atom_vdw[i] = 2.10;
        else if (atom[2*i] == 'Z' && atom[2*i+1] == 'n')
        atom_vdw[i] = 1.39;
        else
        {
            printf("***ERROR: atom2atom_vdw: atom not found -> %c%c\n",\
                    atom[2*i], atom[2*i+1]);
            exit(EXIT_FAILURE);
        }
    }
}


// Precheck if intermolecular distances are too close.
// Doesn't take into account of periodic boundary condition, but is
// much faster than rigourous checking.
// Doesn't check molecule with its own periodic image
//
static int fast_screener_vdw(crystal xtal, float *vdw_matrix)
{
                                    //number of atoms in a molecule
    int N = xtal.num_atoms_in_molecule;
    int m = xtal.Z;                 //number of molecules in a unit cell;
                                // numberof atom in a unit cell = N*m
    int total_atoms = N*m;
    static float mol_len = 0;
    static int first_time = 1;

    if (first_time)
    {
        first_time = 0;
        mol_len = find_mol_len(xtal.Xcord, xtal.Ycord, xtal.Zcord, N);
    }

    float small_number = MAXVDW;

    for(int i = 0; i < total_atoms; i += N)
    {
        float com1[3];
        compute_molecule_COM( xtal, com1, i);
        for(int j = i + N; j < total_atoms; j += N)
        {
            //check if the molecule COM are far. if they are, dont
            //bother checking distances
            if (j == i + N)
            {   float com2[3];
                compute_molecule_COM( xtal, com2, j);
                if( sqrt ((com1[0] - com2[0])*(com1[0] - com2[0])+
                          (com1[1] - com2[1])*(com1[1] - com2[1])+
                          (com1[2] - com2[2])*(com1[2] - com2[2]))
                    >     (mol_len + small_number)                  )
                    continue;
            }

            //molecule COM are close than molecule length
            else
            {
                for(int k = i; k < i + N; k++)
                {
                    for (int z = j; z < j + N; z++ )
                    {
                        if(    (xtal.Xcord[k] - xtal.Xcord[z])*
                               (xtal.Xcord[k] - xtal.Xcord[z])+
                               (xtal.Ycord[k] - xtal.Ycord[z])*
                               (xtal.Ycord[k] - xtal.Ycord[z])+
                               (xtal.Zcord[k] - xtal.Zcord[z])*
                               (xtal.Zcord[k] - xtal.Zcord[z])
                               < *(vdw_matrix + total_atoms*k + z) *
                                 *(vdw_matrix + total_atoms*k + z)  )
                            return 0;
                    }
                }
            }
        }
    }

    return 1;
}


// Core sphere decoder
static int check_pairwise_if_mol_close( float *vdw_matrix,
    int total_atoms, xtal_molecule_pair *mol_pair, float max_dist)
{
    // Unpack
    float *com1 = mol_pair->com1;
    float *com2 = mol_pair->com2;
    float (*L)[3] = mol_pair->L;
    int id1 = mol_pair->index1;
    int id2 = mol_pair->index2;

    float p[3] = {com1[0] - com2[0], com1[1] - com2[1], com1[2] - com2[2]} ;
    int zmax = maxof( fabs(p[2] + max_dist), fabs(p[2] - max_dist) )/ fabs(L[2][2]) + 3;
    //printf("zmax = %d\n", zmax );

    for(int k = zseq_generator(1); k < zmax; k = zseq_generator(0))
    {
        // update max_dist
        float loss1 = (p[2] - k*L[2][2])*(p[2] - k*L[2][2]);
        float max_dist_sqr = max_dist*max_dist - loss1;
        if (max_dist_sqr < 0)
            continue;
        float max_dist_z = sqrt(max_dist_sqr);
        //max_dist_z = max_dist;
        // find ymax
        int ymax =  maxof( fabs(p[2] - L[2][1]*k + max_dist_z),
             fabs(p[2] - L[2][1]*k - max_dist_z) ) / fabs(L[1][1]) + 3;
        //printf("ymax = %d\n", ymax );

        for(int j = yseq_generator(1); j < ymax; j = yseq_generator(0))
        {
            // update max_dist
            float loss2 = (p[1] - L[1][1]*j - L[2][1]*k)*(p[1] - L[1][1]*j - L[2][1]*k);
            float max_dist_sqr =  max_dist_z*max_dist_z - loss2;
            if (max_dist_sqr < 0)
                continue;
            float max_dist_y = sqrt(max_dist_sqr);
            //max_dist_y = max_dist;
            //find xmax
            int xmax = maxof( fabs(p[0] - L[1][0]*j - L[2][0]*k + max_dist_y),
                fabs(p[0] - L[1][0]*j - L[2][0]*k - max_dist_y) )/ fabs(L[0][0]) + 3;
            //printf("xmax = %d\n", xmax );

            for(int i = xseq_generator(1); i < xmax; i = xseq_generator(0))
            {
                // Skip checking the molecule with itself
                if (id1 == id2 && i == 0 && j == 0 && k == 0)
                    continue;

                // Lattice displacement
                float l_disp[3] = { i*L[0][0] + j*L[1][0] + k*L[2][0],
                     j*L[1][1] + k*L[2][1], k*L[2][2] };

                float dist[3] = {p[0] - l_disp[0],
                                 p[1] - l_disp[1],
                                 p[2] - l_disp[2] };

                // Two molecules considered are within the search radius then
                // test if the molecules are too close
                if (NORMSQR(dist) < (max_dist/2 + MAXVDW)*(max_dist/2 + MAXVDW))
                {
                    /*
                    printf("i, j, k = %d, %d, %d dist = %f %f %f\n",i, j, k, dist[0], dist[1], dist[2] );
                    printf("max_dist = %f\n", max_dist);
                    printf("dist = %f\n", NORMSQR(dist));
                    printf("close; ijk = %d %d %d\n", i,j,k );
                    */
                    int ret_val = check_pairwise(vdw_matrix, total_atoms, mol_pair, l_disp);

                    if(!ret_val)
                        return 0;

                }
            }
        }
    }
    return 1;
}

static int check_pairwise(float *vdw_matrix , int total_atoms,
    xtal_molecule_pair *mol_pair, float l_disp[3])
{
    // Unpack
    float *X = mol_pair->X;
    float *Y = mol_pair->Y;
    float *Z = mol_pair->Z;
    int N1 = mol_pair->num_atoms1;
    int N2 = mol_pair->num_atoms2;
    int id1 = mol_pair->index1;
    int id2 = mol_pair->index2;

    //printf("num_atoms1 = %d num_atoms2 = %d\n", N1, N2);

    for(int i1 = 0; i1 < N1; i1++)
        for(int i2 = 0; i2 < N2; i2++)
        {
            // index of an atom in molecule 1 in cord arrays
            // id1 is the index of molecule 1
            int atom1 = i1 + id1;
            int atom2 = i2 + id2;

            float dist[3] = {X[atom1] - X[atom2] - l_disp[0],
                             Y[atom1] - Y[atom2] - l_disp[1],
                             Z[atom1] - Z[atom2] - l_disp[2]
                            };


            if (NORMSQR(dist) < SQR(*(vdw_matrix + atom1*total_atoms + atom2)) )
                return 0;
        }

    return 1;
}


/*
 * Get the length of molecule
 * X, Y, Z are co-ordinates.
 * len is total atoms in the molecule.
 * returns twice the distance of farthest atom from COM
 */
float find_mol_len(float *X, float *Y, float *Z, int len)
{
    float com[3] = {0, 0, 0};
    find_mol_com(X, Y, Z, len, com);
    //printf("com = %f %f %f\n", com[0], com[1], com[2]);
    float first_atom[3] = {X[0], Y[0], Z[0]};
    float max = NORM2(com, first_atom);
    float dist_sqr = 0;

    for(int i = 0; i < len; i++)
    {
        float atom_vec[3] = {X[i], Y[i], Z[i]};
        dist_sqr = NORM2(atom_vec, com);

        if (max < dist_sqr)
            max = dist_sqr;
    }

    return 2*sqrt(max);
}



void find_mol_com(float *X, float *Y, float *Z, int len, float com[3])
{
    com[0] = 0;
    com[1] = 0;
    com[2] = 0;

    for(int i = 0; i < len; i++)
    {
        com[0] += X[i];
        com[1] += Y[i];
        com[2] += Z[i];
    }

    com[0] /= len;
    com[1] /= len;
    com[2] /= len;
}

/* returns the distance between
 * two points under periodic boundary condition. T and T_inv are lattice
 * vector matrix and its inverse. T is in row major form. x1,x2,x3 are
 * the first point and y's are second in cartesian
 * returns the shortest distance as a float.
 * UNUSED(?) BUT USEFUL FUNCTION
 */
float pdist(float T[3][3],
            float T_inv[3][3],
            float x1,
            float x2,
            float x3,
            float y1,
            float y2,
            float y3  )
{
    float p_dist = pdist_appx(T,
            T_inv,
            x1,
            x2,
            x3,
            y1,
            y2,
            y3  );

    float cartesian_distance[3] = {x1 - y1, x2 - y2, x3 - y3};
    float fractional_distance[3];
    for (int i = 0; i < 3; i++)
    {
        float sum = 0;
        for (int j = 0; j < 3; j++)
            sum = sum + T_inv[j][i] * cartesian_distance[j];

        float frac_part = sum - ((long)sum); //if -ve add one
        if (frac_part < 0)
            fractional_distance[i] = frac_part + 1;
        else
            fractional_distance[i] = frac_part;
    }

    //find cartesian distance vector in the bounding box
    //reduced cartesian distance vector
    float red_cart_distance[3];
    for(int i = 0; i < 3; i++)
        red_cart_distance[i] = T[0][i] * fractional_distance[0] +
                               T[1][i] * fractional_distance[1] +
                               T[2][i] * fractional_distance[2];

    //use pdist as radius of sphere are search
    int search_radius = 4;
    int limx = search_radius;
    int limy = search_radius;
    int limz = search_radius;

    float test_dist[3] = {0,0,0};

    for (int i = -limx; i <= limx; i++)
    for (int j = -limy; j <= limy; j++)
    for (int k = -limz; k <= limz; k++)
    {
        test_dist[0] = i*T[0][0] + j*T[1][0] + k*T[2][0] - red_cart_distance[0];
        test_dist[1] = j*T[1][1] + k*T[2][1] - red_cart_distance[1];
        test_dist[2] = k*T[2][2] - red_cart_distance[2];

        if (vector3_norm(test_dist) < p_dist )
            p_dist = vector3_norm(test_dist);
    }

    return p_dist;
}


// For use with some functions
//aproximate version of pdist() used by an earlier version

float pdist_appx(float T[3][3],
            float T_inv[3][3],
            float x1,
            float x2,
            float x3,
            float y1,
            float y2,
            float y3  )
{
    //intialising variables
    float p_dist = 0;

    static float Q[8][3]={ {0,0,0},{1,0,0},{0,1,0},{0,0,1},
                           {1,1,0},{0,1,1},{1,0,1},{1,1,1}  };

    float cartesian_distance[3] = {x1 - y1, x2 - y2, x3 - y3};
    //vector3_subtract(x,y,cartesian_distance);

    float fractional_distance[3];
    for (int i = 0; i < 3; i++)
    {
        float sum = 0;
        for (int j = 0; j < 3; j++)
            sum = sum + T_inv[j][i] * cartesian_distance[j];

        float frac_part = sum - ((long)sum); //if -ve add one
        if (frac_part < 0)
            fractional_distance[i] = frac_part + 1;
        else
            fractional_distance[i] = frac_part;
    }

    for (int z = 0; z < 8; z++)
    {
        float A[3] = {Q[z][0], Q[z][1], Q[z][2]};
        float dist_corner[3];

        dist_corner[0] = fractional_distance[0] - A[0];
        dist_corner[1] = fractional_distance[1] - A[1];
        dist_corner[2] = fractional_distance[2] - A[2];

        //distance vector to 8 corners in cartesian.
        float dist_z = 0;
        for (int i = 0; i < 3; i++)
        {
            float sum = 0;
            for (int j = 0; j < 3 ;j++)
                sum = sum + dist_corner[j] * T[j][i];
            dist_z += sum*sum;
        }
        // length odistnce vector to the zth corner
        // for finding the minimum distance (min dist_z)
        if (z == 0)
            p_dist = dist_z;
        else if (p_dist > dist_z)
            p_dist = dist_z;
    }
    return sqrt(p_dist);
}

// Inline functions
static inline int int_floor(float x)
{
  int i = (int)x; /* truncate */
  return i - ( i > x ); /* convert trunc to floor */
}

static inline int xseq_generator(int reset)
{
    static int count = 0;
    if (reset)
    {
        count = 0;
        return 0;
    }

    if (count > 0)
    {
        count = -count;
    }
    else
    {
        count = -count + 1;
    }

    return count;
}

static inline int yseq_generator(int reset)
{
    static int count = 0;
    if (reset)
    {
        count = 0;
        return 0;
    }

    if (count > 0)
    {
        count = -count;
    }
    else
    {
        count = -count + 1;
    }

    return count;
}

static inline int zseq_generator(int reset)
{
    static int count = 0;
    if (reset)
    {
        count = 0;
        return 0;
    }

    if (count > 0)
    {
        count = -count;
    }
    else
    {
        count = -count + 1;
    }

    return count;
}

static inline int check_lower_triangular(float L[3][3])
{
    if (L[0][1] < CONSTANT_TOL &&
        L[0][2] < CONSTANT_TOL &&
        L[1][2] < CONSTANT_TOL
        )
        return 1;
    else
        return 0;
}

static inline float maxof(float a, float b)
{
    if(a > b)
        return a;
    else
        return b;
}
