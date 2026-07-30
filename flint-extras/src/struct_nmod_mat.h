#ifndef __STRUCT_NMOD_MAT__H
#define __STRUCT_NMOD_MAT__H

#include "pml.h"
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>

/* Toeplitz matrices stored as
 *     a3 a2 a1 a0
 *     a4 a3 a2 a1
 *     a5 a4 a3 a2
 */
typedef struct
{
    ulong nrows,ncols; // Dimensions
    ulong n; // Modulus
    nn_ptr data; // Entries [a_0,...,a_5]
    nmod_t mod;
} struct_nmod_toeplitz_struct;

typedef struct_nmod_toeplitz_struct struct_nmod_toeplitz_t[1];

/* Initialization */
void struct_nmod_toeplitz_init(struct_nmod_toeplitz_t mat,
                                ulong nrows, 
                                ulong ncols,
                                ulong N);

void struct_nmod_toeplitz_set(struct_nmod_toeplitz_t mat,
                               nn_ptr data);

void struct_nmod_toeplitz_init_set(struct_nmod_toeplitz_t mat,
                                    ulong nrows, 
                                    ulong ncols,
                                    ulong N,
                                    nn_ptr data);
    
/* Getters */
ulong struct_nmod_toeplitz_nrows(struct_nmod_toeplitz_t mat);
ulong struct_nmod_toeplitz_ncols(struct_nmod_toeplitz_t mat);

/* Computes the dense representation of a Toeplitz matrix */
void struct_nmod_toeplitz_dense(struct_nmod_toeplitz_t mat, 
                                  nmod_mat_t dense_mat);

/* Computes the matrix-vector product with a Toeplitz matrix */
void struct_nmod_toeplitz_mul_vec(struct_nmod_toeplitz_t mat,
                                    nn_ptr v,
                                    nn_ptr res);

#endif 
