#ifndef __STRUCT_NMOD_MAT__H
#define __STRUCT_NMOD_MAT__H

#include "pml.h"
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>

#ifdef STRUCT_NMOD_MAT_INLINES_C
#define STRUCT_NMOD_MAT_INLINE
#else
#define STRUCT_NMOD_MAT_INLINE static inline
#endif

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

STRUCT_NMOD_MAT_INLINE
ulong *struct_nmod_toeplitz_get_entry(struct_nmod_toeplitz_t mat, ulong i, ulong j){
    return &mat->data[mat->ncols - 1 + i - j];
}

/* Computes the dense representation of a Toeplitz matrix */
void struct_nmod_toeplitz_dense(struct_nmod_toeplitz_t mat, 
                                  nmod_mat_t dense_mat);

/* Computes the sum of two Toeplitz matrices */
void struct_nmod_toeplitz_add(struct_nmod_toeplitz_t a,
                               struct_nmod_toeplitz_t b, 
                               struct_nmod_toeplitz_t res);

/* Computes the matrix-vector product with a Toeplitz matrix */
void struct_nmod_toeplitz_right_mul_vec(struct_nmod_toeplitz_t mat,
                                    nn_ptr v,
                                    nn_ptr res);

/* Computes the matrix-matrix product between a Toeplitz matrix and a dense matrix */
void struct_nmod_toeplitz_right_mul_mat(struct_nmod_toeplitz_t mat,
                                    nmod_mat_t b,
                                    nmod_mat_t res);

/* Computes the vector-matrix product with a Toeplitz matrix */
void struct_nmod_toeplitz_left_mul_vec(struct_nmod_toeplitz_t mat,
                                    nn_ptr v,
                                    nn_ptr res);

/* Computes the matrix-matrix product between a dense matrix and a Toeplitz matrix */
void struct_nmod_toeplitz_left_mul_mat(struct_nmod_toeplitz_t mat,
                                    nmod_mat_t b,
                                    nmod_mat_t res);
#endif 
