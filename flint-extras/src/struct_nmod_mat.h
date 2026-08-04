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

/* ---------------------------
 * Toeplitz matrices stored as
 *     a3 a2 a1 a0
 *     a4 a3 a2 a1
 *     a5 a4 a3 a2
 -------------------------- */
typedef struct
{
    ulong nrows,ncols; // Dimensions
    ulong n; // Modulus TODO Remove this is useless because it is stored in "mod"
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


// TODO ADD GETTERS
/* ---------------------------
 * Lower triangular Toeplitz matrices 
 *     a0 0  0  0
 *     a1 a0 0  0
 *     a2 a1 a0 0
 * -------------------------- */
typedef struct
{
    ulong nrows,ncols; // Dimensions
    nn_ptr data; // Entries [a_0,...,a_5]
    nmod_t mod;
} nmod_mat_lower_toeplitz_struct;

typedef nmod_mat_lower_toeplitz_struct nmod_mat_lower_toeplitz_t[1];

void nmod_mat_lower_toeplitz_init(nmod_mat_lower_toeplitz_t mat,
                                   ulong nrows,
                                   ulong ncols,
                                   ulong N);

void nmod_mat_lower_toeplitz_set(nmod_mat_lower_toeplitz_t mat,
                                  nn_ptr data);

void nmod_mat_lower_toeplitz_init_set(nmod_mat_lower_toeplitz_t mat,
                                       ulong nrows,
                                       ulong ncols,
                                       ulong N,
                                       nn_ptr data);

void nmod_mat_lower_toeplitz_dense(nmod_mat_lower_toeplitz_t mat,
                                    nmod_mat_t res);

// TODO ADD GETTERS
/* ---------------------------
 * Upper triangular Toeplitz matrices 
 *     a3 a2 a1 a0
 *     0  a3 a2 a1
 *     0  0  a3 a2 
 * -------------------------- */
typedef struct
{
    ulong nrows,ncols; // Dimensions
    nn_ptr data; // Entries [a_0,...,a_5]
    nmod_t mod;
} nmod_mat_upper_toeplitz_struct;

typedef nmod_mat_upper_toeplitz_struct nmod_mat_upper_toeplitz_t[1];

void nmod_mat_upper_toeplitz_init(nmod_mat_upper_toeplitz_t mat,
                                   ulong nrows,
                                   ulong ncols,
                                   ulong N);

void nmod_mat_upper_toeplitz_set(nmod_mat_upper_toeplitz_t mat,
                                  nn_ptr data);

void nmod_mat_upper_toeplitz_init_set(nmod_mat_upper_toeplitz_t mat,
                                       ulong nrows,
                                       ulong ncols,
                                       ulong N,
                                       nn_ptr data);

void nmod_mat_upper_toeplitz_dense(nmod_mat_upper_toeplitz_t mat,
                                    nmod_mat_t res);

// This function return the entry (i,j) of the product lmat * umat
ulong nmod_mat_lower_upper_toeplitz_pdt(nmod_mat_lower_toeplitz_t lmat,
                                         nmod_mat_upper_toeplitz_t umat,
                                         ulong i,
                                         ulong j); 

/* ---------------------------
 * Toeplitz-like matrices
 * --------------------------- */
typedef struct
{
    ulong nrows,ncols; // Dimensions
    ulong ngens; // Number of generators
    nmod_mat_lower_toeplitz_t* lower_gens; // Matrices L(g_i)
    nmod_mat_upper_toeplitz_t* upper_gens; // Matrices U(h_i)
    nmod_t mod;
} nmod_toeplitz_like_struct;

typedef nmod_toeplitz_like_struct nmod_toeplitz_like_t[1];

/* Initialization */
void nmod_toeplitz_like_init(nmod_toeplitz_like_t mat,
                             ulong nrows,
                             ulong ncols,
                             ulong N);

void nmod_toeplitz_like_set(nmod_toeplitz_like_t mat,
                             ulong ngens,
                             nmod_mat_lower_toeplitz_t* lower_generators,
                             nmod_mat_upper_toeplitz_t* upper_generators);

// TODO is this relevant (too many arguments ?)
//void nmod_toeplitz_like_init_set();
void nmod_toeplitz_like_dense(nmod_toeplitz_like_t mat, 
                               nmod_mat_t dense_mat);

/* Getters */
// TODO DO THESE INLINE ?
ulong nmod_toeplitz_like_t_nrows(nmod_toeplitz_like_t mat);
ulong nmod_toeplitz_like_t_ncols(nmod_toeplitz_like_t mat);
ulong nmod_toeplitz_like_t_mod(nmod_toeplitz_like_t mat);

ulong nmod_toeplitz_like_get_entry(nmod_toeplitz_like_t mat,
                                    ulong i,
                                    ulong j);

/* Arithmetic */
void nmod_toeplitz_like_add(nmod_toeplitz_like_t a,
                             nmod_toeplitz_like_t b, 
                             nmod_toeplitz_like_t res);


#endif 
