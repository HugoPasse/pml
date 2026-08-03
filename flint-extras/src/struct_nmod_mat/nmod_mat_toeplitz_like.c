#include "struct_nmod_mat.h"

void nmod_toeplitz_like_init(nmod_toeplitz_like_t mat,
                             ulong nrows,
                             ulong ncols,
                             ulong N){

}

void nmod_toeplitz_like_set(nmod_toeplitz_like_t mat,
                             ulong ngens,
                             lower_triangular_nmod_toeplitz_t* lower_generators,
                             upper_triangular_nmod_toeplitz_t* upper_generators){

}

void nmod_toeplitz_like_dense(nmod_toeplitz_like_t mat, 
                               nmod_mat_t dense_mat){

}

ulong nmod_toeplitz_like_t_nrows(nmod_toeplitz_like_t mat){
    return 0;
}

ulong nmod_toeplitz_like_t_ncols(nmod_toeplitz_like_t mat){
    return 0;
}

ulong nmod_toeplitz_like_t_mod(nmod_toeplitz_like_t mat){
    return 0;
}

ulong nmod_toeplitz_like_get_entry(nmod_toeplitz_like_t mat,
                                    ulong i,
                                    ulong j){
    return 0;
}

/* Arithmetic */
void nmod_toeplitz_like_add(nmod_toeplitz_like_t a,
                             nmod_toeplitz_like_t b, 
                             nmod_toeplitz_like_t res){

}

