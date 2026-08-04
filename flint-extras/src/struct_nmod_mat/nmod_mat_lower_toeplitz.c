#include "struct_nmod_mat.h"

void nmod_mat_lower_toeplitz_init(nmod_mat_lower_toeplitz_t mat,
                                   ulong nrows,
                                   ulong ncols,
                                   ulong N){
    mat->nrows = nrows;
    mat->ncols = ncols;
    mat->data = _nmod_vec_init(nrows);
    nmod_init(&mat->mod,N); 
}

void nmod_mat_lower_toeplitz_set(nmod_mat_lower_toeplitz_t mat,
                                  nn_ptr data){
    _nmod_vec_set(mat->data, data, mat->nrows);
}

void nmod_mat_lower_toeplitz_init_set(nmod_mat_lower_toeplitz_t mat,
                                       ulong nrows,
                                       ulong ncols,
                                       ulong N,
                                       nn_ptr data){
    nmod_mat_lower_toeplitz_init(mat,nrows,ncols,N);
    nmod_mat_lower_toeplitz_set(mat,data);
}

void nmod_mat_lower_toeplitz_dense(nmod_mat_lower_toeplitz_t mat,
                                    nmod_mat_t res){
    nmod_mat_init(res, mat->nrows, mat->ncols, mat->mod.n);
    for(int j=0; j<mat->ncols; j++){
        for(int i=0; i<=j; i++){
            int k = mat->ncols - 1 + i - j;
            nmod_mat_set_entry(res, i, j, mat->data[k]);
        }
    }
}
