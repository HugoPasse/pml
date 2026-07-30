#include "struct_nmod_mat.h"

#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>

void struct_nmod_toeplitz_init(struct_nmod_toeplitz_t mat, 
                                ulong nrows, 
                                ulong ncols,
                                ulong N){
    mat->nrows = nrows;
    mat->ncols = ncols;
    mat->data = _nmod_vec_init(nrows+ncols-1); 
     
    mat->mod.n = N;
    mat->mod.ninv = n_preinvert_limb(N);
    mat->mod.norm = flint_clz(N);
}

void struct_nmod_toeplitz_set(struct_nmod_toeplitz_t mat,
                               nn_ptr data){
    _nmod_vec_set(mat->data, data, mat->nrows+mat->ncols-1);
    // Reduce mat->data mod n 
}

void struct_nmod_toeplitz_init_set(struct_nmod_toeplitz_t mat,
                                    ulong nrows,
                                    ulong ncols,
                                    ulong N,
                                    nn_ptr data){
    struct_nmod_toeplitz_init(mat, nrows, ncols, N);
    struct_nmod_toeplitz_set(mat, data);
}

ulong struct_nmod_toeplitz_nrows(struct_nmod_toeplitz_t mat){
    return mat->nrows;
}

ulong struct_nmod_toeplitz_ncols(struct_nmod_toeplitz_t mat){
    return mat->ncols;
}

void struct_nmod_toeplitz_dense(struct_nmod_toeplitz_t mat,
                                  nmod_mat_t dense_mat){

    for(int i=0; i<struct_nmod_toeplitz_nrows(mat); i++){
        for(int j=0; j<struct_nmod_toeplitz_ncols(mat); j++){
            int k = struct_nmod_toeplitz_ncols(mat) - 1 + i - j;
            nmod_mat_set_entry(dense_mat, i, j, mat->data[k]);
        }
    }  
}

void struct_nmod_toeplitz_mul_vec(struct_nmod_toeplitz_t mat,
                                    nn_ptr v,
                                    nn_ptr res){
    _nmod_poly_mulmid(res, mat->data, mat->nrows+mat->ncols-1, v, mat->ncols, mat->ncols-1, mat->nrows+mat->ncols-1, mat->mod);  
}
