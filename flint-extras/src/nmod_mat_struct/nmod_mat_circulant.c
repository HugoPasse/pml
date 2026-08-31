#include "nmod_mat_struct.h"
#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>

void nmod_mat_circulant_init(nmod_mat_circulant_t mat,
                             ulong nrows,
                             ulong ncols,
                             ulong N){
    mat->nrows = nrows; 
    mat->ncols = ncols; 
    mat->data = _nmod_vec_init(ncols);

    mat->mod.n = N;
    mat->mod.ninv = n_preinvert_limb(N);
    mat->mod.norm = flint_clz(N);
}

void nmod_mat_circulant_set(nmod_mat_circulant_t mat,
                            nn_ptr data){
   _nmod_vec_set(mat->data, data, mat->ncols); 
}

void nmod_mat_circulant_init_set(nmod_mat_circulant_t mat,
                             ulong nrows,
                             ulong ncols,
                             nn_ptr data,
                             ulong N){
    nmod_mat_circulant_init(mat, nrows, ncols, N);
    nmod_mat_circulant_set(mat, data);
}

void nmod_mat_circulant_clear(nmod_mat_circulant_t mat){
    if(mat->data){
        flint_free(mat->data);
    }
}

void nmod_mat_circulant_dense(nmod_mat_circulant_t mat,
                                 nmod_mat_t res){
    for(int i=0; i<mat->nrows; i++){
        for(int j=0; j<mat->ncols; j++){
            int k = (mat->nrows * mat->ncols + j-i) % mat->ncols;
            nmod_mat_set_entry(res, i, j, mat->data[k]);
        }
    }
}
