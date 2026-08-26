#include "nmod_mat_struct.h"

void nmod_mat_upper_toeplitz_init(nmod_mat_upper_toeplitz_t mat,
                                   ulong nrows,
                                   ulong ncols,
                                   ulong N){
    mat->nrows = nrows;
    mat->ncols = ncols;
    mat->data = _nmod_vec_init(ncols);
    nmod_init(&mat->mod,N); 
}

void nmod_mat_upper_toeplitz_set(nmod_mat_upper_toeplitz_t mat,
                                  nn_ptr data){
    _nmod_vec_set(mat->data, data, mat->ncols);
}

void nmod_mat_upper_toeplitz_init_set(nmod_mat_upper_toeplitz_t mat,
                                       ulong nrows,
                                       ulong ncols,
                                       ulong N,
                                       nn_ptr data){
    nmod_mat_upper_toeplitz_init(mat,nrows,ncols,N);
    nmod_mat_upper_toeplitz_set(mat,data);
}

void nmod_mat_upper_toeplitz_dense(nmod_mat_upper_toeplitz_t mat,
                                    nmod_mat_t res){
    nmod_mat_init(res, mat->nrows, mat->ncols, mat->mod.n);
    nmod_mat_zero(res);
    for(int i=0; i<mat->nrows; i++){
        for(int j=i; j<mat->ncols; j++){
            int k = mat->ncols - 1 + i - j;
            nmod_mat_set_entry(res, i, j, mat->data[k]);
        }
    }
}

void nmod_mat_upper_toeplitz_as_poly(nmod_mat_upper_toeplitz_t mat,
                                      nmod_poly_t pol){
    nmod_poly_init2(pol, mat->mod.n, mat->ncols);
    for(int i=0; i<mat->ncols; i++){
        // For upper Toeplitz, entries are strored backwards
        nmod_poly_set_coeff_ui(pol, i, mat->data[mat->ncols-1-i]);
    }
}

void nmod_mat_upper_toeplitz_mul_nmod_vec(nmod_mat_upper_toeplitz_t mat,
                                           nn_ptr v,
                                           nn_ptr res){
    nmod_poly_t pol;
    nmod_mat_upper_toeplitz_as_poly(mat, pol);
        
    nmod_poly_t vx;
    nmod_poly_init(vx, mat->mod.n);
    for(int i=0; i<mat->ncols; i++){
        nmod_poly_set_coeff_ui(vx,i,mat->data[i]);
    }
    int k = FLINT_MAX(mat->ncols - mat->nrows, 0);
    nmod_poly_mulhigh(pol,pol,vx,k);
     
    for(int i=0; i<mat->nrows; i++){
        res[i] = nmod_poly_get_coeff_ui(pol,k+i);
    }
    nmod_poly_clear(pol);
    nmod_poly_clear(vx);
}
