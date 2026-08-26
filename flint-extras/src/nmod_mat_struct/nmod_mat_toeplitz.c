#include "nmod_mat_struct.h"

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

// ---------------------------- //
// ----------- SUM ------------ //
// ---------------------------- //

void struct_nmod_toeplitz_add(struct_nmod_toeplitz_t a,
                               struct_nmod_toeplitz_t b, 
                               struct_nmod_toeplitz_t res){
    _nmod_vec_add(res->data, a->data, b->data, b->nrows+b->ncols-1, a->mod);
}

// ---------------------------- //
// ------ MULTIPLICATION ------ //
// ---------------------------- //
void struct_nmod_toeplitz_right_mul_vec(struct_nmod_toeplitz_t mat,
                                    nn_ptr v,
                                    nn_ptr res){
    _nmod_poly_mulmid(res, mat->data, mat->nrows+mat->ncols-1, v, mat->ncols, mat->ncols-1, mat->nrows+mat->ncols-1, mat->mod);  
}

void struct_nmod_toeplitz_right_mul_mat(struct_nmod_toeplitz_t mat,
                                    nmod_mat_t b,
                                    nmod_mat_t res){
    nn_ptr vec_mul = _nmod_vec_init(mat->ncols); 
    nn_ptr vec_res = _nmod_vec_init(mat->nrows); 
    for(int i=0; i<nmod_mat_ncols(b); i++){
        // TODO see if this can be done more efficiently
        for(int j=0; j<nmod_mat_nrows(b); j++){
            vec_mul[j] = nmod_mat_get_entry(b,j,i);
        }

        struct_nmod_toeplitz_right_mul_vec(mat,vec_mul,vec_res);
        // TODO again, see if this can be done more efficiently
        for(int j=0; j<mat->nrows; j++){
            nmod_mat_set_entry(res,j,i,vec_res[j]);
        }
    }
    _nmod_vec_clear(vec_mul);
    _nmod_vec_clear(vec_res);
}

void struct_nmod_toeplitz_left_mul_vec(struct_nmod_toeplitz_t mat,
                                    nn_ptr v,
                                    nn_ptr res){
    //TODO
}

void struct_nmod_toeplitz_left_mul_mat(struct_nmod_toeplitz_t mat,
                                    nmod_mat_t b,
                                    nmod_mat_t res){
    // TODO result matrix indexed in row major order -> copy directly there    
}
