#include "nmod_mat_struct.h"
#include "nmod_poly_mat_approximant.h"

#include <flint/nmod_mat.h>

#include <flint/flint.h>
#include <flint/nmod_types.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

void nmod_mat_toeplitz_init(nmod_mat_toeplitz_t mat, 
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

void nmod_mat_toeplitz_set(nmod_mat_toeplitz_t mat,
                               nn_ptr data){
    _nmod_vec_set(mat->data, data, mat->nrows+mat->ncols-1);
    // Reduce mat->data mod n 
}

void nmod_mat_toeplitz_init_set(nmod_mat_toeplitz_t mat,
                                    ulong nrows,
                                    ulong ncols,
                                    ulong N,
                                    nn_ptr data){
    nmod_mat_toeplitz_init(mat, nrows, ncols, N);
    nmod_mat_toeplitz_set(mat, data);
}


void nmod_mat_toeplitz_randtest(nmod_mat_toeplitz_t mat,
                                  ulong nrows,
                                  ulong ncols,
                                  ulong N,
                                  flint_rand_t state){

    nmod_mat_toeplitz_init(mat, nrows, ncols, N);
    _nmod_vec_randtest(mat->data, state, mat->ncols + mat->nrows - 1, mat->mod);
}

void nmod_mat_toeplitz_clear(nmod_mat_toeplitz_t mat){
    if(mat->data){
        flint_free(mat->data);
    }
}

ulong nmod_mat_toeplitz_nrows(nmod_mat_toeplitz_t mat){
    return mat->nrows;
}

ulong nmod_mat_toeplitz_ncols(nmod_mat_toeplitz_t mat){
    return mat->ncols;
}

void nmod_mat_toeplitz_as_poly(nmod_mat_toeplitz_t mat,
                                nmod_poly_t pol){
    int len = mat->nrows+mat->ncols-1;
    nmod_poly_init2(pol, mat->mod.n, len);
    _nmod_vec_set(pol->coeffs, mat->data, len); 
    pol->length = len;
    _nmod_poly_normalise(pol);
}

void nmod_mat_toeplitz_dense(nmod_mat_toeplitz_t mat,
                                  nmod_mat_t dense_mat){

    for(int i=0; i<nmod_mat_toeplitz_nrows(mat); i++){
        for(int j=0; j<nmod_mat_toeplitz_ncols(mat); j++){
            int k = nmod_mat_toeplitz_ncols(mat) - 1 + i - j;
            nmod_mat_set_entry(dense_mat, i, j, mat->data[k]);
        }
    }  
}

// ---------------------------- //
// ----------- SUM ------------ //
// ---------------------------- //

void nmod_mat_toeplitz_add(nmod_mat_toeplitz_t a,
                               nmod_mat_toeplitz_t b, 
                               nmod_mat_toeplitz_t res){
    _nmod_vec_add(res->data, a->data, b->data, b->nrows+b->ncols-1, a->mod);
}

// ---------------------------- //
// ------ MULTIPLICATION ------ //
// ---------------------------- //
void nmod_mat_toeplitz_right_mul_vec(nmod_mat_toeplitz_t mat,
                                    nn_ptr v,
                                    nn_ptr res){
    _nmod_poly_mulmid(res, mat->data, mat->nrows+mat->ncols-1, v, mat->ncols, mat->ncols-1, mat->nrows+mat->ncols-1, mat->mod);  
}

void nmod_mat_toeplitz_right_mul_mat(nmod_mat_toeplitz_t mat,
                                    nmod_mat_t b,
                                    nmod_mat_t res){
    nn_ptr vec_mul = _nmod_vec_init(mat->ncols); 
    nn_ptr vec_res = _nmod_vec_init(mat->nrows); 
    for(int i=0; i<nmod_mat_ncols(b); i++){
        // TODO see if this can be done more efficiently
        for(int j=0; j<nmod_mat_nrows(b); j++){
            vec_mul[j] = nmod_mat_get_entry(b,j,i);
        }

        nmod_mat_toeplitz_right_mul_vec(mat,vec_mul,vec_res);
        // TODO again, see if this can be done more efficiently
        for(int j=0; j<mat->nrows; j++){
            nmod_mat_set_entry(res,j,i,vec_res[j]);
        }
    }
    _nmod_vec_clear(vec_mul);
    _nmod_vec_clear(vec_res);
}

void nmod_mat_toeplitz_left_mul_vec(nmod_mat_toeplitz_t mat,
                                    nn_ptr v,
                                    nn_ptr res){
    //TODO
}

void nmod_mat_toeplitz_left_mul_mat(nmod_mat_toeplitz_t mat,
                                    nmod_mat_t b,
                                    nmod_mat_t res){
    // TODO result matrix indexed in row major order -> copy directly there    
}


// WARNING: The code here asssumes that the toeplitz matrix is full rank. 
// The issue is that nmod_poly_mat_pmbasis wants appbas allocated with right sizes, 
// i.e. we need a priori knowledge of the rank
// Possible FIX : allocate too many rows for appbas and check a posteriori each row of the basis
void nmod_mat_toeplitz_right_kernel_basis(nmod_mat_toeplitz_t mat,
                                           nmod_mat_t res){
    slong order = mat->nrows + mat->ncols - 1; 
    int r = mat->ncols - mat->nrows; 
    /* Constructing the polynomial matrix [[T(x)],[-1]] */
    // Start with the two entries
    nmod_poly_t pol, minus_one;
    nmod_mat_toeplitz_as_poly(mat, pol);

    nmod_poly_init(minus_one, mat->mod.n);
    nmod_poly_zero(minus_one);
    nmod_poly_set_coeff_ui(minus_one, 0, mat->mod.n-1);
    
    // Set polynomial matrix entries
    nmod_poly_mat_t pmat;
    nmod_poly_mat_init(pmat, 2, 1, mat->mod.n);
    nmod_poly_set(nmod_poly_mat_entry(pmat, 0, 0), pol);
    nmod_poly_set(nmod_poly_mat_entry(pmat, 1, 0), minus_one);
    
    /* We know compute an approximant basis of (pmat,order) in (ncols,ncols)-popov_form */
    nmod_poly_mat_t appbas;
    nmod_poly_mat_init(appbas, r, 2, mat->mod.n);
    slong shifts[2] = {0, 0};
        
    nmod_poly_mat_pmbasis(appbas, shifts, pmat, order);

    for (int i=0; i<r; i++) {
        for (int j=0; j<mat->ncols; j++) {
            nmod_mat_set_entry(res, j, i, nmod_poly_get_coeff_ui(nmod_poly_mat_entry(appbas, i, 0), j));   
        }
    }
}

void nmod_mat_toeplitz_solve_right(nmod_mat_toeplitz_t mat,
                                    nmod_mat_t b,
                                    nmod_mat_t res){

}
