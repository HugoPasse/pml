#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>

#include "nmod_mat_struct.h"

void nmod_toeplitz_like_init(nmod_toeplitz_like_t mat,
                             ulong nrows,
                             ulong ncols,
                             ulong N){
    mat->nrows = nrows;
    mat->ncols = ncols;
    mat->ngens = 0;

    mat->mod.n = N;
    mat->mod.ninv = n_preinvert_limb(N);
    mat->mod.norm = flint_clz(N);
}

// TODO see if generators should be copied or not ?
// Do we want to have nn_ptr as input for lower / upper generators ?
void nmod_toeplitz_like_set(nmod_toeplitz_like_t mat,
                             nn_ptr* G, 
                             nn_ptr* H,
                             int ngens){
    mat->G = flint_malloc(ngens*sizeof(nn_ptr));
    mat->H = flint_malloc(ngens*sizeof(nn_ptr));
    mat->ngens = ngens;
    for (int i=0; i<mat->ngens; i++) {
        (mat->G)[i] = _nmod_vec_init(mat->nrows);
        _nmod_vec_set((mat->G)[i], G[i], mat->nrows);
        (mat->H)[i] = _nmod_vec_init(mat->ncols);
        _nmod_vec_set((mat->G)[i], G[i], mat->nrows);
    }
}

void nmod_toeplitz_like_clear(nmod_toeplitz_like_t mat){
    for (int i=0; i<mat->ngens; i++) {
        _nmod_vec_clear((mat->G)[i]);
        _nmod_vec_clear((mat->H)[i]);
    }
    flint_free(mat->G);
    flint_free(mat->H);
}

void nmod_toeplitz_like_randtest(nmod_toeplitz_like_t mat,
                                  ulong nrows,
                                  ulong ncols,
                                  ulong ngens,
                                  ulong N,
                                  flint_rand_t state){
    nmod_toeplitz_like_init(mat,nrows,ncols,N);
    mat->G = flint_malloc(ngens*sizeof(nn_ptr));
    mat->H = flint_malloc(ngens*sizeof(nn_ptr));
    mat->ngens = ngens;

    for (int i=0; i<ngens; i++){
        (mat->G)[i] = _nmod_vec_init(mat->nrows);
        _nmod_vec_randtest((mat->G)[i], state, mat->nrows, mat->mod); 
        (mat->H)[i] = _nmod_vec_init(mat->ncols);
        _nmod_vec_randtest((mat->H)[i], state, mat->ncols, mat->mod); 
    } 
}

void nmod_toeplitz_like_dense(nmod_toeplitz_like_t mat, 
                               nmod_mat_t dense_mat){
    nmod_mat_zero(dense_mat);
    for(int i=0; i<mat->ngens; i++){
        nmod_mat_lower_toeplitz_t low;
        nmod_mat_lower_toeplitz_init(low, mat->nrows, mat->ncols, mat->mod.n);
        for(int j=0; j<mat->nrows; j++){
            low->data[j] = (mat->G)[i][j];
        }
        nmod_mat_t dense_low;
        nmod_mat_init(dense_low, mat->nrows, mat->ncols, mat->mod.n);
        nmod_mat_lower_toeplitz_dense(low, dense_low);

        nmod_mat_circulant_t circ;
        nmod_mat_t dense_circ;
        nmod_mat_init(dense_circ, mat->ncols, mat->ncols, mat->mod.n);
        nmod_mat_circulant_init(circ, mat->ncols, mat->ncols, mat->mod.n);
        for(int j=0; j<mat->ncols; j++){
            circ->data[j] = (mat->H)[i][j];
        }
        nmod_mat_circulant_dense(circ, dense_circ);

        nmod_mat_t product;
        nmod_mat_init(product, mat->nrows, mat->ncols, mat->mod.n);
        nmod_mat_mul(product, dense_low, dense_circ);
        nmod_mat_add(dense_mat, dense_mat, product);
        
        nmod_mat_lower_toeplitz_clear(low);
        nmod_mat_circulant_clear(circ);
        nmod_mat_clear(product);
        nmod_mat_clear(dense_low);
        nmod_mat_clear(dense_circ);
    } 
    // Remains to multiply by the last mat
    nmod_mat_t final_pdt;
    nmod_mat_init(final_pdt, mat->nrows, mat->nrows, mat->mod.n);
    nmod_mat_zero(final_pdt);
    for(int i=0; i<((mat->nrows-1) / mat->ncols) + 1; i++){
        for(int j=0; j<mat->nrows - i*mat->ncols; j++){
            nmod_mat_set_entry(final_pdt, j+i*mat->ncols, j, 1);
        }
    }
    nmod_mat_mul(dense_mat, final_pdt, dense_mat);
    nmod_mat_clear(final_pdt);
}


ulong nmod_toeplitz_like_t_nrows(nmod_toeplitz_like_t mat){
    return mat->nrows;
}

ulong nmod_toeplitz_like_t_ncols(nmod_toeplitz_like_t mat){
    return mat->ncols;
}

ulong nmod_toeplitz_like_t_mod(nmod_toeplitz_like_t mat){
    return mat->mod.n;
}

ulong nmod_mat_lower_upper_toeplitz_pdt(nmod_mat_lower_toeplitz_t lmat,
                                         nmod_mat_upper_toeplitz_t umat,
                                         ulong i,
                                         ulong j){
    // TODO, not sure if we want to keep this
    return 0;    
}

/* Arithmetic */
void nmod_toeplitz_like_add(nmod_toeplitz_like_t a,
                             nmod_toeplitz_like_t b, 
                             nmod_toeplitz_like_t res){
    int k = a->ngens + b->ngens;
    nmod_mat_init(res->G, a->nrows, k, a->mod.n);
    nmod_mat_init(res->H, a->ncols, k, a->mod.n);
    
    for (int i=0; i<nmod_mat_nrows(a->G); i++) {
        for (int j=0; j<nmod_mat_ncols(a->G); j++) {
            nmod_mat_set_entry(res->G, i, j, nmod_mat_get_entry(a->G, i, j)); 
        }
    }
    for (int i=0; i<nmod_mat_nrows(a->H); i++) {
        for (int j=0; j<nmod_mat_ncols(a->H); j++) {
            nmod_mat_set_entry(res->H, i, j, nmod_mat_get_entry(a->H, i, j)); 
        }
    } 
    for (int i=0; i<nmod_mat_nrows(b->G); i++) {
        for (int j=0; j<nmod_mat_ncols(b->G); j++) {
            nmod_mat_set_entry(res->G, i, j+nmod_mat_ncols(a->G), nmod_mat_get_entry(b->G, i, j)); 
        }
    }
    for (int i=0; i<nmod_mat_nrows(b->H); i++) {
        for (int j=0; j<nmod_mat_ncols(b->H); j++) {
            nmod_mat_set_entry(res->H, i, j+nmod_mat_ncols(a->H), nmod_mat_get_entry(b->H, i, j)); 
        }
    }
}

void nmod_toeplitz_like_mul_nmod_vec_left(nmod_toeplitz_like_t a,
                                           nn_ptr v, 
                                           nn_ptr res){
    // Polynomials associated to the displacement operator
    nn_ptr Q = _nmod_vec_init(a->ncols+1);
    _nmod_vec_zero(Q, a->ncols+1);
    Q[0] = a->mod.n - 1;
    Q[a->ncols] = 1;

    nn_ptr Qtilinv = _nmod_vec_init(a->nrows);
    _nmod_vec_zero(Qtilinv, a->nrows);
    for(int i=0; i*a->ncols < a->nrows; i++){
        Qtilinv[i * a->ncols] = 1;
    } 

    // Reverse the input vector 
    nn_ptr vrev = _nmod_vec_init(a->ncols);
    _nmod_poly_reverse(vrev,v,a->ncols,a->ncols);

    // Setting the result vector to zero
    nn_ptr res_tmp = _nmod_vec_init(a->nrows);
    _nmod_vec_zero(res_tmp, a->nrows);
    _nmod_vec_zero(res, a->nrows);
    nn_ptr tmp = _nmod_vec_init(FLINT_MAX(a->ncols,a->nrows));
    nn_ptr tmp1 = _nmod_vec_init(FLINT_MAX(a->ncols,a->nrows));
    for(int i=0; i<a->ngens; i++){
        _nmod_vec_zero(tmp, FLINT_MAX(a->ncols,a->nrows));
        _nmod_vec_zero(tmp1, FLINT_MAX(a->ncols,a->nrows));
        // First product by M_{hi,Q}
        _nmod_poly_mulmod(tmp,(a->H)[i],a->ncols,vrev,a->ncols,Q,a->ncols+1,a->mod);
        // Product by J_n -> reverse polynomial
        _nmod_poly_reverse(tmp1,tmp,a->ncols,a->ncols);
        // Second product by M_{gi,P}
        _nmod_poly_mullow(tmp, tmp1, a->ncols, (a->G)[i], a->nrows, a->nrows, a->mod);
        // Adding the term to the complete sum
        _nmod_vec_add(res_tmp, res_tmp, tmp, a->nrows, a->mod);
    }
    _nmod_poly_mullow(res, res_tmp, a->nrows, Qtilinv, a->nrows, a->nrows, a->mod);
            
    flint_free(Q);
    flint_free(Qtilinv);
    flint_free(vrev);
    flint_free(tmp);
    flint_free(tmp1);
    flint_free(res_tmp);
}

void nmod_toeplitz_like_mul_nmod_vec_right(nmod_toeplitz_like_t a,
                                           nn_ptr v, 
                                           nn_ptr res){
    // TODO
}

void nmod_toeplitz_like_mul(nmod_toeplitz_like_t a,
                             nmod_toeplitz_like_t b, 
                             nmod_toeplitz_like_t res){
    // TODO
} 

void nmod_toeplitz_like_inv(nmod_toeplitz_like_t a,
                             nmod_toeplitz_like_t res){
    // TODO
}

/* Generators reduction */
void nmod_toeplitz_like_reduce_gens(nmod_toeplitz_like_t mat){
    // TODO
}

// TODO: Qtilinv has an explicit formula, it should be faster to use it
// TODO: Compute the polynomials P,Q,Qtil and Qtilinv at initialization, add as members of struct
//void nmod_toeplitz_like_mul_nmod_vec_left_old(nmod_toeplitz_like_t a,
//                                           nn_ptr v, 
//                                           nn_ptr res){
//    // Polynomials related to the displacement operator
//    nmod_poly_t P,Q,Qtilinv;
//    nmod_poly_init2(P, a->mod.n, a->nrows+1);
//    nmod_poly_init2(Q, a->mod.n, a->ncols+1);
//    nmod_poly_init2(Qtilinv, a->mod.n, a->ncols+1);
//
//    nmod_poly_zero(P); 
//    nmod_poly_set_coeff_ui(P, a->nrows, 1);
//
//    nmod_poly_zero(Q); 
//    nmod_poly_set_coeff_ui(Q, a->ncols, 1);
//    nmod_poly_set_coeff_ui(Q, 0, a->mod.n-1);
//
//    nmod_poly_zero(Qtilinv); 
//    nmod_poly_set_coeff_ui(Qtilinv, a->ncols, a->mod.n-1);
//    nmod_poly_set_coeff_ui(Qtilinv, 0, 1);
//    nmod_poly_invmod(Qtilinv, Qtilinv, P);
//
//    // Multiply by Y_Q
//    // This is taking the reverse of the polynomial of v
//    nmod_poly_t vx, res_poly;
//    nmod_poly_init2(vx, a->mod.n, a->ncols);
//    for(int i=0; i<a->ncols; i++){
//        nmod_poly_set_coeff_ui(vx, i, v[a->ncols-i-1]);
//    }
//
//    nmod_poly_init2(res_poly, a->mod.n, a->nrows);
//    nmod_poly_zero(res_poly); 
//
//    nmod_poly_t tmp, gi;
//    nmod_poly_init(tmp, a->mod.n); // Degree known, to specify 
//    nmod_poly_init(gi, a->mod.n); 
//    for(int i=0; i<a->ngens; i++){
//        // Computing h_i as a polynomial
//        nmod_poly_zero(tmp); 
//        for(int j=0; j<a->ncols; j++){
//            nmod_poly_set_coeff_ui(tmp, j, (a->H)[i][j]);
//        }
//        // First product by M_{hi,Q}
//        nmod_poly_mulmod(tmp,tmp,vx,Q);
//        // Product by J_n -> reverse polynomial
//        nmod_poly_reverse(tmp,tmp,a->ncols);
//        // Computing g_i as a polynomial
//        nmod_poly_zero(gi);
//        for(int j=0; j<a->nrows; j++){
//            nmod_poly_set_coeff_ui(gi, j, (a->G)[i][j]);
//        }
//        // Second product by M_{gi,P}
//        nmod_poly_mulmod(tmp,tmp,gi,P);
//        // Adding the term to the complete sum
//        nmod_poly_add(res_poly, res_poly, tmp);
//    }
//    //nmod_poly_mulmod(res_poly, res_poly, Qtilinv, P);
//    nmod_poly_mullow(res_poly, res_poly, Qtilinv, a->nrows);
//            
//    //_nmod_vec_set(res, res_poly->coeffs, nmod_poly_degree(res_poly));
//    for(int i=0; i<a->nrows; i++){
//       res[i] = nmod_poly_get_coeff_ui(res_poly, i);    
//    }
//    
//    nmod_poly_clear(vx);
//    nmod_poly_clear(tmp);
//    nmod_poly_clear(gi);
//    nmod_poly_clear(res_poly);
//
//    nmod_poly_clear(P);
//    nmod_poly_clear(Q);
//    nmod_poly_clear(Qtilinv);
//}
