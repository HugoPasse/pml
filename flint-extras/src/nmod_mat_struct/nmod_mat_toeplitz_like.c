#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_types.h>
#include <flint/nmod_vec.h>
#include <stdio.h>
#include <time.h>

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
        _nmod_vec_set((mat->H)[i], H[i], mat->ncols);
    }
}

void nmod_toeplitz_like_clear(nmod_toeplitz_like_t mat){
    for (int i=0; i<mat->ngens; i++){
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
    //nmod_mat_init(res->G, a->nrows, k, a->mod.n);
    //nmod_mat_init(res->H, a->ncols, k, a->mod.n);
    nmod_toeplitz_like_init(res, a->nrows, a->ncols, a->mod.n);
    res->ngens = k;
    res->G = flint_malloc(k*sizeof(nn_ptr)); 
    res->H = flint_malloc(k*sizeof(nn_ptr)); 
    
    for(int i=0; i<a->ngens; i++){
        res->G[i] = _nmod_vec_init(a->nrows);
        res->H[i] = _nmod_vec_init(a->ncols);
        _nmod_vec_set(res->G[i], a->G[i], a->nrows);
        _nmod_vec_set(res->H[i], a->H[i], a->ncols);
    }
    
    for(int i=0; i<b->ngens; i++){
        res->G[i+a->ngens] = _nmod_vec_init(b->nrows);
        res->H[i+a->ngens] = _nmod_vec_init(b->ncols);
        _nmod_vec_set(res->G[i+a->ngens], b->G[i], b->nrows);
        _nmod_vec_set(res->H[i+a->ngens], b->H[i], b->ncols);
    }
    

    //for (int i=0; i<nmod_mat_nrows(a->G); i++) {
    //    for (int j=0; j<nmod_mat_ncols(a->G); j++) {
    //        nmod_mat_set_entry(res->G, i, j, nmod_mat_get_entry(a->G, i, j)); 
    //    }
    //}
    //for (int i=0; i<nmod_mat_nrows(a->H); i++) {
    //    for (int j=0; j<nmod_mat_ncols(a->H); j++) {
    //        nmod_mat_set_entry(res->H, i, j, nmod_mat_get_entry(a->H, i, j)); 
    //    }
    //} 
    //for (int i=0; i<nmod_mat_nrows(b->G); i++) {
    //    for (int j=0; j<nmod_mat_ncols(b->G); j++) {
    //        nmod_mat_set_entry(res->G, i, j+nmod_mat_ncols(a->G), nmod_mat_get_entry(b->G, i, j)); 
    //    }
    //}
    //for (int i=0; i<nmod_mat_nrows(b->H); i++) {
    //    for (int j=0; j<nmod_mat_ncols(b->H); j++) {
    //        nmod_mat_set_entry(res->H, i, j+nmod_mat_ncols(a->H), nmod_mat_get_entry(b->H, i, j)); 
    //    }
    //}
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


void nmod_toeplitz_like_mul_nmod_mat_left(nmod_toeplitz_like_t mat,
                                            nmod_mat_t b,
                                            nmod_mat_t res){
    for(int i=0; i<nmod_mat_ncols(b); i++){
        nn_ptr v = _nmod_vec_init(nmod_mat_nrows(b));
        nn_ptr vres = _nmod_vec_init(mat->nrows);
        
        for(int j=0; j<nmod_mat_nrows(b); j++){
            v[j] = nmod_mat_get_entry(b, j, i);
        }
        nmod_toeplitz_like_mul_nmod_vec_left(mat, v, vres);
        
        for(int j=0; j<mat->nrows; j++){
           nmod_mat_set_entry(res, j, i, vres[j]); 
        }
        _nmod_vec_clear(v);
        _nmod_vec_clear(vres);
    }     
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
void _generators_as_nmod_mat(nmod_toeplitz_like_t mat,
                            nmod_mat_t G,
                            nmod_mat_t H){
    nmod_mat_init(G, mat->nrows, mat->ngens, mat->mod.n);  
    nmod_mat_init(H, mat->ncols, mat->ngens, mat->mod.n); 

    for(int i=0; i<mat->ngens; i++){
        for(int j=0; j<mat->nrows; j++){
            nmod_mat_set_entry(G, j, i, mat->G[i][j]);
        }
    }

    for(int i=0; i<mat->ngens; i++){
        for(int j=0; j<mat->ncols; j++){
            nmod_mat_set_entry(H, j, i, mat->H[i][j]);
        }
    }
}

void nmod_toeplitz_like_reduce_gens(nmod_toeplitz_like_t mat){
    // Representing generators as matrices
    nmod_mat_t G;    
    nmod_mat_t H;    
    _generators_as_nmod_mat(mat, G, H);
    
    //nmod_mat_t Ht;    
    //nmod_mat_init(Ht, nmod_mat_ncols(H), nmod_mat_nrows(H), mat->mod.n);
    //nmod_mat_transpose(Ht,H);
    //nmod_mat_t GHt;
    //nmod_mat_init(GHt, mat->nrows, mat->ncols, mat->mod.n);
    //nmod_mat_zero(GHt);
    //nmod_mat_mul(GHt, G, Ht);
    //nmod_mat_print_pretty(GHt);

    // PLU form of G 
    slong *permg = flint_malloc(mat->nrows * sizeof(slong));
    slong rg = nmod_mat_lu(permg, G, 0);
    // Use submatrices to reduce the number of columns with no new allocation
    nmod_mat_t G1, H1;
    nmod_mat_window_init(G1, G, 0, 0, mat->nrows, rg);
    nmod_mat_window_init(H1, H, 0, 0, mat->ncols, rg);
    // Upper triangular matrix
    nmod_mat_t Utg;
    nmod_mat_init(Utg, mat->ngens, rg, mat->mod.n);
    // Fill upper triangular matrix and clear lower triangular matrix
    for(int i=0; i<rg; i++){
        for(int j=i; j<mat->ngens; j++){
            nmod_mat_set_entry(Utg, j, i, nmod_mat_get_entry(G, i, j));
            nmod_mat_set_entry(G, i, j, 0);
        }
        nmod_mat_set_entry(G, i, i, 1);
    }
    // Apply the permutation to the lower triangular part 
    nmod_mat_permute_rows(G1, permg, NULL); 
    nmod_mat_mul(H1, H, Utg);
    
    // --------- Debug ----------------
    //nmod_mat_t H1t;
    //nmod_mat_init(H1t, nmod_mat_ncols(H1), nmod_mat_nrows(H1), mat->mod.n);
    //nmod_mat_transpose(H1t, H1);
    //nmod_mat_zero(GHt);
    //nmod_mat_mul(GHt, G1, H1t);
    //nmod_mat_print_pretty(GHt);

    // PLU form of H 
    slong *permh = flint_malloc(mat->ncols * sizeof(slong));
    slong rh = nmod_mat_lu(permh, H, 0);
    // Use submatrices to reduce the number of columns with no new allocation
    nmod_mat_t G2, H2;
    nmod_mat_window_init(G2, G1, 0, 0, mat->nrows, rh);
    nmod_mat_window_init(H2, H1, 0, 0, mat->ncols, rh);
    // Upper triangular matrix
    nmod_mat_t Uth;
    nmod_mat_init(Uth, rg, rh, mat->mod.n);
    // Fill upper triangular matrix and clear lower triangular matrix
    for(int i=0; i<rh; i++){
        for(int j=i; j<mat->ngens; j++){
            nmod_mat_set_entry(Uth, j, i, nmod_mat_get_entry(H1, i, j));
            nmod_mat_set_entry(H1, i, j, 0);
        }
        nmod_mat_set_entry(H2, i, i, 1);
    }
    // Apply the permutation to the lower triangular part 
    _nmod_vec_print(permh, mat->ncols, mat->mod);
    nmod_mat_print_pretty(H2);
    nmod_mat_permute_rows(H2, permh, NULL); 
    nmod_mat_mul(G2, G1, Uth);

    // Clearing matrix data 
    nmod_toeplitz_like_clear(mat);
    // Convert to nn_ptr*
    mat->ngens = nmod_mat_ncols(H2);
    mat->G = flint_malloc(mat->ngens*sizeof(nn_ptr));
    mat->H = flint_malloc(mat->ngens*sizeof(nn_ptr));
    for(int i=0; i<mat->ngens; i++){
        (mat->G)[i] = _nmod_vec_init(mat->nrows);
        (mat->H)[i] = _nmod_vec_init(mat->ncols);
        for(int j=0; j<mat->nrows; j++){
            (mat->G)[i][j] = nmod_mat_get_entry(G2,j,i);
        }
        for(int j=0; j<mat->ncols; j++){
            (mat->H)[i][j] = nmod_mat_get_entry(H2,j,i);
        }
    } 

    // --------- DEBUG --------- 

    //nmod_mat_t H2t;
    //nmod_mat_init(H2t, nmod_mat_ncols(H2), nmod_mat_nrows(H2), mat->mod.n);
    //nmod_mat_transpose(H2t, H2);
    //nmod_mat_zero(GHt);
    //nmod_mat_mul(GHt, G2, H2t);
    //nmod_mat_print_pretty(GHt);
    
    //nmod_mat_clear(Ht);
    //nmod_mat_clear(H1t);
    //nmod_mat_clear(H2t);
    //nmod_mat_clear(GHt);

    // Free memory
    nmod_mat_window_clear(H2);
    nmod_mat_window_clear(G2);
    nmod_mat_window_clear(H1);
    nmod_mat_window_clear(G1);
    nmod_mat_clear(H);
    nmod_mat_clear(G);
    nmod_mat_clear(Uth);
    nmod_mat_clear(Utg);
    flint_free(permg);
    flint_free(permh);
}
