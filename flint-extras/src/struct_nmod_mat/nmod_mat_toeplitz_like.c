#include <stdlib.h>
#include "struct_nmod_mat.h"

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
                             ulong ngens,
                             nmod_mat_lower_toeplitz_t* lower_generators,
                             nmod_mat_upper_toeplitz_t* upper_generators){
    mat->ngens = ngens;

    mat->lower_gens = lower_generators;
    mat->upper_gens = upper_generators;
}

void nmod_toeplitz_like_set_vec(nmod_toeplitz_like_t mat,
                                 ulong ngens,
                                 ulong rank,
                                 nn_ptr* lower_generators,
                                 nn_ptr* upper_generators){
    mat->ngens = ngens;
    mat->lower_gens = malloc(ngens*sizeof(nmod_mat_lower_toeplitz_t));
    mat->upper_gens = malloc(ngens*sizeof(nmod_mat_upper_toeplitz_t));

    for(int i=0; i<ngens; i++){
        nmod_mat_lower_toeplitz_init_set((mat->lower_gens)[i], mat->nrows, rank, mat->mod.n, lower_generators[i]);
        nmod_mat_upper_toeplitz_init_set((mat->upper_gens)[i], rank, mat->ncols, mat->mod.n, upper_generators[i]);
    }
}

void nmod_toeplitz_like_randtest(nmod_toeplitz_like_t mat,
                                  ulong nrows,
                                  ulong ncols,
                                  ulong ngens,
                                  ulong rank,
                                  ulong N,
                                  flint_rand_t state){

    nmod_toeplitz_like_init(mat,nrows,ncols,N);
    nn_ptr* lower_generators = malloc(ngens * sizeof(nn_ptr));
    nn_ptr* upper_generators = malloc(ngens * sizeof(nn_ptr));
    for(int i=0; i<ngens; i++){
        lower_generators[i] = _nmod_vec_init(nrows+rank-1);
        upper_generators[i] = _nmod_vec_init(ncols+rank-1);

        _nmod_vec_randtest(lower_generators[i], state, nrows+rank-1, mat->mod);
        _nmod_vec_randtest(upper_generators[i], state, ncols+rank-1, mat->mod);
    }
    nmod_toeplitz_like_set_vec(mat,ngens,rank,lower_generators,upper_generators);
}

void nmod_toeplitz_like_dense(nmod_toeplitz_like_t mat, 
                               nmod_mat_t dense_mat){
    nmod_mat_init(dense_mat, mat->nrows, mat->ncols, mat->mod.n);
    nmod_mat_zero(dense_mat);

    for(int i=0; i<mat->ngens; i++){
        nmod_mat_t ldense, udense;
        nmod_mat_init(ldense, mat->nrows, (*(mat->lower_gens+i))->ncols, mat->mod.n);
        nmod_mat_init(udense, (*(mat->upper_gens+i))->nrows, mat->ncols, mat->mod.n);

        nmod_mat_zero(ldense);
        nmod_mat_zero(udense);

        nmod_mat_lower_toeplitz_dense(mat->lower_gens[i], ldense);
        nmod_mat_upper_toeplitz_dense(mat->upper_gens[i], udense);
        
        nmod_mat_addmul(dense_mat, dense_mat, ldense, udense);
    }
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
    dot_params_t params;
    params = _nmod_vec_dot_params(lmat->ncols, lmat->mod);
    int k0,k1,l;
    l  = FLINT_MIN(FLINT_MIN(i,j)+1,lmat->ncols);
    k0 = FLINT_MAX(0, (int)i - (int)lmat->ncols + 1) ;
    k1 = umat->ncols - j - 1; 
    if(j<lmat->ncols-1){
        k0 += FLINT_MAX((int)FLINT_MIN(i,lmat->ncols-1) - (int)j, 0);
    }
    return _nmod_vec_dot_rev(umat->data + k1, lmat->data + k0, l, lmat->mod, params);
}

ulong nmod_toeplitz_like_get_entry(nmod_toeplitz_like_t mat,
                                    ulong i,
                                    ulong j){
    ulong result = 0;
    for(int k=0; k<mat->ngens; k++){
        result += nmod_mat_lower_upper_toeplitz_pdt(*(mat->lower_gens+k), *(mat->upper_gens+k), i, j);
    }
    return result;
}

/* Arithmetic */
void nmod_toeplitz_like_add(nmod_toeplitz_like_t a,
                             nmod_toeplitz_like_t b, 
                             nmod_toeplitz_like_t res){
    res->nrows = a->nrows;
    res->ncols = a->ncols;
    res->ngens = a->ngens + b->ngens;
    nmod_init(&res->mod, a->mod.n);

    for(int i=0; i<a->ngens; i++){
         // Here copy generators
    }
    for(int i=0; i<b->ngens; i++){
        // Here copy generators
    }
    // I think we should here do a copy of the generators
}

