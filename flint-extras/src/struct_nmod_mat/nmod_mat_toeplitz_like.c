#include "struct_nmod_mat.h"

void nmod_toeplitz_like_init(nmod_toeplitz_like_t mat,
                             ulong nrows,
                             ulong ncols,
                             ulong N){

}

void nmod_toeplitz_like_set(nmod_toeplitz_like_t mat,
                             ulong ngens,
                             nmod_mat_lower_toeplitz_t* lower_generators,
                             nmod_mat_upper_toeplitz_t* upper_generators){

}

void nmod_toeplitz_like_dense(nmod_toeplitz_like_t mat, 
                               nmod_mat_t dense_mat){

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
    for(int i=0; i<mat->ngens; i++){
        // TODO dedicated function for upper / lower types to do this
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

