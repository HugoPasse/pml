#include <stdio.h>
#include "struct_nmod_mat.h"

TEST_FUNCTION_START(struct_nmod_toeplitz_mul_vec, state)
{
    long n = 4;
    long m = 3;
    ulong N = 17;
    
    struct_nmod_toeplitz_t mat;
    nn_ptr data = _nmod_vec_init(n+m-1);
    for(int i=0; i<n+m-1; i++){
        data[i] = i;
    }
    
    struct_nmod_toeplitz_init_set(mat, n, m, N, data);        
    
    // debug
    nmod_mat_t dense_mat;
    nmod_mat_init(dense_mat,n,m,N);
    struct_nmod_toeplitz_dense(mat,dense_mat);
    nmod_mat_print_pretty(dense_mat);

    // vector to multiply
    nn_ptr v = _nmod_vec_init(m);
    _nmod_vec_randtest(v,state,m,mat->mod); 
    
    _nmod_vec_print_pretty(v,m,mat->mod);
    // structured part
    nn_ptr res_struct = _nmod_vec_init(n);
    struct_nmod_toeplitz_mul_vec(mat,v,res_struct); 
    
    _nmod_vec_print_pretty(res_struct,n,mat->mod);

    //// dense part
    //nn_ptr res_dense = _nmod_vec_init(n);
    // 
    //
    //nmod_mat_mul_nmod_vec(res_dense,mat,v,m);

    TEST_FUNCTION_END(state);
}


