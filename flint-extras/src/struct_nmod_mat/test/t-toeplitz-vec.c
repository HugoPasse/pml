#include <stdio.h>
#include "struct_nmod_mat.h"

TEST_FUNCTION_START(struct_nmod_toeplitz_mul_vec, state)
{
    for(int i = 0; i < 100 * flint_test_multiplier(); i++)
    {
    
        long n = n_randtest_bits(state,6);
        long m = n_randtest_bits(state,6);
        ulong N = n_randtest_prime(state, 1);
        
        struct_nmod_toeplitz_t mat;
        nn_ptr data = _nmod_vec_init(n+m-1);
        
        nmod_t mod;
        nmod_init(&mod,N);
        _nmod_vec_randtest(data, state, n+m-1, mod); 

        struct_nmod_toeplitz_init_set(mat, n, m, N, data);        
        
        // debug
        nmod_mat_t dense_mat;
        nmod_mat_init(dense_mat,n,m,N);
        struct_nmod_toeplitz_dense(mat,dense_mat);
    
        // vector to multiply
        nn_ptr v = _nmod_vec_init(m);
        _nmod_vec_randtest(v,state,m,mat->mod); 
    
        // structured part
        nn_ptr res_struct = _nmod_vec_init(n);
        struct_nmod_toeplitz_right_mul_vec(mat,v,res_struct); 
        // dense part
        nn_ptr res_dense = _nmod_vec_init(n);
        nmod_mat_mul_nmod_vec(res_dense, dense_mat, v, m);
    
        // Check equality
        int result = _nmod_vec_equal(res_struct, res_dense, n);
        if (!result) TEST_FUNCTION_FAIL("Structured result and dense result do not coincide\n");
    }    
    TEST_FUNCTION_END(state);
}


