#include <stdio.h>
#include "struct_nmod_mat.h"

TEST_FUNCTION_START(struct_nmod_toeplitz_mul_mat, state)
{
    for(int i = 0; i < 100 * flint_test_multiplier(); i++)
    {
    
        long n = n_randtest_bits(state,6);
        long m = n_randtest_bits(state,6);
        long l = n_randtest_bits(state,5);
        ulong N = n_randtest_prime(state, 1);
        
        struct_nmod_toeplitz_t mat;
        nn_ptr data = _nmod_vec_init(n+m-1);
        
        nmod_t mod;
        nmod_init(&mod,N);
        _nmod_vec_randtest(data, state, n+m-1, mod); 

        struct_nmod_toeplitz_init_set(mat, n, m, N, data);        
        
        // vector to multiply
        nmod_mat_t b;
        nmod_mat_init(b,m,l,N);
        nmod_mat_randtest(b, state);

        // structured part
        nmod_mat_t res_struct;
        nmod_mat_init(res_struct,n,l,N);
        struct_nmod_toeplitz_right_mul_mat(mat,b,res_struct); 

        // ---Dense part
        // Init 
        nmod_mat_t dense_mat;
        nmod_mat_init(dense_mat,n,m,N);
        struct_nmod_toeplitz_dense(mat,dense_mat);
        
        // Result
        nmod_mat_t res_dense;
        nmod_mat_init(res_dense,n,l,N);
        nmod_mat_mul(res_dense, dense_mat, b);

        // Check equality
        int result = nmod_mat_equal(res_struct, res_dense);
        if (!result) TEST_FUNCTION_FAIL("Structured result and dense result do not coincide\n");
    }    
    TEST_FUNCTION_END(state);
}


