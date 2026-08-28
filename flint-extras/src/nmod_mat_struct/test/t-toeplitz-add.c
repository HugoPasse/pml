#include <flint/flint.h>
#include <flint/test_helpers.h>
#include "nmod_mat_struct.h"

TEST_FUNCTION_START(nmod_mat_toeplitz_add, state)
{
    for(int i = 0; i < 100 * flint_test_multiplier(); i++)
    {
    
        long n = n_randtest_bits(state,6);
        long m = n_randtest_bits(state,6);
        ulong N = n_randtest_prime(state, 1);
        
        nmod_mat_toeplitz_t mat1, mat2;
        nn_ptr data1 = _nmod_vec_init(n+m-1);
        nn_ptr data2 = _nmod_vec_init(n+m-1);

        nmod_t mod;
        nmod_init(&mod,N);
        _nmod_vec_randtest(data1, state, n+m-1, mod); 
        _nmod_vec_randtest(data2, state, n+m-1, mod); 

        nmod_mat_toeplitz_init_set(mat1, n, m, N, data1);        
        nmod_mat_toeplitz_init_set(mat2, n, m, N, data2);        
        
        // Structured part
        nmod_mat_toeplitz_t res_struct;
        nmod_mat_toeplitz_init(res_struct,n,m,N);
        nmod_mat_toeplitz_add(mat1, mat2, res_struct);

        nmod_mat_t res_struct_dense;
        nmod_mat_init(res_struct_dense,n,m,N);
        nmod_mat_toeplitz_dense(res_struct, res_struct_dense);
        

        // ---Dense part
        // Init 
        nmod_mat_t dense_mat1, dense_mat2;
        nmod_mat_init(dense_mat1,n,m,N);
        nmod_mat_init(dense_mat2,n,m,N);
        nmod_mat_toeplitz_dense(mat1,dense_mat1);
        nmod_mat_toeplitz_dense(mat2,dense_mat2);
        
        // Result
        nmod_mat_t res_dense;
        nmod_mat_init(res_dense,n,m,N);
        nmod_mat_add(res_dense, dense_mat1, dense_mat2);

        // Check equality
        int result = nmod_mat_equal(res_struct_dense, res_dense);
        if (!result) TEST_FUNCTION_FAIL("Structured result and dense result do not coincide\n");
    
        flint_free(data1);
        flint_free(data2);

        nmod_mat_toeplitz_clear(mat1);
        nmod_mat_toeplitz_clear(mat2);
        nmod_mat_toeplitz_clear(res_struct);

        nmod_mat_clear(res_dense);
        nmod_mat_clear(res_struct_dense);
        nmod_mat_clear(dense_mat1);
        nmod_mat_clear(dense_mat2);
    }    
    TEST_FUNCTION_END(state);
}


