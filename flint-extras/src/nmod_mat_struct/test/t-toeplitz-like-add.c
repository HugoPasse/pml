#include <stdio.h>
#include "nmod_mat_struct.h"

TEST_FUNCTION_START(struct_nmod_toeplitz_like_add, state)
{
    for(int k = 0; k < 100 * flint_test_multiplier(); k++)
    {
    
        ulong n = n_randtest_bits(state,6);
        ulong m = n_randtest_bits(state,6);
        ulong k0 = n_randtest_bits(state,4);
        ulong k1 = n_randtest_bits(state,4);

        ulong N = n_randtest_prime(state, 1);
    
        nmod_toeplitz_like_t mat0, mat1, mat_sum;
        nmod_mat_t mat0_dense, mat1_dense, mat_sum_dense, res_dense;

        nmod_toeplitz_like_randtest(mat0, n, m, k0, N, state);
        nmod_toeplitz_like_randtest(mat1, n, m, k1, N, state);
        nmod_toeplitz_like_add(mat0,mat1,mat_sum);
        nmod_toeplitz_like_dense(mat_sum, mat_sum_dense);

        nmod_mat_init(mat0_dense, n, m, N);
        nmod_mat_init(mat1_dense, n, m, N);
        nmod_mat_init(res_dense, n, m, N);
        nmod_toeplitz_like_dense(mat0, mat0_dense);
        nmod_toeplitz_like_dense(mat1, mat1_dense);
        nmod_mat_add(res_dense, mat0_dense, mat1_dense);

        if(!nmod_mat_equal(res_dense, mat_sum_dense)){
            TEST_FUNCTION_FAIL("Sum of Toeplitz like matrices do not coincide");
        }
    }

    TEST_FUNCTION_END(state);
}



