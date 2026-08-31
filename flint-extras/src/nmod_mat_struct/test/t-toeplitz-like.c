#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>
#include <stdio.h>
#include <flint/test_helpers.h>
#include "nmod_mat_struct.h"

TEST_FUNCTION_START(nmod_mat_toeplitz_like, state)
{
    //for(int k = 0; k < 100 * flint_test_multiplier(); k++)
    //{
    
        //long m = n_randtest_bits(state,6);
        //long n = n_randtest_bits(state,6);
        //long k = n_randtest_bits(state,6);
        //ulong N = n_randtest_prime(state, 1);
    
        long m = 15;
        long n = 10;
        long k = 2;
        long N = 101;

        nmod_toeplitz_like_t mat;
        nmod_toeplitz_like_randtest(mat, m, n, k, N, state);
        
        nmod_mat_print_pretty(mat->G);
        nmod_mat_print_pretty(mat->H);

        nmod_mat_t dense_mat;
        nmod_mat_init(dense_mat, m, n, N);
        nmod_toeplitz_like_dense(mat, dense_mat);
        
        flint_printf("Dense representation:\n");
        nmod_mat_print(dense_mat);

        nmod_mat_clear(dense_mat);
        nmod_toeplitz_like_clear(mat);
    //}
    TEST_FUNCTION_END(state);
}


