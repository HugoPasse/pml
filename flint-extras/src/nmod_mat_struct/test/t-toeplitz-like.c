#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>
#include <flint/test_helpers.h>
#include "nmod_mat_struct.h"

TEST_FUNCTION_START(nmod_mat_toeplitz_like, state)
{
    for(int k = 0; k < 100 * flint_test_multiplier(); k++)
    {
    
        long m = n_randtest_bits(state,6);
        long n = n_randtest_bits(state,6);
        long k = n_randtest_bits(state,6);
        ulong N = n_randtest_prime(state, 1);
    
        nmod_toeplitz_like_t mat;
        nmod_toeplitz_like_randtest(mat, m, n, k, N, state);
        
        nmod_mat_t dense_mat;
        nmod_mat_init(dense_mat, m, n, N);
        nmod_toeplitz_like_dense(mat, dense_mat);
        
        nmod_mat_t disp_mat;
        nmod_mat_init(disp_mat, m, n, N);
        nmod_mat_displacement_operator(dense_mat, disp_mat);

        nmod_mat_t GH, Ht;
        nmod_mat_init(Ht, k , n , N);
        nmod_mat_transpose(Ht, mat->H);
        nmod_mat_init(GH, m , n , N);
        nmod_mat_mul(GH, mat->G, Ht);
        
        int result = nmod_mat_equal(GH, disp_mat);
        if (!result) TEST_FUNCTION_FAIL("Structured result and dense result do not coincide\n");

        nmod_mat_clear(dense_mat);
        nmod_mat_clear(disp_mat);
        nmod_mat_clear(GH);
        nmod_mat_clear(Ht);
        
            
        nmod_toeplitz_like_clear(mat);
    }
    TEST_FUNCTION_END(state);
}


