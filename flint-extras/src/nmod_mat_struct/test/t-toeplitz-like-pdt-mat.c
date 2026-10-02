#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_types.h>
#include <flint/nmod_vec.h>
#include <flint/test_helpers.h>
#include <flint/ulong_extras.h>
#include "nmod_mat_struct.h"

TEST_FUNCTION_START(nmod_mat_toeplitz_like_pdt_mat, state)
{
    for(int i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        long n = n_randtest_bits(state, 7);
        long m = n_randtest_bits(state, 7);
        long l = n_randtest_bits(state, 4);

        long k = n_randtest_bits(state, 4);
        ulong N = n_randtest_prime(state, 1);
        
        nmod_toeplitz_like_t mat;
        nmod_toeplitz_like_randtest(mat, n, m, k, N, state);

        nmod_mat_t dense_mat;
        nmod_mat_init(dense_mat, n, m, N);
        nmod_toeplitz_like_dense(mat, dense_mat);

        nmod_mat_t b;
        nmod_mat_init(b, m, l, N);
        nmod_mat_randtest(b, state);

        nmod_mat_t res;
        nmod_mat_init(res, n, l, N);
        nmod_toeplitz_like_mul_nmod_mat_left(mat, b, res);

        nmod_mat_t dense_res;
        nmod_mat_init(dense_res, n, l, N);
        nmod_mat_mul(dense_res, dense_mat, b);

        int result = nmod_mat_equal(res, dense_res);
        if (!result) TEST_FUNCTION_FAIL("Structured result and dense result do not coincide\n");
        
        nmod_toeplitz_like_clear(mat);
        nmod_mat_clear(dense_mat);
        nmod_mat_clear(b);
        nmod_mat_clear(res);
        nmod_mat_clear(dense_res);
    }    
    TEST_FUNCTION_END(state);
}
