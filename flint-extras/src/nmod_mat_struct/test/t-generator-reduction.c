#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>
#include <flint/test_helpers.h>
#include "nmod_mat_struct.h"

TEST_FUNCTION_START(nmod_mat_toeplitz_like_reduce_generators, state)
{
    
    for(int k = 0; k < 100 * flint_test_multiplier(); k++)
    {
        ulong n = n_randtest_bits(state,3);
        ulong m = n_randtest_bits(state,3);
        ulong l = n_randtest_bits(state,2);
        ulong N = n_randtest_prime(state, 1);

        nmod_toeplitz_like_t mat;
        nmod_toeplitz_like_randtest(mat, n, m, l, N, state);
        
        nmod_mat_t densemat0;
        nmod_mat_init(densemat0, n, m, N);
        nmod_toeplitz_like_dense(mat, densemat0);

        flint_printf("Ngens: %d \n", mat->ngens);
        nmod_toeplitz_like_reduce_gens(mat);
        flint_printf("Ngens: %d \n\n", mat->ngens);

        nmod_mat_t densemat1;
        nmod_mat_init(densemat1, n, m, N);
        nmod_toeplitz_like_dense(mat, densemat1);
        
        //nmod_mat_print_pretty(densemat0);
        //nmod_mat_print_pretty(densemat1);

        if(!nmod_mat_equal(densemat0, densemat1)){
            TEST_FUNCTION_FAIL("Compression modified the matrix");
        }
        
        nmod_mat_clear(densemat0);
        nmod_mat_clear(densemat1);
        nmod_toeplitz_like_clear(mat);
    }
    TEST_FUNCTION_END(state);
}
