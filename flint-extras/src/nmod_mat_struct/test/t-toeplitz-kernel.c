#include "nmod_mat_struct.h"

#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>
#include <flint/test_helpers.h>
#include <stdio.h>

TEST_FUNCTION_START(nmod_mat_toeplitz_kernel, state)
{
    for(int i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong nr = n_randtest_bits(state,3);
        ulong nc = nr+4;
        ulong N = n_randtest_prime(state, 1);

        nmod_mat_toeplitz_t mat;
        nmod_mat_toeplitz_randtest(mat,nr,nc, N, state);
        
        nmod_mat_t dense_mat;
        nmod_mat_init(dense_mat, nr, nc, N);
        nmod_mat_toeplitz_dense(mat, dense_mat);
        //nmod_mat_print(dense_mat); flint_printf("\n");

        nmod_mat_t res;
        nmod_mat_toeplitz_right_kernel_basis(mat, res);

        nmod_mat_t pdt;
        nmod_mat_init(pdt, nr, nc-nr, N);
        nmod_mat_mul(pdt, dense_mat, res);
        
        int result = nmod_mat_is_zero(pdt) && (nmod_mat_rank(res) == nmod_mat_ncols(res)) && (nmod_mat_ncols(res) == (mat->ncols - nmod_mat_rank(dense_mat)));
        if (!result) TEST_FUNCTION_FAIL("Structured result and dense result do not coincide\n");

        nmod_mat_clear(dense_mat);
        nmod_mat_clear(res);
        nmod_mat_clear(pdt);

        nmod_mat_toeplitz_clear(mat);
    }
    TEST_FUNCTION_END(state);
}
