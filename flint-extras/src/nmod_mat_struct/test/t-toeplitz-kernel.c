#include "nmod_mat_struct.h"

#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/test_helpers.h>
#include <stdio.h>

TEST_FUNCTION_START(nmod_mat_toeplitz_kernel, state)
{
    ulong N = 101;
    ulong nr = 5;
    ulong nc = 7;

    nmod_mat_toeplitz_t mat;

    nmod_mat_toeplitz_randtest(mat,nr,nc, N, state);
    
    nmod_mat_t dense_mat;
    nmod_mat_init(dense_mat, nr, nc, N);
    nmod_mat_toeplitz_dense(mat, dense_mat);
    nmod_mat_print_pretty(dense_mat);
   
    nmod_mat_t res;
    nmod_mat_init(res, nc, nc-nr, N);
    nmod_mat_zero(res);
    nmod_mat_print_pretty(res);

    printf("ici\n");
    nmod_mat_toeplitz_right_kernel_basis(mat, res);
    nmod_mat_clear(dense_mat);
    nmod_mat_toeplitz_clear(mat);

    nmod_mat_t pdt;
    nmod_mat_init(pdt, nr, nc-nr, N);
    nmod_mat_mul(pdt, dense_mat, res);
    flint_printf("Product: ");
    nmod_mat_print_pretty(pdt);

    TEST_FUNCTION_END(state);
}
