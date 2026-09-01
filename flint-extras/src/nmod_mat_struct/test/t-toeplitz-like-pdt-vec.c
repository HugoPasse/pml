#include <flint/nmod_mat.h>
#include <flint/nmod_types.h>
#include <flint/nmod_vec.h>
#include <flint/test_helpers.h>
#include "nmod_mat_struct.h"

TEST_FUNCTION_START(nmod_mat_toeplitz_like_pdt_vec, state)
{
    for(int i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        long n = 7;
        long m = 7;
        long k = 2;
        ulong N = 101;
        
        nmod_toeplitz_like_t mat;
        nmod_toeplitz_like_randtest(mat, n, m, k, N, state);

        nmod_mat_t dense_mat;
        nmod_mat_init(dense_mat, n, m, N);
        nmod_toeplitz_like_dense(mat, dense_mat);

        nn_ptr v = _nmod_vec_init(m);
        _nmod_vec_randtest(v, state, m, mat->mod);

        nn_ptr res = _nmod_vec_init(n);
        nmod_toeplitz_like_mul_nmod_vec_left(mat, v, res);

        nn_ptr dense_res = _nmod_vec_init(mat->nrows);

        nmod_mat_mul_nmod_vec(dense_res, dense_mat, v, mat->ncols);
        _nmod_vec_reduce(dense_res, dense_res, mat->ncols, mat->mod);
        
        int result = _nmod_vec_equal(res, dense_res, mat->nrows);
        if (!result) TEST_FUNCTION_FAIL("Structured result and dense result do not coincide\n");

        nmod_toeplitz_like_clear(mat);
        nmod_mat_clear(dense_mat);
        _nmod_vec_clear(v);
        _nmod_vec_clear(res);
        _nmod_vec_clear(dense_res);
    }    
    TEST_FUNCTION_END(state);
}


