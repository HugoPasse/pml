#include <stdio.h>
#include "struct_nmod_mat.h"

TEST_FUNCTION_START(struct_nmod_toeplitz_like_pdt_vec, state)
{
    for(int i = 0; i < 100 * flint_test_multiplier(); i++)
    {
    
        long n = 15;
        long m = 7;
        long k = 2;
        ulong N = 101;
        
        nmod_toeplitz_like_t mat;
        nmod_toeplitz_like_randtest(mat, n, m, k, N, state);

        nmod_mat_t dense_mat;
        nmod_toeplitz_like_dense(mat,dense_mat);
        nmod_mat_print_pretty(dense_mat);
        printf("Lower generator");
        _nmod_vec_print_pretty(mat->lower_gens[0]->data, n, mat->mod);
        _nmod_vec_print_pretty(mat->lower_gens[1]->data, n, mat->mod);
        printf("Upper generator");
        _nmod_vec_print_pretty(mat->upper_gens[0]->data, m, mat->mod);
        _nmod_vec_print_pretty(mat->upper_gens[1]->data, m, mat->mod);

        nn_ptr v = _nmod_vec_init(m);
        _nmod_vec_randtest(v, state, m, mat->mod);
        printf("Vector to multiply:\n");
        _nmod_vec_print_pretty(v, m, mat->mod);

        nn_ptr res = _nmod_vec_init(n);
        nmod_toeplitz_like_mul_nmod_vec_left(mat, v, res);

        _nmod_vec_print_pretty(res, mat->nrows, mat->mod);


        nn_ptr dense_res = _nmod_vec_init(mat->nrows);
        nmod_mat_mul_nmod_vec(dense_res, dense_mat, v, mat->ncols);
        _nmod_vec_print_pretty(dense_res, mat->nrows, mat->mod);
        
        printf("---------------------\n");
    }    
    TEST_FUNCTION_END(state);
}


