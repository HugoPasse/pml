#include <stdio.h>
#include "nmod_mat_struct.h"

TEST_FUNCTION_START(nmod_mat_toeplitz_dense, state)
{
    long n = 10;
    long m = 15;
    ulong N = 17;
    
    nmod_mat_toeplitz_t mat;
    nn_ptr data = _nmod_vec_init(n+m-1);
    for(int i=0; i<n+m-1; i++){
        data[i] = i;
    }
    
    nmod_mat_toeplitz_init_set(mat, n, m, N, data);        
    
    nmod_mat_t dense_mat;
    nmod_mat_init(dense_mat,n,m,N);
    nmod_mat_toeplitz_dense(mat,dense_mat);

    nmod_mat_print_pretty(dense_mat);
    
    nmod_mat_clear(dense_mat);

    TEST_FUNCTION_END(state);
}

