#include <stdio.h>
#include "struct_nmod_mat.h"

TEST_FUNCTION_START(struct_nmod_toeplitz_like, state)
{
    for(int k = 0; k < 100 * flint_test_multiplier(); k++)
    {
    
        long n = n_randtest_bits(state,6);
        long m = n_randtest_bits(state,6);
        long k = n_randtest_bits(state,6);
        ulong N = n_randtest_prime(state, 1);
    
        nmod_mat_lower_toeplitz_t lmat;
        nmod_mat_upper_toeplitz_t umat;

        nmod_mat_lower_toeplitz_init(lmat, n, k, N);
        nmod_mat_upper_toeplitz_init(umat, k, m, N);

        nn_ptr ldata, udata;
        ldata = _nmod_vec_init(n);
        udata = _nmod_vec_init(m);
        
        _nmod_vec_randtest(ldata, state, n, lmat->mod);
        _nmod_vec_randtest(udata, state, m, umat->mod);
        
        nmod_mat_lower_toeplitz_set(lmat, ldata);
        nmod_mat_upper_toeplitz_set(umat, udata);

        nmod_mat_t lmat_dense, umat_dense;
        nmod_mat_lower_toeplitz_dense(lmat, lmat_dense);
        nmod_mat_upper_toeplitz_dense(umat, umat_dense);
        
        nmod_mat_t dense_pdt;
        nmod_mat_init(dense_pdt, n, m, N);
        nmod_mat_mul(dense_pdt, lmat_dense, umat_dense);
        
        for(int i=0; i<n; i++){
            for(int j=0; j<m; j++){
                ulong r0 = nmod_mat_lower_upper_toeplitz_pdt(lmat,umat,i,j);
                ulong r1 = nmod_mat_get_entry(dense_pdt,i,j);
                if(!(r0 == r1)){
                    TEST_FUNCTION_FAIL("Incorrect entry at (%d,%d))\n",i,j);
                }
                //printf("i: %d, j: %d, entries:%ld,%ld\n\n",i,j,r, nmod_mat_get_entry(dense_pdt,i,j));
            }
        }
    }
    
    int n = 100;
    int k = 10;
    int m = 70;
    int N = 17;
    nmod_toeplitz_like_t testmat;
    nmod_toeplitz_like_randtest(testmat, n, m, k, N, state);
    
    nmod_mat_t dense_test;
    nmod_mat_init(dense_test, n, m, N);
    
    nmod_toeplitz_like_dense(testmat, dense_test);
    //nmod_mat_print_pretty(dense_test);

    TEST_FUNCTION_END(state);
}


