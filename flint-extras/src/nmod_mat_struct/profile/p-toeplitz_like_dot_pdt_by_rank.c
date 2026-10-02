#include <flint/nmod_mat.h>
#include <stdlib.h>
#include <flint/flint.h>
#include <flint/profiler.h>
#include <flint/ulong_extras.h>
#include <time.h>

#include <flint/nmod_vec.h>

#include "nmod_mat_struct.h"

#define NB_ITER 10

int main(int argc, char** argv){
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125L);

    if(argc != 4){
        flint_printf("Usage: %s nbits nrows ncols\n",argv[0]);
        return 0;
    }

    const slong bits = atoi(argv[1]);
    const slong N = n_randprime(state, bits, 1);
    const slong m = atoi(argv[2]);
    const slong n = atoi(argv[3]);

    flint_printf("Timings with growing displacement rank\n");
    flint_printf("Prime: %u, (%dx%d)-matrices\n",N,m,n);
    
    flint_printf("Disp. Rank |  Time (s) \n");
    for(int i=1; i<NB_ITER; i++){
        int rank = i*10 ; //* FLINT_MIN(n, m) / NB_ITER;
       
        // Setup
        nmod_toeplitz_like_t mat;
        nn_ptr v = _nmod_vec_init(n);
        nn_ptr res = _nmod_vec_init(m);
        nmod_toeplitz_like_randtest(mat, m, n, rank, N,state);
        _nmod_vec_randtest(v, state, n, mat->mod); 
        _nmod_vec_randtest(res, state, m, mat->mod); 
        
        // Timing 
        double FLINT_SET_BUT_UNUSED(tcpu), twall;
        TIMEIT_START
        nmod_toeplitz_like_mul_nmod_vec_left(mat, v, res);      
        TIMEIT_STOP_VALUES(tcpu, twall);

        // Cleaning memory
        nmod_toeplitz_like_clear(mat);
        flint_free(v);
        flint_free(res);

        flint_printf("%10d | %.4e \n", rank, twall);
    }
    return 0;
}
