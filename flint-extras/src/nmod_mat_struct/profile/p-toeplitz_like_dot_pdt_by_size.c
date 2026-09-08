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

    if(argc != 3){
        flint_printf("Usage: %s prime nrows ncols\n",argv[0]);
        return 0;
    }

    const slong N = atoi(argv[1]);
    if(!n_is_prime(N)){
        flint_printf("%d is not prime\n");
        return 0;
    }
    const slong rank = atoi(argv[2]);

    double t;
    clock_t tt;
    long nb_iter;

    flint_printf("Timings with growing displacement rank\n");
    flint_printf("Mat size |  Time (s)  | # Samples\n");
    for(int i=1; i<NB_ITER; i++){
        int m = i*i*rank ; //* FLINT_MIN(n, m) / NB_ITER;
        nb_iter = 0;
        
        t = 0.0;
        while(t < 2){
            nmod_toeplitz_like_t mat;
            nn_ptr v = _nmod_vec_init(m);
            nn_ptr res = _nmod_vec_init(m);
            nmod_toeplitz_like_randtest(mat, m, m, rank, N,state);
            _nmod_vec_randtest(v, state, m, mat->mod); 
            _nmod_vec_randtest(res, state, m, mat->mod); 
            // Structured
            tt = clock();
            nmod_toeplitz_like_mul_nmod_vec_left(mat, v, res);      
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;

            nb_iter += 1;
            // Cleaning memory
            nmod_toeplitz_like_clear(mat);
            flint_free(v);
            flint_free(res);
        }
        t /= nb_iter;
        flint_printf("%8d | %.4e | %d \n", m, t, nb_iter);
    }
    return 0;
}
