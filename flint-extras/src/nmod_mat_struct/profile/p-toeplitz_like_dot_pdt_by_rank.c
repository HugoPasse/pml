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

    //slong primes[59] = {3, 7, 13, 29, 53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317, 196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843, 50331653, 100663319, 201326611, 402653189, 805306457, 1610612741, 3221225473, 6442450967, 12884901893, 25769803799, 51539607599, 103079215111, 206158430209, 412316860441, 824633720837, 1649267441681, 3298534883417, 6597069766657, 13194139533349, 26388279066671, 52776558133303, 105553116266509, 211106232533047, 422212465066001, 844424930132057, 1688849860263953, 3377699720527897, 6755399441055827, 13510798882111519, 27021597764223071, 54043195528445957, 108086391056891941, 216172782113783843, 432345564227567621, 864691128455135281};


    const slong bits = atoi(argv[1]);
    const slong N = n_randprime(state, bits, 1);
    const slong m = atoi(argv[2]);
    const slong n = atoi(argv[3]);

    double t;
    clock_t tt;
    long nb_iter;

    flint_printf("Timings with growing displacement rank\n");
    flint_printf("Disp. Rank |  Time (s)  | # Samples\n");
    for(int i=1; i<NB_ITER; i++){
        int rank = i*10 ; //* FLINT_MIN(n, m) / NB_ITER;
        nb_iter = 0;
        
        t = 0.0;
        while(t < 2){
            nmod_toeplitz_like_t mat;
            nn_ptr v = _nmod_vec_init(n);
            nn_ptr res = _nmod_vec_init(m);
            nmod_toeplitz_like_randtest(mat, m, n, rank, N,state);
            _nmod_vec_randtest(v, state, n, mat->mod); 
            _nmod_vec_randtest(res, state, m, mat->mod); 
            
            // Timing 
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
        flint_printf("%10d | %.4e | %d \n", rank, t, nb_iter);
    }
    return 0;
}
