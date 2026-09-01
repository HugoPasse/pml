#include "nmod_mat_struct.h"
#include <flint/nmod_mat.h>


void nmod_mat_displacement_operator(nmod_mat_t mat,
                                    nmod_mat_t res){
    for (int i=0; i<nmod_mat_nrows(mat)-1; i++) {
        for (int j=0; j<nmod_mat_ncols(mat)-1; j++) {
            nmod_mat_set_entry(res, i+1, j+1, nmod_mat_entry(mat, i, j)); 
        }
    }

    for(int i=0; i<nmod_mat_nrows(mat)-1; i++){
        nmod_mat_set_entry(res, i+1, 0, nmod_mat_entry(mat, i, nmod_mat_ncols(mat)-1)); 
    }

    for (int i=0; i<nmod_mat_ncols(mat); i++) {
        nmod_mat_set_entry(res, 0, i, 0); 
    }
    nmod_mat_sub(res, mat, res);
}
