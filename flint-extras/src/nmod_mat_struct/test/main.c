/* Include functions *********************************************************/
#include <flint/test_helpers.h>
#include "t-reconstruction.c"
#include "t-toeplitz-vec.c"
#include "t-toeplitz-mat.c"
#include "t-toeplitz-add.c"
#include "t-toeplitz-like.c"
#include "t-toeplitz-like-add.c"
#include "t-toeplitz-like-pdt-vec.c"
#include "t-toeplitz-kernel.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    //TEST_FUNCTION(struct_nmod_toeplitz_dense),
    TEST_FUNCTION(nmod_mat_toeplitz_mul_vec),
    TEST_FUNCTION(nmod_mat_toeplitz_mul_mat),
    TEST_FUNCTION(nmod_mat_toeplitz_add),
    TEST_FUNCTION(nmod_mat_toeplitz_kernel),
    TEST_FUNCTION(nmod_mat_toeplitz_like),
    TEST_FUNCTION(nmod_mat_toeplitz_like_add),
    // This does not work yet due to incorrect dense reconstruction
    //TEST_FUNCTION(struct_nmod_toeplitz_like_pdt_vec)
};

/* main function *************************************************************/

TEST_MAIN(tests)
