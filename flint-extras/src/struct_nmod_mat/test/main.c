/* Include functions *********************************************************/
#include <flint/test_helpers.h>
#include "t-reconstruction.c"
#include "t-toeplitz-vec.c"
#include "t-toeplitz-mat.c"
#include "t-toeplitz-add.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    //TEST_FUNCTION(struct_nmod_toeplitz_dense),
    TEST_FUNCTION(struct_nmod_toeplitz_mul_vec),
    TEST_FUNCTION(struct_nmod_toeplitz_mul_mat),
    TEST_FUNCTION(struct_nmod_toeplitz_add),
};

/* main function *************************************************************/

TEST_MAIN(tests)
