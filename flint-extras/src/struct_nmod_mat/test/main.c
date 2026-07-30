/* Include functions *********************************************************/
#include <flint/test_helpers.h>
#include "t-reconstruction.c"
#include "t-toeplitz-vec.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(struct_nmod_toeplitz_dense),
    TEST_FUNCTION(struct_nmod_toeplitz_mul_vec),
};

/* main function *************************************************************/

TEST_MAIN(tests)
