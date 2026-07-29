/* Include functions *********************************************************/
#include <flint/test_helpers.h>
#include "t-reconstruction.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(struct_nmod_toeplitz_dense),
};

/* main function *************************************************************/

TEST_MAIN(tests)
