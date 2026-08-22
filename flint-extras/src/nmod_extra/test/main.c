/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML. It is adapted from the files src/flint.h.in
    in FLINT (GNU LGPL version 3 or later).

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "pml.h"
#if PML_HAVE_AVX2
# include "t-multimod_CRT_reduce.c"
#endif  /* PML_HAVE_AVX2 */
#if PML_HAVE_MACHINE_VECTORS
# include "t-multimod_CRT_CRT.c"
#endif

/* Array of test functions ***************************************************/

#if PML_HAVE_AVX2 || PML_HAVE_MACHINE_VECTORS

test_struct tests[] =
{
#if PML_HAVE_AVX2
    TEST_FUNCTION(nmod_multimod_CRT_reduce),
#endif  /* PML_HAVE_AVX2 */
#if PML_HAVE_MACHINE_VECTORS
    TEST_FUNCTION(nmod_multimod_CRT_CRT),
#endif  /* PML_HAVE_MACHINE_VECTORS */
};

/* main function *************************************************************/

TEST_MAIN(tests)

#else
    int main() {}
#endif

