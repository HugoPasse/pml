/* 
   Copyright (C) 2026 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/


#ifndef MAPML_MATPOLY_DUMMY_EXPORT_H
#define MAPML_MATPOLY_DUMMY_EXPORT_H


#include "mapml_conversion.h"

#ifdef __cplusplus
extern "C" {
#endif



/**********************************************************
 * 
 * maple polynomial round trip (for cheks) 
 * 
 * * Converts a maple vector 
 *   to a maple list polynomial
 *
 *  ALGEB args[1]: polynomial string 
 *        args[2]: modulus 
 * 
 ***********************************************************/


ALGEB polynomial_rt(MKernelVector kv, ALGEB *args); 

/**********************************************************
 * 
 * Maple matrix polynomial round trip (for cheks) 
 * 
 * Converts a maple vector polynomial matrix 
 *   to a maple list polynomial matrix 
 *
 * Internal: mat_poly
 *  
 *  ALGEB args[1]: polynomial matrix: matrix of vectors 
 *        args[2]: modulus 
 * 
 * 
 ***********************************************************/


ALGEB matpoly_rt(MKernelVector kv, ALGEB *args); 



/**********************************************************
 * 
 * Maple matrix polynomial round trip (for cheks) 
 * 
 * Converts a maple vector polynomial matrix 
 *   to a maple list polynomial matrix 
 *
 * Internal: poly_mat
 *  
 *  ALGEB args[1]: polynomial matrix: matrix of vectors 
 *        args[2]: modulus 
 * 
 * 
 ***********************************************************/

ALGEB polymat_rt(MKernelVector kv, ALGEB *args); 


#ifdef __cplusplus
}
#endif

#endif


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
