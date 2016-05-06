//
// File:        io.h
// Package:     STRANDGEN
//              Parallel Infrastructure for Cartesian & Strand Solvers
// Release:     $Name:  $
// Revision:    $Revision:  $
// Modified:    $Date:  $
// Description: Single file for defining io for C++ code.
//              

/*!
 * C++ uses a number of io packages - stdio, iostream, iomanip, etc.
 * Rather than try to figure out which io operation needs which library,
 * the intent of this include is to include all of them, and keep a 
 * running list of which operations are supported by which package.
 */

#ifndef included_io
#define included_io

#ifndef included_STRANDGEN_defs
#define included_STRANDGEN_defs
#include "STRANDGEN_defs.h"
#endif

#ifndef included_stdio
#define included_stdio
#include <stdio.h>
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#ifndef included_iomanip
#define included_iomanip
#include <iomanip>
#endif

#ifndef included_sstream
#define included_sstream
#include <sstream>
#endif

#endif
