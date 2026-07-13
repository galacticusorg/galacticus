// Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
//           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
//    Andrew Benson <abenson@carnegiescience.edu>
//
// This file is part of Galacticus.
//
//    Galacticus is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Galacticus is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

/* Function that returns a pointer to the GSL default random number generator type. */

#include "gsl/gsl_rng.h"

/* Return the GSL default random number generator type. Accessing the `gsl_rng_default` global variable through this
   accessor (rather than binding to it directly from Fortran) ensures that the reference is resolved against the definition in
   the GSL shared library. Binding to the variable directly from Fortran emits it as a "common" symbol which the modern macOS
   linker refuses to reconcile with the definition exported by libgsl. */
const void * GSL_Get_Rng_Default(void)
{
  return gsl_rng_default;
}
