// Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
//           2019
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

//% Implements Fortran-callable wrappers around GSL functions.

#include <stdlib.h>
#include <gsl/gsl_math.h>

gsl_function *gslFunctionConstructor(double (*f) (double x, void * params)) {
  /* Construct a gsl_function object. */
  gsl_function *fGSL;
  fGSL           = (gsl_function *) malloc(sizeof(gsl_function));
  fGSL->function = f;
  fGSL->params   = NULL;
  return fGSL;
}

void gslFunctionDestructor(gsl_function *f) {
  /* Destroy a gsl_function object. */
  free(f);
}
