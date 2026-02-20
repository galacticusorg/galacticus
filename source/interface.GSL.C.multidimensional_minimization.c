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

#include <gsl/gsl_multimin.h>

const gsl_multimin_fdfminimizer_type *gsl_multimin_fdfminimizer_type_get(int i) {
  /* Return a GSL interpolator type */
  const gsl_multimin_fdfminimizer_type *minimizer;
  switch(i) {
  case 1:
    minimizer = gsl_multimin_fdfminimizer_conjugate_fr;
    break;
  case 2:
    minimizer = gsl_multimin_fdfminimizer_conjugate_pr;
    break;
  case 3:
    minimizer = gsl_multimin_fdfminimizer_vector_bfgs2;
    break;
  case 4:
    minimizer = gsl_multimin_fdfminimizer_vector_bfgs;
    break;
  case 5:
    minimizer = gsl_multimin_fdfminimizer_steepest_descent;
    break;
  default:
    minimizer = NULL;
    break;
  }
  return minimizer;
}

const gsl_multimin_fminimizer_type *gsl_multimin_fminimizer_type_get(int i) {
  /* Return a GSL interpolator type */
  const gsl_multimin_fminimizer_type *minimizer;
  switch(i) {
  case 6:
    minimizer = gsl_multimin_fminimizer_nmsimplex2;
    break;
  case 7:
    minimizer = gsl_multimin_fminimizer_nmsimplex2rand;
    break;
  default:
    minimizer = NULL;
    break;
  }
  return minimizer;
}

gsl_multimin_function_fdf *gslMultiminFunctionFdFConstructor(size_t n, double (* f) (const gsl_vector * x, void * params), void (* df) (const gsl_vector * x, void * params, gsl_vector * g), void (* fdf) (const gsl_vector * x, void * params, double * f, gsl_vector * g)) {
  /* Construct a gsl_multimin_function_fdf object. */
  gsl_multimin_function_fdf *fGSL;
  fGSL           = (gsl_multimin_function_fdf *) malloc(sizeof(gsl_multimin_function_fdf));
  fGSL->f      = f;
  fGSL->df     = df;
  fGSL->fdf    = fdf;
  fGSL->n      = n;
  fGSL->params = NULL;
  return fGSL;
}

gsl_multimin_function *gslMultiminFunctionFConstructor(size_t n, double (* f) (const gsl_vector * x, void * params)) {
  /* Construct a gsl_multimin_function object. */
  gsl_multimin_function *fGSL;
  fGSL           = (gsl_multimin_function *) malloc(sizeof(gsl_multimin_function));
  fGSL->f      = f;
  fGSL->n      = n;
  fGSL->params = NULL;
  return fGSL;
}

void gslMultiminFunctionFdFDestructor(gsl_multimin_function_fdf *f) {
  /* Destroy a gsl_multimin_function_fdf object. */
  free(f);
}

void gslMultiminFunctionFDestructor(gsl_multimin_function *f) {
  /* Destroy a gsl_multimin_function object. */
  free(f);
}
