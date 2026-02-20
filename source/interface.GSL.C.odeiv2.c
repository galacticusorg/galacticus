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

#include <gsl/gsl_odeiv2.h>
#include <gslODEInitVal2/gsl_odeiv2.h>
#include <float.h>

const gsl_odeiv2_step_type *gsl_odeiv2_step_type_get(int i) {
  /* Return a GSL ODE stepper type */
  const gsl_odeiv2_step_type *stepper;
  switch(i) {
  case  1:
    stepper = gsl_odeiv2_step_rk2;
    break;
  case  2:
    stepper = gsl_odeiv2_step_rk4;
    break;
  case  3:
    stepper = gsl_odeiv2_step_rkf45;
    break;
  case  4:
    stepper = gsl_odeiv2_step_rkck;
    break;
  case  5:
    stepper = gsl_odeiv2_step_rk8pd;
    break;
  case  6:
    stepper = gsl_odeiv2_step_rk1imp;
    break;
  case  7:
    stepper = gsl_odeiv2_step_rk2imp;
    break;
  case  8:
    stepper = gsl_odeiv2_step_rk4imp;
    break;
  case  9:
    stepper = gsl_odeiv2_step_bsimp;
    break;
  case 10:
    stepper = gsl_odeiv2_step_msadams;
    break;
  case 11:
    stepper = gsl_odeiv2_step_msbdf;
    break;
  case 12:
    stepper = gsl_odeiv2_step_msbdfactive;
    break;
  default:
    stepper = NULL;
    break;
  }
  return stepper;
}

gsl_odeiv2_system *gsl_odeiv2_system_init( 
					  size_t dimension,
					  int (*func)(double t, const double y[], double dydt[], void *p), 
					  int (*jacobian)(double t, const double y[], double * dfdy, double dfdt[], void *p)
					   ) {
  /* Fortran-callable wrapper to initialize a GSL ODEIV2 system object. */
  gsl_odeiv2_system *result;

  result            = (gsl_odeiv2_system *) malloc(sizeof(gsl_odeiv2_system));
  result->function  = func;
  result->jacobian  = jacobian;
  result->params    = NULL;	/* "params" is never used by Galacticus */
  result->dimension = dimension;
  return result;
}

void gsl_odeiv2_system_free(gsl_odeiv2_system *system) {
  /* Fortran-callable wrapper to free GSL ODEIV2 system objects. */
  free(system);
}

double gsl_odeiv2_driver_h (gsl_odeiv2_driver * d)
{
  /* Return the current step size */

  return d->h;
}

void gsl_odeiv2_driver_errors (gsl_odeiv2_driver * d, double yerr[])
{
  /* Return an array of the current errors in the ODE variables */
  int i;
  for (i = 0; i < d->sys->dimension; i++) {
    yerr[i]=d->e->yerr[i];
  }
  return;
}

void gsl_odeiv2_driver_init_errors (gsl_odeiv2_driver * d)
{
  /* Initialize errors */
  int i;
  for (i = 0; i < d->sys->dimension; i++) {
    d->e->yerr[i] = DBL_MAX;
  }
  return;
}
