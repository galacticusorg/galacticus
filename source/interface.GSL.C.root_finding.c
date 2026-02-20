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

#include <gsl/gsl_roots.h>

const gsl_root_fsolver_type *gsl_fsolver_type_get(int i) {
  /* Return a GSL fsolver type */
  const gsl_root_fsolver_type *solver;
  switch(i) {
  case 1:
    solver = gsl_root_fsolver_bisection;
    break;
  case 2:
    solver = gsl_root_fsolver_brent;
    break;
  case 3:
    solver = gsl_root_fsolver_falsepos;
    break;
  default:
    solver = NULL;
    break;
  }
  return solver;
}

const gsl_root_fdfsolver_type *gsl_fdfsolver_type_get(int i) {
  /* Return a GSL fsolver type */
  const gsl_root_fdfsolver_type *solver;
  switch(i) {
  case 4:
    solver = gsl_root_fdfsolver_newton;
    break;
  case 5:
    solver = gsl_root_fdfsolver_secant;
    break;
  case 6:
    solver = gsl_root_fdfsolver_steffenson;
    break;
  default:
    solver = NULL;
    break;
  }
  return solver;
}
