/* ode-initval/odeiv2.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */
/* Modified by Tuomo Keskitalo */

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_types.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>

/* Additional stepper types */

GSL_VAR const gsl_odeiv2_step_type *gsl_odeiv2_step_msbdfactive;

/* LU Decomposition, Gaussian elimination with partial pivoting
 */

int gsl_linalg_LU_decomp_active (gsl_matrix * A, gsl_permutation * p, int *signum);

/* This controller computes errors using different absolute errors for
 * each component
 *
 *    D0 = eps_abs * scale_abs[i] + eps_rel * (a_y |y| + a_dydt h |y'|)
 */

gsl_odeiv2_control *gsl_odeiv2_control_scaled2_new (double eps_abs,
						    double eps_rel, double a_y,
						    double a_dydt,
						    const double scale_abs[],
						    const int is_non_negative[],
						    size_t dim);

int gsl_odeiv2_driver2_apply (gsl_odeiv2_driver * d, double *t,
			      const double t1, double y[],
			      void(*postStep)(double t, double y[], int *s),
			      void(*latentIntegrator)(double *t),
			      void(*stepAnalyzer)(double t, double t1, double y[], double yerr[], double h, int s)
			      );
