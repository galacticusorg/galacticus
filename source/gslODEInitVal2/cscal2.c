/* ode-initval2/cscal.c
 * 
 * Copyright (C) 2002, 2007 Brian Gough
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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>

typedef struct
{
  double eps_abs;
  double eps_rel;
  double a_y;
  double a_dydt;
  double *scale_abs;
  int    *is_non_negative;
}
sc2_control_state_t;

static int
control_set_driver_null (void *vstate, const gsl_odeiv2_driver * d)
{
  /* Dummy set function for those control objects that do not
     need pointer to driver object. */

  return GSL_SUCCESS;
}

static void *
sc2_control_alloc (void)
{
  sc2_control_state_t *s =
    (sc2_control_state_t *) malloc (sizeof (sc2_control_state_t));

  if (s == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for sc2_control_state",
                      GSL_ENOMEM);
    }

  return s;
}

static int
sc2_control_init (void *vstate,
                 double eps_abs, double eps_rel, double a_y, double a_dydt)
{
  sc2_control_state_t *s = (sc2_control_state_t *) vstate;

  if (eps_abs < 0)
    {
      GSL_ERROR ("eps_abs is negative", GSL_EINVAL);
    }
  else if (eps_rel < 0)
    {
      GSL_ERROR ("eps_rel is negative", GSL_EINVAL);
    }
  else if (a_y < 0)
    {
      GSL_ERROR ("a_y is negative", GSL_EINVAL);
    }
  else if (a_dydt < 0)
    {
      GSL_ERROR ("a_dydt is negative", GSL_EINVAL);
    }

  s->eps_rel = eps_rel;
  s->eps_abs = eps_abs;
  s->a_y = a_y;
  s->a_dydt = a_dydt;

  return GSL_SUCCESS;
}

static int
sc2_control_hadjust (void *vstate, size_t dim, unsigned int ord,
                    const double y[], const double yerr[], const double yp[],
                    double *h)
{
  sc2_control_state_t *state = (sc2_control_state_t *) vstate;

  const double eps_abs = state->eps_abs;
  const double eps_rel = state->eps_rel;
  const double a_y = state->a_y;
  const double a_dydt = state->a_dydt;
  const double *scale_abs = state->scale_abs;
  const int    *is_non_negative = state->is_non_negative;

  const double S = 0.9;
  const double h_old = *h;

  double rmax = DBL_MIN;
  size_t i;
  int forbiddenNegatives = 0;

  for (i = 0; i < dim; i++)
    {
      const double D0 =
        eps_rel * (a_y * fabs (y[i]) + a_dydt * fabs (h_old * yp[i]))
        + eps_abs * scale_abs[i];
      const double r = fabs (yerr[i]) / fabs (D0);
      rmax = GSL_MAX_DBL (r, rmax);
      /* Record if we have negative values for any properties that must be non-negative. */
      if (is_non_negative[i] && y[i] < 0.0) {
	forbiddenNegatives = 1;
      }
    }

  if (rmax > 1.1)
    {
      /* decrease step, no more than factor of 5, but a fraction S more
         than scaling suggests (for better accuracy) */
      double r = S / pow (rmax, 1.0 / ord);

      if (r < 0.2)
        r = 0.2;
      
      *h = r * h_old;

      return GSL_ODEIV_HADJ_DEC;
    }
  else if (forbiddenNegatives == 1)
    {
      /* We have negative values for properties that must be non-negative - reduced the step size by a factor 2. */
      *h = 0.5 * h_old;

      return GSL_ODEIV_HADJ_DEC;
    }
  else if (rmax < 0.5)
    {
      /* increase step, no more than factor of 4.9. The upper limit of 4.9 is different from the original cscal implementation
	 choice of 5.0. This is because for an increase of precisely 5.0 the msbdf stepper algorithm can fail (divide by zero
	 error) in cases where it is using 3rd order solutions */
      double r = S / pow (rmax, 1.0 / (ord + 1.0));

      if (r > 4.9)
        r = 4.9;

      if (r < 1.0)              /* don't allow any decrease caused by S<1 */
        r = 1.0;

      *h = r * h_old;

      return GSL_ODEIV_HADJ_INC;
    }
  else
    {
      /* no change */
      return GSL_ODEIV_HADJ_NIL;
    }
}

static int
sc2_control_errlevel (void *vstate, const double y, const double dydt,
                     const double h, const size_t ind, double *errlev)
{
  sc2_control_state_t *state = (sc2_control_state_t *) vstate;

  const double eps_abs = state->eps_abs;
  const double eps_rel = state->eps_rel;
  const double a_y = state->a_y;
  const double a_dydt = state->a_dydt;
  const double *scale_abs = state->scale_abs;

  *errlev = eps_rel * (a_y * fabs (y) + a_dydt * fabs (h * dydt))
    + eps_abs * scale_abs[ind];

  if (*errlev <= 0.0)
    {
      GSL_ERROR_NULL ("errlev <= zero", GSL_ESANITY);
    }

  return GSL_SUCCESS;
}


static void
sc2_control_free (void *vstate)
{
  sc2_control_state_t *state = (sc2_control_state_t *) vstate;
  free (state->scale_abs);
  free (state->is_non_negative);
  free (state);
}

static const gsl_odeiv2_control_type sc2_control_type = { "scaled",      /* name */
  &sc2_control_alloc,
  &sc2_control_init,
  &sc2_control_hadjust,
  &sc2_control_errlevel,
  &control_set_driver_null,
  &sc2_control_free
};

const gsl_odeiv2_control_type *gsl_odeiv2_control_scaled2 = &sc2_control_type;


gsl_odeiv2_control *
gsl_odeiv2_control_scaled2_new (double eps_abs, double eps_rel,
                               double a_y, double a_dydt,
				const double scale_abs[], const int is_non_negative[], size_t dim)
{
  gsl_odeiv2_control *c =
    gsl_odeiv2_control_alloc (gsl_odeiv2_control_scaled2);

  int status = gsl_odeiv2_control_init (c, eps_abs, eps_rel, a_y, a_dydt);

  if (status != GSL_SUCCESS)
    {
      gsl_odeiv2_control_free (c);
      GSL_ERROR_NULL ("error trying to initialize control", status);
    }

  {
    sc2_control_state_t *s = (sc2_control_state_t *) c->state;

    s->scale_abs = (double *) malloc (dim * sizeof (double));

    if (s->scale_abs == 0)
      {
        free (s);
        GSL_ERROR_NULL ("failed to allocate space for scale_abs", GSL_ENOMEM);
      }

    memcpy (s->scale_abs, scale_abs, dim * sizeof (double));

    s->is_non_negative = (int *) malloc (dim * sizeof (int));

    if (s->is_non_negative == 0)
      {
        free (s);
        GSL_ERROR_NULL ("failed to allocate space for is_non_negative", GSL_ENOMEM);
      }

    memcpy (s->is_non_negative, is_non_negative, dim * sizeof (int));
}

  return c;
}
