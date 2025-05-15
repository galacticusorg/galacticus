/* ode-initval2/driver.c
 * 
 * Copyright (C) 2009, 2010 Tuomo Keskitalo
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

/* Driver routine for odeiv2. This is a wrapper for low level GSL
   functions that allows a simple interface to step, control and
   evolve layers.
 */

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gslODEInitVal2/gsl_odeiv2.h>
#include <gsl/gsl_machine.h>

static gsl_odeiv2_driver *
driver_alloc (const gsl_odeiv2_system * sys, const double hstart,
              const gsl_odeiv2_step_type * T)
{
  /* Allocates and initializes an ODE driver system. Step and evolve
     objects are allocated here, but control object is allocated in
     another function.
   */

  gsl_odeiv2_driver *state =
    (gsl_odeiv2_driver *) malloc (sizeof (gsl_odeiv2_driver));

  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for driver state",
                      GSL_ENOMEM);
    }

  if (sys == NULL)
    {
      GSL_ERROR_NULL ("gsl_odeiv2_system must be defined", GSL_EINVAL);
    }

  {
    const size_t dim = sys->dimension;

    if (dim == 0)
      {
        GSL_ERROR_NULL
          ("gsl_odeiv2_system dimension must be a positive integer",
           GSL_EINVAL);
      }

    state->sys = sys;

    state->s = gsl_odeiv2_step_alloc (T, dim);

    if (state->s == NULL)
      {
        free (state);
        GSL_ERROR_NULL ("failed to allocate step object", GSL_ENOMEM);
      }

    state->e = gsl_odeiv2_evolve_alloc (dim);
  }

  if (state->e == NULL)
    {
      gsl_odeiv2_step_free (state->s);
      free (state);
      GSL_ERROR_NULL ("failed to allocate evolve object", GSL_ENOMEM);
    }

  if (hstart > 0.0 || hstart < 0.0)
    {
      state->h = hstart;
    }
  else
    {
      GSL_ERROR_NULL ("invalid hstart", GSL_EINVAL);
    }

  state->h = hstart;
  state->hmin = 0.0;
  state->hmax = GSL_DBL_MAX;
  state->nmax = 0;
  state->n = 0;
  state->c = NULL;

  return state;
}

gsl_odeiv2_driver *
gsl_odeiv2_driver_alloc_scaled2_new (const gsl_odeiv2_system * sys,
                                    const gsl_odeiv2_step_type * T,
                                    const double hstart,
                                    const double epsabs, const double epsrel,
                                    const double a_y, const double a_dydt,
				    const double scale_abs[], const int is_non_negative[])
{
  /* Initializes an ODE driver system with control object of type
     scaled_new. 
   */

  gsl_odeiv2_driver *state = driver_alloc (sys, hstart, T);

  if (state == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate driver object", GSL_ENOMEM);
    }

  if (epsabs >= 0.0 && epsrel >= 0.0)
    {
      state->c = gsl_odeiv2_control_scaled2_new (epsabs, epsrel, a_y, a_dydt,
						 scale_abs, is_non_negative, sys->dimension);

      if (state->c == NULL)
        {
          gsl_odeiv2_driver_free (state);
          GSL_ERROR_NULL ("failed to allocate control object", GSL_ENOMEM);
        }
    }
  else
    {
      gsl_odeiv2_driver_free (state);
      GSL_ERROR_NULL ("epsabs and epsrel must be positive", GSL_EINVAL);
    }

  /* Distribute pointer to driver object */

  gsl_odeiv2_step_set_driver (state->s, state);
  gsl_odeiv2_evolve_set_driver (state->e, state);
  gsl_odeiv2_control_set_driver (state->c, state);

  return state;
}

int
gsl_odeiv2_driver2_apply (gsl_odeiv2_driver * d, double *t,
			  const double t1, double y[],
			  void(*postStep)(double t, double y[], int *s),
			  void(*latentIntegrator)(double *t),
			  void(*stepAnalyzer)(double t, double t1, double y[], double yerr[], double h, int s)
			  )
{
  /* Main driver function that evolves the system from t to t1. In
     beginning vector y contains the values of dependent variables at
     t. This function returns values at t=t1 in y. In case of
     unrecoverable error, y and t contains the values after the last
     successful step.
   */

  /* Modified by Andrew Benson to allow call to a function to integrate latent variables after each time step */
  
  int sign = 0;
  d->n = 0;

  /* Determine integration direction sign */

  if (d->h > 0.0)
    {
      sign = 1;
    }
  else
    {
      sign = -1;
    }

  /* Check that t, t1 and step direction are sensible */

  if (sign * (t1 - *t) < 0.0)
    {
      GSL_ERROR_NULL
        ("integration limits and/or step direction not consistent",
         GSL_EINVAL);
    }

  /* Evolution loop */

  while (sign * (t1 - *t) > 0.0)
    {
      int s = gsl_odeiv2_evolve_apply (d->e, d->c, d->s, d->sys,
                                       t, t1, &(d->h), y);

      if (stepAnalyzer != NULL ) 
	{
	  stepAnalyzer(*t,t1,y,d->e->yerr,d->e->last_step,s);
	}

      if (s != GSL_SUCCESS)
        {
          return s;
        }

      /* Call function to perform post-step behavior. */
      if (postStep != NULL)
	{
	  int s = GSL_SUCCESS;
	  postStep(*t,y,&s);
	  if (s != 0)
	    {
	      /* Post-step processing returned failure, which indicates that the state of the system was changed. We must reset
		 the evolver (specifically setting d->e->count = 0) to avoid re-using the dydt's computed at the end of the
		 current step at the start of the next step, as they could be invalid due to this change of state. */
	      gsl_odeiv2_evolve_reset(d->e);
	    }
	}
      
      /* Call function to perform integration of latent variables */
      if (latentIntegrator != NULL)
	{
	  latentIntegrator(t); 
	}
      
      /* Check for maximum allowed steps */

      if ((d->nmax > 0) && (d->n > d->nmax))
        {
          return GSL_EMAXITER;
        }

      /* Set step size if maximum size is exceeded */

      if (fabs (d->h) > d->hmax)
        {
          d->h = sign * d->hmax;
        }

      /* Check for too small step size */

      if (fabs (d->h) < d->hmin)
        {
          return GSL_ENOPROG;
        }

      d->n++;
    }

  return GSL_SUCCESS;
}
