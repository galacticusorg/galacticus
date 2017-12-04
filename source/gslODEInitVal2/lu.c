/* linalg/lu.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009 Gerard Jungman, Brian Gough
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

#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_linalg.h>

#define REAL double

/* Factorise a general N x N matrix A into,
 *
 *   P A = L U
 *
 * where P is a permutation matrix, L is unit lower triangular and U
 * is upper triangular.
 *
 * L is stored in the strict lower triangular part of the input
 * matrix. The diagonal elements of L are unity and are not stored.
 *
 * U is stored in the diagonal and upper triangular part of the
 * input matrix.  
 * 
 * P is stored in the permutation p. Column j of P is column k of the
 * identity matrix, where k = permutation->data[j]
 *
 * signum gives the sign of the permutation, (-1)^n, where n is the
 * number of interchanges in the permutation. 
 *
 * See Golub & Van Loan, Matrix Computations, Algorithm 3.4.1 (Gauss
 * Elimination with Partial Pivoting).
 *
 * Tests for zero entries, aij, and skips the final loop over k in such cases.
 */

int
gsl_linalg_LU_decomp_active (gsl_matrix * A, gsl_permutation * p, int *signum)
{
  if (A->size1 != A->size2)
    {
      GSL_ERROR ("LU decomposition requires square matrix", GSL_ENOTSQR);
    }
  else if (p->size != A->size1)
    {
      GSL_ERROR ("permutation length must match matrix size", GSL_EBADLEN);
    }
  else
    {
      const size_t N = A->size1;
      size_t i, j, k;

      *signum = 1;
      gsl_permutation_init (p);

      for (j = 0; j < N - 1; j++)
        {
          /* Find maximum in the j-th column */

          REAL ajj, max = fabs (gsl_matrix_get (A, j, j));
          size_t i_pivot = j;

          for (i = j + 1; i < N; i++)
            {
              REAL aij = fabs (gsl_matrix_get (A, i, j));

              if (aij > max)
                {
                  max = aij;
                  i_pivot = i;
                }
            }

          if (i_pivot != j)
            {
              gsl_matrix_swap_rows (A, j, i_pivot);
              gsl_permutation_swap (p, j, i_pivot);
              *signum = -(*signum);
            }

          ajj = gsl_matrix_get (A, j, j);

          if (ajj != 0.0)
            {
              for (i = j + 1; i < N; i++)
                {
                  REAL aij = gsl_matrix_get (A, i, j) / ajj;

		  /* If aij is zero then there is no modification to any matrix element, and we can skip the following */
		  if (aij != 0.0)
		    {
		      
		      gsl_matrix_set (A, i, j, aij);
		      
		      for (k = j + 1; k < N; k++)
			{
			  REAL aik = gsl_matrix_get (A, i, k);
			  REAL ajk = gsl_matrix_get (A, j, k);
			  gsl_matrix_set (A, i, k, aik - aij * ajk);
			}
		      
		    }
		  
                }
            }
        }
      
      return GSL_SUCCESS;
    }
}
