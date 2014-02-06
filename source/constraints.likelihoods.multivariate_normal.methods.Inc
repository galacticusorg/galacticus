!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!% Defines methods for a multivariate normal likelihood function.

function likelihoodMultivariateNormalConstructor(means,covariance)
  !% Constructor for multivariate-normal likelihood class.
  implicit none
  type            (likelihoodMultivariateNormal)                                :: likelihoodMultivariateNormalConstructor
  double precision                              , intent(in   ), dimension(:  ) :: means
  double precision                              , intent(in   ), dimension(:,:) :: covariance

  ! Store means and covariance.
  likelihoodMultivariateNormalConstructor%means     =means
  likelihoodMultivariateNormalConstructor%covariance=covariance
  ! Find the inverse of the covariance matrix.
  likelihoodMultivariateNormalConstructor%inverseCovariance=likelihoodMultivariateNormalConstructor%covariance%invert()
  return
end function likelihoodMultivariateNormalConstructor

double precision function likelihoodMultivariateNormalEvaluate(self,simulationState,temperature)
  !% Return the log-likelihood for a multivariate-normal likelihood function.
  use Constraints_State
  implicit none
  class           (likelihoodMultivariateNormal), intent(inout) :: self
  class           (state                       ), intent(in   ) :: simulationState
  double precision                              , intent(in   ) :: temperature
  type            (vector                      )                :: stateVector, difference

  stateVector                         =simulationState%get()
  difference                          =stateVector-self%means
  likelihoodMultivariateNormalEvaluate=-0.5d0*(difference*(self%inverseCovariance*difference))
  return
end function likelihoodMultivariateNormalEvaluate
