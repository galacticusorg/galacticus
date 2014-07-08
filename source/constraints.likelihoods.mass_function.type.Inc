
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

  !% Defines a likelihood function for fitting of mass functions.

  type, public, extends(likelihood) :: likelihoodMassFunction
     private
     logical                                               :: useSurveyLimits       , modelSurfaceBrightness
     double precision                                      :: haloMassMinimum       , haloMassMaximum       , &
          &                                                   redshiftMinimum       , redshiftMaximum       , &
          &                                                   logHaloMassMinimum    , logHaloMassMaximum    , &
          &                                                   surfaceBrightnessLimit
     double precision        , dimension(:  ), allocatable :: mass                  , massFunctionObserved  , &
          &                                                   massMinimum           , massMaximum           , &
          &                                                   massFunction
     double precision        , dimension(:,:), allocatable :: covarianceMatrix
     type            (vector)                              :: means
     type            (matrix)                              :: covariance            , inverseCovariance
   contains
     procedure :: evaluate        => likelihoodMassFunctionEvaluate
     procedure :: functionChanged => likelihoodMassFunctionFunctionChanged
  end type likelihoodMassFunction
  
  interface likelihoodMassFunction
     !% Constructors for the mass function likelihood object.
     module procedure likelihoodMassFunctionConstructor
  end interface likelihoodMassFunction
