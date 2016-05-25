!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Defines a likelihood function for fitting of halo mass functions.

  type, public, extends(likelihood) :: likelihoodHaloMassFunction
     private
     double precision        , dimension(:  ), allocatable :: mass            , massFunction     , &
          &                                                   massMinimum     , massMaximum
     double precision        , dimension(:,:), allocatable :: covarianceMatrix
     double precision                                      :: time            , massParticle
     type            (vector)                              :: means
     type            (matrix)                              :: covariance      , inverseCovariance
     integer                                               :: errorModel
   contains
     procedure :: evaluate        => likelihoodHaloMassFunctionEvaluate
     procedure :: functionChanged => likelihoodHaloMassFunctionFunctionChanged
  end type likelihoodHaloMassFunction
  
  interface likelihoodHaloMassFunction
     !% Constructors for the halo mass function likelihood object.
     module procedure likelihoodHaloMassFunctionConstructor
  end interface likelihoodHaloMassFunction

  ! Mass function error model enumeration.
  integer, parameter :: haloMassFunctionErrorModelNone                =0
  integer, parameter :: haloMassFunctionErrorModelPowerLaw            =1
  integer, parameter :: haloMassFunctionErrorModelSphericalOverdensity=2
  
