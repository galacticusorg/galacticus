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

  !% Defines a likelihood function for fitting of SEDs to broad-band photometry.

  type, public, extends(likelihood) :: likelihoodSEDFit
     private
     integer                                     :: photometryCount         , dustType
     integer         , allocatable, dimension(:) :: filterIndex             , luminosityIndex         , &
          &                                         postprocessingChainIndex
     double precision, allocatable, dimension(:) :: magnitude               , error                   , &
          &                                         redshift                , massToLightRatio        , &
          &                                         age                     , wavelengthEffective
   contains
     procedure :: evaluate        => likelihoodSEDFitEvaluate
     procedure :: functionChanged => likelihoodSEDFitFunctionChanged
  end type likelihoodSEDFit
  
  interface likelihoodSEDFit
     !% Constructors for the SED fit likelihood object.
     module procedure likelihoodSEDFitConstructor
  end interface likelihoodSEDFit

  integer, parameter :: likelihoodSEDFitDustTypeNull           =0
  integer, parameter :: likelihoodSEDFitDustTypeCharlotFall2000=1
