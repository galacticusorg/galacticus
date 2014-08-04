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

  !% Defines a likelihood function which implements ``posterior as prior''.

  type, public, extends(likelihood) :: likelihoodPosteriorPrior
     private
     class           (likelihood      ), pointer                     :: wrappedLikelihoodFunction
     integer                           , allocatable, dimension(:  ) :: neighborIndices
     double precision                  , allocatable, dimension(:  ) :: neighborDistances        , rootVariance
     double precision                  , allocatable, dimension(:,:) :: states
     class           (nearestNeighbors), allocatable                 :: searcher
     integer                                                         :: dimCount                 , convergedStateCount  , &
          &                                                             neighborCount
     double precision                                                :: tolerance                , logPriorNormalization
     logical                                                         :: initialized
   contains
     final     ::                    likelihoodPosteriorPriorDestructor
     procedure :: evaluate        => likelihoodPosteriorPriorEvaluate
     procedure :: functionChanged => likelihoodPosteriorPriorFunctionChanged
     procedure :: initialize      => likelihoodPosteriorPriorInitialize
  end type likelihoodPosteriorPrior
  
  interface likelihoodPosteriorPrior
     !% Constructors for the ``posterior as prior'' likelihood object.
     module procedure likelihoodPosteriorPriorConstructor
  end interface likelihoodPosteriorPrior
