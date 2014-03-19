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

  !% Defines a Gaussian regression likelihood function.

  type, public, extends(likelihood) :: likelihoodGaussianRegression
     private
     class           (likelihood      ), pointer                     :: simulatorLikelihoodFunction
     integer                                                         :: accumulatedStateCount      , emulatorRebuildCount       , &
          &                                                             polynomialOrder            , polynomialCoefficientCount , &
          &                                                             reportCount                , simulationCount            , &
          &                                                             evaluationCount            , emulatorCheckCount         , &
          &                                                             emulatorFailCount
     integer         (kind=FGSL_Size_T)                              :: regressionMatrixSize
     double precision                  , allocatable, dimension(:  ) :: simulatorLikelihood        , polynomialCoefficient      , &
          &                                                             likelihoodSums             , coefficients               , &
          &                                                             stateOffset                , weight                     , &
          &                                                             likelihoodResiduals        , stateScales
     double precision                  , allocatable, dimension(:,:) :: simulationState            , stateSums                  , &
          &                                                             regressionMatrix           , statesCombined
     logical                                                         :: initialized                , regressionMatrixIsSingular
     type            (fgsl_matrix      )                             :: regressionMatrixLU
     type            (fgsl_permutation )                             :: permutations
     double precision                                                :: C0                         , C1                         , &
          &                                                             CR                         , sigmaBuffer                , &
          &                                                             logLikelihoodBuffer        , logLikelihoodErrorTolerance
     logical                                                         :: emulateOutliers
   contains
     final     ::                    likelihoodGaussianRegressionDestructor
     procedure :: evaluate        => likelihoodGaussianRegressionEvaluate
     procedure :: functionChanged => likelihoodGaussianRegressionFunctionChanged
  end type likelihoodGaussianRegression
  
  interface likelihoodGaussianRegression
     !% Constructors for the Gaussian regression likelihood object.
     module procedure likelihoodGaussianRegressionConstructor
  end interface likelihoodGaussianRegression

  type polynomialIterator
     !% An object used for iterating over coefficients of polynomials.
     private
     integer                                    :: order  , rank, orderCurrent, stateCurrent, count
     integer        , allocatable, dimension(:) :: indices
   contains
     !@ <objectMethods>
     !@   <object>polynomialIterator</object>
     !@   <objectMethod>
     !@     <method>reset</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Reset the iterator object to the start of its sequence.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>iterate</method>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@     <description>Move to the next iteration of polynomial coefficient indices. Returns true if successful. If no more iterations are available, returns false.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>index</method>
     !@     <type>\intone</type>
     !@     <arguments>\intone\ i\argin</arguments>
     !@     <description>Return the $i^{\rm th}$ index of the polynomial coefficient.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>currentOrder</method>
     !@     <type>\intone</type>
     !@     <arguments></arguments>
     !@     <description>Return the current order of the polynomial coefficient.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>counter</method>
     !@     <type>\intone</type>
     !@     <arguments></arguments>
     !@     <description>Return an incremental counter (i.e. begins at $0$ and increases by $1$ on each iteration).</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: index        => polynomialIteratorIndex
     procedure :: currentOrder => polynomialIteratorCurrentOrder
     procedure :: counter      => polynomialIteratorCounter
     procedure :: iterate      => polynomialIteratorIterate
     procedure :: reset        => polynomialIteratorReset
  end type polynomialIterator
  
  interface polynomialIterator
     module procedure polynomialIteratorConstructor
  end interface polynomialIterator

  ! Workspaces used when fitting the semi-variogram.
  double precision, allocatable, dimension(:) :: separationsNormalized  , semiVariancesNormalized, separationsLimited
  double precision                            :: separationNormalization, semiVarianceNormalization
  
