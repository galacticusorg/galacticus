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

  !% Defines a convergence class that GelmanRubin converges (for testing purposes).

  type, public, extends(convergence) :: convergenceGelmanRubin
     !% Defines a convergence class using the Gelman-Rubin statistic.
     private
     double precision                            :: thresholdHatR             , outlierSignificance , &
          &                                         outlierLogLikelihoodOffset
     integer                                     :: burnCount                 , testCount           , &
          &                                         stepCount                 , outlierCountMaximum , &
          &                                         reportCount               , estimateCount       , &
          &                                         logFileUnit               , convergedAtStepCount
     logical                                     :: converged
     double precision, allocatable, dimension(:) :: correctedHatR
     logical         , allocatable, dimension(:) :: chainMask
   contains
     final     ::                   convergenceGelmanRubinDestructor
     procedure :: isConverged     => convergenceGelmanRubinIsConverged
     procedure :: convergedAtStep => convergenceGelmanRubinConvergedAtStep
     procedure :: reset           => convergenceGelmanRubinReset
     procedure :: logReport       => convergenceGelmanRubinLogReport
     procedure :: chainIsOutlier  => convergenceGelmanRubinChainIsOutlier
 end type convergenceGelmanRubin
  
  interface convergenceGelmanRubin
     !% Constructors for the ``GelmanRubin'' convergence class.
     module procedure convergenceGelmanRubinConstructor
  end interface convergenceGelmanRubin
