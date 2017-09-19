!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
!!    Andrew Benson <abenson@carnegiescience.edu>
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

  !% Defines a likelihood function for fitting of halo spin distributions.

  type, public, extends(likelihood) :: likelihoodHaloSpinDistribution
     private
     double precision        , dimension(:  ), allocatable :: spin                , distribution                      , &
          &                                                   spinMinimum         , spinMaximum                       , &
          &                                                   distributionError
     double precision                                      :: time                , massParticle                      , &
          &                                                   massHaloMinimum     , energyEstimateParticleCountMaximum
     integer                                               :: particleCountMinimum, distributionType
   contains
     procedure :: evaluate        => likelihoodHaloSpinDistributionEvaluate
     procedure :: functionChanged => likelihoodHaloSpinDistributionFunctionChanged
  end type likelihoodHaloSpinDistribution
  
  interface likelihoodHaloSpinDistribution
     !% Constructors for the halo spin distribution likelihood object.
     module procedure likelihoodHaloSpinDistributionConstructor
  end interface likelihoodHaloSpinDistribution
  
  !# <enumeration>
  !#  <name>distributionType</name>
  !#  <description>Used to specify the type of intrinsic spin distribution.</description>
  !#  <entry label="logNormal" />
  !#  <entry label="bett2007"  />
  !# </enumeration>
