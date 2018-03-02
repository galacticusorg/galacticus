!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% Defines a distribution for cosmological density fields with normally distributed amplitude, but conditioned to have not
  !% collapse on a larger scale \citep[][eqn.~9]{mo_analytic_1996}.

  type, public, extends(distribution1D) :: distributionPeakBackground
     double precision                                   :: varianceBackground, thresholdCollapse, &
          &                                                normalization
     type            (table1DLinearLinear)              :: cdf
     class           (table1D            ), allocatable :: cdfInverse   
   contains
     procedure :: density    => distributionPeakBackgroundDensity
     procedure :: cumulative => distributionPeakBackgroundCumulative
     procedure :: inverse    => distributionPeakBackgroundInverse
     procedure :: minimum    => distributionPeakBackgroundMinimum
     procedure :: maximum    => distributionPeakBackgroundMaximum
  end type distributionPeakBackground

  interface distributionPeakBackground
     !% Constructors for the {\normalfont \ttfamily peakBackground} distribution object.
     module procedure distributionPeakBackgroundConstructor
  end interface distributionPeakBackground
