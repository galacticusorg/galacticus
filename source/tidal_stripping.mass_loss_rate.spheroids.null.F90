!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implementation of a null tidal stripping of spheroids class.

  !# <tidalStrippingSpheroids name="tidalStrippingSpheroidsNull">
  !#  <description>A null model of tidal stripping in galactic spheroids.</description>
  !# </tidalStrippingSpheroids>
  type, extends(tidalStrippingSpheroidsClass) :: tidalStrippingSpheroidsNull
     !% Implementation of a null model of tidal stripping of galactic spheroids.
     private
   contains
     procedure :: rateMassLoss => nullRateMassLoss
  end type tidalStrippingSpheroidsNull

  interface tidalStrippingSpheroidsNull
     !% Constructors for the {\normalfont \ttfamily null} model of tidal stripping of spheroids class.
     module procedure nullConstructorParameters
  end interface tidalStrippingSpheroidsNull

contains

  function nullConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily null} timescale for star formation feedback in spheroids class which takes a
    !% parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(tidalStrippingSpheroidsNull)                :: self
    type(inputParameters            ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    self=tidalStrippingSpheroidsNull()
    return
  end function nullConstructorParameters

  double precision function nullRateMassLoss(self,node)
    !% Returns a zero mass loss rate due to tidal stripping.
    implicit none
    class(tidalStrippingSpheroidsNull), intent(inout) :: self
    type (treeNode                   ), intent(inout) :: node
    !GCC$ attributes unused :: self, node

    nullRateMassLoss=0.0d0
    return
  end function nullRateMassLoss
