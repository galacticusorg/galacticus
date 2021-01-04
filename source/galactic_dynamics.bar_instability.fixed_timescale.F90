!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

  !% Implementation of a simple model for galactic disk bar instability in which the timescale is fixed.

  !# <galacticDynamicsBarInstability name="galacticDynamicsBarInstabilityFixedTimescale">
  !#  <description>A simple model for galactic disk bar instability in which the timescale is fixed.</description>
  !# </galacticDynamicsBarInstability>
  type, extends(galacticDynamicsBarInstabilityClass) :: galacticDynamicsBarInstabilityFixedTimescale
     !% Implementation of a simple model for galactic disk bar instability in which the timescale is fixed.
     private
     double precision :: timescale_
   contains
     procedure :: timescale => fixedTimescaleTimescale
  end type galacticDynamicsBarInstabilityFixedTimescale

  interface galacticDynamicsBarInstabilityFixedTimescale
     !% Constructors for the {\normalfont \ttfamily fixedTimescale} model for galactic disk bar instability class.
     module procedure fixedTimescaleConstructorParameters
     module procedure fixedTimescaleConstructorInternal
  end interface galacticDynamicsBarInstabilityFixedTimescale

contains

  function fixedTimescaleConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily fixedTimescale} model for galactic disk bar instability class which takes a
    !% parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticDynamicsBarInstabilityFixedTimescale)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    double precision                                                              :: timescale

    !# <inputParameter>
    !#   <name>timescale</name>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The timescale for bar instability.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    self=galacticDynamicsBarInstabilityFixedTimescale(timescale)
    !# <inputParametersValidate source="parameters"/>
    return
  end function fixedTimescaleConstructorParameters

  function fixedTimescaleConstructorInternal(timescale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily fixedTimescale} model for galactic disk bar instability class.
    implicit none
    type            (galacticDynamicsBarInstabilityFixedTimescale)                :: self
    double precision                                              , intent(in   ) :: timescale_
    !# <constructorAssign variables="timescale_"/>

    return
  end function fixedTimescaleConstructorInternal

  subroutine fixedTimescaleTimescale(self,node,timescale,externalDrivingSpecificTorque)
    !% Assume a constant timescale for depletion of a disk to a pseudo-bulge via bar instability.
    implicit none
    class           (galacticDynamicsBarInstabilityFixedTimescale), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    double precision                                              , intent(  out) :: externalDrivingSpecificTorque, timescale
    !$GLC attributes unused :: node

    timescale                    =self%timescale_
    externalDrivingSpecificTorque=0.0d0
    return
  end subroutine fixedTimescaleTimescale
