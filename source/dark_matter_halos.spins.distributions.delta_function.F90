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

  !% An implementation of the dark matter halo spin distribution which assumes a
  !% $\delta$-function distribution.

  !# <haloSpinDistribution name="haloSpinDistributionDeltaFunction">
  !#  <description>A $\delta$-function halo spin distribution.</description>
  !# </haloSpinDistribution>
  type, extends(haloSpinDistributionClass) :: haloSpinDistributionDeltaFunction
     !% A dark matter halo spin distribution concentration class which assumes a
     !% $\delta$-function distribution.
     private
     double precision :: spin
   contains
     final     ::                 deltaFunctionDestructor
     procedure :: sample       => deltaFunctionSample
     procedure :: distribution => deltaFunctionDistribution
  end type haloSpinDistributionDeltaFunction

  interface haloSpinDistributionDeltaFunction
     !% Constructors for the {\normalfont \ttfamily deltaFunction} dark matter halo spin
     !% distribution class.
     module procedure deltaFunctionConstructorParameters
     module procedure deltaFunctionConstructorInternal
  end interface haloSpinDistributionDeltaFunction

contains

  function deltaFunctionConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily deltaFunction} dark matter halo spin
    !% distribution class which takes a parameter list as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(haloSpinDistributionDeltaFunction)                :: deltaFunctionConstructorParameters
    type(inputParameters                  ), intent(inout) :: parameters

    ! Check and read parameters.
    !# <inputParameter>
    !#   <name>spin</name>
    !#   <source>parameters</source>
    !#   <variable>deltaFunctionConstructorParameters%spin</variable>
    !#   <defaultValue>0.03687d0</defaultValue>
    !#   <defaultSource>\citep{bett_spin_2007}</defaultSource>
    !#   <description>The fixed value of spin in a $\delta$-function spin distribution.</description>
    !# </inputParameter>
    !# <inputParametersValidate source="parameters"/>
    return
  end function deltaFunctionConstructorParameters

  function deltaFunctionConstructorInternal(spin)
    !% Internal constructor for the {\normalfont \ttfamily deltaFunction} dark matter halo spin
    !% distribution class.
    implicit none
    type(haloSpinDistributionDeltaFunction) :: deltaFunctionConstructorInternal
    double precision, intent(in   ) :: spin

    deltaFunctionConstructorInternal%spin=spin
    return
  end function deltaFunctionConstructorInternal

  subroutine deltaFunctionDestructor(self)
    !% Destructor for the {\normalfont \ttfamily deltaFunction} dark matter halo spin
    !% distribution class.
    implicit none
    type(haloSpinDistributionDeltaFunction), intent(inout) :: self
    !$GLC attributes unused :: self

    ! Nothing to do.
    return
  end subroutine deltaFunctionDestructor

  double precision function deltaFunctionSample(self,node)
    !% Sample from a $\delta$-function spin parameter distribution for the given {\normalfont
    !% \ttfamily node}.
    implicit none
    class(haloSpinDistributionDeltaFunction), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node
    !$GLC attributes unused :: node

    deltaFunctionSample=self%spin
    return
  end function deltaFunctionSample

  double precision function deltaFunctionDistribution(self,node)
    !% Return the spin parameter distribution for the given {\normalfont \ttfamily node}.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(haloSpinDistributionDeltaFunction), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    deltaFunctionDistribution=0.0d0
    call Galacticus_Error_Report('distribution function can not be evaluated'//{introspection:location})
    return
  end function deltaFunctionDistribution
