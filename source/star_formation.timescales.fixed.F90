!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

  !!{
  Implementation of a fixed timescale for star formation.
  !!}

  !![
  <starFormationTimescale name="starFormationTimescaleFixed">
   <description>
    A star formation timescale class which fixed timescale for star formation {\normalfont \ttfamily [timescale]} (in Gyr).
   </description>
  </starFormationTimescale>
  !!]
  type, extends(starFormationTimescaleClass) :: starFormationTimescaleFixed
     !!{
     Implementation of a fixed timescale for star formation.
     !!}
     private
     double precision :: timescale_
   contains
     procedure :: timescale => fixedTimescale
  end type starFormationTimescaleFixed

  interface starFormationTimescaleFixed
     !!{
     Constructors for the \refClass{starFormationTimescaleFixed} timescale for star formation.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface starFormationTimescaleFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{starFormationTimescaleFixed} timescale for star formation class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationTimescaleFixed)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    double precision                                             :: timescale

    !![
    <inputParameter>
      <name>timescale</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The timescale for star formation in the fixed timescale model.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=starFormationTimescaleFixed(timescale)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(timescale) result(self)
    !!{
    Internal constructor for the \refClass{starFormationTimescaleFixed} timescale for star formation class.
    !!}
    implicit none
    type            (starFormationTimescaleFixed)                :: self
    double precision                             , intent(in   ) :: timescale

    self%timescale_=timescale
    return
  end function fixedConstructorInternal

  double precision function fixedTimescale(self,component)
    !!{
    Returns the timescale (in Gyr) for star formation in the given {\normalfont \ttfamily component}, assuming a fixed
    timescale.
    !!}
    implicit none
    class(starFormationTimescaleFixed), intent(inout) :: self
    class(nodeComponent              ), intent(inout) :: component
    !$GLC attributes unused :: component

    fixedTimescale=self%timescale_
    return
  end function fixedTimescale
