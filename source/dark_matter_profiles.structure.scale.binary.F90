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
  An implementation of dark matter halo profile scale radii which switches between two methods based on a filter.
  !!}

  use :: Galactic_Filters, only : galacticFilter, galacticFilterClass

  !![
  <darkMatterProfileScaleRadius name="darkMatterProfileScaleRadiusBinary">
   <description>A dark matter halo profile scale radii class which switches between two methods based on a filter.</description>
  </darkMatterProfileScaleRadius>
  !!]
  type, extends(darkMatterProfileScaleRadiusClass) :: darkMatterProfileScaleRadiusBinary
     !!{
     A dark matter halo profile scale radii class which switches between two methods based on a filter.
     !!}
     private
     class(darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadiusAccept_ => null(), darkMatterProfileScaleRadiusReject_ => null()
     class(galacticFilterClass              ), pointer :: galacticFilter_                     => null()
   contains
     final     ::           binaryDestructor
     procedure :: radius => binaryRadius
  end type darkMatterProfileScaleRadiusBinary

  interface darkMatterProfileScaleRadiusBinary
     !!{
     Constructors for the \refClass{darkMatterProfileScaleRadiusBinary} dark matter halo profile concentration class.
     !!}
     module procedure binaryConstructorParameters
     module procedure binaryConstructorInternal
  end interface darkMatterProfileScaleRadiusBinary

contains

  function binaryConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily binary} dark matter halo profile concentration class.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileScaleRadiusBinary)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(darkMatterProfileScaleRadiusClass ), pointer       :: darkMatterProfileScaleRadiusAccept_, darkMatterProfileScaleRadiusReject_
     class(galacticFilterClass              ), pointer       :: galacticFilter_

    if (.not.parameters%isPresent('darkMatterProfileScaleRadiusAccept')) call Error_Report('a scale radius class for accepted nodes must be given'//{introspection:location})
    if (.not.parameters%isPresent('darkMatterProfileScaleRadiusReject')) call Error_Report('a scale radius class for rejected nodes must be given'//{introspection:location})
    !![
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadiusAccept_" parameterName="darkMatterProfileScaleRadiusAccept" source="parameters"/>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadiusReject_" parameterName="darkMatterProfileScaleRadiusReject" source="parameters"/>
    <objectBuilder class="galacticFilter"               name="galacticFilter_"                                                                        source="parameters"/>
    !!]
    self=darkMatterProfileScaleRadiusBinary(darkMatterProfileScaleRadiusAccept_,darkMatterProfileScaleRadiusReject_,galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileScaleRadiusAccept_"/>
    <objectDestructor name="darkMatterProfileScaleRadiusReject_"/>
    <objectDestructor name="galacticFilter_"                    />
    !!]
    return
  end function binaryConstructorParameters

  function binaryConstructorInternal(darkMatterProfileScaleRadiusAccept_,darkMatterProfileScaleRadiusReject_,galacticFilter_) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileScaleRadiusBinary} dark matter halo profile concentration class.
    !!}
    implicit none
    type (darkMatterProfileScaleRadiusBinary)                        :: self
    class(darkMatterProfileScaleRadiusClass ), intent(in   ), target :: darkMatterProfileScaleRadiusAccept_, darkMatterProfileScaleRadiusReject_
    class(galacticFilterClass               ), intent(in   ), target :: galacticFilter_
    !![
    <constructorAssign variables="*darkMatterProfileScaleRadiusAccept_, *darkMatterProfileScaleRadiusReject_, *galacticFilter_"/>
    !!]

    return
  end function binaryConstructorInternal

  subroutine binaryDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileScaleRadiusBinary} dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileScaleRadiusBinary), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileScaleRadiusAccept_"/>
    <objectDestructor name="self%darkMatterProfileScaleRadiusReject_"/>
    <objectDestructor name="self%galacticFilter_"                    />
    !!]
    return
  end subroutine binaryDestructor

  double precision function binaryRadius(self,node)
    !!{
    Return the scale radius of the dark matter halo profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class (darkMatterProfileScaleRadiusBinary), intent(inout), target :: self
    type  (treeNode                          ), intent(inout), target :: node

    if (self%galacticFilter_%passes(node)) then
       binaryRadius=self%darkMatterProfileScaleRadiusAccept_%radius(node)
    else
       binaryRadius=self%darkMatterProfileScaleRadiusReject_%radius(node)
    end if
    return
  end function binaryRadius
