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
  Implements a node property extractor which reports if a node is considered to be physically plausible.
  !!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorIsPhysicallyPlausible">
   <description>A node property extractor class which reports if a node is considered to be physically plausible.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorIsPhysicallyPlausible
     !!{
     A node property extractor class which reports if a node is considered to be physically plausible.
      !!}
     private
   contains
     procedure :: extract     => isPhysicallyPlausibleExtract
     procedure :: name        => isPhysicallyPlausibleName
     procedure :: description => isPhysicallyPlausibleDescription
  end type nodePropertyExtractorIsPhysicallyPlausible

  interface nodePropertyExtractorIsPhysicallyPlausible
     !!{
     Constructors for the \refClass{nodePropertyExtractorIsPhysicallyPlausible} output analysis class.
     !!}
     module procedure isPhysicallyPlausibleConstructorParameters
  end interface nodePropertyExtractorIsPhysicallyPlausible

contains

  function isPhysicallyPlausibleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorIsPhysicallyPlausible} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorIsPhysicallyPlausible)                :: self
    type(inputParameters                           ), intent(inout) :: parameters

    self=nodePropertyExtractorIsPhysicallyPlausible()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function isPhysicallyPlausibleConstructorParameters

  function isPhysicallyPlausibleExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily isPhysicallyPlausible} node property extractor.
    !!}
    implicit none
    integer         (kind_int8                                 )                          :: isPhysicallyPlausibleExtract
    class           (nodePropertyExtractorIsPhysicallyPlausible), intent(inout)           :: self
    type            (treeNode                                  ), intent(inout), target   :: node
    double precision                                            , intent(in   )           :: time
    type            (multiCounter                              ), intent(inout), optional :: instance
    !$GLC attributes unused :: self, instance, time

    if (node%isPhysicallyPlausible) then
       isPhysicallyPlausibleExtract=1
    else
       isPhysicallyPlausibleExtract=0
    end if
    return
  end function isPhysicallyPlausibleExtract

  function isPhysicallyPlausibleName(self)
    !!{
    Return the name of the isPhysicallyPlausible property.
    !!}
    implicit none
    type (varying_string                            )                :: isPhysicallyPlausibleName
    class(nodePropertyExtractorIsPhysicallyPlausible), intent(inout) :: self
    !$GLC attributes unused :: self

    isPhysicallyPlausibleName=var_str('nodeIsPhysicallyPlausible')
    return
  end function isPhysicallyPlausibleName

  function isPhysicallyPlausibleDescription(self)
    !!{
    Return a description of the isPhysicallyPlausible property.
    !!}
    implicit none
    type (varying_string                            )                :: isPhysicallyPlausibleDescription
    class(nodePropertyExtractorIsPhysicallyPlausible), intent(inout) :: self
    !$GLC attributes unused :: self

    isPhysicallyPlausibleDescription=var_str('Indicates if the node considered to be physically plausible (0|1).')
    return
  end function isPhysicallyPlausibleDescription
