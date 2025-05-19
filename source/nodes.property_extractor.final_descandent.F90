!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Implements an ISM mass output analysis property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorFinalDescendant">
   <description>An ISM mass output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorFinalDescendant
     !!{
     A stellar mass output analysis class.
     !!}
     private
   contains
     procedure :: extract     => finalDescendantExtract
     procedure :: name        => finalDescendantName
     procedure :: description => finalDescendantDescription
  end type nodePropertyExtractorFinalDescendant

  interface nodePropertyExtractorFinalDescendant
     !!{
     Constructors for the {\normalfont \ttfamily finalDescendant} output analysis class.
     !!}
     module procedure finalDescendantConstructorParameters
  end interface nodePropertyExtractorFinalDescendant

contains

  function finalDescendantConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily finalDescendant} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorFinalDescendant)                :: self
    type(inputParameters                     ), intent(inout) :: parameters

    self=nodePropertyExtractorFinalDescendant()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function finalDescendantConstructorParameters

  function finalDescendantExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily finalDescendant} node property extractor.
    !!}
    implicit none
    integer         (kind_int8                           )                          :: finalDescendantExtract
    class           (nodePropertyExtractorFinalDescendant), intent(inout)           :: self
    type            (treeNode                            ), intent(inout), target   :: node
    double precision                                      , intent(in   )           :: time
    type            (multiCounter                        ), intent(inout), optional :: instance
    type            (treeNode                            ), pointer                 :: nodeDescendant
    !$GLC attributes unused :: self, instance, time

    nodeDescendant => node
    do while (associated(nodeDescendant%parent))
       if (nodeDescendant%isSatellite()) then
          nodeDescendant => nodeDescendant%mergesWith()
       else
          nodeDescendant => nodeDescendant%parent
       end if
    end do
    finalDescendantExtract=nodeDescendant%index()
    return
  end function finalDescendantExtract


  function finalDescendantName(self)
    !!{
    Return the name of the finalDescendant property.
    !!}
    implicit none
    type (varying_string                      )                :: finalDescendantName
    class(nodePropertyExtractorFinalDescendant), intent(inout) :: self
    !$GLC attributes unused :: self

    finalDescendantName=var_str('finalDescendantIndex')
    return
  end function finalDescendantName

  function finalDescendantDescription(self)
    !!{
    Return a description of the finalDescendant property.
    !!}
    implicit none
    type (varying_string                      )                :: finalDescendantDescription
    class(nodePropertyExtractorFinalDescendant), intent(inout) :: self
    !$GLC attributes unused :: self

    finalDescendantDescription=var_str('ID of the node which this node will have descended into at the base of the tree.')
    return
  end function finalDescendantDescription



