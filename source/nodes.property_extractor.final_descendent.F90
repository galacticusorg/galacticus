!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a module which implements an ISM mass output analysis property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorFinalDescendent">
   <description>An ISM mass output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorFinalDescendent
     !!{
     A stelalr mass output analysis class.
     !!}
     private
   contains
     procedure :: extract     => finalDescendentExtract
     procedure :: name        => finalDescendentName
     procedure :: description => finalDescendentDescription
  end type nodePropertyExtractorFinalDescendent

  interface nodePropertyExtractorFinalDescendent
     !!{
     Constructors for the ``finalDescendent'' output analysis class.
     !!}
     module procedure finalDescendentConstructorParameters
  end interface nodePropertyExtractorFinalDescendent

contains

  function finalDescendentConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily finalDescendent} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorFinalDescendent)                :: self
    type(inputParameters                     ), intent(inout) :: parameters

    self=nodePropertyExtractorFinalDescendent()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function finalDescendentConstructorParameters

  function finalDescendentExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily finalDescendent} node property extractor.
    !!}
    implicit none
    integer         (kind_int8                           )                          :: finalDescendentExtract
    class           (nodePropertyExtractorFinalDescendent), intent(inout)           :: self
    type            (treeNode                            ), intent(inout), target   :: node
    double precision                                      , intent(in   )           :: time
    type            (multiCounter                        ), intent(inout), optional :: instance
    type            (treeNode                            ), pointer                 :: nodeDescendent
    !$GLC attributes unused :: self, instance, time

    nodeDescendent => node
    do while (associated(nodeDescendent%parent))
       if (nodeDescendent%isSatellite()) then
          nodeDescendent => nodeDescendent%mergesWith()
       else
          nodeDescendent => nodeDescendent%parent
       end if
    end do
    finalDescendentExtract=nodeDescendent%index()
    return
  end function finalDescendentExtract


  function finalDescendentName(self)
    !!{
    Return the name of the finalDescendent property.
    !!}
    implicit none
    type (varying_string                      )                :: finalDescendentName
    class(nodePropertyExtractorFinalDescendent), intent(inout) :: self
    !$GLC attributes unused :: self

    finalDescendentName=var_str('finalDescendentIndex')
    return
  end function finalDescendentName

  function finalDescendentDescription(self)
    !!{
    Return a description of the finalDescendent property.
    !!}
    implicit none
    type (varying_string                      )                :: finalDescendentDescription
    class(nodePropertyExtractorFinalDescendent), intent(inout) :: self
    !$GLC attributes unused :: self

    finalDescendentDescription=var_str('ID of the node which this node will have descended into at the base of the tree.')
    return
  end function finalDescendentDescription



