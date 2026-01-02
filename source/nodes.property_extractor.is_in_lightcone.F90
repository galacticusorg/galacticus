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
Implements a property extractor for ``in lightcone'' status.
!!}

  use :: Geometry_Lightcones, only : geometryLightconeClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorIsInLightcone">
   <description>An ``in lightcone'' status property extractor.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorIsInLightcone
     !!{
     An ``in lightcone'' status property extractor.
     !!}
     private
     class(geometryLightconeClass), pointer :: geometryLightcone_ => null()
   contains
     final     ::                isInLightconeDestructor
     procedure :: extract     => isInLightconeExtract
     procedure :: name        => isInLightconeName
     procedure :: description => isInLightconeDescription
  end type nodePropertyExtractorIsInLightcone

  interface nodePropertyExtractorIsInLightcone
     !!{
     Constructors for the \refClass{nodePropertyExtractorIsInLightcone} output analysis class.
     !!}
     module procedure isInLightconeConstructorParameters
     module procedure isInLightconeConstructorInternal
  end interface nodePropertyExtractorIsInLightcone

contains

  function isInLightconeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorIsInLightcone} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorIsInLightcone)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(geometryLightconeClass            ), pointer       :: geometryLightcone_

    !![
    <objectBuilder class="geometryLightcone" name="geometryLightcone_" source="parameters"/>
    !!]
    self=nodePropertyExtractorIsInLightcone(geometryLightcone_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="geometryLightcone_"/>
    !!]
    return
  end function isInLightconeConstructorParameters

  function isInLightconeConstructorInternal(geometryLightcone_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorIsInLightcone} node property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorIsInLightcone)                        :: self
    class(geometryLightconeClass            ), intent(in   ), target :: geometryLightcone_
    !![
    <constructorAssign variables="*geometryLightcone_"/>
    !!]
    
    return
  end function isInLightconeConstructorInternal

  subroutine isInLightconeDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorIsInLightcone} node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorIsInLightcone), intent(inout) :: self

    !![
    <objectDestructor name="self%geometryLightcone_"/>
    !!]
    return
  end subroutine isInLightconeDestructor

  function isInLightconeExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily isInLightcone} node property extractor.
    !!}
    implicit none
    integer         (kind_int8                         )                          :: isInLightconeExtract
    class           (nodePropertyExtractorIsInLightcone), intent(inout)           :: self
    type            (treeNode                          ), intent(inout), target   :: node
    double precision                                    , intent(in   )           :: time
    type            (multiCounter                      ), intent(inout), optional :: instance
    !$GLC attributes unused :: instance, time

    if (self%geometryLightcone_%isInLightcone(node,atPresentEpoch=.true.)) then
       isInLightconeExtract=1_c_size_t
    else
       isInLightconeExtract=0_c_size_t
    end if
    return
  end function isInLightconeExtract


  function isInLightconeName(self)
    !!{
    Return the name of the ``in lightcone'' property.
    !!}
    implicit none
    type (varying_string                    )                :: isInLightconeName
    class(nodePropertyExtractorIsInLightcone), intent(inout) :: self
    !$GLC attributes unused :: self

    isInLightconeName=var_str('isInLightcone')
    return
  end function isInLightconeName

  function isInLightconeDescription(self)
    !!{
    Return a description of the ``in lightcone'' property.
    !!}
    implicit none
    type (varying_string                     )                :: isInLightconeDescription
    class(nodePropertyExtractorIsInLightcone), intent(inout) :: self
    !$GLC attributes unused :: self

    isInLightconeDescription=var_str('Node is within the lightcone (1 if true, 0 otherwise).')
    return
  end function isInLightconeDescription


