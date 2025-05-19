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

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassBlackHoles">
   <description>
     A node property extractor which extracts a list of all super-massive black hole masses.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorMassBlackHoles
     !!{
     A property extractor which extracts a list of all super-massive black hole masses.
     !!}
     private
   contains
     procedure :: elementCount => massBlackHolesElementCount
     procedure :: extract      => massBlackHolesExtract
     procedure :: names        => massBlackHolesNames
     procedure :: descriptions => massBlackHolesDescriptions
     procedure :: unitsInSI    => massBlackHolesUnitsInSI
  end type nodePropertyExtractorMassBlackHoles

  interface nodePropertyExtractorMassBlackHoles
     !!{
     Constructors for the {\normalfont \ttfamily massBlackHoles} output extractor class.
     !!}
     module procedure massBlackHolesConstructorParameters
  end interface nodePropertyExtractorMassBlackHoles

contains

  function massBlackHolesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily massBlackHoles} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorMassBlackHoles)                :: self
    type(inputParameters                    ), intent(inout) :: parameters

    self=nodePropertyExtractorMassBlackHoles()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massBlackHolesConstructorParameters

  integer function massBlackHolesElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorMassBlackHoles), intent(inout) :: self

    massBlackHolesElementCount=1
    return
  end function massBlackHolesElementCount

  function massBlackHolesExtract(self,node,instance) result(mass)
    !!{
    Implement an output extractor for the masses of all supermassive black holes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    double precision                                     , dimension(:,:), allocatable :: mass
    class           (nodePropertyExtractorMassBlackHoles), intent(inout)               :: self
    type            (treeNode                           ), intent(inout)               :: node
    type            (multiCounter                       ), intent(inout) , optional    :: instance
    class           (nodeComponentBlackHole             )                , pointer     :: blackHole
    integer                                                                            :: i        , countBlackHoles
    !$GLC attributes unused :: instance

    countBlackHoles=node%blackHoleCount()
    allocate(mass(countBlackHoles,1))
    do i=1,countBlackHoles
       blackHole      => node     %blackHole(instance=i)
       mass     (i,1) =  blackHole%mass     (          )
    end do
    return
  end function massBlackHolesExtract

  subroutine massBlackHolesNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily massBlackHoles} properties.
    !!}
    implicit none
    class(nodePropertyExtractorMassBlackHoles), intent(inout)                             :: self
    type (varying_string                     ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self

    allocate(names(1))
    names(1)=var_str('massBlackHoles')
    return
  end subroutine massBlackHolesNames

  subroutine massBlackHolesDescriptions(self,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily massBlackHoles} properties.
    !!}
    implicit none
    class(nodePropertyExtractorMassBlackHoles), intent(inout)                             :: self
    type (varying_string                     ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(1))
    descriptions(1)=var_str('Masses of super-massive black holes in this galaxy.')
    return
  end subroutine massBlackHolesDescriptions

  function massBlackHolesUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily massBlackHoles} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    double precision                                     , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorMassBlackHoles), intent(inout)              :: self
    !$GLC attributes unused :: self

    allocate(unitsInSI(1))
    unitsInSI(1)=massSolar
    return
  end function massBlackHolesUnitsInSI
