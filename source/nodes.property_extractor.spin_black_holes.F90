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

  !![
  <nodePropertyExtractor name="nodePropertyExtractorSpinBlackHoles">
   <description>
     A node property extractor which extracts a list of all super-massive black hole spins.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorSpinBlackHoles
     !!{
     A property extractor which extracts a list of all super-massive black hole spins.
     !!}
     private
   contains
     procedure :: elementCount => spinBlackHolesElementCount
     procedure :: extract      => spinBlackHolesExtract
     procedure :: names        => spinBlackHolesNames
     procedure :: descriptions => spinBlackHolesDescriptions
     procedure :: unitsInSI    => spinBlackHolesUnitsInSI
  end type nodePropertyExtractorSpinBlackHoles

  interface nodePropertyExtractorSpinBlackHoles
     !!{
     Constructors for the \refClass{nodePropertyExtractorSpinBlackHoles} output extractor class.
     !!}
     module procedure spinBlackHolesConstructorParameters
  end interface nodePropertyExtractorSpinBlackHoles

contains

  function spinBlackHolesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorSpinBlackHoles} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorSpinBlackHoles)                :: self
    type(inputParameters                    ), intent(inout) :: parameters

    self=nodePropertyExtractorSpinBlackHoles()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function spinBlackHolesConstructorParameters

  integer function spinBlackHolesElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorSpinBlackHoles), intent(inout) :: self

    spinBlackHolesElementCount=1
    return
  end function spinBlackHolesElementCount

  function spinBlackHolesExtract(self,node,instance) result(spin)
    !!{
    Implement an output extractor for the spins of all supermassive black holes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    double precision                                     , dimension(:,:), allocatable :: spin
    class           (nodePropertyExtractorSpinBlackHoles), intent(inout)               :: self
    type            (treeNode                           ), intent(inout)               :: node
    type            (multiCounter                       ), intent(inout) , optional    :: instance
    class           (nodeComponentBlackHole             )                , pointer     :: blackHole
    integer                                                                            :: i        , countBlackHoles
    !$GLC attributes unused :: instance

    countBlackHoles=node%blackHoleCount()
    allocate(spin(countBlackHoles,1))
    do i=1,countBlackHoles
       blackHole      => node     %blackHole(instance=i)
       spin     (i,1) =  blackHole%spin     (          )
    end do
    return
  end function spinBlackHolesExtract

  subroutine spinBlackHolesNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily spinBlackHoles} properties.
    !!}
    implicit none
    class(nodePropertyExtractorSpinBlackHoles), intent(inout)                             :: self
    type (varying_string                     ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self

    allocate(names(1))
    names(1)=var_str('spinBlackHoles')
    return
  end subroutine spinBlackHolesNames

  subroutine spinBlackHolesDescriptions(self,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily spinBlackHoles} properties.
    !!}
    implicit none
    class(nodePropertyExtractorSpinBlackHoles), intent(inout)                             :: self
    type (varying_string                     ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(1))
    descriptions(1)=var_str('Spins of super-massive black holes in this galaxy.')
    return
  end subroutine spinBlackHolesDescriptions

  function spinBlackHolesUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily spinBlackHoles} properties in the SI system.
    !!}
    implicit none
    double precision                                     , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorSpinBlackHoles), intent(inout)              :: self
    !$GLC attributes unused :: self

    allocate(unitsInSI(1))
    unitsInSI(1)=1.0d0
    return
  end function spinBlackHolesUnitsInSI
