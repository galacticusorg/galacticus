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
  <nodePropertyExtractor name="nodePropertyExtractorMassAccretionHistory">
   <description>
     A node property extractor which extracts the mass accretion history for each node.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorMassAccretionHistory
     !!{
     A property extractor which extracts the mass accretion history for each node.
     !!}
     private
     integer :: massAccretionHistoryTimeID, massAccretionHistoryMassID
   contains
     procedure :: elementCount => massAccretionHistoryElementCount
     procedure :: extract      => massAccretionHistoryExtract
     procedure :: names        => massAccretionHistoryNames
     procedure :: descriptions => massAccretionHistoryDescriptions
     procedure :: unitsInSI    => massAccretionHistoryUnitsInSI
  end type nodePropertyExtractorMassAccretionHistory

  interface nodePropertyExtractorMassAccretionHistory
     !!{
     Constructors for the {\normalfont \ttfamily massAccretionHistory} output extractor class.
     !!}
     module procedure massAccretionHistoryConstructorParameters
     module procedure massAccretionHistoryConstructorInternal
  end interface nodePropertyExtractorMassAccretionHistory

contains

  function massAccretionHistoryConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily massAccretionHistory} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorMassAccretionHistory)                :: self
    type(inputParameters                          ), intent(inout) :: parameters

    self=nodePropertyExtractorMassAccretionHistory()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massAccretionHistoryConstructorParameters

  function massAccretionHistoryConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily massAccretionHistory} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorMassAccretionHistory) :: self
    
    !![
    <addMetaProperty component="basic" name="massAccretionHistoryTime" id="self%massAccretionHistoryTimeID" rank="1" isCreator="no"/>
    <addMetaProperty component="basic" name="massAccretionHistoryMass" id="self%massAccretionHistoryMassID" rank="1" isCreator="no"/>
    !!]
    return
  end function massAccretionHistoryConstructorInternal

  integer function massAccretionHistoryElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorMassAccretionHistory), intent(inout) :: self

    massAccretionHistoryElementCount=2
    return
  end function massAccretionHistoryElementCount

  function massAccretionHistoryExtract(self,node,instance) result(massAccretionHistory)
    !!{
    Implement a massAccretionHistory output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    double precision                                           , dimension(:,:), allocatable :: massAccretionHistory
    class           (nodePropertyExtractorMassAccretionHistory), intent(inout)               :: self
    type            (treeNode                                 ), intent(inout)               :: node
    type            (multiCounter                             ), intent(inout) , optional    :: instance
    class           (nodeComponentBasic                       )                , pointer     :: basic
    double precision                                           , dimension(:  ), allocatable :: times           , masses
    !$GLC attributes unused :: instance
    !$GLC attributes initialized :: times, masses
    
    basic  => node %basic                    (                               )
    times  =  basic%floatRank1MetaPropertyGet(self%massAccretionHistoryTimeID)
    masses =  basic%floatRank1MetaPropertyGet(self%massAccretionHistoryMassID)
    allocate(massAccretionHistory(size(times),2))
    massAccretionHistory(:,1)=times
    massAccretionHistory(:,2)=masses
    return
  end function massAccretionHistoryExtract
  
  subroutine massAccretionHistoryNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily massAccretionHistory} properties.
    !!}
    implicit none
    class(nodePropertyExtractorMassAccretionHistory), intent(inout)                             :: self
    type (varying_string                           ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self

    allocate(names(2))
    names(1)=var_str('haloAccretionHistoryTime')
    names(2)=var_str('haloAccretionHistoryMass')
    return
  end subroutine massAccretionHistoryNames

  subroutine massAccretionHistoryDescriptions(self,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily massAccretionHistory} properties.
    !!}
    implicit none
    class(nodePropertyExtractorMassAccretionHistory), intent(inout)                             :: self
    type (varying_string                           ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(2))
    descriptions(1)=var_str('The time at which the DMO mass of the halo is tabulated.')
    descriptions(2)=var_str('The DMO mass of the halo as a function of time.'         )
    return
  end subroutine massAccretionHistoryDescriptions

  function massAccretionHistoryUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily massAccretionHistory} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, gigaYear
    implicit none
    double precision                                           , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorMassAccretionHistory), intent(inout)              :: self
    !$GLC attributes unused :: self

    allocate(unitsInSI(2))
    unitsInSI(1)=gigaYear
    unitsInSI(2)=massSolar    
    return
  end function massAccretionHistoryUnitsInSI
