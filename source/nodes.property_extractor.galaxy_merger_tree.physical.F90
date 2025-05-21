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
  <nodePropertyExtractor name="nodePropertyExtractorGalaxyMergerTreePhysical">
   <description>
     A node property extractor which extracts the physical properties of galaxy merger trees.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorGalaxyMergerTreePhysical
     !!{
     A property extractor which extracts the physical properties of galaxy merger trees.
     !!}
     private
     integer :: timeID, propertyID
   contains
     procedure :: elementCount => galaxyMergerTreePhysicalElementCount
     procedure :: extract      => galaxyMergerTreePhysicalExtract
     procedure :: names        => galaxyMergerTreePhysicalNames
     procedure :: descriptions => galaxyMergerTreePhysicalDescriptions
     procedure :: unitsInSI    => galaxyMergerTreePhysicalUnitsInSI
  end type nodePropertyExtractorGalaxyMergerTreePhysical

  interface nodePropertyExtractorGalaxyMergerTreePhysical
     !!{
     Constructors for the \refClass{nodePropertyExtractorGalaxyMergerTreePhysical} output extractor class.
     !!}
     module procedure galaxyMergerTreePhysicalConstructorParameters
     module procedure galaxyMergerTreePhysicalConstructorInternal
  end interface nodePropertyExtractorGalaxyMergerTreePhysical

contains

  function galaxyMergerTreePhysicalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorGalaxyMergerTreePhysical} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorGalaxyMergerTreePhysical)                :: self
    type(inputParameters                              ), intent(inout) :: parameters

    self=nodePropertyExtractorGalaxyMergerTreePhysical()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function galaxyMergerTreePhysicalConstructorParameters

  function galaxyMergerTreePhysicalConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorGalaxyMergerTreePhysical} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorGalaxyMergerTreePhysical) :: self
    
    !![
    <addMetaProperty component="basic" name="galaxyMergerTreeProperty" id="self%propertyID" rank="1" isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="basic" name="galaxyMergerTreeTime"     id="self%timeID"     rank="1" isEvolvable="no" isCreator="no"/>
      !!]
    return
  end function galaxyMergerTreePhysicalConstructorInternal

  integer function galaxyMergerTreePhysicalElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    use :: Node_Property_Extractor_Galaxy_Merger_Trees, only : nodePropertyExtractorGalaxyMergerTreeCount
    implicit none
    class(nodePropertyExtractorGalaxyMergerTreePhysical), intent(inout) :: self

    galaxyMergerTreePhysicalElementCount=1+nodePropertyExtractorGalaxyMergerTreeCount
    return
  end function galaxyMergerTreePhysicalElementCount

  function galaxyMergerTreePhysicalExtract(self,node,instance) result(galaxyMergerTree)
    !!{
    Implement a galaxyMergerTreePhysical output extractor.
    !!}
    use :: Galacticus_Nodes                           , only : nodeComponentBasic
    use :: Node_Property_Extractor_Galaxy_Merger_Trees, only : nodePropertyExtractorGalaxyMergerTreeCount
    implicit none
    double precision                                               , dimension(:,:), allocatable :: galaxyMergerTree
    class           (nodePropertyExtractorGalaxyMergerTreePhysical), intent(inout)               :: self
    type            (treeNode                                     ), intent(inout)               :: node
    type            (multiCounter                                 ), intent(inout) , optional    :: instance
    class           (nodeComponentBasic                           )                , pointer     :: basic
    double precision                                               , dimension(:  ), allocatable :: times           , properties
    integer                                                                                      :: i               , j
    !$GLC attributes unused :: instance
    !$GLC attributes initialized :: times, properties
    
    basic      => node %basic                    (               )
    times      =  basic%floatRank1MetaPropertyGet(self%    timeID)
    properties =  basic%floatRank1MetaPropertyGet(self%propertyID)
    allocate(galaxyMergerTree(size(times),1+nodePropertyExtractorGalaxyMergerTreeCount))
    galaxyMergerTree(:,1)=times
    do i=1,nodePropertyExtractorGalaxyMergerTreeCount
       do j=1,size(times)
          galaxyMergerTree(j,1+i)=properties((j-1)*nodePropertyExtractorGalaxyMergerTreeCount+i)
       end do
    end do
    return
  end function galaxyMergerTreePhysicalExtract
  
  subroutine galaxyMergerTreePhysicalNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily galaxyMergerTreePhysical} properties.
    !!}
    use :: Error                                      , only : Error_Report
    use :: Node_Property_Extractor_Galaxy_Merger_Trees, only : nodePropertyExtractorGalaxyMergerTreeCount, nodePropertyExtractorGalaxyMergerTree_
    use :: String_Handling                            , only : String_Upper_Case_First
    implicit none
    class  (nodePropertyExtractorGalaxyMergerTreePhysical), intent(inout)                            :: self
    type   (varying_string                               ), intent(inout), dimension(:), allocatable :: names
    class  (*                                            )                             , pointer     :: extractor_
    integer                                                                                          :: i
    !$GLC attributes unused :: self

    allocate(names(1+nodePropertyExtractorGalaxyMergerTreeCount))
    names(1)=var_str('galaxyMergerTreeTime'    )
    do i=1,nodePropertyExtractorGalaxyMergerTreeCount
       extractor_ => nodePropertyExtractorGalaxyMergerTree_(i)%extractor_ 
       select type (extractor_)
       class is (nodePropertyExtractorScalar)
          names(i+1)=var_str('galaxyMergerTree')//String_Upper_Case_First(char(extractor_%name()))
       class default
          names(i+1)=var_str('')
          call Error_Report("unexpected class"//{introspection:location})
       end select
    end do
    return
  end subroutine galaxyMergerTreePhysicalNames

  subroutine galaxyMergerTreePhysicalDescriptions(self,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily galaxyMergerTreePhysical} properties.
    !!}
    use :: Error                                      , only : Error_Report
    use :: Node_Property_Extractor_Galaxy_Merger_Trees, only : nodePropertyExtractorGalaxyMergerTreeCount, nodePropertyExtractorGalaxyMergerTree_
    implicit none
    class  (nodePropertyExtractorGalaxyMergerTreePhysical), intent(inout)                             :: self
    type   (varying_string                               ), intent(inout), dimension(:) , allocatable :: descriptions
    integer                                                                                           :: i
    class  (*                                            )                              , pointer     :: extractor_
    !$GLC attributes unused :: self

    allocate(descriptions(1+nodePropertyExtractorGalaxyMergerTreeCount))
    descriptions(1)=var_str('Sample times in the galaxy merger tree.')
    do i=1,nodePropertyExtractorGalaxyMergerTreeCount
       extractor_ => nodePropertyExtractorGalaxyMergerTree_(i)%extractor_ 
       select type (extractor_)
       class is (nodePropertyExtractorScalar)
          descriptions(i+1)=var_str('Property in the galaxy merger tree: ')//extractor_%description()
       class default
          descriptions(i+1)=var_str('')
          call Error_Report("unexpected class"//{introspection:location})
       end select
    end do
    return
  end subroutine galaxyMergerTreePhysicalDescriptions

  function galaxyMergerTreePhysicalUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily galaxyMergerTreePhysical} properties in the SI system.
    !!}
    use :: Error                                      , only : Error_Report
    use :: Numerical_Constants_Astronomical           , only : gigaYear
    use :: Node_Property_Extractor_Galaxy_Merger_Trees, only : nodePropertyExtractorGalaxyMergerTreeCount, nodePropertyExtractorGalaxyMergerTree_
    implicit none
    double precision                                               , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorGalaxyMergerTreePhysical), intent(inout)              :: self
    integer                                                                                     :: i
    class           (*                                            )               , pointer     :: extractor_
    !$GLC attributes unused :: self

    allocate(unitsInSI(1+nodePropertyExtractorGalaxyMergerTreeCount))
    unitsInSI(1)=gigaYear
    do i=1,nodePropertyExtractorGalaxyMergerTreeCount
       extractor_ => nodePropertyExtractorGalaxyMergerTree_(i)%extractor_ 
       select type (extractor_)
       class is (nodePropertyExtractorScalar)
          unitsInSI(1+i)=extractor_%unitsInSI()
       class default
          unitsInSI(1+i)=0.0d0
          call Error_Report("unexpected class"//{introspection:location})
       end select
    end do
    return
  end function galaxyMergerTreePhysicalUnitsInSI
