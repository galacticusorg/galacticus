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
  Implements a property extractor class for the star formation history of a component.
  !!}
  
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType
  use :: Star_Formation_Histories  , only : starFormationHistoryClass
  use :: Output_Times              , only : outputTimesClass
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorStarFormationHistoryMass">
    <description>A property extractor class for the star formation history of a component.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList2D) :: nodePropertyExtractorStarFormationHistoryMass
     !!{
     A property extractor class for the star formation history of a component.
     !!}
     private
     class(starFormationHistoryClass   ), pointer :: starFormationHistory_ => null()
     class(outputTimesClass            ), pointer :: outputTimes_          => null()
     type (enumerationComponentTypeType)          :: component
   contains
     final     ::                 starFormationHistoryMassDestructor
     procedure :: elementCount => starFormationHistoryMassElementCount
     procedure :: extract      => starFormationHistoryMassExtract
     procedure :: names        => starFormationHistoryMassNames
     procedure :: descriptions => starFormationHistoryMassDescriptions
     procedure :: unitsInSI    => starFormationHistoryMassUnitsInSI
     procedure :: metaData     => starFormationHistoryMassMetaData
  end type nodePropertyExtractorStarFormationHistoryMass
  
  interface nodePropertyExtractorStarFormationHistoryMass
     !!{
     Constructors for the {\normalfont \ttfamily starFormationHistoryMass} output analysis class.
     !!}
     module procedure starFormationHistoryMassConstructorParameters
     module procedure starFormationHistoryMassConstructorInternal
  end interface nodePropertyExtractorStarFormationHistoryMass
      
contains

  function starFormationHistoryMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormationHistoryMass} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode
    implicit none
    type (nodePropertyExtractorStarFormationHistoryMass)                :: self
    type (inputParameters                              ), intent(inout) :: parameters
    class(starFormationHistoryClass                    ), pointer       :: starFormationHistory_
    class(outputTimesClass                             ), pointer       :: outputTimes_
    type (varying_string                               )                :: component
    
    !![
    <inputParameter>
      <name>component</name>
      <source>parameters</source>
      <description>The component from which to extract star formation history.</description>
    </inputParameter>
    <objectBuilder class="starFormationHistory" name="starFormationHistory_" source="parameters"/>
    <objectBuilder class="outputTimes"          name="outputTimes_"          source="parameters"/>
    !!]
    self=nodePropertyExtractorStarFormationHistoryMass(enumerationComponentTypeEncode(char(component),includesPrefix=.false.),starFormationHistory_,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationHistory_"/>
    <objectDestructor name="outputTimes_"         />
    !!]
    return
  end function starFormationHistoryMassConstructorParameters

  function starFormationHistoryMassConstructorInternal(component,starFormationHistory_,outputTimes_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily starFormationHistoryMass} property extractor class.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDisk, componentTypeSpheroid, componentTypeNuclearStarCluster, componentTypeAll
    use :: Error                     , only : Error_Report
    implicit none
    type (nodePropertyExtractorStarFormationHistoryMass)                        :: self
    class(starFormationHistoryClass                    ), intent(in   ), target :: starFormationHistory_
    class(outputTimesClass                             ), intent(in   ), target :: outputTimes_
    type (enumerationComponentTypeType                 ), intent(in   )         :: component
    !![
    <constructorAssign variables="component, *starFormationHistory_, *outputTimes_"/>
    !!]
    
    if     (                                                                                                                                 &
         &   component /= componentTypeDisk                                                                                                  &
         &  .and.                                                                                                                            &
         &   component /= componentTypeSpheroid                                                                                              &
         &  .and.                                                                                                                            &
         &   component /= componentTypeNuclearStarCluster                                                                                    &
         &  .and.                                                                                                                            &
         &   component /= componentTypeAll                                                                                                   &
         & ) call Error_Report("only 'disk', 'spheroid', 'nuclearStarCluster' and 'all' components are supported"//{introspection:location})    
    return
  end function starFormationHistoryMassConstructorInternal

  subroutine starFormationHistoryMassDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily sed} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorStarFormationHistoryMass), intent(inout) :: self
    
    !![
    <objectDestructor name="self%starFormationHistory_"/>
    <objectDestructor name="self%outputTimes_"         />
    !!]
    return
  end subroutine starFormationHistoryMassDestructor

  integer function starFormationHistoryMassElementCount(self)
    !!{
    Return the number of elements in the {\normalfont \ttfamily starFormationHistoryMass} property extractors.
    !!}
    implicit none
    class(nodePropertyExtractorStarFormationHistoryMass), intent(inout) :: self

    starFormationHistoryMassElementCount=1
    return
  end function starFormationHistoryMassElementCount

  function starFormationHistoryMassExtract(self,node,instance)
    !!{
    Implement a {\normalfont \ttfamily starFormationHistoryMass} property extractor.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentDisk, nodeComponentSpheroid, nodeComponentNSC
    use :: Galactic_Structure_Options, only : componentTypeDisk, componentTypeSpheroid, componentTypeNuclearStarCluster, componentTypeAll
    use :: Histories                 , only : history
    implicit none
    double precision                                               , dimension(:,:,:  ), allocatable :: starFormationHistoryMassExtract
    class           (nodePropertyExtractorStarFormationHistoryMass), intent(inout)                   :: self
    type            (treeNode                                     ), intent(inout)                   :: node
    type            (multiCounter                                 ), intent(inout)     , optional    :: instance
    double precision                                               , dimension(:,:  )  , allocatable :: masses
    class           (nodeComponentDisk                            )                    , pointer     :: disk
    class           (nodeComponentSpheroid                        )                    , pointer     :: spheroid
    class           (nodeComponentNSC                             )                    , pointer     :: nuclearStarCluster
    type            (history                                      )                                  :: starFormationHistory           , starFormationHistoryDisk              , &
         &                                                                                              starFormationHistorySpheroid   , starFormationHistoryNuclearStarCluster
    !$GLC attributes unused :: instance

    ! Get the relevant star formation history.
    select case (self%component%ID)
    case (componentTypeDisk              %ID)
       disk                 => node              %disk                ()
       starFormationHistory =  disk              %starFormationHistory()
    case (componentTypeSpheroid          %ID)
       spheroid             => node              %spheroid            ()
       starFormationHistory =  spheroid          %starFormationHistory()
    case (componentTypeNuclearStarCluster%ID)
       nuclearStarCluster   => node              %NSC                 ()
       starFormationHistory =  nuclearStarCluster%starFormationHistory()
    case (componentTypeAll     %ID)
       spheroid                               => node              %spheroid            ()
       disk                                   => node              %disk                ()
       nuclearStarCluster                     => node              %NSC                 ()
       starFormationHistoryDisk               =  disk              %starFormationHistory()
       starFormationHistorySpheroid           =  spheroid          %starFormationHistory()
       starFormationHistoryNuclearStarCluster =  nuclearStarCluster%starFormationHistory()
       if      (     starFormationHistoryDisk%exists()) then
          if      (     starFormationHistorySpheroid%exists() .and.      starFormationHistoryNuclearStarCluster%exists() ) then
             starFormationHistory= starFormationHistoryDisk               &
                  &               +starFormationHistorySpheroid           &
                  &               +starFormationHistoryNuclearStarCluster

          else if (.not.starFormationHistorySpheroid%exists() .and.      starFormationHistoryNuclearStarCluster%exists() ) then
             starFormationHistory= starFormationHistoryDisk               &
                  &               +starFormationHistoryNuclearStarCluster
          else if (     starFormationHistorySpheroid%exists() .and. .not.starFormationHistoryNuclearStarCluster%exists() ) then
             starFormationHistory= starFormationHistoryDisk               &
                  &               +starFormationHistorySpheroid
          else 
             starFormationHistory= starFormationHistoryDisk
          end if
       else if (.not.starFormationHistoryDisk%exists()) then
          if      (     starFormationHistorySpheroid%exists() .and.      starFormationHistoryNuclearStarCluster%exists() ) then
             starFormationHistory= starFormationHistorySpheroid           &
                  &               +starFormationHistoryNuclearStarCluster
          else if (.not.starFormationHistorySpheroid%exists() .and.      starFormationHistoryNuclearStarCluster%exists() ) then
             starFormationHistory= starFormationHistoryNuclearStarCluster
          else if (     starFormationHistorySpheroid%exists() .and. .not.starFormationHistoryNuclearStarCluster%exists() ) then
             starFormationHistory= starFormationHistorySpheroid
          end if
       end if
    end select
    if (starFormationHistory%exists()) then
       masses=self%starFormationHistory_%masses(node,starFormationHistory,allowTruncation=.true.)
       allocate(starFormationHistoryMassExtract(size(masses,dim=1),size(masses,dim=2),1))
       starFormationHistoryMassExtract(:,:,1)=masses
    else
       allocate(starFormationHistoryMassExtract(0                 ,0                 ,1))
    end if
    return
  end function starFormationHistoryMassExtract

  subroutine starFormationHistoryMassNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily starFormationHistoryMass} properties.
    !!}
    use :: Galactic_Structure_Options, only : enumerationComponentTypeDecode
    implicit none
    class(nodePropertyExtractorStarFormationHistoryMass), intent(inout)                             :: self
    type(varying_string                                ), intent(inout), dimension(:) , allocatable :: names

    allocate(names(1))
    names(1)=enumerationComponentTypeDecode(self%component,includePrefix=.false.)//"StarFormationHistoryMass"
    return
  end subroutine starFormationHistoryMassNames

  subroutine starFormationHistoryMassDescriptions(self,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily starFormationHistoryMass} property.
    !!}
    use :: Galactic_Structure_Options, only : enumerationComponentTypeDecode
    implicit none
    class(nodePropertyExtractorStarFormationHistoryMass), intent(inout)                             :: self
    type (varying_string                               ), intent(inout), dimension(:) , allocatable :: descriptions

    allocate(descriptions(1))
    descriptions(1)="Star formation history for the "//enumerationComponentTypeDecode(self%component,includePrefix=.false.)//" [M☉ Gyr¯¹]."
    return
  end subroutine starFormationHistoryMassDescriptions

  function starFormationHistoryMassUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily starFormationHistoryMass} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, gigaYear
    implicit none
    double precision                                               , allocatable  , dimension(:) :: starFormationHistoryMassUnitsInSI
    class           (nodePropertyExtractorStarFormationHistoryMass), intent(inout)               :: self

    allocate(starFormationHistoryMassUnitsInSI(1))
    starFormationHistoryMassUnitsInSI(1)=massSolar/gigaYear
    return
  end function starFormationHistoryMassUnitsInSI
  
  subroutine starFormationHistoryMassMetaData(self,node,time,metaDataRank0,metaDataRank1)
    !!{
    Return metadata associated with the {\normalfont \ttfamily starFormationHistoryMass} properties.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBasic
    use :: Star_Formation_Histories, only : starFormationHistoryAgesFixedPerOutput
    implicit none
    class           (nodePropertyExtractorStarFormationHistoryMass), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    double precision                                               , intent(in   ) :: time
    type            (doubleHash                                   ), intent(inout) :: metaDataRank0
    type            (rank1DoubleHash                              ), intent(inout) :: metaDataRank1
    integer         (c_size_t                                     )                :: indexOutput
    !$GLC attributes unused :: metaDataRank0

    call    metaDataRank1%set('metallicity',self%starFormationHistory_%metallicityBoundaries(                                              ))
    if (self%starFormationHistory_%ageDistribution() == starFormationHistoryAgesFixedPerOutput) then
       indexOutput =  self%outputTimes_%index(time,findClosest=.true.)
       call metaDataRank1%set('time'       ,self%starFormationHistory_%times                (indexOutput=indexOutput,allowTruncation=.true.))
    end if
    return
  end subroutine starFormationHistoryMassMetaData
