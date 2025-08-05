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
  Implements a node operator class that counts the number of recent major mergers between nodes prior to each output time.
  !!}

  use :: Output_Times           , only : outputTimesClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <enumeration>
   <name>intervalType</name>
   <description>Options for ``recent'' major merger interval types.</description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>private</visibility>
   <entry label="absolute" />
   <entry label="dynamical"/>
  </enumeration>
  !!]

  !![
  <nodeOperator name="nodeOperatorNodeMajorMergerRecentCount">
   <description>A node operator class that counts the number of recent major mergers between nodes prior to each output time.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorNodeMajorMergerRecentCount
     !!{
     A node operator class that counts the number of recent major mergers between nodes prior to each output time.
     !!}
     private
     class           (outputTimesClass           ), pointer :: outputTimes_                 => null()
     class           (darkMatterHaloScaleClass   ), pointer :: darkMatterHaloScale_         => null()
     integer                                                :: nodeMajorMergerRecentCountID
     type            (enumerationIntervalTypeType)          :: intervalType
     double precision                                       :: intervalRecent                        , massRatioMajor
     logical                                                :: intervalFromInfall
  contains
     final     ::                   nodeMajorMergerRecentCountDestructor
     procedure :: nodesMerge     => nodeMajorMergerRecentCountNodesMerge
     procedure :: nodeInitialize => nodeMajorMergerRecentCountNodeInitialize
     procedure :: nodePromote    => nodeMajorMergerRecentCountNodePromote
  end type nodeOperatorNodeMajorMergerRecentCount
  
  interface nodeOperatorNodeMajorMergerRecentCount
     !!{
     Constructors for the \refClass{nodeOperatorNodeMajorMergerRecentCount} node operator class.
     !!}
     module procedure nodeMajorMergerRecentCountConstructorParameters
     module procedure nodeMajorMergerRecentCountConstructorInternal
  end interface nodeOperatorNodeMajorMergerRecentCount
  
contains

  function nodeMajorMergerRecentCountConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorNodeMajorMergerRecentCount} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorNodeMajorMergerRecentCount)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (outputTimesClass                      ), pointer       :: outputTimes_
    class           (darkMatterHaloScaleClass              ), pointer       :: darkMatterHaloScale_
    double precision                                                        :: intervalRecent      , massRatioMajor
    type            (varying_string                        )                :: intervalType
    logical                                                                 :: intervalFromInfall

    !![
    <inputParameter>
      <name>massRatioMajor</name>
      <defaultValue>0.25d0</defaultValue>
      <description>
	The mass ratio ($M_2/M_1$ where $M_2 &lt; M_1$) of merging halos above which the merger should be considered to be
	``major''.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>intervalRecent</name>
      <defaultValue>2.0d0</defaultValue>
      <description>
	The time interval used to define ``recent'' mergers. This parameter is in units of Gyr if {\normalfont \ttfamily
	[intervalType]}$=${\normalfont \ttfamily absolute}, or in units of the halo dynamical time if {\normalfont \ttfamily
	[intervalType]}$=${\normalfont \ttfamily dynamical}.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>intervalType</name>
      <defaultValue>var_str('dynamical')</defaultValue>
      <description>
	Specifies the units for the {\normalfont \ttfamily [intervalRecent]} parameter. If set to {\normalfont \ttfamily absolute}
	then {\normalfont \ttfamily [intervalRecent]} is given in Gyr, while if set to {\normalfont \ttfamily dynamical}
	{\normalfont \ttfamily [intervalRecent]} is given in units of the halo dynamical time.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>intervalFromInfall</name>
      <defaultValue>.false.</defaultValue>
      <description>
	Specifies whether ``recent'' for satellite galaxies is measured from the current time, or from the time at which they were
	last isolated.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="outputTimes"         name="outputTimes_"         source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodeOperatorNodeMajorMergerRecentCount(massRatioMajor,intervalRecent,enumerationIntervalTypeEncode(char(intervalType),includesPrefix=.false.),intervalFromInfall,outputTimes_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"        />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function nodeMajorMergerRecentCountConstructorParameters

  function nodeMajorMergerRecentCountConstructorInternal(massRatioMajor,intervalRecent,intervalType,intervalFromInfall,outputTimes_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorNodeMajorMergerRecentCount} node operator class.
    !!}
    use :: Galacticus_Nodes, only : defaultBasicComponent
    implicit none
    type            (nodeOperatorNodeMajorMergerRecentCount)                        :: self
    double precision                                        , intent(in   )         :: intervalRecent     , massRatioMajor
    type            (enumerationIntervalTypeType           ), intent(in   )         :: intervalType
    logical                                                 , intent(in   )         :: intervalFromInfall
    class           (outputTimesClass                      ), intent(in   ), target :: outputTimes_
    class           (darkMatterHaloScaleClass              ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="massRatioMajor, intervalRecent, intervalType, intervalFromInfall, *outputTimes_, *darkMatterHaloScale_"/>
    !!]
    
    !![
    <addMetaProperty component="basic" name="nodeMajorMergerRecentCount" id="self%nodeMajorMergerRecentCountID" type="integer" rank="1" isCreator="yes"/>
    !!]
    return
  end function nodeMajorMergerRecentCountConstructorInternal

  subroutine nodeMajorMergerRecentCountDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorNodeMajorMergerRecentCount} node operator class.
    !!}
    implicit none
    type(nodeOperatorNodeMajorMergerRecentCount), intent(inout) :: self
     
    !![
    <objectDestructor name="self%outputTimes_"        />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine nodeMajorMergerRecentCountDestructor

  subroutine nodeMajorMergerRecentCountNodeInitialize(self,node)
    !!{
    Record counts of galaxy-galaxy major mergers.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorNodeMajorMergerRecentCount), intent(inout), target  :: self
    type (treeNode                              ), intent(inout), target  :: node
    class(nodeComponentBasic                    )               , pointer :: basic

    basic => node%basic()
    call basic%integerRank1MetaPropertySet(self%nodeMajorMergerRecentCountID,spread(0,1,self%outputTimes_%count()))
    return
  end subroutine nodeMajorMergerRecentCountNodeInitialize

  subroutine nodeMajorMergerRecentCountNodePromote(self,node)
    !!{
    Record counts of galaxy-galaxy major mergers.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorNodeMajorMergerRecentCount), intent(inout) :: self
    type (treeNode                              ), intent(inout) :: node
    class(nodeComponentBasic                    ), pointer       :: basic, basicParent
    
    basic       => node       %basic()
    basicParent => node%parent%basic()
    call basic      %integerRank1MetaPropertySet(                                                                            &
         &                                                                                self%nodeMajorMergerRecentCountID, &
         &                                       +basic      %integerRank1MetaPropertyGet(self%nodeMajorMergerRecentCountID) &
         &                                       +basicParent%integerRank1MetaPropertyGet(self%nodeMajorMergerRecentCountID) &
         &                                      )
    call basicParent%integerRank1MetaPropertySet(                                                                            &
         &                                                                                self%nodeMajorMergerRecentCountID, &
         &                                       spread(0,1,self%outputTimes_%count())                                       &
         &                                      )
    return
  end subroutine nodeMajorMergerRecentCountNodePromote
  
  subroutine nodeMajorMergerRecentCountNodesMerge(self,node)
    !!{
    Record counts of galaxy-galaxy major mergers.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (nodeOperatorNodeMajorMergerRecentCount), intent(inout)              :: self
    type            (treeNode                              ), intent(inout)              :: node
    class           (nodeComponentBasic                    ), pointer                    :: basic                , basicParent  , &
         &                                                                                  basicDescendantParent
    type            (treeNode                              ), pointer                    :: nodeDescendant
    integer                                                 , dimension(:) , allocatable :: countMajorMergers
    integer         (c_size_t                              )                             :: i
    double precision                                                                     :: timeBase             , intervalRecent

    ! Return immediately if this is not a major merger.
    basic       => node       %basic()
    basicParent => node%parent%basic()
    if (basic%mass() < self%massRatioMajor*basicParent%mass()) return 
    ! Determine the interval.
    select case (self%intervalType%ID)
    case (intervalTypeAbsolute %ID)
       intervalRecent=self%intervalRecent
    case (intervalTypeDynamical%ID)
       intervalRecent=self%intervalRecent*self%darkMatterHaloScale_%timescaleDynamical(node)
    case default
       intervalRecent=0.0d0
       call Error_Report('unrecognized recent time interval type'//{introspection:location})
    end select    
    ! Get current count of mergers.
    countMajorMergers=basicParent%integerRank1MetaPropertyGet(self%nodeMajorMergerRecentCountID)
    ! Check each output time.
    do i=1_c_size_t,self%outputTimes_%count()
       ! Determine the base time to measure the interval from.
       if (self%intervalFromInfall) then
          if (node%parent%isSatellite()) then
             timeBase=basicParent%timeLastIsolated()
          else
             timeBase=self%outputTimes_%time(i)
             nodeDescendant => node%parent
             do while (associated(nodeDescendant))
                if (nodeDescendant%isPrimaryProgenitor()) then
                   nodeDescendant => nodeDescendant%parent
                else
                   if (associated(nodeDescendant%parent)) then
                      basicDescendantParent => nodeDescendant%parent%basic()
                      timeBase=min(timeBase,basicDescendantParent%time())
                   end if
                   exit
                end if
             end do
          end if
       else
          timeBase=self%outputTimes_%time(i)
       end if
       ! Check for a recent merger.
       if     (                                             &
            &   basic%time() <= timeBase                    &
            &  .and.                                        &
            &   basic%time() >  timeBase-intervalRecent     &
            & ) countMajorMergers(i)=countMajorMergers(i)+1
    end do
    ! Update the count of recent mergers.
    call basicParent%integerRank1MetaPropertySet(self%nodeMajorMergerRecentCountID,countMajorMergers)
    return
  end subroutine nodeMajorMergerRecentCountNodesMerge
