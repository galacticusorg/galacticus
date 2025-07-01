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
  Implements a star formation histories class which records star formation in logarithmically-sized time
  bins of fixed age and split by metallicity.
  !!}

  use :: Geometry_Lightcones, only : geometryLightcone
  use :: Cosmology_Functions, only : cosmologyFunctions
  
  !![
  <starFormationHistory name="starFormationHistoryFixedAges" recursive="yes">
   <description>
     A star formation histories class which records star formation in logarithmically-sized time bins of fixed age and split by
     metallicity. The minimum age is specified via the {\normalfont \ttfamily [ageMinimum]} parameter (the maximum age is always
     the age of the universe), with the number of ages specified via {\normalfont \ttfamily [countAges]}. (One additional bin, at
     age zero, is always added.) This class is intended for use with lightcone output where the lightcone crossing times for each
     node can be computed in advance. One star formation history is computed for each lightcone crossing.
     
     The time associated with each bin is the maximum time for which star formation will be accumulated to the bin, with the
     minimum time corresponding to the value associated with the previous bin (or $t=0$ for the first bin).
     
     The metallicity bins are arranged logarithmically in metallicity with {\normalfont \ttfamily [countMetallicities]} bins
     between {\normalfont \ttfamily [metallicityMinimum]} and {\normalfont \ttfamily [metallicityMaximum]} (specified in Solar
     units). Note that the metallicity associated with each bin is the maximum metallicity for that bin, with the minimum
     metallicity corresponding to the value associated with the previous bin (or zero metallicity for the first bin). Note that a
     final bin, extending to infinite metallicity, is always added automatically. If {\normalfont \ttfamily
     [countMetallicities]}$=0$ is set, then the star formation history is not split by metallicity (i.e. a single metallicity bin
     encompassing all metallicities from zero to infinity is used). Alternatively, specific metallicity bin boundaries can be set
     via the {\normalfont \ttfamily [metallicityBoundaries]} parameter---a final boundary corresponding to infinity is always
     added automatically.
   </description>
   <deepCopy>
     <ignore variables="recursiveSelf"/>
   </deepCopy>
   <stateStorable>
     <exclude variables="recursiveSelf"/>
   </stateStorable>
  </starFormationHistory>
  !!]
  type, extends(starFormationHistoryClass) :: starFormationHistoryFixedAges
     !!{
     A star formation histories class which records star formation split by metallicity.
     !!}
     private
     logical                                                                    :: isRecursive                  , parentDeferred
     class           (starFormationHistoryFixedAges), pointer                   :: recursiveSelf       => null()
     class           (geometryLightconeClass       ), pointer                   :: geometryLightcone_  => null()
     class           (cosmologyFunctionsClass      ), pointer                   :: cosmologyFunctions_ => null()
     double precision                                                           :: ageMinimum                   , ageMaximum        , &
          &                                                                        metallicityMaximum           , metallicityMinimum
     integer         (c_size_t                     )                            :: countAges                    , countMetallicities
     integer                                                                    :: timesCrossingID              , countRetain       , &
          &                                                                        createdInID
     integer         (kind_int8                    )                            :: uniqueIDPrevious
     double precision                                                           :: timePrevious
     double precision                               , allocatable, dimension(:) :: metallicityTable
   contains
     final     ::                          fixedAgesDestructor
     procedure :: create                => fixedAgesCreate
     procedure :: update                => fixedAgesUpdate
     procedure :: move                  => fixedAgesMove
     procedure :: rate                  => fixedAgesRate
     procedure :: scales                => fixedAgesScales
     procedure :: times                 => fixedAgesTimes
     procedure :: timeNext              => fixedAgesTimeNext
     procedure :: masses                => fixedAgesMasses
     procedure :: metallicityBoundaries => fixedAgesMetallicityBoundaries
     procedure :: ageDistribution       => fixedAgesAgeDistribution
     procedure :: descriptor            => fixedAgesDescriptor
     procedure :: deepCopy              => fixedAgesDeepCopy
     procedure :: deepCopyReset         => fixedAgesDeepCopyReset
     procedure :: deepCopyFinalize      => fixedAgesDeepCopyFinalize
  end type starFormationHistoryFixedAges

  interface starFormationHistoryFixedAges
     !!{
     Constructors for the \refClass{starFormationHistoryFixedAges} star formation history class.
     !!}
     module procedure fixedAgesConstructorParameters
     module procedure fixedAgesConstructorInternal
  end interface starFormationHistoryFixedAges

  ! Effective infinite metallicity.
  double precision, parameter :: metallicityInfinite=huge(1.0d0)

contains

  recursive function fixedAgesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{starFormationHistoryFixedAges} star formation history class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationHistoryFixedAges), target                     :: self
    type            (inputParameters              ), intent(inout), target      :: parameters
    class           (geometryLightconeClass       )               , pointer     :: geometryLightcone_
    class           (cosmologyFunctionsClass      )               , pointer     :: cosmologyFunctions_
    double precision                               , dimension(:) , allocatable :: metallicityBoundaries
    double precision                                                            :: metallicityMinimum   , metallicityMaximum, &
         &                                                                         ageMinimum
    integer         (c_size_t)                                                  :: countMetallicities   , countAges

    !![
    <objectBuilder class="geometryLightcone"  name="geometryLightcone_"  source="parameters"/>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <inputParameter>
      <name>ageMinimum</name>
      <defaultValue>0.01d0</defaultValue>
      <description>The minimum age to use in tabulations of star formation histories [Gyr].</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>countAges</name>
      <defaultValue>10_c_size_t</defaultValue>
      <description>The maximum number of ages to track in any star formation history.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    if (parameters%isPresent('metallicityBoundaries')) then
       countMetallicities=parameters%count('metallicityBoundaries')
       allocate(metallicityBoundaries(countMetallicities+1))
       !![
       <inputParameter>
         <name>metallicityBoundaries</name>
         <description>The metallicities corresponding to boundaries between metallicity bins to use when tabulating star formation histories.</description>
         <source>parameters</source>
         <variable>metallicityBoundaries(1:size(metallicityBoundaries)-1)</variable>
         <type>real</type>
         <cardinality>0..*</cardinality>
       </inputParameter>
       !!]
       metallicityBoundaries(size(metallicityBoundaries))=metallicityInfinite
    else
       !![
       <inputParameter>
         <name>countMetallicities</name>
         <defaultValue>10_c_size_t</defaultValue>
         <description>The number of bins in metallicity to use when tabulating star formation histories.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>metallicityMinimum</name>
         <defaultValue>1.0d-4</defaultValue>
         <description>The upper limit to the metallicity in the lowest metallicity bin when tabulating star formation histories [Solar units].</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>metallicityMaximum</name>
         <defaultValue>1.0d+1</defaultValue>
         <description>The upper limit to the metallicity in the highest metallicity bin when tabulating star formation histories [Solar units].</description>
         <source>parameters</source>
       </inputParameter>
       !!]
    end if
    !![
    <conditionalCall>
     <call>self=starFormationHistoryFixedAges(cosmologyFunctions_,geometryLightcone_,ageMinimum,countAges{conditions})</call>
     <argument name="metallicityBoundaries" value="metallicityBoundaries" condition="     parameters%isPresent('metallicityBoundaries')"/>
     <argument name="countMetallicities"    value="countMetallicities"    condition=".not.parameters%isPresent('metallicityBoundaries')"/>
     <argument name="metallicityMinimum"    value="metallicityMinimum"    condition=".not.parameters%isPresent('metallicityBoundaries')"/>
     <argument name="metallicityMaximum"    value="metallicityMaximum"    condition=".not.parameters%isPresent('metallicityBoundaries')"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="geometryLightcone_" />
    !!]
    return
  end function fixedAgesConstructorParameters

  recursive function fixedAgesConstructorInternal(cosmologyFunctions_,geometryLightcone_,ageMinimum,countAges,metallicityBoundaries,countMetallicities,metallicityMinimum,metallicityMaximum) result(self)
    !!{
    Internal constructor for the \refClass{starFormationHistoryFixedAges} star formation history class.
    !!}
    use :: Error                     , only : Error_Report
    use :: Galactic_Structure_Options, only : componentTypeMax, componentTypeMin
    use :: Numerical_Ranges          , only : Make_Range      , rangeTypeLogarithmic
    implicit none
    type            (starFormationHistoryFixedAges)                                        :: self
    double precision                               , intent(in   ), dimension(:), optional :: metallicityBoundaries
    double precision                               , intent(in   )              , optional :: metallicityMinimum   , metallicityMaximum
    double precision                               , intent(in   )                         :: ageMinimum
    integer         (c_size_t                     ), intent(in   )                         :: countAges
    integer         (c_size_t                     ), intent(in   )              , optional :: countMetallicities
    class           (cosmologyFunctionsClass      ), intent(in   ), target                 :: cosmologyFunctions_
    class           (geometryLightconeClass       ), intent(in   ), target                 :: geometryLightcone_
    
    !![
    <constructorAssign variables="ageMinimum, countAges, metallicityMinimum, metallicityMaximum, countMetallicities, *cosmologyFunctions_, *geometryLightcone_"/>
    !!]

    ! Validate metallicity argument and construct the table of metallicities.
    if (present(metallicityBoundaries)) then
       if     (                                &
            &   present(countMetallicities   ) &
            &  .or.                            &
            &   present(metallicityMinimum   ) &
            &  .or.                            &
            &   present(metallicityMaximum   ) &
            & ) call Error_Report('specify either a list of metallicity boundaries, or a range, not both'//{introspection:location})
       allocate(self%metallicityTable(size(metallicityBoundaries)))
       self%metallicityTable  =     metallicityBoundaries
       self%countMetallicities=size(metallicityBoundaries)-1
    else
       if     (                                &
            &   present(metallicityBoundaries) &
            & ) call Error_Report('specify either a list of metallicity boundaries, or a range, not both'//{introspection:location})
       if     (                                     &
            &   .not.present(countMetallicities   ) &
            &  .or.                                 &
            &   .not.present(metallicityMinimum   ) &
            &  .or.                                 &
            &   .not.present(metallicityMaximum   ) &
            & ) call Error_Report('metallicity range is incompletely specified'//{introspection:location})
       ! Construct a table of metallicities at which to tabulate. Add an extra bin since we want to catch all metallicities,
       ! including those below and above the maximum. A single bin is not allowed, but zero bins implies that no metallicity
       ! resolution is required.
       select case (countMetallicities)
       case (:-1,1)
          call Error_Report('number of bins must be 0, or greater than 1'//{introspection:location})
       case default
          allocate(self%metallicityTable(countMetallicities+1))
          if (countMetallicities > 1) self%metallicityTable(1:countMetallicities)=Make_Range(metallicityMinimum,metallicityMaximum,int(countMetallicities),rangeType=rangeTypeLogarithmic)
          self%metallicityTable(countMetallicities+1)=metallicityInfinite
       end select
    end if
    ! Create a meta-property to store lightcone crossing times.
    !![
    <addMetaProperty component="basic" name="starFormationHistoriesFixedAgeTimesCrossing" id="self%timesCrossingID"                    rank="1" isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="basic" name="starFormationHistoriesFixedAgeCreatedIn"     id="self%createdInID"     type="longInteger" rank="0"                  isCreator="yes"/>
    !!]
    ! Initialize state.
    self%uniqueIDPrevious =-huge(0_kind_int8)
    self%timePrevious     =-huge(0.0d0      )
    ! Set the maximum age.
    self%ageMaximum=self%cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0)
    ! Set recursive properties.
    self%parentDeferred=.false.
    self%isRecursive   =.false.
    return
  end function fixedAgesConstructorInternal

  subroutine fixedAgesDestructor(self)
    !!{
    Destructor for the \refClass{starFormationHistoryFixedAges} star formation histories class.
    !!}
    implicit none
    type(starFormationHistoryFixedAges), intent(inout) :: self

    !![
    <objectDestructor name="self%geometryLightcone_" />
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine fixedAgesDestructor

  subroutine fixedAgesCreate(self,node,historyStarFormation,timeBegin,timeEnd)
    !!{
    Create the history required for storing star formation history.
    !!}
    use :: Display             , only : displayMessage    , displayIndent      , displayUnindent
    use :: Galacticus_Nodes    , only : nodeComponentBasic
    use :: Histories           , only : history
    use :: Error               , only : Error_Report
    use :: ISO_Varying_String  , only : varying_string    , var_str
    use :: String_Handling     , only : operator(//)
    use :: Numerical_Ranges    , only : Make_Range        , rangeTypeLogarithmic
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    class           (starFormationHistoryFixedAges), intent(inout)              :: self
    type            (treeNode                     ), intent(inout), target      :: node
    type            (history                      ), intent(inout)              :: historyStarFormation
    double precision                               , intent(in   )              :: timeBegin
    double precision                               , intent(in   ), optional    :: timeEnd
    double precision                               , dimension(:) , allocatable :: timesNodeCrossing           , timesNodeCrossingPrevious        , &
         &                                                                         timesNodeCrossingTmp
    double precision                               , dimension(0)               :: timesNodeCrossingNull
    class           (nodeComponentBasic           ), pointer                    :: basic
    double precision                               , parameter                  :: toleranceRelative    =1.0d-6, toleranceAbsolute        =1.0d-9
    double precision                                                            :: timeNodeStart               , timeNodeEnd                      , &
         &                                                                         timeNodeCrossing
    type            (varying_string               )                             :: message
    character       (len=16                       )                             :: label
    integer                                                                     :: i
    !$GLC attributes unused :: timeBegin, timeEnd

    ! Call the recursive copy if necessary.
    if (self%isRecursive) then
       call self%recursiveSelf%create(node,historyStarFormation,timeBegin,timeEnd)
       return
    end if
    ! Find lightcone crossing times.
    basic                     =>     node %basic()
    timeNodeStart             =  max(basic%time (),self%geometryLightcone_%timeMinimum())
    timeNodeEnd               =                    self%geometryLightcone_%timeMaximum()
    timeNodeCrossing          =  self %geometryLightcone_%timeLightconeCrossing    (node                ,timeNodeStart,timeNodeEnd,timesNodeCrossing)
    timesNodeCrossingPrevious =  basic                   %floatRank1MetaPropertyGet(self%timesCrossingID                                            )
    if (size(timesNodeCrossingPrevious) > 0) then
       ! Rounding errors can lead to tiny shifts in crossing times which can (very occasionally) lead to a previously-found
       ! crossing time being missed if the node is now very close to that crossing time. Check for such occurrences here and add
       ! the crossing time back into the list.
       if (size(timesNodeCrossing) == size(timesNodeCrossingPrevious)-1) then
          if (Values_Agree(basic%time(),timesNodeCrossingPrevious(1),absTol=toleranceAbsolute) .and. basic%time() <= timesNodeCrossingPrevious(1)) then
             call move_alloc(timesNodeCrossing,timesNodeCrossingTmp)
             allocate(timesNodeCrossing(size(timesNodeCrossingTmp)+1))
             timesNodeCrossing(1)=timesNodeCrossingPrevious(1)
             if (size(timesNodeCrossingTmp) > 0) timesNodeCrossing(2:size(timesNodeCrossing))=timesNodeCrossingTmp
             deallocate(timesNodeCrossingTmp)
          end if
       end if
       ! Check for cases where the first crossing time is equal to (or very close to) the current time. These crossing have
       ! already been processed and so can be ignored.
       if (size(timesNodeCrossing) == size(timesNodeCrossingPrevious)+1) then
          if (Values_Agree(timesNodeCrossing(1),basic%time(),absTol=toleranceAbsolute)) then
             ! The first crossing time is equal to the current time.
             call move_alloc(timesNodeCrossing,timesNodeCrossingTmp)
             allocate(timesNodeCrossing(size(timesNodeCrossingTmp)-1))
             if (size(timesNodeCrossing) > 0) timesNodeCrossing=timesNodeCrossingTmp(2:size(timesNodeCrossingTmp))
             deallocate(timesNodeCrossingTmp)
          end if
       end if
       ! Validate consistency in the lightcone crossing times.
       if (size(timesNodeCrossing) /= size(timesNodeCrossingPrevious)) then
          write (label,'(e16.10)') basic%time()
          call displayIndent(var_str("number of crossing times has changed for node ")//node%index()//' at time '//trim(adjustl(label))//' Gyr')
          call displayMessage("times (new | old) are:")
          do i=1,max(size(timesNodeCrossing),size(timesNodeCrossingPrevious))
             if (i <= size(timesNodeCrossing)) then
                write (label,'(e16.10)') timesNodeCrossing        (i)
                message=         label
             else
                message=         "              "
             end if
             message=message//" | "
             if (i <= size(timesNodeCrossingPrevious)) then
                write (label,'(e16.10)') timesNodeCrossingPrevious(i)
                message=message//label
             else
                message=message//"              "
             end if
             call displayMessage(message)
          end do
          call displayUnindent("")
          call Error_Report("invalid crossing times"//{introspection:location})
       end if
       if (.not.all(Values_Agree(timesNodeCrossing,basic%floatRank1MetaPropertyGet(self%timesCrossingID),relTol=toleranceRelative))) then
          call displayIndent(var_str("crossing times have changed for node ")//node%index())
          call displayMessage("times (new | old | difference) are:")
          do i=1,size(timesNodeCrossing)
             write (label,'(e16.10)') timesNodeCrossing(i)
             message=         label//" | "
             write (label,'(e16.10)')                          timesNodeCrossingPrevious(i)
             message=message//label//" | "
             write (label,'(e16.10)') abs(timesNodeCrossing(i)-timesNodeCrossingPrevious(i))
             message=message//label
             call displayMessage(message)
          end do
          call displayUnindent("")
          call Error_Report("invalid crossing times"//{introspection:location})
       end if
    end if
    ! Create the star formation histories as needed.
    call basic%longIntegerRank0MetaPropertySet(self%createdInID,node%index())
    if (allocated(timesNodeCrossing)) then
       call basic               %floatRank1MetaPropertySet(self%timesCrossingID,timesNodeCrossing)
       call historyStarFormation%create                   (size(timesNodeCrossing)*int(self%countMetallicities+1),int(self%countAges+1))
       ! Here we set the times relative to t=0. This ensures that times are increasing which is expected by history objects.
       historyStarFormation%time(1:self%countAges  )=-Make_Range(self%ageMaximum,self%ageMinimum,int(self%countAges),rangeTypeLogarithmic)
       historyStarFormation%time(  self%countAges+1)=+0.0d0       
       deallocate(timesNodeCrossing)
    else
       ! No lightcone crossings for this node - store a null set of crossing times.
       call basic               %floatRank1MetaPropertySet(self%timesCrossingID,timesNodeCrossingNull)
       call historyStarFormation%create                   (                   0,                    0)
    end if
    return
  end subroutine fixedAgesCreate

  subroutine fixedAgesRate(self,node,historyStarFormation,abundancesFuel,rateStarFormation)
    !!{
    Set the rate the star formation history for {\normalfont \ttfamily node}.
    !!}
    use :: Abundances_Structure, only : abundances        , metallicityTypeLinearByMassSolar
    use :: Arrays_Search       , only : searchArray
    use :: Galacticus_Nodes    , only : nodeComponentBasic, treeNode
    use :: Error               , only : Error_Report
    use :: String_Handling     , only : operator(//)
    use :: ISO_Varying_String  , only : operator(//)      , var_str
    implicit none
    class           (starFormationHistoryFixedAges), intent(inout)               :: self
    type            (treeNode                     ), intent(inout)               :: node
    type            (history                      ), intent(inout)               :: historyStarFormation
    type            (abundances                   ), intent(in   )               :: abundancesFuel
    double precision                               , intent(in   )               :: rateStarFormation
    class           (nodeComponentBasic           ), pointer                     :: basic
    double precision                               , allocatable  , dimension(:) :: timesCrossing       , times
    integer         (c_size_t                     )                              :: iHistory            , iMetallicity, &
         &                                                                          countHistories      , i
    double precision                                                             :: fuelMetallicity     , time

    ! Call the recursive copy if necessary.
    if (self%isRecursive) then
       call self%recursiveSelf%rate(node,historyStarFormation,abundancesFuel,rateStarFormation)
       return
    end if
    ! Check if history exists.
    if (historyStarFormation%exists()) then
       basic => node %basic()
       time  =  basic%time ()
       ! Iterate over all star formation histories being computed.
       countHistories=size(historyStarFormation%data,dim=2)/int(self%countMetallicities+1)
       timesCrossing =basic%floatRank1MetaPropertyGet(self%timesCrossingID)
       if (countHistories /= size(timesCrossing)) then
          block
            type     (varying_string) :: message
            character(len=16        ) :: label
            write (label,'(e16.10)') time
            message=var_str("inconsistent number of star formation histories (")//countHistories//") and lightcone crossing times ("//size(timesCrossing)//") in node "//node%index()//" at time "//trim(adjustl(label))//" Gyr (star formation history was created in progenitor "//basic%longIntegerRank0MetaPropertyGet(self%createdInID)//")"
            call Error_Report(message//{introspection:location})
          end block
       end if
       if (countHistories > 0) then
          allocate(times(self%countAges+1))
          do i=1,countHistories
             ! Find the times corresponding to this star formation history. Note that the times stored in the history object are
             ! relative to t=0, so we increment them by the actual crossing time.
             times=timesCrossing(i)+historyStarFormation%time
             ! Find the point in the table at which to accumulate the star formation rate.
             if (time <= times(1)) then
                iHistory=1
             else
                iHistory=searchArray(times,time)+1
             end if
             ! Find the metallicity bin to accumulate to.
             if    (self%countMetallicities    == 0              ) then
                iMetallicity   =                                                  +1
             else
                fuelMetallicity=abundancesFuel%metallicity(metallicityType=metallicityTypeLinearByMassSolar)
                if (self%metallicityTable  (1) >  fuelMetallicity) then
                   iMetallicity=                                                  +1
                else
                   iMetallicity=searchArray(self%metallicityTable,fuelMetallicity)+1
                end if
             end if
             ! Accumulate to the appropriate time.
             historyStarFormation%data(iHistory,(i-1)*(self%countMetallicities+1)+iMetallicity)=rateStarFormation
          end do
       end if
    end if
    return
  end subroutine fixedAgesRate

  subroutine fixedAgesUpdate(self,node,indexOutput,historyStarFormation)
    !!{
    Update the star formation history after outputting.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (starFormationHistoryFixedAges), intent(inout)              :: self
    type            (treeNode                     ), intent(inout), target      :: node
    integer         (c_size_t                     ), intent(in   )              :: indexOutput
    type            (history                      ), intent(inout)              :: historyStarFormation
    double precision                               , dimension(:) , allocatable :: timesNodeCrossing          , timesNodeCrossingFuture
    class           (nodeComponentBasic           )               , pointer     :: basic
    double precision                               , parameter                  :: toleranceRelative   =1.0d-6
    type            (history                      )                             :: newHistory
    !$GLC attributes unused :: node

    ! Call the recursive copy if necessary.
    if (self%isRecursive) then
       call self%recursiveSelf%update(node,indexOutput,historyStarFormation)
       return
    end if
    ! Determine if this is a new update.
    basic => node%basic()
    if     (                                           &
         &   node %uniqueID() /= self%uniqueIDPrevious &
         &  .or.                                       &
         &   basic%time    () /= self%timePrevious     &
         & ) then
       self%uniqueIDPrevious=node %uniqueID()
       self%timePrevious    =basic%time    ()
       ! Trim crossing times that are no longer needed. Add some tolerance here as crossing times can differ by tiny amounts due
       ! to numerical precision differences between when these are created for the SFH and when differential evolution finds
       ! lightcone crossings.
       timesNodeCrossing      =basic%floatRank1MetaPropertyGet(self%timesCrossingID)
       timesNodeCrossingFuture=pack(timesNodeCrossing,timesNodeCrossing > basic%time()*(1.0d0+toleranceRelative))
       self%countRetain       =size(timesNodeCrossingFuture)
       call basic%floatRank1MetaPropertySet(self%timesCrossingID,timesNodeCrossingFuture)
    end if
    if (historyStarFormation%exists()) then
       if (self%countRetain > 0) then
          ! Retain the required number of histories, discarding the rest,
          call newHistory%create(self%countRetain*int(self%countMetallicities+1),int(self%countAges+1))
          newHistory%rangeType=historyStarFormation%rangeType
          newHistory%time     =historyStarFormation%time
          newHistory%data     =historyStarFormation%data(:,size(historyStarFormation%data,dim=2)-self%countRetain*int(self%countMetallicities+1)+1:size(historyStarFormation%data,dim=2))      
          call historyStarFormation%destroy()
          historyStarFormation=newHistory
          call newHistory          %destroy()
       else
          ! No histories to be retained - simply destroy the history.
          call historyStarFormation%destroy()
       end if
    end if
    return
  end subroutine fixedAgesUpdate

  subroutine fixedAgesMove(self,node1,node2,starFormationHistory1,starFormationHistory2)
    !!{
    Move one history into another.
    !!}
    use :: Galacticus_Nodes    , only : nodeComponentBasic
    use :: ISO_Varying_String  , only : varying_string    , var_str      , operator(//)
    use :: String_Handling     , only : operator(//)
    use :: Numerical_Comparison, only : Values_Agree
    use :: Error               , only : Error_Report
    use :: Display             , only : displayMessage    , displayIndent, displayUnindent
    implicit none
    class           (starFormationHistoryFixedAges), intent(inout)              :: self
    type            (treeNode                     ), intent(inout)              :: node1                       , node2
    type            (history                      ), intent(inout)              :: starFormationHistory1       , starFormationHistory2
    class           (nodeComponentBasic           )               , pointer     :: basic1                      , basic2
    double precision                               , dimension(:) , allocatable :: timesNodeCrossing1          , timesNodeCrossing2
    double precision                               , parameter                  :: toleranceRelative    =1.0d-6
    type            (varying_string               )                             :: message
    character       (len=12                       )                             :: label
    integer                                                                     :: i

    ! Call the recursive copy if necessary.
    if (self%isRecursive) then
       call self%recursiveSelf%move(node1,node2,starFormationHistory1,starFormationHistory2)
       return
    end if
    ! Extract and validate crossing times.
    basic1             => node1 %basic                    (                    )
    basic2             => node2 %basic                    (                    )
    timesNodeCrossing1 =  basic1%floatRank1MetaPropertyGet(self%timesCrossingID)
    timesNodeCrossing2 =  basic2%floatRank1MetaPropertyGet(self%timesCrossingID)
    ! If there are no histories to move, we can simply return.
    if (size(timesNodeCrossing2) == 0) return
    ! Validate histories.
    if (size(timesNodeCrossing1) == 0) then
       ! If crossing times exist in the moved star formation history, but not in the history to which it is moved, we must
       ! establish them now.
       if (size(timesNodeCrossing2) > 0) call basic1%floatRank1MetaPropertySet(self%timesCrossingID,timesNodeCrossing2)
    else
       ! Check that crossing times agree.
       if (size(timesNodeCrossing1) /= size(timesNodeCrossing2)) then
          message=var_str('number of crossing times differs in merging nodes ')//node1%index()//" ("//size(timesNodeCrossing1)//") and "//node2%index()//" ("//size(timesNodeCrossing2)//")"
          call displayIndent(message)
          write (label,'(e12.6)') basic1%time()
          message='nodes exist at times '//trim(adjustl(label))//' Gyr and '
          write (label,'(e12.6)') basic2%time()
          message=message//trim(adjustl(label))//' Gyr'
          call displayMessage(message)
          call displayMessage(var_str('star formation histories created in progenitors ')//basic1%longIntegerRank0MetaPropertyGet(self%createdInID)//' and '//basic2%longIntegerRank0MetaPropertyGet(self%createdInID))
          call displayMessage("times (target | mergee) are:")
          do i=1,max(size(timesNodeCrossing1),size(timesNodeCrossing2))
             if (i <= size(timesNodeCrossing1)) then
                write (label,'(e12.6)') timesNodeCrossing1(i)
             else
                label=repeat(" ",12)
             end if
             message=         label//" | "
             if (i <= size(timesNodeCrossing2)) then
                write (label,'(e12.6)') timesNodeCrossing2(i)
             else
                label=repeat(" ",12)
             end if
             message=message//label
             call displayMessage(message)
          end do
          call displayUnindent("")
          call Error_Report("invalid crossing times"//{introspection:location})
       else if (.not.all(Values_Agree(timesNodeCrossing1,timesNodeCrossing2,relTol=toleranceRelative))) then
          call displayIndent(var_str('crossing times differ in merging nodes ')//node1%index()//" and "//node2%index())
          write (label,'(e12.6)') basic1%time()
          message='nodes exist at times '//trim(adjustl(label))//' Gyr and '
          write (label,'(e12.6)') basic2%time()
          message=message//trim(adjustl(label))//' Gyr'
          call displayMessage(message)
          call displayMessage(var_str('star formation histories created in progenitors ')//basic1%longIntegerRank0MetaPropertyGet(self%createdInID)//' and '//basic2%longIntegerRank0MetaPropertyGet(self%createdInID))
          call displayMessage("times (target | mergee | difference) are:")
          do i=1,size(timesNodeCrossing1)
             write (label,'(e12.6)') timesNodeCrossing1(i)
             message=         label//" | "
             write (label,'(e12.6)')                           timesNodeCrossing2(i)
             message=message//label//" | "
             write (label,'(e12.6)') abs(timesNodeCrossing1(i)-timesNodeCrossing2(i))
             message=message//label
             call displayMessage(message)
          end do
          call displayUnindent("")
          call Error_Report("invalid crossing times"//{introspection:location})
       end if
    end if
    call starFormationHistory1%increment(starFormationHistory2,autoExtend=.true.)
    call starFormationHistory2%reset    (                                       )
    return
  end subroutine fixedAgesMove
  
  subroutine fixedAgesScales(self,historyStarFormation,node,massStellar,massGas,abundancesStellar)
    !!{
    Set the scalings for error control on the absolute values of star formation histories.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (starFormationHistoryFixedAges), intent(inout)               :: self
    double precision                               , intent(in   )               :: massStellar               , massGas
    type            (abundances                   ), intent(in   )               :: abundancesStellar
    type            (history                      ), intent(inout)               :: historyStarFormation
    type            (treeNode                     ), intent(inout)               :: node
    class           (nodeComponentBasic           ), pointer                     :: basic
    double precision                               , parameter                   :: massMinimum         =1.0d0
    double precision                               , allocatable  , dimension(:) :: timeSteps                 , timesCrossing
    integer         (c_size_t                     )                              :: i                         , j
    !$GLC attributes unused :: abundancesStellar

    ! Call the recursive copy if needed.
    if (self%isRecursive) then
       call self%recursiveSelf%scales(historyStarFormation,node,massStellar,massGas,abundancesStellar)
       return
    end if
    if (.not.historyStarFormation%exists()) return
    ! Get the set of crossing times for this node.
    basic         => node %basic                    (                    )
    timesCrossing =  basic%floatRank1MetaPropertyGet(self%timesCrossingID)
    ! Compute suitable scales for all star formation history bins.
    call historyStarFormation%timeSteps(timeSteps)
    if (size(timesCrossing) > 0) then
       do i=1,size(timesCrossing)
          do j=1,self%countMetallicities+1
             ! The scale is set to a representative stellar mass scale multiplied by the fraction of the total history time in
             ! each time bin.
             historyStarFormation%data(:,(i-1)*(self%countMetallicities+1)+j)=+max(massStellar+massGas,massMinimum) &
                  &                                                           *timeSteps                            &
                  &                                                           /timesCrossing(i)
          end do
       end do
    end if
    deallocate(timeSteps)
    return
  end subroutine fixedAgesScales
  
  function fixedAgesMetallicityBoundaries(self)
    !!{
    Return the boundaries of the metallicities used in this tabulation.
    !!}
    implicit none
    double precision                               , allocatable  , dimension(:) :: fixedAgesMetallicityBoundaries
    class           (starFormationHistoryFixedAges), intent(inout)               :: self

    ! Call the recursive copy if necessary.
    if (self%isRecursive) then
       fixedAgesMetallicityBoundaries=self%recursiveSelf%metallicityBoundaries()
       return
    end if
    ! Return our table of metallicities.
    allocate(fixedAgesMetallicityBoundaries(0:size(self%metallicityTable)-1))
    fixedAgesMetallicityBoundaries(0:size(self%metallicityTable)-1)=self%metallicityTable(1:size(self%metallicityTable))
    return
  end function fixedAgesMetallicityBoundaries

  function fixedAgesTimes(self,node,indexOutput,starFormationHistory,allowTruncation,timeStart) result(times)
    !!{
    Return the times used in this tabulation.
    !!}
    use :: Galacticus_Nodes    , only : nodeComponentBasic
    use :: Error               , only : Error_Report
    use :: Numerical_Comparison, only : Values_Agree
    use :: ISO_Varying_String  , only : varying_string    , assignment(=), operator(//)
    use :: String_Handling     , only : operator(//)
    implicit none
    double precision                               , allocatable  , dimension(:) :: times
    class           (starFormationHistoryFixedAges), intent(inout)               :: self
    type            (treeNode                     ), intent(inout), optional     :: node
    integer         (c_size_t                     ), intent(in   ), optional     :: indexOutput
    type            (history                      ), intent(in   ), optional     :: starFormationHistory
    logical                                        , intent(in   ), optional     :: allowTruncation
    double precision                               , intent(  out), optional     :: timeStart
    double precision                               , allocatable  , dimension(:) :: timesCrossing       , times_
    class           (nodeComponentBasic           ), pointer                     :: basic
    type            (varying_string               )                              :: message
    character       (len=16                       )                              :: label
    !![
    <optionalArgument name="allowTruncation" defaultsTo=".false."/>
    !!]
    
    ! Call the recursive copy if necessary.
    if (self%isRecursive) then
       times=self%recursiveSelf%times(node,indexOutput,starFormationHistory,allowTruncation,timeStart)
       return
    end if
    ! Check that the current time matches the next tabulated time.
    basic         => node %basic                    (                    )
    timesCrossing =  basic%floatRank1MetaPropertyGet(self%timesCrossingID)
    if (.not.Values_Agree(basic%time(),timesCrossing(1),relTol=1.0d-6)) then
       write (label,'(e16.10)') basic%time         ( )
       message="time ("//label//") "
       write (label,'(e16.10)')       timesCrossing(1)
       message=message//"does not match expected time ("//label//") for node "//node%index()//" (star formation history created in progenitor "//basic%longIntegerRank0MetaPropertyGet(self%createdInID)//")"
       call Error_Report(message//{introspection:location})
    end if
    ! Set the times for this output. Note that the times stored in the history object are relative to t=0, so we increment them by
    ! the actual crossing time.
    if (allowTruncation_) then
       ! Truncation is allowed - return only the masses for times that are greater than zero. (Since we use a fixed grid of ages
       ! it is possible to have negative times.)
       allocate(times_(self%countAges+1))
       times_=timesCrossing(1)+starFormationHistory%time
       times=pack(times_,times_ > 0.0d0)
    else
       ! Truncation is now allowed - return all times.
       allocate(times(self%countAges+1))
       times=timesCrossing(1)+starFormationHistory%time
    end if
    ! Set the start time.
    if (present(timeStart))                 &
         & timeStart=+     timesCrossing(1) &
         &           -self%ageMaximum
    return
  end function fixedAgesTimes

  double precision function fixedAgesTimeNext(self,node,starFormationHistory) result(timeNext)
    !!{
    Return the next tabulation time boundary across all histories.
    !!}
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    use            :: Galacticus_Nodes, only : nodeComponentBasic
    use            :: Arrays_Search   , only : searchArray
    implicit none
    class           (starFormationHistoryFixedAges), intent(inout)               :: self
    type            (treeNode                     ), intent(inout)               :: node
    type            (history                      ), intent(in   )               :: starFormationHistory
    double precision                               , allocatable  , dimension(:) :: timesCrossing
    class           (nodeComponentBasic           ), pointer                     :: basic
    integer         (c_size_t                     )                              :: i                   , j

    ! Call the recursive copy if necessary.
    if (self%isRecursive) then
       timeNext=self%recursiveSelf%timeNext(node,starFormationHistory)
       return
    end if
    ! Find the next boundary time across all histories.
    basic         => node %basic                    (                    )
    timesCrossing =  basic%floatRank1MetaPropertyGet(self%timesCrossingID)
    timeNext      =  huge(0.0d0)
    do i=1,size(timesCrossing)
       j=searchArray(starFormationHistory%time,basic%time()-timesCrossing(i))
       if     (                                                                 &
            &   basic%time() >= starFormationHistory%time(j+1)+timesCrossing(i) &
            &  .and.                                                            &
            &   j+1          <  size(starFormationHistory%time)                 &
            & ) j=j+1
       timeNext=min(timeNext,starFormationHistory%time(j+1)+timesCrossing(i))
    end do
    return
  end function fixedAgesTimeNext

  function fixedAgesMasses(self,node,starFormationHistory,allowTruncation) result(masses)
    !!{
    Return the times used in this tabulation.
    !!}
    use :: Galacticus_Nodes    , only : nodeComponentBasic
    use :: Error               , only : Error_Report
    use :: Numerical_Comparison, only : Values_Agree
    use :: ISO_Varying_String  , only : varying_string    , assignment(=), operator(//)
    use :: String_Handling     , only : operator(//)
    implicit none
    double precision                               , allocatable  , dimension(:,:) :: masses
    class           (starFormationHistoryFixedAges), intent(inout)                 :: self
    type            (treeNode                     ), intent(inout)                 :: node
    type            (history                      ), intent(in   )                 :: starFormationHistory
    logical                                        , intent(in   ), optional       :: allowTruncation
    double precision                               , allocatable  , dimension(:  ) :: timesCrossing       , times_
    class           (nodeComponentBasic           ), pointer                       :: basic
    type            (varying_string               )                                :: message
    character       (len=16                       )                                :: label
    !![
    <optionalArgument name="allowTruncation" defaultsTo=".false."/>
    !!]
    
    ! Call the recursive copy if necessary.
    if (self%isRecursive) then
       masses=self%recursiveSelf%masses(node,starFormationHistory,allowTruncation)
       return
    end if
    ! Check that the current time matches the next tabulated time.
    basic         => node %basic                    (                    )
    timesCrossing =  basic%floatRank1MetaPropertyGet(self%timesCrossingID)
    if (.not.Values_Agree(basic%time(),timesCrossing(1),relTol=1.0d-6)) then
       write (label,'(e16.10)') basic%time         ( )
       message="time ("//label//") "
       write (label,'(e16.10)')       timesCrossing(1)
       message=message//"does not match expected time ("//label//") for node "//node%index()
       call Error_Report(message//{introspection:location})
    end if
    ! Set the times for this output. Note that the times stored in the history object are relative to t=0, so we increment them by
    ! the actual crossing time.
    if (allowTruncation_) then
       ! Truncation is allowed - return only the times that are greater than zero. (Since we use a fixed grid of ages it is
       ! possible to have negative times.)
       allocate(times_(self%countAges+1))
       times_=timesCrossing(1)+starFormationHistory%time
       masses=reshape(pack(starFormationHistory%data(:,1:self%countMetallicities+1),spread(times_ > 0.0d0,dim=2,nCopies=self%countMetallicities+1)),[count(times_ > 0.0d0,kind=c_size_t),self%countMetallicities+1_c_size_t])
    else
       ! Truncation is not allowed - return all masses.
       allocate(masses(self%countMetallicities+1,self%countAges+1))
       masses=starFormationHistory%data(:,1:self%countMetallicities+1)
    end if
    return
  end function fixedAgesMasses

  function fixedAgesAgeDistribution(self) result(ageDistribution)
    !!{
    Indicate that star formation history ages are fixed globally.
    !!}
    implicit none
    type (enumerationStarFormationHistoryAgesType)                :: ageDistribution
    class(starFormationHistoryFixedAges          ), intent(inout) :: self

    ageDistribution=starFormationHistoryAgesFixed
    return
  end function fixedAgesAgeDistribution

  subroutine fixedAgesDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters  , only : inputParameters
    use :: ISO_Varying_String, only : assignment(=)  , char, operator(//)
    implicit    none
    class    (starFormationHistoryFixedAges), intent(inout)           :: self
    type     (inputParameters              ), intent(inout)           :: descriptor
    logical                                 , intent(in   ), optional :: includeClass              , includeFileModificationTimes
    character(len=18                       )                          :: parameterLabel
    type     (inputParameters              )                          :: parameters
    integer                                                           :: i
    type     (varying_string               )                          :: metallicityBoundariesLabel

    ! Call the recursive copy if necessary.
    if (self%isRecursive) then
       call self%recursiveSelf%descriptor(descriptor,includeClass,includeFileModificationTimes)
       return
    end if
    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('starFormationHistory','fixedAges')
    parameters=descriptor%subparameters('starFormationHistory')
    write (parameterLabel,'(e17.10)') self%ageMinimum
    call parameters%addParameter('ageMinimum'       ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%ageMaximum
    call parameters%addParameter('ageMaximum'       ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(i17)   ') self%countAges
    call parameters%addParameter('countAges'        ,trim(adjustl(parameterLabel)))
    metallicityBoundariesLabel=""
    do i=1,size(self%metallicityTable)
       write (parameterLabel,'(e17.10)') self%metallicityTable(i)
       metallicityBoundariesLabel=metallicityBoundariesLabel//trim(adjustl(parameterLabel))//" "
    end do
    call parameters%addParameter('metallicityBoundaries',char(metallicityBoundariesLabel))
    call self%geometryLightcone_%descriptor(parameters,includeClass,includeFileModificationTimes)
    return
  end subroutine fixedAgesDescriptor

  subroutine fixedAgesDeepCopyReset(self)
    !!{
    Perform a deep copy reset of the object.
    !!}
    implicit none
    class(starFormationHistoryFixedAges), intent(inout) :: self

    self                           %   copiedSelf => null()
    if (.not.self%isRecursive) self%recursiveSelf => null()
    if (associated(self%geometryLightcone_)) call self%geometryLightcone_%deepCopyReset()
    return
  end subroutine fixedAgesDeepCopyReset
  
  subroutine fixedAgesDeepCopyFinalize(self)
    !!{
    Finalize a deep reset of the object.
    !!}
    implicit none
    class(starFormationHistoryFixedAges), intent(inout) :: self

    if (self%isRecursive) call fixedAgesFindParent(self)
    if (associated(self%geometryLightcone_)) call self%geometryLightcone_%deepCopyFinalize()
    return
  end subroutine fixedAgesDeepCopyFinalize
  
  subroutine fixedAgesDeepCopy(self,destination)
    !!{
    Perform a deep copy of the object.
    !!}
    use :: Error             , only : Error_Report
#ifdef OBJECTDEBUG
    use :: Display           , only : displayMessage            , verbosityLevelSilent
    use :: MPI_Utilities     , only : mpiSelf
    use :: Function_Classes  , only : debugReporting
    use :: ISO_Varying_String, only : operator(//)              , var_str
    use :: String_Handling   , only : operator(//)
#endif
    implicit none
    class(starFormationHistoryFixedAges), intent(inout), target :: self
    class(starFormationHistoryClass    ), intent(inout)         :: destination

    call self%starFormationHistoryClass%deepCopy(destination)
    select type (destination)
    type is (starFormationHistoryFixedAges)
       destination%isRecursive                    =self%isRecursive
       if (self%isRecursive) then
          if (associated(self%recursiveSelf%recursiveSelf)) then
             ! If the parent self's recursiveSelf pointer is set, it indicates that it was deep-copied, and the pointer points to
             ! that copy. In that case we set the parent self of our destination to that copy.
             destination%recursiveSelf  => self%recursiveSelf%recursiveSelf
          else
            ! The parent self does not appear to have been deep-copied yet. Retain the same parent self pointer in our copy, but
             ! indicate that we need to look for the new parent later.
             destination%recursiveSelf  => self%recursiveSelf
             destination%parentDeferred =  .true.
          end if
       else
          ! This is a parent of a recursively-constructed object. Record the location of our copy so that it can be used to set
          ! the parent in deep copies of the child object.
          call fixedAgesDeepCopyAssign(self,destination)
          destination%recursiveSelf             => null()
          destination%ageMinimum                =  self%ageMinimum           
          destination%ageMaximum                =  self%ageMaximum           
          destination%metallicityMaximum        =  self%metallicityMaximum        
          destination%metallicityMinimum        =  self%metallicityMinimum        
          destination%countAges                 =  self%countAges            
          destination%countMetallicities        =  self%countMetallicities        
          destination%metallicityTable          =  self%metallicityTable          
          destination%timesCrossingID           =  self%timesCrossingID          
          destination%createdInID               =  self%createdInID          
          destination%parentDeferred            = .false.
          nullify(destination%geometryLightcone_)
          if (associated(self%geometryLightcone_)) then
             if (associated(self%geometryLightcone_%copiedSelf)) then
                select type(s => self%geometryLightcone_%copiedSelf)
                   class is (geometryLightconeClass)
                   destination%geometryLightcone_ => s
                   class default
                   call Error_Report('copiedSelf has incorrect type'//{introspection:location})
                end select
                call self%geometryLightcone_%copiedSelf%referenceCountIncrement()
             else
                allocate(destination%geometryLightcone_,mold=self%geometryLightcone_)
                call self%geometryLightcone_%deepCopy(destination%geometryLightcone_)
                self%geometryLightcone_%copiedSelf => destination%geometryLightcone_
                call destination%geometryLightcone_%autoHook()
             end if
#ifdef OBJECTDEBUG
             if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): darkmatterprofiledmo_ : [destination] : ')//loc(destination)//' : '//loc(destination%geometryLightcone_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
          end if
          nullify(destination%cosmologyFunctions_)
          if (associated(self%cosmologyFunctions_)) then
             if (associated(self%cosmologyFunctions_%copiedSelf)) then
                select type(s => self%cosmologyFunctions_%copiedSelf)
                   class is (cosmologyFunctionsClass)
                   destination%cosmologyFunctions_ => s
                   class default
                   call Error_Report('copiedSelf has incorrect type'//{introspection:location})
                end select
                call self%cosmologyFunctions_%copiedSelf%referenceCountIncrement()
             else
                allocate(destination%cosmologyFunctions_,mold=self%cosmologyFunctions_)
                call self%cosmologyFunctions_%deepCopy(destination%cosmologyFunctions_)
                self%cosmologyFunctions_%copiedSelf => destination%cosmologyFunctions_
                call destination%cosmologyFunctions_%autoHook()
             end if
#ifdef OBJECTDEBUG
             if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): darkmatterprofiledmo_ : [destination] : ')//loc(destination)//' : '//loc(destination%cosmologyFunctions_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
          end if
       end if
    class default
       call Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine fixedAgesDeepCopy

  subroutine fixedAgesDeepCopyAssign(self,destination)
    !!{
    Perform pointer assignment during a deep copy of the object.
    !!}
    implicit none
    class(starFormationHistoryFixedAges), intent(inout)         :: self
    class(starFormationHistoryClass    ), intent(inout), target :: destination

    select type (destination)
    type is (starFormationHistoryFixedAges)
       self%recursiveSelf => destination
    end select
    return
  end subroutine fixedAgesDeepCopyAssign

  subroutine fixedAgesFindParent(self)
    !!{
    Find the deep-copied parent of a recursive child.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(starFormationHistoryFixedAges), intent(inout) :: self

    if (self%parentDeferred) then
       if (associated(self%recursiveSelf%recursiveSelf)) then
          self%recursiveSelf => self%recursiveSelf%recursiveSelf
       else
         call Error_Report("recursive child's parent was not copied"//{introspection:location})
       end if
       self%parentDeferred=.false.
    end if
    return
  end subroutine fixedAgesFindParent
