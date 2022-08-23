!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which implements a star formation histories class which records \emph{in situ} star formation.
!!}

  use :: Output_Times, only : outputTimes, outputTimesClass

  !![
  <starFormationHistory name="starFormationHistoryInSitu">
   <description>
    A star formation histories class which records \emph{in situ} star formation. The star formation history is tabulated on a
    grid of time and is split between in-situ and accreted star formation. The time grid is the same as (and controlled by the
    same parameters) are for the {\normalfont \ttfamily metallicitySplit} method. Output follows the conventional format, with
    2D star formation history datasets to represent the history as a function of time and origin. The first element in the
    origin dimension records in-situ star formation, while the second element records total star formation.
   </description>
  </starFormationHistory>
  !!]
  type, extends(starFormationHistoryClass) :: starFormationHistoryInSitu
     !!{
     A star formation histories class which records \emph{in situ} star formation.
     !!}
     private
     class           (outputTimesClass), pointer :: outputTimes_ => null()
     double precision                            :: timeStep              , timeStepFine, &
          &                                         timeFine
   contains
     !![
     <methods>
       <method description="Make the star formation history." method="make" />
     </methods>
     !!]
     final     ::                          inSituDestructor
     procedure :: create                => inSituCreate
     procedure :: rate                  => inSituRate
     procedure :: output                => inSituOutput
     procedure :: scales                => inSituScales
     procedure :: make                  => inSituMake
     procedure :: autoHook              => inSituAutoHook
     procedure :: metallicityBoundaries => inSituMetallicityBoundaries
  end type starFormationHistoryInSitu

  interface starFormationHistoryInSitu
     !!{
     Constructors for the ``inSitu'' star formation history class.
     !!}
     module procedure inSituConstructorParameters
     module procedure inSituConstructorInternal
  end interface starFormationHistoryInSitu

  ! Type used to store timestep range information.
  type inSituTimeStepRange
     private
     integer                                        :: count
     double precision                               :: timeBegin, timeEnd
     type            (inSituTimeStepRange), pointer :: next
  end type inSituTimeStepRange

contains

  function inSituConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``inSitu'' star formation history class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationHistoryInSitu)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (outputTimesClass          ), pointer       :: outputTimes_
    double precision                                            :: timeStep    , timeStepFine, &
         &                                                         timeFine

    !![
    <inputParameter>
      <name>timeStep</name>
      <defaultValue>0.1d0</defaultValue>
      <description>The time step to use in tabulations of star formation histories [Gyr].</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>timeStepFine</name>
      <defaultValue>0.01d0</defaultValue>
      <description>The fine time step to use in tabulations of star formation histories [Gyr].</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>timeFine</name>
      <defaultValue>0.1d0</defaultValue>
      <description>The period prior to each output for which the fine time step is used in tabulations of star formation histories [Gyr].</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    self=starFormationHistoryInSitu(timeStep,timeStepFine,timeFine,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function inSituConstructorParameters

  function inSituConstructorInternal(timeStep,timeStepFine,timeFine,outputTimes_) result(self)
    !!{
    Internal constructor for the ``inSitu'' star formation history class.
    !!}
    implicit none
    type            (starFormationHistoryInSitu)                        :: self
    double precision                            , intent(in   )         :: timeStep    , timeStepFine, &
         &                                                                 timeFine
    class           (outputTimesClass          ), intent(in   ), target :: outputTimes_
    !![
    <constructorAssign variables="timeStep, timeStepFine, timeFine, *outputTimes_"/>
    !!]

    return
  end function inSituConstructorInternal

  subroutine inSituAutoHook(self)
    !!{
    Attach to the satellite merging event hook.
    !!}
    use :: Events_Hooks, only : openMPThreadBindingAllLevels, satelliteMergerEvent
    implicit none
    class(starFormationHistoryInSitu), intent(inout) :: self

    call satelliteMergerEvent%attach(self,inSituSatelliteMerger,openMPThreadBindingAllLevels)
    return
  end subroutine inSituAutoHook

  subroutine inSituDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily inSitu} star formation histories class.
    !!}
    implicit none
    type(starFormationHistoryInSitu), intent(inout) :: self

    !![
    <objectDestructor name="self%outputTimes_"/>
    !!]
    return
  end subroutine inSituDestructor

  subroutine inSituCreate(self,node,historyStarFormation,timeBegin)
    !!{
    Create the history required for storing star formation history.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (starFormationHistoryInSitu), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    type            (history                   ), intent(inout) :: historyStarFormation
    double precision                            , intent(in   ) :: timeBegin
    class           (nodeComponentBasic        ), pointer       :: basic
    double precision                                            :: timeBeginActual     , timeEnd

    ! Find the start and end times for this history.
    basic           =>               node %basic()
    timeBeginActual =  min(timeBegin,basic%time ())
    timeEnd         =  self%outputTimes_%timeNext(timeBegin)
    call self%make(historyStarFormation,timeBeginActual,timeEnd)
    return
  end subroutine inSituCreate

  subroutine inSituRate(self,node,historyStarFormation,abundancesFuel,rateStarFormation)
    !!{
    Set the rate the star formation history for {\normalfont \ttfamily node}.
    !!}
    use :: Arrays_Search   , only : searchArray
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (starFormationHistoryInSitu), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    type            (history                   ), intent(inout) :: historyStarFormation
    type            (abundances                ), intent(in   ) :: abundancesFuel
    double precision                            , intent(in   ) :: rateStarFormation
    class           (nodeComponentBasic        ), pointer       :: basic
    integer                                                     :: historyCount
    integer         (c_size_t                  )                :: iHistory
    double precision                                            :: timeNode
    !$GLC attributes unused :: self, historyStarFormation, abundancesFuel

    basic                                 =>              node                %basic()
    timeNode                              =               basic               %time ()
    historyCount                          =          size(historyStarFormation%time            )
    iHistory                              =  searchArray(historyStarFormation%time   ,timeNode)+1
    historyStarFormation%data(iHistory,:) =  rateStarFormation
    return
  end subroutine inSituRate

  subroutine inSituOutput(self,node,nodePassesFilter,historyStarFormation,indexOutput,indexTree,componentType,treeLock)
    !!{
    Output the star formation history for {\normalfont \ttfamily node}.
    !!}
    use :: Output_HDF5               , only : outputFile
    use :: Galacticus_Nodes          , only : mergerTree                    , nodeComponentBasic, treeNode
    use :: Galactic_Structure_Options, only : enumerationComponentTypeDecode
    use :: HDF5_Access               , only : hdf5Access
    use :: IO_HDF5                   , only : hdf5Object
    use :: String_Handling           , only : operator(//)
    implicit none
    class           (starFormationHistoryInSitu  ), intent(inout)         :: self
    type            (treeNode                    ), intent(inout), target :: node
    logical                                       , intent(in   )         :: nodePassesFilter
    type            (history                     ), intent(inout)         :: historyStarFormation
    integer         (c_size_t                    ), intent(in   )         :: indexOutput
    integer         (kind=kind_int8              ), intent(in   )         :: indexTree
    type            (enumerationComponentTypeType), intent(in   )         :: componentType
    type            (ompLock                     ), intent(inout)         :: treeLock
    class           (nodeComponentBasic          ), pointer               :: basicParent
    type            (treeNode                    ), pointer               :: nodeParent
    double precision                                                      :: timeBegin           , timeEnd
    type            (varying_string              )                        :: groupName
    type            (hdf5Object                  )                        :: historyGroup        , outputGroup, &
         &                                                                   treeGroup
    type            (history                     )                        :: newHistory
    !$GLC attributes unused :: treeLock

    if (.not.historyStarFormation%exists()) return
    if (nodePassesFilter) then
       !$ call hdf5Access%set()
       historyGroup=outputFile%openGroup("starFormationHistories","Star formation history data.")
       groupName   ="Output"
       groupName   =groupName//indexOutput
       outputGroup =historyGroup%openGroup(char(groupName),"Star formation histories for all trees at each output.")
       groupName   ="mergerTree"
       groupName   =groupName//indexTree
       treeGroup   =outputGroup %openGroup(char(groupName),"Star formation histories for each tree."               )
       groupName   =enumerationComponentTypeDecode(componentType,includePrefix=.false.)//"Time"
       groupname   =groupName           //node%index()
       call treeGroup%writeDataset(historyStarFormation%time,char(groupName),"Star formation history times of the "         //char(enumerationComponentTypeDecode(componentType,includePrefix=.false.))//" component.")
       groupName   =enumerationComponentTypeDecode(componentType,includePrefix=.false.)//"SFH"
       groupname   =groupName//node%index()
       call treeGroup%writeDataset(historyStarFormation%data,char(groupName),"Star formation history stellar masses of the "//char(enumerationComponentTypeDecode(componentType,includePrefix=.false.))//" component.")
       call treeGroup   %close()
       call outputGroup %close()
       call historyGroup%close()
       !$ call hdf5Access%unset()
    end if
    timeBegin=historyStarFormation%time(1)
    if (indexOutput < self%outputTimes_%count()) then
       timeEnd=self%outputTimes_%time(indexOutput+1)
    else
       nodeParent => node
       do while (associated(nodeParent%parent))
          nodeParent => nodeParent%parent
       end do
       basicParent => nodeParent %basic()
       timeEnd     =  basicParent%time ()
    end if
    call self%make(newHistory,timeBegin,timeEnd,historyStarFormation%time)
    newHistory%data(1:size(historyStarFormation%time),:)=historyStarFormation%data(:,:)
    call historyStarFormation%destroy()
    historyStarFormation=newHistory
    call newHistory%destroy(recordMemory=.false.)
    return
  end subroutine inSituOutput

  subroutine inSituScales(self,historyStarFormation,massStellar,abundancesStellar)
    !!{
    Set the scalings for error control on the absolute values of star formation histories.
    !!}
    use :: Memory_Management, only : deallocateArray
    implicit none
    class           (starFormationHistoryInSitu), intent(inout)               :: self
    double precision                            , intent(in   )               :: massStellar
    type            (abundances                ), intent(in   )               :: abundancesStellar
    type            (history                   ), intent(inout)               :: historyStarFormation
    double precision                            , allocatable  , dimension(:) :: timeSteps
    double precision                            , parameter                   :: massStellarMinimum  =1.0d0
    integer                                                                   :: i
    !$GLC attributes unused :: self, abundancesStellar

    if (.not.historyStarFormation%exists()) return
    call historyStarFormation%timeSteps(timeSteps)
    forall(i=1:2)
       historyStarFormation%data(:,i)=max(massStellar,massStellarMinimum)/timeSteps
    end forall
    call deallocateArray(timeSteps)
    return
  end subroutine inSituScales

  subroutine inSituMake(self,historyStarFormation,timeBegin,timeEnd,timesCurrent)
    !!{
    Create the history required for storing star formation history.
    !!}
    use :: Error           , only : Error_Report
    use :: Numerical_Ranges, only : Make_Range  , rangeTypeLinear
    implicit none
    class           (starFormationHistoryInSitu), intent(inout)                         :: self
    type            (history                   ), intent(inout)                         :: historyStarFormation
    double precision                            , intent(in   )                         :: timeBegin           , timeEnd
    double precision                            , intent(in   ), dimension(:), optional :: timesCurrent
    type            (inSituTimeStepRange       ), pointer                               :: timeStepFirst       , timeStepNext , &
         &                                                                                 timeStepCurrent
    integer                                                                             :: countTimeCoarse     , countTimeFine, &
         &                                                                                 countTime
    logical                                                                             :: timeStepFirstFound
    double precision                                                                    :: timeCoarseBegin     , timeCoarseEnd, &
         &                                                                                 timeFineBegin       , timeNext     , &
         &                                                                                 timeNow

    ! Exit with a null history if it would contain no time.
    if (timeEnd <= timeBegin) then
       call historyStarFormation%destroy()
       return
    end if
    ! If we have a set of times tabulated already, do some sanity checks.
    if (present(timesCurrent)) then
       ! Complain if the beginning time is before the given list of times.
       if (timeBegin < timesCurrent(1                 )) call Error_Report('requested begin time is before currently tabulated times'//{introspection:location})
       ! Complain if the end time is less than the maximum tabulated time.
       if (timeEnd   < timesCurrent(size(timesCurrent))) call Error_Report('requested end time is within currently tabulated times'  //{introspection:location})
    end if

    ! Step through time, creating a set of timesteps as needed.
    if (present(timesCurrent)) then
       timeNow         =  timesCurrent(size(timesCurrent))
    else
       timeNow         =  timeBegin
    end if
    countTime          =  0
    timeStepFirstFound =  .false.
    timeStepCurrent    => null()
    do while (timeNow < timeEnd)
       ! Get the time of the next output
       timeNext=self%outputTimes_%timeNext(timeNow)
       ! Unphysical (negative) value indicates no next output.
       if (timeNext < 0.0d0 .or. timeNext > timeEnd) timeNext=timeEnd
       ! Construct coarse and fine timesteps for this output, recording the parameters of each range.
       ! Determine the number of fine timestep bins required and the time at which we begin using fine timesteps.
       if (self%timeFine > 0.0d0) then
          countTimeFine  =int(min(timeNext-timeNow,self%timeFine)/self%timeStepFine)+1
          timeFineBegin  =timeNext-self%timeStepFine*dble(countTimeFine-1)
          timeCoarseBegin=timeNow      +self%timeStep
          timeCoarseEnd  =timeFineBegin-self%timeStepFine
       else
          countTimeFine  =0
          timeFineBegin  =timeNext
          timeCoarseBegin=timeNow      +self%timeStep
          timeCoarseEnd  =timeNext
       end if
       ! Determine the number of coarse time bins required for this history.
       if (timeCoarseEnd > timeCoarseBegin) then
          countTimeCoarse=max(int((timeCoarseEnd-timeCoarseBegin)/self%timeStep)+1,2)
       else if (countTimeFine == 0) then
          countTimeCoarse=2
          timeCoarseBegin=(timeCoarseEnd-timeNow)/3.0d0+timeNow
       else
          countTimeCoarse=0
       end if
       ! Create the time steps.
       if (timeStepFirstFound) then
          allocate(timeStepCurrent%next)
          timeStepCurrent => timeStepCurrent%next
       else
          allocate(timeStepFirst)
          timeStepCurrent => timeStepFirst
          if (countTimeCoarse > 0) then
             countTimeCoarse=countTimeCoarse+1
             timeCoarseBegin=max(timeCoarseBegin-self%timeStep    ,0.0d0)
          else
             countTimeFine  =countTimeFine  +1
             timeFineBegin  =max(timeFineBegin  -self%timeStepFine,0.0d0)
          end if
          timeStepFirstFound=.true.
       end if
       if (countTimeCoarse > 0) then
          timeStepCurrent%count    =  countTimeCoarse
          timeStepCurrent%timeBegin=  timeCoarseBegin
          timeStepCurrent%timeEnd  =  timeCoarseEnd
          allocate(timeStepCurrent%next)
          timeStepCurrent          => timeStepCurrent%next
       end if
       if (countTimeFine > 0) then
          timeStepCurrent%count    =  countTimeFine
          timeStepCurrent%timeBegin=  timeFineBegin
          timeStepCurrent%timeEnd  =  timeNext
       end if
       timeStepCurrent%next => null()
       ! Increment the total number of steps required.
       countTime=countTime+countTimeFine+countTimeCoarse
       ! Increment the time.
       timeNow=timeNext
    end do
    ! Shift the end point for the final step to the overall end time.
    if (timeStepFirstFound) timeStepCurrent%timeEnd=timeNext
    ! Copy in existing times if necessary.
    if (present(timesCurrent)) then
       countTime=countTime+size(timesCurrent)
       if (timeStepFirstFound) countTime=countTime-1
    end if
    call historyStarFormation%create(2,countTime)
    countTime=0
    if (present(timesCurrent)) then
       if (timeStepFirstFound) then
          historyStarFormation%time(countTime+1:countTime+size(timesCurrent)-1)=timesCurrent(1:size(timesCurrent)-1)
          countTime=size(timesCurrent)-1
       else
          historyStarFormation%time(countTime+1:countTime+size(timesCurrent)  )=timesCurrent(1:size(timesCurrent)  )
          countTime=size(timesCurrent)
       end if
    end if
    ! Create new times if necessary.
    if (timeStepFirstFound) then
       timeStepCurrent => timeStepFirst
       do while (associated(timeStepCurrent))
          ! Populate the time array.
          if      (timeStepCurrent%count == 1) then
             historyStarFormation%time(countTime+1                                )=                                     timeStepCurrent%timeEnd
          else if (timeStepCurrent%count >  1) then
             historyStarFormation%time(countTime+1:countTime+timeStepCurrent%count)=Make_Range(timeStepCurrent%timeBegin,timeStepCurrent%timeEnd,timeStepCurrent%count,rangeTypeLinear)
          end if
          countTime=countTime+timeStepCurrent%count
          ! Jump to the next time step.
          timeStepNext => timeStepCurrent%next
          deallocate(timeStepCurrent)
          timeStepCurrent => timeStepNext
       end do
    end if
    return
  end subroutine inSituMake

  subroutine inSituSatelliteMerger(self,node)
    !!{
    Zero any in-situ star formation history for galaxy about to merge.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid, treeNode
    implicit none
    class(starFormationHistoryInSitu), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node
    class(nodeComponentDisk         ), pointer       :: disk
    class(nodeComponentSpheroid     ), pointer       :: spheroid
    type (history                   )                :: historyStarFormationDisk, historyStarFormationSpheroid

    select type (self)
    class is (starFormationHistoryInSitu)
       disk                                   => node    %disk                ()
       spheroid                               => node    %spheroid            ()
       historyStarFormationDisk               =  disk    %starFormationHistory()
       historyStarFormationSpheroid           =  spheroid%starFormationHistory()
       if (historyStarFormationDisk    %exists()) then
          historyStarFormationDisk    %data(:,1)=0.0d0
          call disk    %starFormationHistorySet(    historyStarFormationDisk)
       end if
       if (historyStarFormationSpheroid%exists()) then
          historyStarFormationSpheroid%data(:,1)=0.0d0
          call spheroid%starFormationHistorySet(historyStarFormationSpheroid)
       end if
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine inSituSatelliteMerger

  function inSituMetallicityBoundaries(self)
    !!{
    Return the boundaries of the metallicities used in this tabulation.
    !!}
    implicit none
    double precision                            , allocatable  , dimension(:) :: inSituMetallicityBoundaries
    class           (starFormationHistoryInSitu), intent(inout)               :: self

    allocate(inSituMetallicityBoundaries(0:1))
    inSituMetallicityBoundaries(0:1)=[0.0d0,huge(0.0d0)]
    return
  end function inSituMetallicityBoundaries
