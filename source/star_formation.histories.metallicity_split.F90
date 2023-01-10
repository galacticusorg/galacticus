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
Contains a module which implements a star formation histories class which records star formation split by metallicity.
!!}

  use :: Output_Times, only : outputTimes, outputTimesClass

  !![
  <starFormationHistory name="starFormationHistoryMetallicitySplit">
   <description>
    A star formation histories class which records star formation split by metallicity. The star formation history is tabulated
    on a grid of time and metallicity. The binning in time is chosen such that bins are at most of size {\normalfont \ttfamily
    [timeStep]} between the time at which each galaxy formed and the final output time, and at most of size
    {\normalfont \ttfamily [timeStepFine]} in the period {\normalfont \ttfamily
    [timeFine]} prior to each output time (all times specified in Gyr). The allows fine binning of recent
    star formation just prior to each output. Usually, the metallicity binning is arranged logarithmically in metallicity with
    {\normalfont \ttfamily [countMetallicities]} bins between {\normalfont \ttfamily
    [metallicityMinimum]} and {\normalfont \ttfamily [metallicityMaximum]} (specified
    in Solar units). Note that the metallicity associated with each bin is the minimum metallicity for that bin (the maximum
    being the metallicity value associated with the next bin, except for the final bin which extends to infinite
    metallicity). If {\normalfont \ttfamily [countMetallicities]}$=0$ is set, then the star formation history
    is not split by metallicity (i.e. a single metallicity bin encompassing all metallicities from zero to infinity is
    used). Alternatively, specific metallicity bin boundaries can be set via the {\normalfont \ttfamily
    [metallicityBoundaries]} parameter---a final boundary corresponding to infinity is always added
    automatically.
  
    Star formation histories are output as follows:
    \begin{verbatim}
    HDF5 "galacticus.hdf5" {
    GROUP "starFormationHistories" {
       COMMENT "Star formation history data."
       DATASET "metallicities" {
        COMMENT "Metallicities at which star formation histories are tabulated"
          DATATYPE  H5T_IEEE_F64LE
          DATASPACE  SIMPLE { ( [countMetallicities] ) / ( [countMetallicities] ) }
       }
       GROUP "Output1" {
          COMMENT "Star formation histories for all trees at each out"
          GROUP "mergerTree1" {
             COMMENT "Star formation histories for each tree."
             DATASET "diskSFH&lt;nodeID&gt;" {
             COMMENT "Star formation history stellar masses of the disk "
                DATATYPE  H5T_IEEE_F64LE
                DATASPACE  SIMPLE { }
             }
             DATASET "diskTime&lt;nodeID&gt;" {
             COMMENT "Star formation history times of the disk component"
                DATATYPE  H5T_IEEE_F64LE
                DATASPACE  SIMPLE { }
             }
             DATASET "spheroidSFH&lt;nodeID&gt;" {
             COMMENT "Star formation history stellar masses of the spher"
                DATATYPE  H5T_IEEE_F64LE
                DATASPACE  SIMPLE { }
             }
             DATASET "spheroidTime&lt;nodeID&gt;" {
             COMMENT "Star formation history times of the spheroid compo"
                DATATYPE  H5T_IEEE_F64LE
                DATASPACE  SIMPLE { }
             }
          }
          GROUP "mergerTree2" {
          .
          .
          .
          }
       }
       GROUP "Output1" {
       .
       .
       .
       }
    }
    }
    \end{verbatim}
    where {\normalfont \ttfamily nodeID} is the index of the relevant node. 
   </description>
  </starFormationHistory>
  !!]
  type, extends(starFormationHistoryClass) :: starFormationHistoryMetallicitySplit
     !!{
     A star formation histories class which records star formation split by metallicity.
     !!}
     private
     class           (outputTimesClass), pointer                   :: outputTimes_            => null()
     integer                                                       :: countMetallicities
     double precision                                              :: timeStep                         , timeStepFine          , &
          &                                                           timeFine                         , metallicityMaximum    , &
          &                                                           metallicityMinimum
     double precision                  , allocatable, dimension(:) :: metallicityTable                 , metallicityBoundaries_
     logical                                                       :: metallicityTableWritten
   contains
     !![
     <methods>
       <method description="Make the star formation history." method="make" />
     </methods>
     !!]
     final     ::                          metallicitySplitDestructor
     procedure :: create                => metallicitySplitCreate
     procedure :: rate                  => metallicitySplitRate
     procedure :: output                => metallicitySplitOutput
     procedure :: scales                => metallicitySplitScales
     procedure :: make                  => metallicitySplitMake
     procedure :: metallicityBoundaries => metallicitySplitMetallicityBoundaries
  end type starFormationHistoryMetallicitySplit

  interface starFormationHistoryMetallicitySplit
     !!{
     Constructors for the ``metallicitySplit'' star formation history class.
     !!}
     module procedure metallicitySplitConstructorParameters
     module procedure metallicitySplitConstructorInternal
  end interface starFormationHistoryMetallicitySplit

  ! Type used to store timestep range information.
  type metallicitySplitTimeStepRange
     private
     integer                                                  :: count
     double precision                                         :: timeBegin          , timeEnd
     type            (metallicitySplitTimeStepRange), pointer :: next      => null()
  end type metallicitySplitTimeStepRange

  ! Effective infinite metallicity.
  double precision, parameter :: metallicitySplitMetallicityInfinite=huge(1.0d0)

contains

  function metallicitySplitConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``metallicitySplit'' star formation history class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationHistoryMetallicitySplit)                              :: self
    type            (inputParameters                     ), intent(inout)               :: parameters
    class           (outputTimesClass                    ), pointer                     :: outputTimes_
    double precision                                      , allocatable  , dimension(:) :: metallicityBoundaries
    double precision                                                                    :: timeStep             , timeStepFine      , &
         &                                                                                 timeFine             , metallicityMinimum, &
         &                                                                                 metallicityMaximum
    integer                                                                             :: countMetallicities

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
    !!]
    if (parameters%isPresent('metallicityBoundaries')) then
       self%countMetallicities=parameters%count('metallicityBoundaries')
       allocate(metallicityBoundaries(self%countMetallicities))
       !![
       <inputParameter>
         <name>metallicityBoundaries</name>
         <description>The metallicities corresponding to boundaries between metallicity bins to use when tabulating star formation histories.</description>
         <source>parameters</source>
         <variable>metallicityBoundaries</variable>
         <type>real</type>
         <cardinality>0..*</cardinality>
       </inputParameter>
       !!]
    else
       !![
       <inputParameter>
         <name>countMetallicities</name>
         <defaultValue>10</defaultValue>
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
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    <conditionalCall>
     <call>self=starFormationHistoryMetallicitySplit(outputTimes_,timeStep,timeStepFine,timeFine{conditions})</call>
     <argument name="metallicityBoundaries" value="metallicityBoundaries" condition="     parameters%isPresent('metallicityBoundaries')"/>
     <argument name="countMetallicities"    value="countMetallicities"    condition=".not.parameters%isPresent('metallicityBoundaries')"/>
     <argument name="metallicityMinimum"    value="metallicityMinimum"    condition=".not.parameters%isPresent('metallicityBoundaries')"/>
     <argument name="metallicityMaximum"    value="metallicityMaximum"    condition=".not.parameters%isPresent('metallicityBoundaries')"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function metallicitySplitConstructorParameters

  function metallicitySplitConstructorInternal(outputTimes_,timeStep,timeStepFine,timeFine,metallicityBoundaries,countMetallicities,metallicityMinimum,metallicityMaximum) result(self)
    !!{
    Internal constructor for the ``metallicitySplit'' star formation history class.
    !!}
    use :: Error           , only : Error_Report
    use :: Numerical_Ranges, only : Make_Range  , rangeTypeLogarithmic
    implicit none
    type            (starFormationHistoryMetallicitySplit)                                        :: self
    double precision                                      , intent(in   ), dimension(:), optional :: metallicityBoundaries
    double precision                                      , intent(in   )              , optional :: metallicityMinimum   , metallicityMaximum
    double precision                                      , intent(in   )                         :: timeStep             , timeStepFine, &
         &                                                                                           timeFine
    integer                                               , intent(in   )              , optional :: countMetallicities
    class           (outputTimesClass                    ), intent(in   ), target                 :: outputTimes_
    !![
    <constructorAssign variables="timeStep, timeStepFine, timeFine, metallicityMinimum, metallicityMaximum, countMetallicities, *outputTimes_"/>
    !!]

    if (present(metallicityBoundaries)) then
       if     (                                &
            &   present(countMetallicities   ) &
            &  .or.                            &
            &   present(metallicityMinimum   ) &
            &  .or.                            &
            &   present(metallicityMaximum   ) &
            & ) call Error_Report('specify either a list of metallicity boundaries, or a range, not both'//{introspection:location})
       allocate(self%metallicityTable(size(metallicityBoundaries)+1))
       self%metallicityTable      (1:size(metallicityBoundaries)  )=metallicityBoundaries
       self%metallicityTable      (  size(metallicityBoundaries)+1)=metallicitySplitMetallicityInfinite
       self%metallicityBoundaries_                                 =metallicityBoundaries
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
          if (countMetallicities > 1) self%metallicityTable(1:countMetallicities)=Make_Range(metallicityMinimum,metallicityMaximum,countMetallicities,rangeType=rangeTypeLogarithmic)
          self%metallicityTable(countMetallicities+1)=metallicitySplitMetallicityInfinite
       end select
       self%metallicityBoundaries_=self%metallicityTable(1:countMetallicities)
    end if
    self%metallicityTableWritten=.false.
    return
  end function metallicitySplitConstructorInternal

  subroutine metallicitySplitDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily metallicitySplit} star formation histories class.
    !!}
    implicit none
    type(starFormationHistoryMetallicitySplit), intent(inout) :: self

    !![
    <objectDestructor name="self%outputTimes_"/>
    !!]
    return
  end subroutine metallicitySplitDestructor

  subroutine metallicitySplitCreate(self,node,historyStarFormation,timeBegin)
    !!{
    Create the history required for storing star formation history.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (starFormationHistoryMetallicitySplit), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    type            (history                             ), intent(inout) :: historyStarFormation
    double precision                                      , intent(in   ) :: timeBegin
    class           (nodeComponentBasic                  ), pointer       :: basic
    double precision                                                      :: timeBeginActual     , timeEnd

    basic           => node%basic()
    timeBeginActual =  min(timeBegin,basic%time())
    timeEnd         =  self%outputTimes_%timeNext(basic%time())
    call self%make(historyStarFormation,timeBeginActual,timeEnd)
    return
  end subroutine metallicitySplitCreate

  subroutine metallicitySplitRate(self,node,historyStarFormation,abundancesFuel,rateStarFormation)
    !!{
    Set the rate the star formation history for {\normalfont \ttfamily node}.
    !!}
    use :: Abundances_Structure, only : abundances        , metallicityTypeLinearByMassSolar
    use :: Arrays_Search       , only : searchArray
    use :: Galacticus_Nodes    , only : nodeComponentBasic, treeNode
    use :: Error               , only : Error_Report
    implicit none
    class           (starFormationHistoryMetallicitySplit), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    type            (history                             ), intent(inout) :: historyStarFormation
    type            (abundances                          ), intent(in   ) :: abundancesFuel
    double precision                                      , intent(in   ) :: rateStarFormation
    class           (nodeComponentBasic                  ), pointer       :: basic
    integer                                                               :: historyCount
    integer         (c_size_t                            )                :: iHistory            , iMetallicity
    double precision                                                      :: fuelMetallicity     , timeNode

    ! Check if history exists.
    if (historyStarFormation%exists()) then
       basic        =>      node                %basic()
       timeNode     =       basic               %time ()
       historyCount =  size(historyStarFormation%time   )
       ! Find the point in the table at which to accumulate the star formation rate.
       iHistory=searchArray(historyStarFormation%time,timeNode)+1
       ! Find the metallicity bin to accumulate to.
       fuelMetallicity=abundancesFuel%metallicity(metallicityType=metallicityTypeLinearByMassSolar)
       if (fuelMetallicity < self%metallicityTable(1) .or. self%countMetallicities == 0) then
          iMetallicity=                                                   +1
       else
          iMetallicity=searchArray(self%metallicityTable,fuelMetallicity)+1
       end if
       ! Accumulate to the appropriate time.
       historyStarFormation%data(iHistory,iMetallicity)=rateStarFormation
    else
       ! No history exists - this is acceptable only if the star formation rate is zero.
       if (rateStarFormation > 0.0d0) call Error_Report('non-zero star formation rate, but star formation history is uninitialized'//{introspection:location})
    end if
    return
  end subroutine metallicitySplitRate

  subroutine metallicitySplitOutput(self,node,nodePassesFilter,historyStarFormation,indexOutput,indexTree,componentType,treeLock)
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
    class           (starFormationHistoryMetallicitySplit), intent(inout)         :: self
    type            (treeNode                            ), intent(inout), target :: node
    logical                                               , intent(in   )         :: nodePassesFilter
    type            (history                             ), intent(inout)         :: historyStarFormation
    integer         (c_size_t                            ), intent(in   )         :: indexOutput
    integer         (kind=kind_int8                      ), intent(in   )         :: indexTree
    type            (enumerationComponentTypeType        ), intent(in   )         :: componentType
    type            (ompLock                             ), intent(inout)         :: treeLock
    class           (nodeComponentBasic                  ), pointer               :: basicParent
    type            (treeNode                            ), pointer               :: nodeParent
    double precision                                                              :: timeBegin           , timeEnd
    type            (varying_string                      )                        :: groupName
    type            (hdf5Object                          )                        :: historyGroup        , outputGroup, &
         &                                                                           treeGroup
    type            (history                             )                        :: newHistory
    !$GLC attributes unused :: treeLock

 
    if (.not.historyStarFormation%exists()) return
    if (nodePassesFilter) then
       !$ call hdf5Access%set()
       if (.not.self%metallicityTableWritten) then
          historyGroup=outputFile%openGroup("starFormationHistories","Star formation history data.")
          call historyGroup%writeDataset(self%metallicityTable,"metallicities","Metallicities at which star formation histories are tabulated.")
          call historyGroup%close       (                                                                                                      )
          self%metallicityTableWritten=.true.
       end if
       historyGroup=outputFile          %openGroup("starFormationHistories","Star formation history data."                          )
       groupName=var_str("Output"    )//indexOutput
       outputGroup =historyGroup        %openGroup(char(groupName)         ,"Star formation histories for all trees at each output.")
       groupName=var_str("mergerTree")//indexTree
       treeGroup    =outputGroup        %openGroup(char(groupName)         ,"Star formation histories for each tree."               )
       groupName=enumerationComponentTypeDecode(componentType,includePrefix=.false.)//"Time"//node%index()
       call treeGroup%writeDataset(historyStarFormation%time,char(groupName),"Star formation history times of the "         //char(enumerationComponentTypeDecode(componentType,includePrefix=.false.))//" component.")
       groupName=enumerationComponentTypeDecode(componentType,includePrefix=.false.)//"SFH" //node%index()
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
       basicParent => nodeParent%basic()
       timeEnd=basicParent%time()
    end if
    call self%make(newHistory,timeBegin,timeEnd,historyStarFormation%time)
    newHistory%data(1:size(historyStarFormation%time),:)=historyStarFormation%data(:,:)
    call historyStarFormation%destroy()
    historyStarFormation=newHistory
    call newHistory%destroy()
    return
  end subroutine metallicitySplitOutput

  subroutine metallicitySplitScales(self,historyStarFormation,massStellar,abundancesStellar)
    !!{
    Set the scalings for error control on the absolute values of star formation histories.
    !!}
    implicit none
    class           (starFormationHistoryMetallicitySplit), intent(inout)               :: self
    double precision                                      , intent(in   )               :: massStellar
    type            (abundances                          ), intent(in   )               :: abundancesStellar
    type            (history                             ), intent(inout)               :: historyStarFormation
    double precision                                      , parameter                   :: massStellarMinimum  =1.0d0
    double precision                                      , allocatable  , dimension(:) :: timeSteps
    integer                                                                             :: iMetallicity
    !$GLC attributes unused :: abundancesStellar

    if (.not.historyStarFormation%exists()) return
    call historyStarFormation%timeSteps(timeSteps)
    forall(iMetallicity=1:self%countMetallicities+1)
       historyStarFormation%data(:,iMetallicity)=max(massStellar,massStellarMinimum)/timeSteps
    end forall
    deallocate(timeSteps)
    return
  end subroutine metallicitySplitScales

  subroutine metallicitySplitMake(self,historyStarFormation,timeBegin,timeEnd,timesCurrent)
    !!{
    Create the history required for storing star formation history.
    !!}
    use :: Error           , only : Error_Report
    use :: Numerical_Ranges, only : Make_Range  , rangeTypeLinear
    implicit none
    class           (starFormationHistoryMetallicitySplit), intent(inout)                         :: self
    type            (history                             ), intent(inout)                         :: historyStarFormation
    double precision                                      , intent(in   )                         :: timeBegin           , timeEnd
    double precision                                      , intent(in   ), dimension(:), optional :: timesCurrent
    type            (metallicitySplitTimeStepRange       ), pointer                               :: timeStepFirst       , timeStepNext , &
         &                                                                                           timeStepNow
    integer                                                                                       :: countTimeCoarse     , countTimeFine, &
         &                                                                                           countTime
    logical                                                                                       :: timeStepFirstFound
    double precision                                                                              :: timeCoarseBegin     , timeCoarseEnd, &
         &                                                                                           timeFineBegin       , timeNext     , &
         &                                                                                           timeNow

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
    timeStepFirstFound = .false.
    timeStepNow        => null()
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
          allocate(timeStepNow%next)
          timeStepNow => timeStepNow%next
       else
          allocate(timeStepFirst)
          timeStepNow => timeStepFirst
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
          timeStepNow%count    =  countTimeCoarse
          timeStepNow%timeBegin=  timeCoarseBegin
          timeStepNow%timeEnd  =  timeCoarseEnd
          allocate(timeStepNow%next)
          timeStepNow          => timeStepNow%next
       end if
       if (countTimeFine > 0) then
          timeStepNow%count    =  countTimeFine
          timeStepNow%timeBegin=  timeFineBegin
          timeStepNow%timeEnd  =  timeNext
       end if
       timeStepNow%next => null()
       ! Increment the total number of steps required.
       countTime=countTime+countTimeFine+countTimeCoarse
       ! Increment the time.
       timeNow=timeNext
    end do
    ! Shift the end point for the final step to the overall end time.
    if (timeStepFirstFound) timeStepNow%timeEnd=timeNext
    ! Copy in existing times if necessary.
    if (present(timesCurrent)) then
       countTime=countTime+size(timesCurrent)
       if (timeStepFirstFound) countTime=countTime-1
    end if
    call historyStarFormation%create(self%countMetallicities+1,countTime)
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
       timeStepNow => timeStepFirst
       do while (associated(timeStepNow))
          ! Populate the time array.
          if      (timeStepNow%count == 1) then
             historyStarFormation%time(countTime+1                            )=                                 timeStepNow%timeEnd
          else if (timeStepNow%count >  1) then
             historyStarFormation%time(countTime+1:countTime+timeStepNow%count)=Make_Range(timeStepNow%timeBegin,timeStepNow%timeEnd,timeStepNow%count,rangeTypeLinear)
          end if
          countTime=countTime+timeStepNow%count
          ! Jump to the next time step.
          timeStepNext => timeStepNow%next
          deallocate(timeStepNow)
          timeStepNow => timeStepNext
       end do
    end if
    return
  end subroutine metallicitySplitMake

  function metallicitySplitMetallicityBoundaries(self)
    !!{
    Return the boundaries of the metallicities used in this tabulation.
    !!}
    implicit none
    double precision                                      , allocatable  , dimension(:) :: metallicitySplitMetallicityBoundaries
    class           (starFormationHistoryMetallicitySplit), intent(inout)               :: self

    allocate(metallicitySplitMetallicityBoundaries(0:size(self%metallicityTable)-1))
    metallicitySplitMetallicityBoundaries(0:size(self%metallicityTable)-1)=self%metallicityTable(1:size(self%metallicityTable))
    return
  end function metallicitySplitMetallicityBoundaries
