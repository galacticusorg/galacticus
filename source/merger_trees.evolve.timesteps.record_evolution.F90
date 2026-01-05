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
Implements a merger tree evolution timestepping class which limits the step to the next epoch at which to record evolution of the
main branch galaxy.
!!}

  use :: Cosmology_Functions    , only : cosmologyFunctions, cosmologyFunctionsClass
  use :: Numerical_Interpolation, only : interpolator
  use :: Output_Times           , only : outputTimes       , outputTimesClass

  !![
  <mergerTreeEvolveTimestep name="mergerTreeEvolveTimestepRecordEvolution">
   <description>
    A merger tree evolution timestepping class which enforces that
    \begin{equation}
     \Delta t \le t_{\mathrm{record},i} - t
    \end{equation}
    where $t$ is the current time, $t_{\mathrm{record},i}$ is the $i^\mathrm{th}$ time at which the evolution of main branch galaxies
    is to be output and $i$ is chosen to be the smallest $i$ such that $t_{\mathrm{record},i} &gt; t$. If there is no $i$ for which
    $t_{\mathrm{record},i} &gt; t$ this criterion is not applied. If this criterion is the limiting criterion for $\Delta t$ then the
    properties of the galaxy will be recorded at the end of the timestep.
  
    Timesteps are logarithmically spaced in cosmic time between {\normalfont \ttfamily [timeBegin]} and \newline {\normalfont
    \ttfamily [timeEnd]}, with the total number of timesteps specified by {\normalfont \ttfamily [countSteps]}.
  
    This recorded evolution will be written to the group {\normalfont \ttfamily mainProgenitorEvolution} in the \glc\ output
    file. Within that group two datasets, {\normalfont \ttfamily time} and {\normalfont \ttfamily expansionFactor}, give the
    times and expansion factors at which evolution was recorded. Then for each merger tree two datasets, {\normalfont \ttfamily
    stellarMass&lt;N&gt;} and {\normalfont \ttfamily totalMass&lt;N&gt;} (where {\normalfont \ttfamily &lt;N&gt;} is the merger tree index), give
    the stellar and total baryonic mass of the main branch progenitor at each timestep.
   </description>
  </mergerTreeEvolveTimestep>
  !!]
  type, extends(mergerTreeEvolveTimestepClass) :: mergerTreeEvolveTimestepRecordEvolution
     !!{
     Implementation of a merger tree evolution timestepping class which limits the step to the next epoch at which to record
     evolution of the main branch galaxy.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctions_    => null()
     class           (outputTimesClass       ), pointer                   :: outputTimes_           => null()
     logical                                                              :: oneTimeDatasetsWritten
     integer                                                              :: countSteps
     double precision                                                     :: timeBegin                       , timeEnd
     double precision                         , allocatable, dimension(:) :: expansionFactor                 , massStellar, &
          &                                                                  time                            , massTotal
     type            (interpolator           )                            :: interpolator_
   contains
     !![
     <methods>
       <method description="Reset the record of galaxy evolution." method="reset" />
     </methods>
     !!]
     final     ::                 recordEvolutionDestructor
     procedure :: timeEvolveTo => recordEvolutionTimeEvolveTo
     procedure :: autoHook     => recordEvolutionAutoHook
     procedure :: reset        => recordEvolutionReset
  end type mergerTreeEvolveTimestepRecordEvolution

  interface mergerTreeEvolveTimestepRecordEvolution
     !!{
     Constructors for the \refClass{mergerTreeEvolveTimestepRecordEvolution} merger tree evolution timestep class.
     !!}
     module procedure recordEvolutionConstructorParameters
     module procedure recordEvolutionConstructorInternal
  end interface mergerTreeEvolveTimestepRecordEvolution

contains

  function recordEvolutionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveTimestepRecordEvolution} merger tree evolution timestep class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeEvolveTimestepRecordEvolution)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    class           (outputTimesClass                       ), pointer       :: outputTimes_
    double precision                                                         :: timeBegin          , timeEnd, &
         &                                                                      ageUniverse
    integer                                                                  :: countSteps

    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="outputTimes"        name="outputTimes_"        source="parameters"/>
    !!]
    ageUniverse=cosmologyFunctions_%cosmicTime(1.0d0)
    !![
    <inputParameter>
      <name>timeBegin</name>
      <defaultValue>0.05d0*ageUniverse</defaultValue>
      <description>The earliest time at which to tabulate the evolution of main branch progenitor galaxies (in Gyr).</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>timeEnd</name>
      <defaultValue>ageUniverse</defaultValue>
      <description>The latest time at which to tabulate the evolution of main branch progenitor galaxies (in Gyr).</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>countSteps</name>
      <defaultValue>100</defaultValue>
      <description>The number of steps (spaced logarithmically in cosmic time) at which to tabulate the evolution of main branch progenitor galaxies.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=mergerTreeEvolveTimestepRecordEvolution(timeBegin,timeEnd,countSteps,cosmologyFunctions_,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="outputTimes_"       />
    !!]
    return
  end function recordEvolutionConstructorParameters

  function recordEvolutionConstructorInternal(timeBegin,timeEnd,countSteps,cosmologyFunctions_,outputTimes_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeEvolveTimestepRecordEvolution} merger tree evolution timestep class.
    !!}
    use, intrinsic :: ISO_C_Binding    , only : c_size_t
    use            :: Numerical_Ranges , only : Make_Range   , rangeTypeLogarithmic
    implicit none
    type            (mergerTreeEvolveTimestepRecordEvolution)                        :: self
    class           (cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    class           (outputTimesClass                       ), intent(in   ), target :: outputTimes_
    double precision                                         , intent(in   )         :: timeBegin          , timeEnd
    integer                                                  , intent(in   )         :: countSteps
    integer         (c_size_t                               )                        :: timeIndex
    !![
    <constructorAssign variables="timeBegin, timeEnd, countSteps, *cosmologyFunctions_, *outputTimes_"/>
    !!]

    allocate(self%time           (self%countSteps))
    allocate(self%expansionFactor(self%countSteps))
    allocate(self%massStellar    (self%countSteps))
    allocate(self%massTotal      (self%countSteps))
    self%time=Make_Range(self%timeBegin,self%timeEnd,self%countSteps,rangeTypeLogarithmic)
    do timeIndex=1,self%countSteps
       self%expansionFactor(timeIndex)=self%cosmologyFunctions_%expansionFactor(self%time(timeIndex))
    end do
    call self%reset()
    self%oneTimeDatasetsWritten=.false.
    self%interpolator_         =interpolator(self%time)
    return
  end function recordEvolutionConstructorInternal

  subroutine recordEvolutionAutoHook(self)
    !!{
    Create a hook to the merger tree extra output event to allow us to write out our data.
    !!}
    use :: Events_Hooks, only : mergerTreeExtraOutputEvent, openMPThreadBindingAtLevel
    implicit none
    class(mergerTreeEvolveTimestepRecordEvolution), intent(inout) :: self

    call mergerTreeExtraOutputEvent%attach(self,recordEvolutionOutput,label='mergerTreeEvolveTimestepRecordEvolution',openMPThreadBinding=openMPThreadBindingAtLevel)
    return
  end subroutine recordEvolutionAutoHook

  subroutine recordEvolutionDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeEvolveTimestepRecordEvolution} merger tree evolution timestep class.
    !!}
    use :: Events_Hooks, only : mergerTreeExtraOutputEvent
    implicit none
    type(mergerTreeEvolveTimestepRecordEvolution), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    <objectDestructor name="self%outputTimes_"       />
    !!]
    if (mergerTreeExtraOutputEvent%isAttached(self,recordEvolutionOutput)) call mergerTreeExtraOutputEvent%detach(self,recordEvolutionOutput)
    return
  end subroutine recordEvolutionDestructor

  double precision function recordEvolutionTimeEvolveTo(self,timeEnd,node,task,taskSelf,report,lockNode,lockType)
    !!{
    Determines the timestep to go to the next tabulation point for galaxy evolution storage.
    !!}
    use            :: Evolve_To_Time_Reports , only : Evolve_To_Time_Report
    use            :: Galacticus_Nodes       , only : nodeComponentBasic   , treeNode
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: ISO_Varying_String     , only : varying_string
    implicit none
    class           (mergerTreeEvolveTimestepRecordEvolution), intent(inout), target            :: self
    double precision                                         , intent(in   )                    :: timeEnd
    type            (treeNode                               ), intent(inout), target            :: node
    procedure       (timestepTask                           ), intent(  out), pointer           :: task
    class           (*                                      ), intent(  out), pointer           :: taskSelf
    logical                                                  , intent(in   )                    :: report
    type            (treeNode                               ), intent(  out), pointer, optional :: lockNode
    type            (varying_string                         ), intent(  out)         , optional :: lockType
    class           (nodeComponentBasic                     )               , pointer           :: basic
    integer         (c_size_t                               )                                   :: indexTime
    double precision                                                                            :: time
    !$GLC attributes unused :: timeEnd

    recordEvolutionTimeEvolveTo =  huge(0.0d0)
    if (present(lockNode)) lockNode => null()
    if (present(lockType)) lockType =  ""
    task                            => null()
    taskSelf                        => null()
    if (node%isOnMainBranch()) then
       basic     => node %basic()
       time      =  basic%time ()
       indexTime =  self%interpolator_%locate(time)
       if (time < self%time(indexTime+1)) then
          recordEvolutionTimeEvolveTo     =  self%time(indexTime+1)
          if (present(lockNode)) lockNode => node
          if (present(lockType)) lockType =  "record evolution"
          task                            => recordEvolutionStore
          taskSelf                        => self
       end if
    end if
    if (report) call Evolve_To_Time_Report("record evolution: ",recordEvolutionTimeEvolveTo)
    return
  end function recordEvolutionTimeEvolveTo

  subroutine recordEvolutionStore(self,tree,node,deadlockStatus)
    !!{
    Store properties of the main progenitor galaxy.
    !!}
    use            :: Galactic_Structure_Options, only : massTypeGalactic     , massTypeStellar
    use            :: Error                     , only : Error_Report
    use            :: Galacticus_Nodes          , only : mergerTree           , nodeComponentBasic, treeNode
    use            :: Mass_Distributions        , only : massDistributionClass
    use, intrinsic :: ISO_C_Binding             , only : c_size_t
    implicit none
    class           (*                            ), intent(inout)          :: self
    type            (mergerTree                   ), intent(in   )          :: tree
    type            (treeNode                     ), intent(inout), pointer :: node
    type            (enumerationDeadlockStatusType), intent(inout)          :: deadlockStatus
    class           (nodeComponentBasic           )               , pointer :: basic
    class           (massDistributionClass        )               , pointer :: massDistributionGalactic, massDistributionStellar
    integer         (c_size_t                     )                         :: indexTime
    double precision                                                        :: time
    !$GLC attributes unused :: deadlockStatus, tree

    select type (self)
    class is (mergerTreeEvolveTimestepRecordEvolution)
       basic => node %basic()
       time  =  basic%time ()
       if (time == self%time(self%countSteps)) then
          indexTime=self%countSteps
       else
          indexTime=self%interpolator_%locate(time)
       end if
       massDistributionStellar             => node                    %massDistribution(massType=massTypeStellar )
       massDistributionGalactic            => node                    %massDistribution(massType=massTypeGalactic)
       self%massStellar        (indexTime) =  massDistributionStellar %massTotal       (                         )
       self%massTotal          (indexTime) =  massDistributionGalactic%massTotal       (                         )
       !![
       <objectDestructor name="massDistributionStellar" />
       <objectDestructor name="massDistributionGalactic"/>
       !!]
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine recordEvolutionStore

  subroutine recordEvolutionOutput(self,node,iOutput,treeIndex,nodePassesFilter,treeLock)
    !!{
    Store main branch evolution to the output file.
    !!}
    use            :: Error                           , only : Error_Report
    use            :: Output_HDF5                     , only : outputFile
    use            :: HDF5_Access                     , only : hdf5Access
    use            :: IO_HDF5                         , only : hdf5Object
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: ISO_Varying_String              , only : var_str              , varying_string
    use            :: Kind_Numbers                    , only : kind_int8
    use            :: Numerical_Constants_Astronomical, only : gigaYear             , massSolar
    use            :: String_Handling                 , only : operator(//)
    use            :: Locks                           , only : ompLock
    implicit none
    class  (*             ), intent(inout) :: self
    type   (treeNode      ), intent(inout) :: node
    integer(c_size_t      ), intent(in   ) :: iOutput
    integer(kind=kind_int8), intent(in   ) :: treeIndex
    logical                , intent(in   ) :: nodePassesFilter
    type   (ompLock       ), intent(inout) :: treeLock
    type   (varying_string)                :: datasetName
    type   (hdf5Object    )                :: outputGroup     , dataset
    !$GLC attributes unused :: treeLock

    select type (self)
    class is (mergerTreeEvolveTimestepRecordEvolution)
       if (nodePassesFilter.and.iOutput == self%outputTimes_%count().and.node%isOnMainBranch()) then
          !$ call hdf5Access%set()
          outputGroup=outputFile%openGroup("mainProgenitorEvolution","Evolution data of main progenitors.")
          if (.not.self%oneTimeDatasetsWritten) then
             call outputGroup%writeDataset  (self%time           ,"time"           ,"The time of the main progenitor."            ,datasetReturned=dataset)
             call dataset    %writeAttribute(gigaYear            ,"unitsInSI"                                                                             )
             call dataset    %close         (                                                                                                             )
             call outputGroup%writeDataset  (self%expansionFactor,"expansionFactor","The expansion factor of the main progenitor."                        )
             self%oneTimeDatasetsWritten=.true.
          end if
          datasetName=var_str("stellarMass")//treeIndex
          call outputGroup%writeDataset  (self%massStellar,char(datasetName),"The stellar mass of the main progenitor."           ,datasetReturned=dataset)
          call dataset    %writeAttribute(massSolar       ,"unitsInSI"                                                                                    )
          call dataset    %close         (                                                                                                                )
          datasetName=var_str("totalMass"  )//treeIndex
          call outputGroup%writeDataset  (self%massTotal  ,char(datasetName),"The total baryonic mass of the main progenitor."    ,datasetReturned=dataset)
          call dataset    %writeAttribute(massSolar       ,"unitsInSI"                                                                                    )
          call dataset    %close         (                                                                                                                )
          call outputGroup%close         (                                                                                                                )
          !$ call hdf5Access%unset()
          call    self      %reset()
       end if
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine recordEvolutionOutput

  subroutine recordEvolutionReset(self)
    !!{
    Resets recorded datasets to zero.
    !!}
    implicit none
    class(mergerTreeEvolveTimestepRecordEvolution), intent(inout) :: self

    self%massStellar=0.0d0
    self%massTotal  =0.0d0
    return
  end subroutine recordEvolutionReset
