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
  Implements a merger tree evolve profiler that collects simple data.
  !!}

  use :: Hashes, only : integerHash
  
  !![
  <mergerTreeEvolveProfiler name="mergerTreeEvolveProfilerSimple">
   <description>
    A merger tree evolve profiler that collects simple data. Each step taken by the ODE evolver is then analyzed. First, a
    record of the size of the time step taken is recorded. Second, the property which is currently limiting the time step size
    (i.e. that which has the largest error over the step as judged using the same heuristics as the ODE solver uses to
    determine step size) is determined and a record of this is kept.
  
    At the end of a run the accumulated data is written to the \glc\ output file, into a group named {\normalfont \ttfamily
    metaData/evolverProfiler}. A histogram of time step sizes is written to {\normalfont \ttfamily timeStepCount} with bins
    specified in {\normalfont \ttfamily timeStep}---these bins can be adjusted using {\normalfont \ttfamily [timeStepMinimum]},
    {\normalfont \ttfamily [timeStepMaximum]} and {\normalfont \ttfamily [timeStepPointsPerDecade]}. A histogram of which
    properties limited step size is written to {\normalfont \ttfamily propertyHitCount} with the associated property names
    written to {\normalfont \ttfamily [propertyNames]}.
   </description>
  </mergerTreeEvolveProfiler>
  !!]
  type, extends(mergerTreeEvolveProfilerClass) :: mergerTreeEvolveProfilerSimple
     !!{
     A merger tree evolve profiler that collects simple data.
     !!}
     private
     double precision                                         :: timeStepMaximum                   , timeStepMinimum           , &
          &                                                      timeStepSmallest
     integer                                                  :: timeStepPointsPerDecade
     double precision             , allocatable, dimension(:) :: timeStep                          , timeCPU                   , &
          &                                                      timeCPUInterrupted
     integer         (c_size_t   ), allocatable, dimension(:) :: timeStepCount                     , evaluationCount           , &
          &                                                      timeStepCountInterrupted          , evaluationCountInterrupted
     type            (integerHash)                            :: propertyHits
     type            (treeNode   ), pointer                   :: node                     => null()
   contains
     final     ::                   simpleDestructor
     procedure :: stepDescriptor => simpleStepDescriptor
     procedure :: profile        => simpleProfile
  end type mergerTreeEvolveProfilerSimple

  interface mergerTreeEvolveProfilerSimple
     !!{
     Constructors for the {\normalfont \ttfamily simple} merger tree evolve profiler class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface mergerTreeEvolveProfilerSimple

contains
  
  function simpleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily simple} merger tree evolve profiler class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (mergerTreeEvolveProfilerSimple)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    double precision                                                :: timeStepMinimum        , timeStepMaximum
    integer                                                         :: timeStepPointsPerDecade

    !![
    <inputParameter>
      <name>timeStepMinimum</name>
      <defaultValue>1.0d-6</defaultValue>
      <description>The smallest timestep to use in profiling ODE solver steps.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>timeStepMaximum</name>
      <defaultValue>1.0d+1</defaultValue>
      <description>The largest timestep to use in profiling ODE solver steps.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>timeStepPointsPerDecade</name>
      <defaultValue>3</defaultValue>
      <description>The number of bins per decade of timestep to use when profiling ODE solver steps.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=mergerTreeEvolveProfilerSimple(timeStepMinimum,timeStepMaximum,timeStepPointsPerDecade)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function simpleConstructorParameters

  subroutine simpleStepDescriptor(self,descriptor)
    !!{
    Set the descriptor for the current step.
    !!}
    implicit none
    class(mergerTreeEvolveProfilerSimple), intent(inout) :: self
    type (varying_string                ), intent(in   ) :: descriptor
    !$GLC attributes unused :: self, descriptor
    
    return
  end subroutine simpleStepDescriptor
  
  function simpleConstructorInternal(timeStepMinimum,timeStepMaximum,timeStepPointsPerDecade) result(self)
    use :: Numerical_Ranges, only : Make_Range, rangeTypeLogarithmic
    implicit none
    type            (mergerTreeEvolveProfilerSimple)                :: self
    double precision                                , intent(in   ) :: timeStepMinimum        , timeStepMaximum
    integer                                         , intent(in   ) :: timeStepPointsPerDecade
    integer                                                         :: timeStepPoints
    !![
    <constructorAssign variables="timeStepMinimum, timeStepMaximum, timeStepPointsPerDecade"/>
    !!]
    
    timeStepPoints=int(log10(timeStepMaximum/timeStepMinimum)*dble(timeStepPointsPerDecade))+1
    allocate(self%timeStep                  (timeStepPoints))
    allocate(self%timeStepCount             (timeStepPoints))
    allocate(self%timeStepCountInterrupted  (timeStepPoints))
    allocate(self%evaluationCount           (timeStepPoints))
    allocate(self%evaluationCountInterrupted(timeStepPoints))
    allocate(self%timeCPU                   (timeStepPoints))
    allocate(self%timeCPUInterrupted        (timeStepPoints))
    self%timeStep                  =Make_Range(timeStepMinimum,timeStepMaximum,timeStepPoints,rangeTypeLogarithmic)
    self%timeStepCount             =0_c_size_t
    self%evaluationCount           =0_c_size_t
    self%timeCPU                   =0.0d0
    self%timeStepCountInterrupted  =0_c_size_t
    self%evaluationCountInterrupted=0_c_size_t
    self%timeCPUInterrupted        =0_c_size_t
    self%timeStepSmallest          =huge(0.0d0)
    call self%propertyHits%initialize()
    return
  end function simpleConstructorInternal

  subroutine simpleDestructor(self)
    !!{
    Output collected meta-data on tree evolution.
    !!}
    use :: Display                         , only : displayIndent , displayUnindent
    use :: Output_HDF5                     , only : outputFile
    use :: HDF5_Access                     , only : hdf5Access
    use :: IO_HDF5                         , only : hdf5Object
    use :: Numerical_Constants_Astronomical, only : gigaYear
    use :: ISO_Varying_String              , only : assignment(=)
    implicit none
    type            (mergerTreeEvolveProfilerSimple), intent(inout)               :: self
    type            (varying_string                ), allocatable  , dimension(:) :: propertyNames
    integer                                         , allocatable  , dimension(:) :: propertyHitCount
    double precision                                , allocatable  , dimension(:) :: timeCPUPrevious
    integer         (c_size_t                      ), allocatable  , dimension(:) :: timeStepCountPrevious, evaluationCountPrevious
    type            (hdf5Object                    )                              :: metaDataDataset      , metaDataGroup          , &
         &                                                                           profilerDataGroup
    integer                                                                       :: i                    , hitCount
    type            (varying_string                )                              :: message
    character       (len=12                        )                              :: label
    
    if (.not.allocated(self%timestepCount     )) return
    if (           all(self%timestepCount == 0)) return
    ! Store meta-data to the output file.
    !$ call hdf5Access%set  ()
    metaDataGroup    =outputFile   %openGroup('metaData'       ,'Galacticus meta data.'     )
    profilerDataGroup=metaDataGroup%openGroup('evolverProfiler','Meta-data on tree evolver.')
    if (profilerDataGroup%hasDataset('timeStepCount'             )) then
       call profilerDataGroup%readDataset('timeStepCount'             ,  timeStepCountPrevious)
       self%  timeStepCount           =self%  timeStepCount           +  timeStepCountPrevious
       deallocate(  timeStepCountPrevious)
    end if
    if (profilerDataGroup%hasDataset('evaluationCount'           )) then
       call profilerDataGroup%readDataset('evaluationCount'           ,evaluationCountPrevious)
       self%evaluationCount           =self%evaluationCount           +evaluationCountPrevious
       deallocate(evaluationCountPrevious)
    end if
    if (profilerDataGroup%hasDataset('timeCPU'                   )) then
       call profilerDataGroup%readDataset('timeCPU'                   ,  timeCPUPrevious      )
       self%  timeCPU                 =self%  timeCPU                 +  timeCPUPrevious
       deallocate(  timeCPUPrevious      )
    end if
    if (profilerDataGroup%hasDataset('timeStepCountInterrupted'  )) then
       call profilerDataGroup%readDataset('timeStepCountInterrupted'  ,  timeStepCountPrevious)
       self%  timeStepCountInterrupted=self%  timeStepCountInterrupted+  timeStepCountPrevious
       deallocate(  timeStepCountPrevious)
    end if
    if (profilerDataGroup%hasDataset('evaluationCountInterrupted')) then
       call profilerDataGroup%readDataset('evaluationCountInterrupted',evaluationCountPrevious)
       self%evaluationCountInterrupted=self%evaluationCountInterrupted+evaluationCountPrevious
       deallocate(evaluationCountPrevious)
    end if
    if (profilerDataGroup%hasDataset('timeCPUInterrupted'        )) then
       call profilerDataGroup%readDataset('timeCPUInterrupted'        ,  timeCPUPrevious      )
       self%  timeCPUInterrupted      =self%        timeCPUInterrupted+        timeCPUPrevious
       deallocate(  timeCPUPrevious      )
    end if
    if (profilerDataGroup%hasDataset('propertyNames')) then
       call profilerDataGroup%readDataset("propertyNames"   ,propertyNames   )
       call profilerDataGroup%readDataset("propertyHitCount",propertyHitCount)
       do i=1,size(propertyNames)
          if (self%propertyHits%exists(propertyNames(i))) then
             hitCount=self%propertyHits%value(propertyNames(i))+propertyHitCount(i)
          else
             hitCount=                                          propertyHitCount(i)
          end if
          call self%propertyHits%set(propertyNames(i),hitCount)
       end do
       deallocate(propertyNames   )
       deallocate(propertyHitCount)
       call profilerDataGroup%remove("propertyNames"   )
       call profilerDataGroup%remove("propertyHitCount")
    end if
    call self%propertyHits%keys  (propertyNames   )
    call self%propertyHits%values(propertyHitCount)
    call profilerDataGroup%writeDataset  (self%timeStep                  ,"timeStep"                  ,"Timestep [Gyr]"                             ,datasetReturned=metaDataDataset)
    call metaDataDataset  %writeAttribute(     gigaYear                  ,"unitsInSI"                                                                                               )
    call metaDataDataset  %close         (                                                                                                                                          )
    call profilerDataGroup%writeDataset  (self%  timeStepCount           ,"timeStepCount"             ,"Timestep histogram []"                                                      )
    call profilerDataGroup%writeDataset  (self%evaluationCount           ,"evaluationCount"           ,"Evaluations at this timestep []"                                            )
    call profilerDataGroup%writeDataset  (self%  timeCPU                 ,"timeCPU"                   ,"CPU time histogram [s]"                                                     )
    call profilerDataGroup%writeDataset  (self%  timeStepCountInterrupted,"timeStepCountInterrupted"  ,"Interrupted timestep histogram []"                                          )
    call profilerDataGroup%writeDataset  (self%evaluationCountInterrupted,"evaluationCountInterrupted","Interrupted evaluations at this timestep []"                                )
    call profilerDataGroup%writeDataset  (self%        timeCPUInterrupted,"timeCPUInterrupted"        ,"Interrupted CPU time histogram [s]"                                         )
    call profilerDataGroup%writeDataset  (     propertyNames             ,"propertyNames"             ,"Property names"                                                             )
    call profilerDataGroup%writeDataset  (     propertyHitCount          ,"propertyHitCount"          ,"Property hit count"                                                         )
    call profilerDataGroup%close         (                                                                                                                                          )
    call metaDataGroup    %close         (                                                                                                                                          )
    !$ call hdf5Access%unset()
    ! Report on the node causing the smallest timestep.
    if (associated(self%node)) then
       write (label,'(e12.6)') self%timeStepSmallest
       message='Snapshot of node causing the smallest timestep of '//trim(adjustl(label))//' Gyr'
       call displayIndent(message)
       call self%node%serializeASCII()
       call self%node%destroy       ()
       deallocate(self%node)
       call displayUnindent('')       
    end if
    return
  end subroutine simpleDestructor

  subroutine simpleProfile(self,node,time,timeStart,timeEnd,timestep,countEvaluations,interrupted,propertyIndex,propertyName,propertyValue,propertyRate,propertyScale,propertyError,timeCPU)
    !!{
    Profile the differential evolution step.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: Arrays_Search, only : searchArray
    implicit none
    class           (mergerTreeEvolveProfilerSimple), intent(inout)               :: self
    type            (treeNode                      ), intent(in   )               :: node
    double precision                                , intent(in   )               :: time            , timeStep     , &
         &                                                                           timeStart       , timeEnd      , &
         &                                                                           timeCPU
    integer         (c_size_t                      ), intent(in   )               :: countEvaluations
    logical                                         , intent(in   )               :: interrupted
    integer         (c_size_t                      ), intent(in   )               :: propertyIndex
    type            (varying_string                ), intent(in   )               :: propertyName
    double precision                                , intent(in   ), dimension(:) :: propertyValue   , propertyScale, &
         &                                                                           propertyError   , propertyRate
    integer                                                                       :: hitCount
    integer         (c_size_t                      )                              :: i
    !$GLC attributes unused :: time, timeStart, timeEnd, propertyIndex, propertyValue, propertyScale, propertyError, propertyRate

    ! Accumulate timestep.
    i                                    =searchArray(self%timeStep,timeStep)
    self   %  timeStepCount           (i)=self%  timeStepCount           (i)+1_c_size_t
    self   %evaluationCount           (i)=self%evaluationCount           (i)+countEvaluations
    self   %  timeCPU                 (i)=self%  timeCPU                 (i)+timeCPU
    if (interrupted) then
       self%  timeStepCountInterrupted(i)=self%  timeStepCountInterrupted(i)+1_c_size_t
       self%evaluationCountInterrupted(i)=self%evaluationCountInterrupted(i)+countEvaluations
       self   %  timeCPUInterrupted   (i)=self%  timeCPUInterrupted      (i)+timeCPU
    end if
    ! Accumulate property hit counts.
    if (self%propertyHits%exists(propertyName)) then
       hitCount=self%propertyHits%value(propertyName)+1
    else
       hitCount=1
    end if
    call self%propertyHits%set(propertyName,hitCount)
    ! Store a copy of the node corresponding to the smallest seen timestep.
    if (timestep < self%timeStepSmallest) then
       if (associated(self%node)) then
          call self%node%destroy()
          deallocate(self%node)
       end if
       allocate(self%node)
       call node%copyNodeTo(self%node)
       self%node%hostTree    => null()
       self%timeStepSmallest =  timestep
    end if
    return
  end subroutine simpleProfile
