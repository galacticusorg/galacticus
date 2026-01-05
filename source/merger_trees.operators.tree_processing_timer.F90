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
  Implements a merger tree operator class which records and outputs tree processing time information.
  !!}

  use, intrinsic :: ISO_C_Binding, only : c_size_t
  use            :: Kind_Numbers , only : kind_int8

  !![
  <mergerTreeOperator name="mergerTreeOperatorTreeProcessingTimer">
   <description>
    A merger tree operator class which records and outputs tree processing time information. Tree timing data to be recorded
    and output to the {\normalfont \ttfamily metaData/treeTiming} group. Three datasets are written to this group:
    \begin{description}
     \item[{\normalfont \ttfamily treeMasses}] Gives the base node masses of the recorded trees (in units of $M_\odot$);
     \item[{\normalfont \ttfamily treeConstructTimes}] Gives the time (in seconds) taken to construct each merger tree;
     \item[{\normalfont \ttfamily treeEvolveTimes}] Gives the time (in seconds) taken to evolve each merger tree.
    \end{description}
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorTreeProcessingTimer
     !!{
     A merger tree operator class which records and outputs tree processing time information.
     !!}
     private
     logical                                                :: collectMemoryUsageData=.false.
     double precision                                       :: timePostEvolution             , timePreConstruction, &
          &                                                    timePreEvolution              , mass
     integer         (kind_int8)                            :: treeID
     integer         (c_size_t )                            :: memoryUsagePeak               , countNodes
     real                                                   :: time
     integer                                                :: countTrees
     double precision           , allocatable, dimension(:) :: timesConstruct                , timesEvolve        , &
          &                                                    masses
     integer         (kind_int8), allocatable, dimension(:) :: treeIDs
     integer         (c_size_t ), allocatable, dimension(:) :: memoryUsagesPeak              , countsNodes
   contains
     final     ::                           treeProcessingTimerDestructor
     procedure :: operatePreConstruction => treeProcessingTimerOperatePreConstruction
     procedure :: operatePreEvolution    => treeProcessingTimerOperatePreEvolution
     procedure :: operatePostEvolution   => treeProcessingTimerOperatePostEvolution
     procedure :: finalize               => treeProcessingTimerFinalize
     procedure :: autoHook               => treeProcessingTimerAutoHook
  end type mergerTreeOperatorTreeProcessingTimer

  interface mergerTreeOperatorTreeProcessingTimer
     !!{
     Constructors for the \refClass{mergerTreeOperatorTreeProcessingTimer} merger tree operator class.
     !!}
     module procedure treeProcessingTimerConstructorParameters
     module procedure treeProcessingTimerConstructorInternal
  end interface mergerTreeOperatorTreeProcessingTimer

  ! Array size increment.
  integer, parameter :: countIncrement=100

contains

  function treeProcessingTimerConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeOperatorTreeProcessingTimer} merger tree operator class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (mergerTreeOperatorTreeProcessingTimer)                :: self
    type   (inputParameters                      ), intent(inout) :: parameters
    logical                                                       :: collectMemoryUsageData
    
    !![
    <inputParameter>
      <name>collectMemoryUsageData</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not to collect and output data on the memory used while processing trees.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=treeProcessingTimerConstructorInternal(collectMemoryUsageData)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function treeProcessingTimerConstructorParameters

  function treeProcessingTimerConstructorInternal(collectMemoryUsageData) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeOperatorTreeProcessingTimer} merger tree operator class.
    !!}
    implicit none
    type (mergerTreeOperatorTreeProcessingTimer)                :: self
    logical                                     , intent(in   ) :: collectMemoryUsageData
    !![
    <constructorAssign variables="collectMemoryUsageData"/>
    !!]

    self%countTrees=0
    return
  end function treeProcessingTimerConstructorInternal

  subroutine treeProcessingTimerAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : postEvolveEvent, openMPThreadBindingAtLevel
    implicit none
    class(mergerTreeOperatorTreeProcessingTimer), intent(inout) :: self

    if (self%collectMemoryUsageData) call postEvolveEvent%attach(self,treeProcessingTimerPostEvolve,openMPThreadBindingAtLevel,label='mergerTreeOperatorTreeProcessingTimer')
    return
  end subroutine treeProcessingTimerAutoHook

  subroutine treeProcessingTimerDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeOperatorTreeProcessingTimer} merger tree operator class.
    !!}
    use :: Events_Hooks, only : postEvolveEvent
    implicit none
    type(mergerTreeOperatorTreeProcessingTimer), intent(inout) :: self

    if (self%collectMemoryUsageData .and. postEvolveEvent%isAttached(self,treeProcessingTimerPostEvolve)) &
         & call postEvolveEvent%detach(self,treeProcessingTimerPostEvolve)
    return
  end subroutine treeProcessingTimerDestructor

  subroutine treeProcessingTimerOperatePreConstruction(self)
    !!{
    Record the CPU time prior to construction of a tree.
    !!}
    !$ use :: OMP_Lib, only : OMP_Get_WTime, OMP_In_Parallel
    implicit none
    class(mergerTreeOperatorTreeProcessingTimer), intent(inout) :: self
    real                                                        :: time
        
    ! Record the CPU time prior to construction.
    !$ if (OMP_In_Parallel()) then
    !$ self%timePreConstruction=OMP_Get_WTime()
    !$ else
    call CPU_Time(time)
    self%timePreConstruction=dble(time)
    !$ end if
    self%timePreEvolution =-1.0
    self%timePostEvolution=-1.0
    self%memoryUsagePeak  = 0_c_size_t
    return
  end subroutine treeProcessingTimerOperatePreConstruction

  subroutine treeProcessingTimerOperatePreEvolution(self,tree)
    !!{
    Record the CPU time prior to evolving {\normalfont \ttfamily tree}.
    !!}
    use    :: Galacticus_Nodes   , only : mergerTree              , nodeComponentBasic, treeNode
    use    :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    !$ use :: OMP_Lib            , only : OMP_Get_WTime           , OMP_In_Parallel
    implicit none
    class(mergerTreeOperatorTreeProcessingTimer), intent(inout), target  :: self
    type (mergerTree                           ), intent(inout), target  :: tree
    type (treeNode                             )               , pointer :: node
    class(nodeComponentBasic                   )               , pointer :: basic
    type (treeNode                             )               , pointer :: nodeWork
    type (mergerTreeWalkerAllNodes             )                         :: treeWalker
    real                                                                 :: time
    
    ! Record the mass and ID of the tree.
    node        => tree %nodeBase
    basic       => node %basic   ()
    self%mass   =  basic%mass    ()
    self%treeID =  tree %index
    ! Count nodes in the tree.
    self%countNodes=0_c_size_t
    treeWalker     =mergerTreeWalkerAllNodes(tree,spanForest=.true.)
    do while (treeWalker%next(nodeWork))
       self%countNodes=self%countNodes+1_c_size_t
    end do
    ! Record the CPU time.
    !$ if (OMP_In_Parallel()) then
    !$ self%timePreEvolution=OMP_Get_WTime()
    !$ else
    call CPU_Time(time)
    self%timePreEvolution=dble(time)
    !$ end if
    return
  end subroutine treeProcessingTimerOperatePreEvolution
  
  subroutine treeProcessingTimerOperatePostEvolution(self)
    !!{
    Record the CPU time after evolving a tree.
    !!}
    !$ use :: OMP_Lib          , only : OMP_Get_WTime, OMP_In_Parallel
    implicit none
    class           (mergerTreeOperatorTreeProcessingTimer), intent(inout)               :: self
    double precision                                       , allocatable  , dimension(:) :: timesConstructTemporary  , timesEvolveTemporary, &
         &                                                                                  massesTemporary
    integer         (kind_int8                            ), allocatable  , dimension(:) :: treeIDsTemporary
    integer         (c_size_t                             ), allocatable  , dimension(:) :: memoryUsagesPeakTemporary, countsNodesTemporary
    real                                                                                 :: time
    
    ! Record the final CPU time.
    !$ if (omp_in_parallel()) then
    !$ self%timePostEvolution=OMP_Get_WTime()
    !$ else
    call CPU_Time(time)
    self%timePostEvolution=dble(time)
    !$ end if
    ! Check that the tree was actually processed. If no pre-evolution time was recorded then no tree existed to be
    ! processed. In that case we ignore the results.
    if (self%timePreEvolution > 0.0d0) then
       ! Ensure that record arrays are sufficiently sized.
       if (.not.allocated(self%masses)) then
          allocate(self%masses          (countIncrement))
          allocate(self%timesConstruct  (countIncrement))
          allocate(self%timesEvolve     (countIncrement))
          allocate(self%treeIDs         (countIncrement))
          allocate(self%memoryUsagesPeak(countIncrement))
          allocate(self%countsNodes     (countIncrement))
       else if (self%countTrees >= size(self%masses)) then
          call Move_Alloc   (self%masses          ,massesTemporary          )
          call Move_Alloc   (self%timesConstruct  ,timesConstructTemporary  )
          call Move_Alloc   (self%timesEvolve     ,timesEvolveTemporary     )
          call Move_Alloc   (self%treeIDs         ,treeIDsTemporary         )
          call Move_Alloc   (self%memoryUsagesPeak,memoryUsagesPeakTemporary)
          call Move_Alloc   (self%countsNodes     ,countsNodesTemporary     )
          allocate(self%masses          (size(massesTemporary)+countIncrement))
          allocate(self%timesConstruct  (size(massesTemporary)+countIncrement))
          allocate(self%timesEvolve     (size(massesTemporary)+countIncrement))
          allocate(self%treeIDs         (size(massesTemporary)+countIncrement))
          allocate(self%memoryUsagesPeak(size(massesTemporary)+countIncrement))
          allocate(self%countsNodes     (size(massesTemporary)+countIncrement))
          self%masses          (1:size(massesTemporary))=massesTemporary
          self%timesConstruct  (1:size(massesTemporary))=timesConstructTemporary
          self%timesEvolve     (1:size(massesTemporary))=timesEvolveTemporary
          self%treeIDs         (1:size(massesTemporary))=treeIDsTemporary
          self%memoryUsagesPeak(1:size(massesTemporary))=memoryUsagesPeakTemporary
          self%countsNodes     (1:size(massesTemporary))=countsNodesTemporary
          deallocate(massesTemporary          )
          deallocate(timesConstructTemporary  )
          deallocate(timesEvolveTemporary     )
          deallocate(treeIDsTemporary         )
          deallocate(memoryUsagesPeakTemporary)
          deallocate(countsNodesTemporary     )
       end if
       ! Store the timing data.
       self%countTrees                       =self%countTrees       +1
       self%masses          (self%countTrees)=self%mass
       self%timesConstruct  (self%countTrees)=self%timePreEvolution -self%timePreConstruction
       self%timesEvolve     (self%countTrees)=self%timePostEvolution-self%timePreEvolution
       self%treeIDs         (self%countTrees)=self%treeID
       self%memoryUsagesPeak(self%countTrees)=self%memoryUsagePeak
       self%countsNodes     (self%countTrees)=self%countNodes
    end if
    return
  end subroutine treeProcessingTimerOperatePostEvolution
  
  subroutine treeProcessingTimerPostEvolve(self,node)
    !!{
    Record memory usage.
    !!}
    use :: Galacticus_Nodes   , only : mergerTree              , treeNode
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    implicit none
    class  (*                       ), intent(inout)          :: self
    type   (treeNode                ), intent(inout), target  :: node
    type   (mergerTree              )               , pointer :: tree       , treeWork
    type   (treeNode                )               , pointer :: nodeWork
    integer(c_size_t                )                         :: memoryUsage
    type   (mergerTreeWalkerAllNodes)                         :: treeWalker

    select type (self)
       class is (mergerTreeOperatorTreeProcessingTimer)
       tree        => node%hostTree%firstTree
       treeWork    => tree
       treeWalker  =  mergerTreeWalkerAllNodes(tree,spanForest=.true.)
       memoryUsage =  0_c_size_t
       do while (associated(treeWork))
          memoryUsage =  memoryUsage+treeWork%sizeof  ()
          treeWork    =>             treeWork%nextTree
       end do
       do while (treeWalker%next(nodeWork))
          memoryUsage=memoryUsage+nodeWork%sizeOf()
       end do
       self%memoryUsagePeak=max(self%memoryUsagePeak,memoryUsage)
    end select
    return
  end subroutine treeProcessingTimerPostEvolve

  subroutine treeProcessingTimerFinalize(self)
    !!{
    Outputs collected meta-data on tree processing times.
    !!}
    use :: Output_HDF5                     , only : outputFile
    use :: HDF5_Access                     , only : hdf5Access
    use :: IO_HDF5                         , only : hdf5Object
    use :: HDF5                            , only : hsize_t
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class  (mergerTreeOperatorTreeProcessingTimer), intent(inout) :: self
    integer(hsize_t                              ), parameter     :: chunkSize      =100_hsize_t
    type   (hdf5Object                           )                :: metaDataDataset            , metaDataGroup, timingDataGroup
    logical                                                       :: preExists

    if (self%countTrees > 0) then
       !$ call hdf5Access%set()
       metaDataGroup  =outputFile     %openGroup ('metaData'  ,'Galacticus meta data.'    )
       timingDataGroup=metaDataGroup  %openGroup ('treeTiming','Meta-data on tree timing.')
       preExists      =timingDataGroup%hasDataset('treeID'                                )
       call timingDataGroup       %writeDataset  (self%treeIDs         (1:self%countTrees),"treeID"         ,"Tree ID"                    ,chunkSize=chunkSize                                ,appendTo=.true.)
       call timingDataGroup       %writeDataset  (self%masses          (1:self%countTrees),"treeMass"       ,"Tree mass [M⊙]"             ,chunkSize=chunkSize,datasetReturned=metaDataDataset,appendTo=.true.)
       if (.not.preExists) &
            & call metaDataDataset%writeAttribute(massSolar                               ,"unitsInSI"                                                                                                        )
       call timingDataGroup       %writeDataset  (self%timesConstruct  (1:self%countTrees),"timeConstruct"  ,"Tree construction time [s]" ,chunkSize=chunkSize,datasetReturned=metaDataDataset,appendTo=.true.)
       if (.not.preExists) &
            & call metaDataDataset%writeAttribute(1.0d0                                   ,"unitsInSI"                                                                                                        )
       call timingDataGroup       %writeDataset  (self%timesEvolve     (1:self%countTrees),"timeEvolve"     ,"Tree evolution time [s]"    ,chunkSize=chunkSize,datasetReturned=metaDataDataset,appendTo=.true.)
       if (.not.preExists) &
            & call metaDataDataset%writeAttribute(1.0d0                                   ,"unitsInSI"                                                                                                        )
       call timingDataGroup       %writeDataset  (self%countsNodes     (1:self%countTrees),"countNodes"     ,"Number of nodes in the tree",chunkSize=chunkSize                                ,appendTo=.true.)
       if (self%collectMemoryUsageData) &
            & call timingDataGroup%writeDataset  (self%memoryUsagesPeak(1:self%countTrees),"memoryUsagePeak","Peak memory usage [bytes]"  ,chunkSize=chunkSize                                ,appendTo=.true.)
       !$ call hdf5Access%unset()
    end if
    return
  end subroutine treeProcessingTimerFinalize
