!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which records and outputs timing data for processing trees.

module Galacticus_Meta_Tree_Timing
  !% Records and outputs timing data for processing trees.
  implicit none
  private
  public :: Meta_Tree_Timing_Pre_Construction, Meta_Tree_Timing_Pre_Evolve, Meta_Tree_Timing_Post_Evolve, Meta_Tree_Timing_Output

  ! Flag indicating if the module is initialized.
  logical                                     :: metaTimingDataInitialized=.false.

  ! Flag indicating if timing data is to be collected.
  logical                                     :: metaCollectTimingData

  ! Record of processing times.
  real                                        :: timePostEvolution                , timePreConstruction, &
       &                                         timePreEvolution
  double precision                            :: treeMass
  !$omp threadprivate(timePreConstruction,timePreEvolution,timePostEvolution,treeMass)
  ! Arrays for storing timing.
  integer         , parameter                 :: treeArrayIncreaseSize    =100
  integer                                     :: treesRecordedCount       =0
  double precision, allocatable, dimension(:) :: treeConstructTimes               , treeEvolveTimes    , &
       &                                         treeMasses

contains

  subroutine Meta_Tree_Timing_Initialize()
    !% Initialize the tree timing meta-data module.
    use Input_Parameters
    implicit none

    ! Check if module is initialized.
    if (.not.metaTimingDataInitialized) then
       !$omp critical (Meta_Tree_Timing_Pre_Construct_Initialize)
       if (.not.metaTimingDataInitialized) then
          ! Get parameter specifying if output is required.
          !@ <inputParameter>
          !@   <name>metaCollectTimingData</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not collect and output data on the time spent processing trees.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('metaCollectTimingData',metaCollectTimingData,defaultValue=.false.)
          ! Flag that module is initialized.
          metaTimingDataInitialized=.true.
       end if
       !$omp end critical (Meta_Tree_Timing_Pre_Construct_Initialize)
    end if
    return
  end subroutine Meta_Tree_Timing_Initialize

  !# <mergerTreePreTreeConstructionTask>
  !#   <unitName>Meta_Tree_Timing_Pre_Construction</unitName>
  !# </mergerTreePreTreeConstructionTask>
  subroutine Meta_Tree_Timing_Pre_Construction()
    !% Record the CPU time prior to construction of a tree.
    implicit none

    ! Ensure the module is initialized.
    call Meta_Tree_Timing_Initialize()

    if (metaCollectTimingData) then
       ! Record the CPU time prior to construction.
       call CPU_Time(timePreConstruction)
       timePreEvolution =-1.0
       timePostEvolution=-1.0
    end if

    return
  end subroutine Meta_Tree_Timing_Pre_Construction

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Meta_Tree_Timing_Pre_Evolve</unitName>
  !# </mergerTreePreEvolveTask>
  subroutine Meta_Tree_Timing_Pre_Evolve(thisTree)
    !% Record the CPU time prior to evolving {\tt thisTree}.
    use Galacticus_Nodes
    implicit none
    type (mergerTree        ), intent(in   ) :: thisTree
    type (treeNode          ), pointer       :: thisNode
    class(nodeComponentBasic), pointer       :: thisBasicComponent

    ! Ensure the module is initialized.
    call Meta_Tree_Timing_Initialize()

    if (metaCollectTimingData) then
       ! Record the CPU time.
       call CPU_Time(timePreEvolution)
       ! Record the mass of the tree.
       thisNode           => thisTree          %baseNode
       thisBasicComponent => thisNode          %basic()
       treeMass           =  thisBasicComponent%mass ()
    end if

    return
  end subroutine Meta_Tree_Timing_Pre_Evolve

  !# <mergerTreePostEvolveTask>
  !#   <unitName>Meta_Tree_Timing_Post_Evolve</unitName>
  !# </mergerTreePostEvolveTask>
  subroutine Meta_Tree_Timing_Post_Evolve()
    !% Record the CPU time after evolving a tree.
    use Memory_Management
    implicit none
    double precision, allocatable, dimension(:) :: treeConstructTimesTemporary, treeEvolveTimesTemporary, &
         &                                         treeMassesTemporary

    ! Ensure the module is initialized.
    call Meta_Tree_Timing_Initialize()

    if (metaCollectTimingData) then
       ! Record the final CPU time.
       call CPU_Time(timePostEvolution)
       !$omp critical (Meta_Tree_Timing_Pre_Construct_Record)
       ! Ensure that record arrays are sufficiently sized.
       if (.not.allocated(treeMasses)) then
          call Alloc_Array(treeMasses        ,[                 treeArrayIncreaseSize])
          call Alloc_Array(treeConstructTimes,[                 treeArrayIncreaseSize])
          call Alloc_Array(treeEvolveTimes   ,[                 treeArrayIncreaseSize])
       else if (treesRecordedCount >= size(treeMasses)) then
          call Move_Alloc(treeMasses        ,treeMassesTemporary        )
          call Move_Alloc(treeConstructTimes,treeConstructTimesTemporary)
          call Move_Alloc(treeEvolveTimes   ,treeEvolveTimesTemporary   )
          call Alloc_Array(treeMasses        ,[size(treeMassesTemporary)+treeArrayIncreaseSize])
          call Alloc_Array(treeConstructTimes,[size(treeMassesTemporary)+treeArrayIncreaseSize])
          call Alloc_Array(treeEvolveTimes   ,[size(treeMassesTemporary)+treeArrayIncreaseSize])
          treeMasses        (1:size(treeMassesTemporary))=treeMassesTemporary
          treeConstructTimes(1:size(treeMassesTemporary))=treeConstructTimesTemporary
          treeEvolveTimes   (1:size(treeMassesTemporary))=treeEvolveTimesTemporary
          call Dealloc_Array(treeMassesTemporary        )
          call Dealloc_Array(treeConstructTimesTemporary)
          call Dealloc_Array(treeEvolveTimesTemporary   )
       end if
       ! Store the timing data.
       treesRecordedCount=treesRecordedCount+1
       treeMasses        (treesRecordedCount)=treeMass
       treeConstructTimes(treesRecordedCount)=timePreEvolution -timePreConstruction
       treeEvolveTimes   (treesRecordedCount)=timePostEvolution-timePreEvolution
       !$omp end critical (Meta_Tree_Timing_Pre_Construct_Record)
    end if

    return
  end subroutine Meta_Tree_Timing_Post_Evolve

  !# <hdfPreCloseTask>
  !#  <unitName>Meta_Tree_Timing_Output</unitName>
  !# </hdfPreCloseTask>
  subroutine Meta_Tree_Timing_Output
    !% Outputs collected meta-data on tree evolution.
    use Galacticus_HDF5
    use Numerical_Constants_Astronomical
    implicit none
    type(hdf5Object) :: metaDataDataset, metaDataGroup, timingDataGroup

    ! Ensure the module is initialized.
    call Meta_Tree_Timing_Initialize()

    ! Output tree evolution meta-data if any was collected.
    if (metaCollectTimingData) then

       ! Open output groups.
       metaDataGroup  =galacticusOutputFile%openGroup('metaData'  ,'Galacticus meta data.'    )
       timingDataGroup=metaDataGroup       %openGroup('treeTiming','Meta-data on tree timing.')

       ! Write timing data.
       call timingDataGroup%writeDataset(treeMasses        (1:treesRecordedCount),"treeMasses"        ,"Tree mass [M⊙]"       ,datasetReturned=metaDataDataset)
       call metaDataDataset%writeAttribute(massSolar,"unitsInSI")
       call metaDataDataset%close()
       call timingDataGroup%writeDataset(treeConstructTimes(1:treesRecordedCount),"treeConstructTimes","Tree construction time [s]",datasetReturned=metaDataDataset)
       call metaDataDataset%writeAttribute(1.0d0    ,"unitsInSI")
       call metaDataDataset%close()
       call timingDataGroup%writeDataset(treeEvolveTimes   (1:treesRecordedCount),"treeEvolveTimes"   ,"Tree evolution time [s]"   ,datasetReturned=metaDataDataset)
       call metaDataDataset%writeAttribute(1.0d0    ,"unitsInSI")
       call metaDataDataset%close()

       ! Close output groups.
       call timingDataGroup%close()
       call metaDataGroup  %close()

    end if

    return
  end subroutine Meta_Tree_Timing_Output

end module Galacticus_Meta_Tree_Timing
