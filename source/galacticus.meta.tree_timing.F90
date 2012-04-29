!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which records and outputs timing data for processing trees.

module Galacticus_Meta_Tree_Timing
  !% Records and outputs timing data for processing trees.
  use Kind_Numbers
  implicit none
  private
  public :: Meta_Tree_Timing_Pre_Construction, Meta_Tree_Timing_Pre_Evolve, Meta_Tree_Timing_Post_Evolve, Meta_Tree_Timing_Output

  ! Flag indicating if the module is initialized.
  logical                                     :: metaTimingDataInitialized=.false.
  
  ! Flag indicating if timing data is to be collected.
  logical                                     :: metaCollectTimingData
  
  ! Record of processing times.
  integer(kind=kind_int8)                     :: treeIndex
  real                                        :: timePreConstruction,timePreEvolution,timePostEvolution
  double precision                            :: treeMass
  !$omp threadprivate(treeIndex,timePreConstruction,timePreEvolution,timePostEvolution,treeMass)

  ! Arrays for storing timing.
  integer,          parameter                 :: treeArrayIncreaseSize=100
  integer                                     :: treesRecordedCount   =  0
  double precision, allocatable, dimension(:) :: treeMasses,treeConstructTimes,treeEvolveTimes

contains

  subroutine Meta_Tree_Timing_Initialize()
    !% Initialize the tree timing meta-data module.
    use Input_Parameters
    implicit none

    ! Check if module is initialized.
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
    use Merger_Trees
    use Tree_Nodes
    implicit none
    type(mergerTree), intent(in) :: thisTree
    type(treeNode),   pointer    :: thisNode

    ! Ensure the module is initialized.
    call Meta_Tree_Timing_Initialize()

    if (metaCollectTimingData) then
       ! Record the CPU time.
       call CPU_Time(timePreEvolution)
       ! Record the mass of the tree.
       thisNode => thisTree%baseNode
       treeMass=Tree_Node_Mass(thisNode)
    end if

    return
  end subroutine Meta_Tree_Timing_Pre_Evolve
  
  !# <mergerTreePostEvolveTask>
  !#   <unitName>Meta_Tree_Timing_Post_Evolve</unitName>
  !# </mergerTreePostEvolveTask>
  subroutine Meta_Tree_Timing_Post_Evolve()
    !% Record the CPU time after evolving a tree.
    use Merger_Trees
    use Memory_Management
    implicit none
    double precision, allocatable, dimension(:) :: treeMassesTemporary,treeConstructTimesTemporary,treeEvolveTimesTemporary

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
    use ISO_Varying_String
    use IO_HDF5
    use Galacticus_HDF5
    use Numerical_Constants_Astronomical
    implicit none
    type(hdf5Object) :: metaDataGroup,timingDataGroup,metaDataDataset

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
