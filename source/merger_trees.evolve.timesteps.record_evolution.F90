!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements a time-stepping criterion for merger tree evolution which permits evolution of the main
!% branch galaxy to be stored.

module Merger_Tree_Timesteps_Record_Evolution
  !% Implements a time-stepping criterion for merger tree evolution which permits evolution of the main
  !% branch galaxy to be stored.
  use FGSL
  private
  public :: Merger_Tree_Timestep_Record_Evolution, Merger_Tree_Record_Evolution_Output

  ! Variable inidicating if module is initialized and active.
  logical          :: timestepRecordEvolutionInitialized=.false.

  ! Variable indicating if one-time datasets have been written yet.
  logical          :: oneTimeDatasetsWritten=.false.

  ! Variable indicating if evolution should be recorded.
  logical          :: timestepRecordEvolution

  ! Variables which control the distribution of timesteps.
  integer          :: timestepRecordEvolutionSteps
  double precision :: timestepRecordEvolutionBegin,timestepRecordEvolutionEnd

  ! Storage arrays.
  double precision, dimension(:), allocatable :: evolutionTime,evolutionExpansion,evolutionStellarMass,evolutionTotalMass

  ! Interpolation variables.
  type(fgsl_interp_accel)                     :: interpolationAccelerator
  !$omp threadprivate(interpolationAccelerator)

contains

  !# <timeStepsTask>
  !#  <unitName>Merger_Tree_Timestep_Record_Evolution</unitName>
  !# </timeStepsTask>
  subroutine Merger_Tree_Timestep_Record_Evolution(thisNode,timeStep,End_Of_Timestep_Task)
    !% Determines the timestep to go to the next tabulation point for galaxy evolution storage.
    use Tree_Nodes
    use Input_Parameters
    use Cosmology_Functions
    use Memory_Management
    use Numerical_Ranges
    use Numerical_Interpolation
    use Merger_Trees_Evolve_Timesteps_Template
    implicit none
    type(treeNode),                           intent(inout), pointer :: thisNode
    procedure(End_Of_Timestep_Task_Template), intent(inout), pointer :: End_Of_Timestep_Task
    double precision,                         intent(inout)          :: timeStep
    integer                                                          :: timeIndex
    double precision                                                 :: time,ourTimeStep
    
    !$omp critical (timestepRecordEvolutionInitialize)
    if (.not.timestepRecordEvolutionInitialized) then
       ! Get module parameters.
       !@ <inputParameter>
       !@   <name>timestepRecordEvolution</name>
       !@   <defaultValue>false</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not the evolution of the main branch galaxy should be recorded.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('timestepRecordEvolution',timestepRecordEvolution,defaultValue=.false.)
       if (timestepRecordEvolution) then
          ! Get time at present day.
          time=Cosmology_Age(aExpansion=0.999d0)
          ! Get module parameters.
          !@ <inputParameter>
          !@   <name>timestepRecordEvolutionBegin</name>
          !@   <defaultValue>5\% of the age of the Universe</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The earliest time at which to tabulate the evolution of main branch progenitor galaxies (in Gyr).
          !@   </description>
          !@ </inputParameter>
          call Get_Input_Parameter('timestepRecordEvolutionBegin',timestepRecordEvolutionBegin,defaultValue=0.05d0*time)
          !@ <inputParameter>
          !@   <name>timestepRecordEvolutionEnd</name>
          !@   <defaultValue>The age of the Universe</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The latest time at which to tabulate the evolution of main branch progenitor galaxies (in Gyr).
          !@   </description>
          !@ </inputParameter>
          call Get_Input_Parameter('timestepRecordEvolutionEnd'  ,timestepRecordEvolutionEnd  ,defaultValue=       time)
          !@ <inputParameter>
          !@   <name>timestepRecordEvolutionSteps</name>
          !@   <defaultValue>30</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The number of steps (spaced logarithmically in cosmic time) at which to tabulate the evolution of main branch progenitor galaxies.
          !@   </description>
          !@ </inputParameter>
          call Get_Input_Parameter('timestepRecordEvolutionSteps',timestepRecordEvolutionSteps,defaultValue=100        )
          ! Allocate storage arrays.
          call Alloc_Array(evolutionTime       ,[timestepRecordEvolutionSteps])
          call Alloc_Array(evolutionExpansion  ,[timestepRecordEvolutionSteps])
          call Alloc_Array(evolutionStellarMass,[timestepRecordEvolutionSteps])
          call Alloc_Array(evolutionTotalMass  ,[timestepRecordEvolutionSteps])
          ! Initialize arrays.
          evolutionTime=Make_Range(timestepRecordEvolutionBegin,timestepRecordEvolutionEnd,timestepRecordEvolutionSteps,rangeTypeLogarithmic)
          do timeIndex=1,timestepRecordEvolutionSteps
             evolutionExpansion(timeIndex)=Expansion_Factor(evolutionTime(timeIndex))
          end do
          call Reset_Records()
       end if
       timestepRecordEvolutionInitialized=.true.
    end if
    !$omp end critical (timestepRecordEvolutionInitialize)
       
    ! Adjust timestep if applicable.
    if (timestepRecordEvolution.and.thisNode%isOnMainBranch()) then
       ! Get current cosmic time.
       time=Tree_Node_Time(thisNode)
       
       ! Determine how long until next available timestep.
       timeIndex=Interpolate_Locate(timestepRecordEvolutionSteps,evolutionTime,interpolationAccelerator,time)
       if (time < evolutionTime(timeIndex+1)) then
          ! Find next time for storage.
          ourTimeStep=evolutionTime(timeIndex+1)-time
          
          ! Set return value if our timestep is smaller than current one.
          if (ourTimeStep <= timeStep) then
             timeStep=ourTimeStep
             End_Of_Timestep_Task => Merger_Tree_Record_Evolution_Store
          end if
       end if
    end if
    return
  end subroutine Merger_Tree_Timestep_Record_Evolution

  subroutine Merger_Tree_Record_Evolution_Store(thisTree,thisNode)
    !% Store properties of the main progenitor galaxy.
    use Merger_Trees
    use Tree_Nodes
    use Numerical_Interpolation
    use Galactic_Structure_Options
    use Galactic_Structure_Enclosed_Masses
    implicit none
    type(mergerTree), intent(in)             :: thisTree
    type(treeNode),   intent(inout), pointer :: thisNode
    integer                                  :: timeIndex
    double precision                         :: time

    ! Get current cosmic time.
    time=Tree_Node_Time(thisNode)

    ! Determine how long until next available timestep.
    if (time == evolutionTime(timestepRecordEvolutionSteps)) then
       timeIndex=timestepRecordEvolutionSteps
    else
       timeIndex=Interpolate_Locate(timestepRecordEvolutionSteps,evolutionTime,interpolationAccelerator,time)
    end if

    ! Accumulate the properties.
    evolutionStellarMass(timeIndex)=Galactic_Structure_Enclosed_Mass(thisNode,massType=massTypeStellar )
    evolutionTotalMass  (timeIndex)=Galactic_Structure_Enclosed_Mass(thisNode,massType=massTypeGalactic)

    return
  end subroutine Merger_Tree_Record_Evolution_Store

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Merger_Tree_Record_Evolution_Output</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Merger_Tree_Record_Evolution_Output(thisNode,iOutput,treeIndex,nodePassesFilter)
    !% Store Fourier-space halo profiles to the output file.
    use Tree_Nodes
    use IO_HDF5
    use Galacticus_HDF5
    use Galacticus_Output_Times
    use ISO_Varying_String
    use String_Handling
    use Kind_Numbers
    use Numerical_Constants_Astronomical
    implicit none
    type(treeNode),          intent(inout), pointer :: thisNode
    integer,                 intent(in)             :: iOutput
    integer(kind=kind_int8), intent(in)             :: treeIndex
    logical,                 intent(in)             :: nodePassesFilter
    type(varying_string)                            :: datasetName
    type(hdf5Object)                                :: outputGroup,thisDataset
 
    ! If halo model output was requested, output the Fourier-space halo profiles.
    if (nodePassesFilter.and.timestepRecordEvolution.and.iOutput == Galacticus_Output_Time_Count().and.thisNode%isOnMainBranch())&
         & then
       ! Create a group for the profile datasets.
       outputGroup=IO_HDF5_Open_Group(galacticusOutputFile,"mainProgenitorEvolution","Evolution data of main progenitors.")

       ! Write one time datasets if necessary.
       if (.not.oneTimeDatasetsWritten) then
          ! Write the datasets.
          call outputGroup%writeDataset(evolutionTime     ,"time"           ,"The time of the main progenitor."            ,datasetReturned=thisDataset)
          call thisDataset%writeAttribute(gigaYear,"unitsInSI")
          call thisDataset%close()
          call outputGroup%writeDataset(evolutionExpansion,"expansionFactor","The expansion factor of the main progenitor."                            )
          ! Record that these datasets have been written.
          oneTimeDatasetsWritten=.true.
       end if
       
       ! Write datasets to the group.
       datasetName="stellarMass"
       datasetName=datasetName//treeIndex
       call outputGroup%writeDataset(evolutionStellarMass,char(datasetName),"The stellar mass of the main progenitor."       ,datasetReturned=thisDataset)
       call thisDataset%writeAttribute(massSolar,"unitsInSI")
       call thisDataset%close()
       datasetName="totalMass"
       datasetName=datasetName//treeIndex
       call outputGroup%writeDataset(evolutionTotalMass  ,char(datasetName),"The total baryonic mass of the main progenitor.",datasetReturned=thisDataset)
       call thisDataset%writeAttribute(massSolar,"unitsInSI")
       call thisDataset%close()
       ! Close the output group.
       call outputGroup%close()
       ! Reset the recorded masses to zero.
       call Reset_Records()
    end if
    return
  end subroutine Merger_Tree_Record_Evolution_Output
  
  subroutine Reset_Records()
    !% Resets recorded datasets to zero.
    implicit none
    
    ! Reset all recorded values to zero.
    evolutionStellarMass=0.0d0
    evolutionTotalMass  =0.0d0
    return
  end subroutine Reset_Records

end module Merger_Tree_Timesteps_Record_Evolution
