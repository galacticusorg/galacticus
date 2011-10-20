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


!% Contains a module which implements building of merger trees using the algorithm of \cite{cole_hierarchical_2000}.

module Merger_Tree_Build_Cole2000
  !% Implements building of merger trees using the algorithm of \cite{cole_hierarchical_2000}.
  use Merger_Trees
  use FGSL
  implicit none
  private
  public :: Merger_Tree_Build_Cole2000_Initialize, Merger_Tree_Build_Cole2000_Snapshot, Merger_Tree_Build_Cole2000_State_Store,&
       & Merger_Tree_Build_Cole2000_State_Retrieve
  
  ! Variables controlling merger tree accuracy.
  double precision :: mergerTreeBuildCole2000MergeProbability,mergerTreeBuildCole2000AccretionLimit&
       &,mergerTreeBuildCole2000MassResolution

  ! Random number sequence variables
  type(fgsl_rng)   :: pseudoSequenceObject,clonedPseudoSequenceObject
  logical          :: reset=.true.,resetSnapshot
  !$omp threadprivate(pseudoSequenceObject,reset,clonedPseudoSequenceObject,resetSnapshot)

contains

  !# <mergerTreeBuildMethod>
  !#  <unitName>Merger_Tree_Build_Cole2000_Initialize</unitName>
  !# </mergerTreeBuildMethod>
  subroutine Merger_Tree_Build_Cole2000_Initialize(mergerTreeBuildMethod,Merger_Tree_Build)
    !% Initializes the \cite{cole_hierarchical_2000} merger tree building module.
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type(varying_string),          intent(in)    :: mergerTreeBuildMethod
    procedure(),          pointer, intent(inout) :: Merger_Tree_Build

    ! Check if our method is to be used.    
    if (mergerTreeBuildMethod == 'Cole2000') then
       ! Assign pointer to our merger tree building subroutine.
       Merger_Tree_Build => Merger_Tree_Build_Do_Cole2000
       ! Read parameters controlling tree accuracy.
       !@ <inputParameter>
       !@   <name>mergerTreeBuildCole2000MergeProbability</name>
       !@   <defaultValue>0.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The largest probability of branching allowed in a timestep in merger trees built by the \cite{cole_hierarchical_2000} method.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildCole2000MergeProbability',mergerTreeBuildCole2000MergeProbability,defaultValue=1.0d-1)
       !@ <inputParameter>
       !@   <name>mergerTreeBuildCole2000AccretionLimit</name>
       !@   <defaultValue>0.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The largest fractional mass change due to subresolution accretion allowed in a timestep in merger trees built by the
       !@     \cite{cole_hierarchical_2000} method.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildCole2000AccretionLimit'  ,mergerTreeBuildCole2000AccretionLimit  ,defaultValue=1.0d-1)
       !@ <inputParameter>
       !@   <name>mergerTreeBuildCole2000MassResolution</name>
       !@   <defaultValue>$5\times 10^9$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The minimum mass (in units of $M_\odot$) of halos to be resolved in merger trees built using the
       !@     \cite{cole_hierarchical_2000} method.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildCole2000MassResolution'  ,mergerTreeBuildCole2000MassResolution  ,defaultValue&
            &=5.0d9 )
    end if
    return
  end subroutine Merger_Tree_Build_Cole2000_Initialize

  subroutine Merger_Tree_Build_Do_Cole2000(thisTree)
    !% Build a merger tree.
    use Tree_Nodes
    use Critical_Overdensity
    use Merger_Tree_Branching
    use Pseudo_Random
    use Kind_Numbers
    implicit none
    type(mergerTree),       intent(inout) :: thisTree
    type(treeNode),         pointer       :: thisNode,newNode1,newNode2
    integer(kind=kind_int8)               :: nodeIndex
    double precision                      :: branchingProbability,accretionFraction,deltaCritical,collapseTime,uniformRandom &
         &,deltaW ,nodeMass1,nodeMass2,time,deltaCritical1,deltaCritical2,baseNodeTime
    logical                               :: doBranch

    nodeIndex=1                   ! Initialize the node index counter to unity.
    thisNode => thisTree%baseNode ! Point to the base node.
    ! Convert time for base node to critical overdensity (which we use as a time coordinate in this module).
    baseNodeTime=Tree_Node_Time(thisNode)
    deltaCritical=Critical_Overdensity_for_Collapse(time=Tree_Node_Time(thisNode),mass=Tree_Node_Mass(thisNode))
    call Tree_Node_Time_Set(thisNode,deltaCritical)
    ! Begin tree build loop.
    do while (associated(thisNode))
       ! If halo is above the resolution limit, then evolve it.
       if (Tree_Node_Mass(thisNode)>mergerTreeBuildCole2000MassResolution) then
          ! Find branching probability.
          branchingProbability=Tree_Branching_Probability(Tree_Node_Mass(thisNode),Tree_Node_Time(thisNode)&
               &,mergerTreeBuildCole2000MassResolution)
          ! Find accretion rate.
          accretionFraction=Tree_Subresolution_Fraction(Tree_Node_Mass(thisNode),Tree_Node_Time(thisNode)&
               &,mergerTreeBuildCole2000MassResolution)

          ! A negative accretion fraction indicates that the node is so close to the resolution limit that
          ! an accretion rate cannot be determined (given available numerical accuracy). In such cases we
          ! consider the node to have reached the end of its resolved evolution and so walk to the next node.
          if (accretionFraction < 0.0d0) then
             call thisNode%walkTreeConstruction(thisNode)
             cycle
          end if

          ! Finding maximum allowed step in w.
          deltaW=min(mergerTreeBuildCole2000AccretionLimit/accretionFraction,Tree_Maximum_Step(Tree_Node_Mass(thisNode)&
               &,Tree_Node_Time(thisNode),mergerTreeBuildCole2000MassResolution))
          if (branchingProbability > 0.0d0) deltaW=min(deltaW,mergerTreeBuildCole2000MergeProbability/branchingProbability)
          ! Scale values to the determined timestep.
          branchingProbability=branchingProbability*deltaW
          accretionFraction   =accretionFraction   *deltaW
          ! Decide if a branching occurs.
          if (branchingProbability > 0.0d0) then
             uniformRandom=Pseudo_Random_Get(pseudoSequenceObject,reset)
             doBranch=(uniformRandom <= branchingProbability)
          else
             doBranch=.false.
          end if
          ! Determine the critical overdensity for collapse for the new halo(s).
          deltaCritical=Tree_Node_Time(thisNode)+deltaW
          ! Create new nodes.
          select case (doBranch)
          case (.true.)
             ! Branching occurs - create two progenitors.
             nodeIndex=nodeIndex+1
             call thisTree%createNode(newNode1,nodeIndex)
             ! Compute mass of one of the new nodes. First convert the realized probability back to a rate.
             branchingProbability=uniformRandom/deltaW
             nodeMass1=Tree_Branch_Mass(Tree_Node_Mass(thisNode),Tree_Node_Time(thisNode),mergerTreeBuildCole2000MassResolution&
                  &,branchingProbability)
             ! Compute the time corresponding to this branching event.
             time=Time_of_Collapse(criticalOverdensity=deltaCritical,mass=Tree_Node_Mass(thisNode))
             ! Set properties of first new node.
             deltaCritical1=Critical_Overdensity_for_Collapse(time=time,mass=nodeMass1)
             call Tree_Node_Mass_Set(newNode1,nodeMass1     )
             call Tree_Node_Time_Set(newNode1,deltaCritical1)
             ! Create second progenitor.
             nodeIndex=nodeIndex+1
             call thisTree%createNode(newNode2,nodeIndex)
             ! Compute mass of second new node.
             nodeMass2=Tree_Node_Mass(thisNode)*(1.0d0-accretionFraction)-nodeMass1
             ! Set properties of second new node.
             deltaCritical2=Critical_Overdensity_for_Collapse(time=time,mass=nodeMass2)
             call Tree_Node_Mass_Set(newNode2,nodeMass2     )
             call Tree_Node_Time_Set(newNode2,deltaCritical2)
             ! Create links from old to new nodes and vice-versa. (Ensure that child node is the more massive progenitor.)
             if (nodeMass2 > nodeMass1) then
                thisNode%childNode   => newNode2
                newNode2%siblingNode => newNode1
             else
                thisNode%childNode   => newNode1
                newNode1%siblingNode => newNode2
             end if
             newNode1%parentNode  => thisNode
             newNode2%parentNode  => thisNode
          case (.false.)
             ! No branching occurs - create one progenitor.
             nodeIndex=nodeIndex+1
             call thisTree%createNode(newNode1,nodeIndex)
             ! Compute new mass accounting for sub-resolution accretion.
             nodeMass1=Tree_Node_Mass(thisNode)*(1.0d0-accretionFraction)
             ! Compute the time corresponding to this branching event.
             time=Time_of_Collapse(criticalOverdensity=deltaCritical,mass=Tree_Node_Mass(thisNode))
             ! Set properties of the new node.
             deltaCritical1=Critical_Overdensity_for_Collapse(time=time,mass=nodeMass1)
             call Tree_Node_Mass_Set(newNode1,nodeMass1     )
             call Tree_Node_Time_Set(newNode1,deltaCritical1)
             ! Create links from old to new node and vice-versa.
             thisNode%childNode  => newNode1
             newNode1%parentNode => thisNode
          end select
       end if

       ! Walk to the next node.
       ! <gfortan 4.6> explicitly specify the target as thisNode since we can't use the "_Same_Node" tree walking procedures.
       call thisNode%walkTreeConstruction(thisNode)

    end do

    ! Walk the tree and convert w to time.
    thisNode => thisTree%baseNode
    do while (associated(thisNode))
       collapseTime=Time_of_Collapse(criticalOverdensity=Tree_Node_Time(thisNode),mass=Tree_Node_Mass(thisNode))
       call Tree_Node_Time_Set(thisNode,collapseTime)
       ! <gfortan 4.6> explicitly specify the target as thisNode since we can't use the "_Same_Node" tree walking procedures.
       call thisNode%walkTree(thisNode)
    end do
    call Tree_Node_Time_Set(thisTree%baseNode,baseNodeTime)

    return
  end subroutine Merger_Tree_Build_Do_Cole2000

  !# <galacticusStateSnapshotTask>
  !#  <unitName>Merger_Tree_Build_Cole2000_Snapshot</unitName>
  !# </galacticusStateSnapshotTask>
  subroutine Merger_Tree_Build_Cole2000_Snapshot
    !% Store a snapshot of the random number generator internal state.
    implicit none

    if (.not.reset) clonedPseudoSequenceObject=FGSL_Rng_Clone(pseudoSequenceObject)
    resetSnapshot=reset
    return
  end subroutine Merger_Tree_Build_Cole2000_Snapshot
  
  !# <galacticusStateStoreTask>
  !#  <unitName>Merger_Tree_Build_Cole2000_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Merger_Tree_Build_Cole2000_State_Store(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    use Pseudo_Random
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) resetSnapshot
    if (.not.resetSnapshot) call Pseudo_Random_Store(clonedPseudoSequenceObject,fgslStateFile)
    return
  end subroutine Merger_Tree_Build_Cole2000_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Merger_Tree_Build_Cole2000_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Merger_Tree_Build_Cole2000_State_Retrieve(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    use Pseudo_Random
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    read (stateFile) reset
    if (.not.reset) call Pseudo_Random_Retrieve(pseudoSequenceObject,fgslStateFile)
    return
  end subroutine Merger_Tree_Build_Cole2000_State_Retrieve
  
end module Merger_Tree_Build_Cole2000
