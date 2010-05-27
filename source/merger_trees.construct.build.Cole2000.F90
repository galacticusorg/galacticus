!% Contains a module which implements building of merger trees using the algorithm of \cite{cole_hierarchical_2000}.

module Merger_Tree_Build_Cole2000
  !% Implements building of merger trees using the algorithm of \cite{cole_hierarchical_2000}.
  use Merger_Trees
  use FGSL
  private
  public :: Merger_Tree_Build_Cole2000_Initialize, Merger_Tree_Build_Cole2000_Snapshot, Merger_Tree_Build_Cole2000_State_Store,&
       & Merger_Tree_Build_Cole2000_State_Retrieve
  
  ! Variables controlling merger tree accuracy.
  double precision :: mergerTreeBuildCole2000MergeProbability,mergerTreeBuildCole2000AccretionLimit&
       &,mergerTreeBuildCole2000MassResolution

  ! Random number sequence variables
  type(fgsl_rng)   :: pseudoSequenceObject,clonedPseudoSequenceObject
  logical          :: reset=.true.,resetSnapshot
  !$omp threadprivate(pseudoSequenceObject,reset)

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
    if (mergerTreeBuildMethod.eq.'Cole2000') then
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
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildCole2000MassResolution'  ,mergerTreeBuildCole2000MassResolution  ,defaultValue&
            &=5.0d9 )
    end if
    return
  end subroutine Merger_Tree_Build_Cole2000_Initialize

  subroutine Merger_Tree_Build_Do_Cole2000(thisTree)
    !% Build a merger tree.
    use Tree_Nodes
    use Tree_Node_Methods
    use Critical_Overdensity
    use Merger_Tree_Branching
    use Pseudo_Random
    implicit none
    type(mergerTree), intent(inout) :: thisTree
    type(treeNode),   pointer       :: thisNode,newNode1,newNode2
    integer                         :: nodeIndex
    double precision                :: branchingProbability,accretionFraction,deltaCritical,collapseTime,uniformRandom ,deltaW &
         &,nodeMass1,nodeMass2
    logical                         :: doBranch

    nodeIndex=1                   ! Initialize the node index counter to unity.
    thisNode => thisTree%baseNode ! Point to the base node.
    ! Convert time for base node to critical overdensity (which we use as a time coordinate in this module).
    deltaCritical=Critical_Overdensity_for_Collapse(Tree_Node_Time(thisNode))
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
             ! Set properties of first new node.
             call Tree_Node_Mass_Set(newNode1,nodeMass1)
             call Tree_Node_Time_Set(newNode1,deltaCritical)
             ! Create second progenitor.
             nodeIndex=nodeIndex+1
             call thisTree%createNode(newNode2,nodeIndex)
             ! Compute mass of second new node.
             nodeMass2=Tree_Node_Mass(thisNode)*(1.0d0-accretionFraction)-nodeMass1
             ! Set properties of second new node.
             call Tree_Node_Mass_Set(newNode2,nodeMass2)
             call Tree_Node_Time_Set(newNode2,deltaCritical)
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
             ! Set properties of the new node.
             call Tree_Node_Mass_Set(newNode1,nodeMass1)
             call Tree_Node_Time_Set(newNode1,deltaCritical)
             ! Create links from old to new node and vice-versa.
             thisNode%childNode  => newNode1
             newNode1%parentNode => thisNode
          end select
       end if

       ! Walk to the next node.
       call thisNode%walkTreeConstruction()

    end do

    ! Walk the tree and convert w to time.
    thisNode => thisTree%baseNode
    do while (associated(thisNode))
       collapseTime=Time_of_Collapse(Tree_Node_Time(thisNode))
       call Tree_Node_Time_Set(thisNode,collapseTime)
       call thisNode%walkTree()
    end do

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
