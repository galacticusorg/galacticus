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

!% Contains a module which implements building of merger trees using the algorithm of \cite{cole_hierarchical_2000}.

module Merger_Tree_Build_Cole2000
  !% Implements building of merger trees using the algorithm of \cite{cole_hierarchical_2000}.
  use FGSL
  implicit none
  private
  public :: Merger_Tree_Build_Cole2000_Initialize, Merger_Tree_Build_Cole2000_Snapshot, Merger_Tree_Build_Cole2000_State_Store,&
       & Merger_Tree_Build_Cole2000_State_Retrieve

  ! Variables controlling merger tree accuracy.
  double precision           :: mergerTreeBuildCole2000AccretionLimit         , mergerTreeBuildCole2000EarliestTime  , &
       &                        mergerTreeBuildCole2000HighestRedshift        , mergerTreeBuildCole2000MassResolution, &
       &                        mergerTreeBuildCole2000MergeProbability

  ! Random number sequence variables
  type            (fgsl_rng) :: clonedPseudoSequenceObject                    , pseudoSequenceObject
  logical                    :: reset                                  =.true., resetSnapshot
  !$omp threadprivate(pseudoSequenceObject,reset,clonedPseudoSequenceObject,resetSnapshot)
  ! Variables used in integrands.
  double precision           :: currentTime
  !$omp threadprivate(currentTime)

contains

  !# <mergerTreeBuildMethod>
  !#  <unitName>Merger_Tree_Build_Cole2000_Initialize</unitName>
  !# </mergerTreeBuildMethod>
  subroutine Merger_Tree_Build_Cole2000_Initialize(mergerTreeBuildMethod,Merger_Tree_Build)
    !% Initializes the \cite{cole_hierarchical_2000} merger tree building module.
    use Input_Parameters
    use ISO_Varying_String
    use Cosmology_Functions
    implicit none
    type     (varying_string               ), intent(in   )          :: mergerTreeBuildMethod
    procedure(Merger_Tree_Build_Do_Cole2000), intent(inout), pointer :: Merger_Tree_Build

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
       !@ <inputParameter>
       !@   <name>mergerTreeBuildCole2000HighestRedshift</name>
       !@   <defaultValue>$10^5$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The highest redshift to which merger trees will be built in the \cite{cole_hierarchical_2000} method.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildCole2000HighestRedshift'  ,mergerTreeBuildCole2000HighestRedshift  ,defaultValue&
            &=1.0d5 )
       mergerTreeBuildCole2000EarliestTime=Cosmology_Age(Expansion_Factor_From_Redshift(mergerTreeBuildCole2000HighestRedshift))
    end if
    return
  end subroutine Merger_Tree_Build_Cole2000_Initialize

  subroutine Merger_Tree_Build_Do_Cole2000(thisTree)
    !% Build a merger tree.
    use Galacticus_Nodes
    use Halo_Mass_Function
    use Critical_Overdensity
    use Merger_Tree_Branching
    use Pseudo_Random
    use Kind_Numbers
    implicit none
    type            (mergerTree        ), intent(inout), target :: thisTree
    type            (treeNode          ), pointer               :: newNode1          , newNode2          , thisNode
    class           (nodeComponentBasic), pointer               :: newBasicComponent1, newBasicComponent2, thisBasicComponent
    integer         (kind=kind_int8    )                        :: nodeIndex
    double precision                                            :: accretionFraction , baseNodeTime      , branchingProbability, &
         &                                                         collapseTime      , deltaCritical     , deltaCritical1      , &
         &                                                         deltaCritical2    , deltaW            , nodeMass1           , &
         &                                                         nodeMass2         , time              , uniformRandom
    logical                                                     :: doBranch

    nodeIndex          =  1                 ! Initialize the node index counter to unity.
    thisNode           => thisTree%baseNode ! Point to the base node.
    thisBasicComponent => thisNode%basic()  ! Get the basic component of the node.

    ! Convert time for base node to critical overdensity (which we use as a time coordinate in this module).
    baseNodeTime=thisBasicComponent%time()
    deltaCritical=Critical_Overdensity_for_Collapse(time=thisBasicComponent%time(),mass=thisBasicComponent%mass())
    call thisBasicComponent%timeSet(deltaCritical)
    ! Begin tree build loop.
    do while (associated(thisNode))
       ! Get the basic component of the node.
       thisBasicComponent => thisNode%basic()

       ! If halo is above the resolution limit, then evolve it.
       if     (                                                                                                                                       &
            &   thisBasicComponent%mass() > mergerTreeBuildCole2000MassResolution                                                                     &
            &  .and.                                                                                                                                  &
            &   Time_of_Collapse(criticalOverdensity=thisBasicComponent%time(),mass=thisBasicComponent%mass()) > mergerTreeBuildCole2000EarliestTime  &
            & ) then

          ! Find branching probability.
          branchingProbability=Tree_Branching_Probability (thisBasicComponent%mass(),thisBasicComponent%time()&
               &,mergerTreeBuildCole2000MassResolution)
          ! Find accretion rate.
          accretionFraction   =Tree_Subresolution_Fraction(thisBasicComponent%mass(),thisBasicComponent%time()&
               &,mergerTreeBuildCole2000MassResolution)

          ! A negative accretion fraction indicates that the node is so close to the resolution limit that
          ! an accretion rate cannot be determined (given available numerical accuracy). In such cases we
          ! consider the node to have reached the end of its resolved evolution and so walk to the next node.
          if (accretionFraction < 0.0d0) then
             call thisNode%walkTreeUnderConstruction(thisNode)
             cycle
          end if

          ! Finding maximum allowed step in w.
          deltaW=min(mergerTreeBuildCole2000AccretionLimit/accretionFraction,Tree_Maximum_Step(thisBasicComponent%mass()&
               &,thisBasicComponent%time(),mergerTreeBuildCole2000MassResolution))
          deltaW=Tree_Maximum_Step(thisBasicComponent%mass(),thisBasicComponent%time(),mergerTreeBuildCole2000MassResolution)
          if (accretionFraction    > 0.0d0) deltaW=min(deltaW,mergerTreeBuildCole2000AccretionLimit  /accretionFraction   )
          if (branchingProbability > 0.0d0) deltaW=min(deltaW,mergerTreeBuildCole2000MergeProbability/branchingProbability)
          ! Scale values to the determined timestep.
          branchingProbability=branchingProbability*deltaW
          accretionFraction   =accretionFraction   *deltaW
          ! Accretion fraction must be less than unity. Reduce timestep (and branching probability and accretion fraction) by
          ! factors of two until this condition is satisfied.
          do while (accretionFraction >= 1.0d0)
             deltaW              =deltaW              *0.5d0
             branchingProbability=branchingProbability*0.5d0
             accretionFraction   =accretionFraction   *0.5d0
          end do

          ! Decide if a branching occurs.
          if (branchingProbability > 0.0d0) then
             uniformRandom=Pseudo_Random_Get(pseudoSequenceObject,reset)
             doBranch=(uniformRandom <= branchingProbability)
          else
             doBranch=.false.
          end if
          ! Determine the critical overdensity for collapse for the new halo(s).
          deltaCritical=thisBasicComponent%time()+deltaW
          ! Create new nodes.
          select case (doBranch)
          case (.true.)
             ! Branching occurs - create two progenitors.
             nodeIndex=nodeIndex+1
             newNode1 => treeNode(nodeIndex,thisTree)
             newBasicComponent1 => newNode1%basic(autoCreate=.true.)
             ! Compute mass of one of the new nodes. First convert the realized probability back to a rate.
             branchingProbability=uniformRandom/deltaW
             nodeMass1=Tree_Branch_Mass(thisBasicComponent%mass(),thisBasicComponent%time(),mergerTreeBuildCole2000MassResolution&
                  &,branchingProbability)

             ! Compute the time corresponding to this branching event.
             time=Time_of_Collapse(criticalOverdensity=deltaCritical,mass=thisBasicComponent%mass())
             ! Set properties of first new node.
             deltaCritical1=Critical_Overdensity_for_Collapse(time=time,mass=nodeMass1)
             call newBasicComponent1%massSet(nodeMass1     )
             call newBasicComponent1%timeSet(deltaCritical1)
             ! Create second progenitor.
             nodeIndex=nodeIndex+1
             newNode2 => treeNode(nodeIndex,thisTree)
             newBasicComponent2 => newNode2%basic(autoCreate=.true.)
             ! Compute mass of second new node.
             nodeMass2=thisBasicComponent%mass()*(1.0d0-accretionFraction)-nodeMass1
             ! Set properties of second new node.
             deltaCritical2=Critical_Overdensity_for_Collapse(time=time,mass=nodeMass2)
             call newBasicComponent2%massSet(nodeMass2     )
             call newBasicComponent2%timeSet(deltaCritical2)
             ! Create links from old to new nodes and vice-versa. (Ensure that child node is the more massive progenitor.)
             if (nodeMass2 > nodeMass1) then
                thisNode%firstChild => newNode2
                newNode2%sibling    => newNode1
             else
                thisNode%firstChild => newNode1
                newNode1%sibling    => newNode2
             end if
             newNode1%parent        => thisNode
             newNode2%parent        => thisNode
          case (.false.)
             ! No branching occurs - create one progenitor.
             nodeIndex=nodeIndex+1
             newNode1 => treeNode(nodeIndex,thisTree)
             newBasicComponent1 => newNode1%basic(autoCreate=.true.)
             ! Compute new mass accounting for sub-resolution accretion.
             nodeMass1=thisBasicComponent%mass()*(1.0d0-accretionFraction)
             ! Compute the time corresponding to this branching event.
             time=Time_of_Collapse(criticalOverdensity=deltaCritical,mass=thisBasicComponent%mass())
             ! Set properties of the new node.
             deltaCritical1=Critical_Overdensity_for_Collapse(time=time,mass=nodeMass1)
             call newBasicComponent1%massSet(nodeMass1     )
             call newBasicComponent1%timeSet(deltaCritical1)
             ! Create links from old to new node and vice-versa.
             thisNode%firstChild => newNode1
             newNode1%parent     => thisNode
          end select
       end if

       ! Walk to the next node.
       ! <gfortan 4.6> explicitly specify the target as thisNode since we can't use the "_Same_Node" tree walking procedures.
       call thisNode%walkTreeUnderConstruction(thisNode)

    end do

    ! Walk the tree and convert w to time.
    thisNode => thisTree%baseNode
    do while (associated(thisNode))
       ! Get the basic component of the node.
       thisBasicComponent => thisNode%basic()
       ! Compute the collapse time.
       collapseTime=Time_of_Collapse(criticalOverdensity=thisBasicComponent%time(),mass=thisBasicComponent%mass())
       call thisBasicComponent%timeSet(collapseTime)
       ! <gfortan 4.6> explicitly specify the target as thisNode since we can't use the "_Same_Node" tree walking procedures.
       call thisNode%walkTree(thisNode)
    end do
    thisBasicComponent => thisTree%baseNode%basic()
    call thisBasicComponent%timeSet(baseNodeTime)

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
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

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
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    read (stateFile) reset
    if (.not.reset) call Pseudo_Random_Retrieve(pseudoSequenceObject,fgslStateFile)
    return
  end subroutine Merger_Tree_Build_Cole2000_State_Retrieve

end module Merger_Tree_Build_Cole2000
