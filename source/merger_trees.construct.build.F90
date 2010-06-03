!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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






!% Contains a module which implements building of merger trees after drawing masses at random from a mass function.

module Merger_Tree_Build
  !% Implements building of merger trees after drawing masses at random from a mass function.
  use Merger_Trees
  use ISO_Varying_String
  private
  public :: Merger_Tree_Build_Initialize

  ! Variables giving the mass range and sampling frequency for mass function sampling.
  double precision     :: mergerTreeBuildHaloMassMinimum,mergerTreeBuildHaloMassMaximum,mergerTreeBuildTreesBaseRedshift &
       &,mergerTreeBuildTreesBaseTime,mergerTreeBuildTreesHaloMassExponent
  integer              :: mergerTreeBuildTreesPerDecade,mergerTreeBuildTreesBeginAtTree
  type(varying_string) :: mergerTreeBuildTreesHaloMassDistribution

  ! Array of halo masses to use.
  integer                                     :: treeCount,nextTreeIndex
  double precision, allocatable, dimension(:) :: treeHaloMass,treeWeight

  ! Name of merger tree builder method.
  type(varying_string) :: mergerTreeBuildMethod
  ! Pointer to the subroutine that builds the merger tree.
  procedure(Merger_Tree_Builder_Template), pointer :: Merger_Tree_Builder => null()
  interface Merger_Tree_Builder_Template
     subroutine Merger_Tree_Builder_Template(thisTree)
       import mergerTree
       type(mergerTree), intent(inout) :: thisTree
     end subroutine Merger_Tree_Builder_Template
  end interface
  
contains

  !# <mergerTreeConstructMethod>
  !#  <unitName>Merger_Tree_Build_Initialize</unitName>
  !# </mergerTreeConstructMethod>
  subroutine Merger_Tree_Build_Initialize(mergerTreeConstructMethod,Merger_Tree_Construct)
    !% Initializes the merger tree building module.
    use Input_Parameters
    use Memory_Management
    use Cosmology_Functions
    use FGSL
    use Quasi_Random
    use Halo_Mass_Function
    use Sort
    use Galacticus_Error
    use Numerical_Ranges
    !# <include directive="mergerTreeBuildMethod" type="moduleUse">
    include 'merger_trees.build.modules.inc'
    !# </include>
    implicit none
    type(varying_string),          intent(in)    :: mergerTreeConstructMethod
    procedure(),          pointer, intent(inout) :: Merger_Tree_Construct
    integer                                      :: iTree
    type(fgsl_qrng)                              :: quasiSequenceObject
    logical                                      :: quasiSequenceReset=.true.
    double precision                             :: expansionFactor,massMinimum,massMaximum

    ! Check if our method is to be used.
    if (mergerTreeConstructMethod.eq.'build') then
       ! Assign pointer to our merger tree construction subroutine.
       Merger_Tree_Construct => Merger_Tree_Build_Do
       ! Read parameters for halo mass sampling.
       !@ <inputParameter>
       !@   <name>mergerTreeBuildHaloMassMinimum</name>
       !@   <defaultValue>$10^{10}$</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The minimum mass of merger tree base halos to consider when building merger trees, in units of $M_\odot$.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildHaloMassMinimum'  ,mergerTreeBuildHaloMassMinimum  ,defaultValue=1.0d10)
       !@ <inputParameter>
       !@   <name>mergerTreeBuildHaloMassMaximum</name>
       !@   <defaultValue>$10^{15}$</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum mass of merger tree base halos to consider when building merger trees, in units of $M_\odot$.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildHaloMassMaximum'  ,mergerTreeBuildHaloMassMaximum  ,defaultValue=1.0d15)
       !@ <inputParameter>
       !@   <name>mergerTreeBuildTreesPerDecade</name>
       !@   <defaultValue>10</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The number of merger trees to build per decade of base halo mass.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildTreesPerDecade'   ,mergerTreeBuildTreesPerDecade   ,defaultValue=10    )
       !@ <inputParameter>
       !@   <name>mergerTreeBuildTreesBaseRedshift</name>
       !@   <defaultValue>0</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The redshift at which to plant the base node when building merger trees.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildTreesBaseRedshift',mergerTreeBuildTreesBaseRedshift,defaultValue=0.0d0 )
       !@ <inputParameter>
       !@   <name>mergerTreeBuildTreesBeginAtTree</name>
       !@   <defaultValue>1</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The index (in order of increasing base halo mass) of the tree at which to begin when building merger trees.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildTreesBeginAtTree' ,mergerTreeBuildTreesBeginAtTree ,defaultValue=1     )
       nextTreeIndex=mergerTreeBuildTreesBeginAtTree
       !@ <inputParameter>
       !@   <name>mergerTreeBuildTreesHaloMassDistribution</name>
       !@   <defaultValue>uniform</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The method to be used to construct a distribution of base halo masses.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildTreesHaloMassDistribution',mergerTreeBuildTreesHaloMassDistribution,defaultValue="uniform")
       !@ <inputParameter>
       !@   <name>mergerTreeBuildTreesHaloMassExponent</name>
       !@   <defaultValue>1</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Halo masses will be (pseudo-)uniformly distributed in $[\log(M)]^{1/(1+\alpha)}$ where $\alpha=${\tt mergerTreeBuildTreesHaloMassExponent}.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildTreesHaloMassExponent',mergerTreeBuildTreesHaloMassExponent,defaultValue=1.0d0)
       ! Generate a randomly sampled set of halo masses.
       treeCount=max(2,int(dlog10(mergerTreeBuildHaloMassMaximum/mergerTreeBuildHaloMassMinimum)*dble(mergerTreeBuildTreesPerDecade)))
       call Alloc_Array(treeHaloMass,treeCount,'treeHaloMass')
       call Alloc_Array(treeWeight  ,treeCount,'treeWeight'  )

      ! Create a distribution of halo masses.
       select case (char(mergerTreeBuildTreesHaloMassDistribution))
       case ("quasi")
          ! Use a quasi-random sequence to generate halo masses.
          do iTree=1,treeCount
             treeHaloMass(iTree)=Quasi_Random_Get(quasiSequenceObject,reset=quasiSequenceReset)
          end do
          call Quasi_Random_Free(quasiSequenceObject)
          call Sort_Do(treeHaloMass)
       case ("uniform")
          ! Use a uniform distribution in logarithm of halo mass.
          treeHaloMass=Make_Range(0.0d0,1.0d0,treeCount,rangeType=rangeTypeLinear)
       case default
          call Galacticus_Error_Report('Merger_Tree_Build_Initialize','unknown halo mass distribution option')
       end select
       treeHaloMass=dexp((treeHaloMass**(1.0d0+mergerTreeBuildTreesHaloMassExponent))*dlog(mergerTreeBuildHaloMassMaximum&
            &/mergerTreeBuildHaloMassMinimum) +dlog(mergerTreeBuildHaloMassMinimum))

       ! Find the cosmic time at which the trees are based.
       expansionFactor=Expansion_Factor_from_Redshift(mergerTreeBuildTreesBaseRedshift)
       mergerTreeBuildTreesBaseTime=Cosmology_Age(expansionFactor)       
       ! Compute the weight (number of trees per unit volume) for each tree.
       do iTree=1,treeCount
          ! Get the minimum mass of the interval occupied by this tree.
          if (iTree==1) then
             MassMinimum=mergerTreeBuildHaloMassMinimum
          else
             MassMinimum=dsqrt(treeHaloMass(iTree)*treeHaloMass(iTree-1))
          end if
          ! Get the maximum mass of the interval occupied by this tree.
          if (iTree==treeCount) then
             MassMaximum=mergerTreeBuildHaloMassMaximum
          else
             MassMaximum=dsqrt(treeHaloMass(iTree)*treeHaloMass(iTree+1))
          end if
          ! Get the integral of the halo mass function over this range.
          treeWeight(iTree)=Halo_Mass_Function_Integrated(mergerTreeBuildTreesBaseTime,MassMinimum,MassMaximum)
       end do
       ! Determine which tree builder to use.
       !@ <inputParameter>
       !@   <name>mergerTreeBuildMethod</name>
       !@   <defaultValue>Cole2000</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used to build merger trees.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildMethod',mergerTreeBuildMethod,defaultValue='Cole2000')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="mergerTreeBuildMethod" type="code" action="subroutine">
       !#  <subroutineArgs>mergerTreeBuildMethod,Merger_Tree_Builder</subroutineArgs>
       include 'merger_trees.build.inc'
       !# </include>
       if (.not.associated(Merger_Tree_Builder)) call Galacticus_Error_Report('Merger_Tree_Build','method '&
            &//char(mergerTreeBuildMethod)//' is unrecognized')
    end if
    return
  end subroutine Merger_Tree_Build_Initialize

  subroutine Merger_Tree_Build_Do(thisTree)
    !% Build a merger tree.
    use Tree_Node_Methods
    use Galacticus_State
    implicit none
    type(mergerTree), intent(inout) :: thisTree

    ! Get a base halo mass and initialize. Do this within an OpenMP critical section so that threads don't try to get the same
    ! tree.
    !$omp critical (Merger_Tree_Build_Do)
    if (nextTreeIndex<=treeCount) then
       ! Retrieve stored internal state if possible.
       call Galacticus_State_Retrieve
       ! Take a snapshot of the internal state and store it.
       call Galacticus_State_Snapshot
       call Galacticus_State_Store
       ! Give the tree an index.
       thisTree%index=nextTreeIndex
       ! Create the base node.
       call thisTree%createNode(thisTree%baseNode,1)
       ! Assign a weight to the tree.
       thisTree%volumeWeight=treeWeight(nextTreeIndex)
       ! Assign a mass to it.
       call Tree_Node_Mass_Set(thisTree%baseNode,treeHaloMass(nextTreeIndex))
       ! Assign a time.
       call Tree_Node_Time_Set(thisTree%baseNode,mergerTreeBuildTreesBaseTime)
       ! Increment the tree index counter.
       nextTreeIndex=nextTreeIndex+1
    end if
    !$omp end critical (Merger_Tree_Build_Do)
    ! If we got a tree, we can now process it (in paralell if running under OpenMP).
    if (associated(thisTree%baseNode)) then
       ! Call routine to actually build the tree.
       call Merger_Tree_Builder(thisTree)
    end if
    return
  end subroutine Merger_Tree_Build_Do

end module Merger_Tree_Build
