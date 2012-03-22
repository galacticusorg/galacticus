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


!% Contains a module which implements building of merger trees after drawing masses at random from a mass function.

module Merger_Tree_Build
  !% Implements building of merger trees after drawing masses at random from a mass function.
  use Merger_Trees
  use ISO_Varying_String
  implicit none
  private
  public :: Merger_Tree_Build_Initialize

  ! Variables giving the mass range and sampling frequency for mass function sampling.
  double precision     :: mergerTreeBuildHaloMassMinimum,mergerTreeBuildHaloMassMaximum,mergerTreeBuildTreesBaseRedshift &
       &,mergerTreeBuildTreesBaseTime,mergerTreeBuildTreesHaloMassExponent
  integer              :: mergerTreeBuildTreesPerDecade,mergerTreeBuildTreesBeginAtTree
  type(varying_string) :: mergerTreeBuildTreesHaloMassDistribution,mergerTreeBuildTreeMassesFile

  ! Direction in which to process trees.
  logical              :: mergerTreeBuildTreesProcessDescending

  ! Array of halo masses to use.
  integer                                     :: treeCount,nextTreeIndex
  double precision, allocatable, dimension(:) :: treeHaloMass,treeWeight

  ! Name of merger tree builder method.
  type(varying_string) :: mergerTreeBuildMethod
  ! Pointer to the subroutine that builds the merger tree.
  procedure(Merger_Tree_Builder_Template), pointer :: Merger_Tree_Builder => null()
  abstract interface
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
    use ISO_Varying_String
    use FoX_dom
    !# <include directive="mergerTreeBuildMethod" type="moduleUse">
    include 'merger_trees.build.modules.inc'
    !# </include>
    implicit none
    type(varying_string),          intent(in)    :: mergerTreeConstructMethod
    procedure(),          pointer, intent(inout) :: Merger_Tree_Construct
    type(Node),           pointer                :: doc,thisTree
    type(NodeList),       pointer                :: rootMassList
    integer                                      :: iTree,ioErr
    type(fgsl_qrng)                              :: quasiSequenceObject
    logical                                      :: quasiSequenceReset=.true.
    double precision                             :: expansionFactor,massMinimum,massMaximum

    ! Check if our method is to be used.
    if (mergerTreeConstructMethod == 'build') then
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
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildHaloMassMinimum'  ,mergerTreeBuildHaloMassMinimum  ,defaultValue=1.0d10)
       !@ <inputParameter>
       !@   <name>mergerTreeBuildHaloMassMaximum</name>
       !@   <defaultValue>$10^{15}$</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum mass of merger tree base halos to consider when building merger trees, in units of $M_\odot$.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildHaloMassMaximum'  ,mergerTreeBuildHaloMassMaximum  ,defaultValue=1.0d15)
       !@ <inputParameter>
       !@   <name>mergerTreeBuildTreesPerDecade</name>
       !@   <defaultValue>10</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The number of merger trees to build per decade of base halo mass.
       !@   </description>
       !@   <type>integer</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildTreesPerDecade'   ,mergerTreeBuildTreesPerDecade   ,defaultValue=10    )
       !@ <inputParameter>
       !@   <name>mergerTreeBuildTreesBaseRedshift</name>
       !@   <defaultValue>0</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The redshift at which to plant the base node when building merger trees.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildTreesBaseRedshift',mergerTreeBuildTreesBaseRedshift,defaultValue=0.0d0 )
       !@ <inputParameter>
       !@   <name>mergerTreeBuildTreesBeginAtTree</name>
       !@   <defaultValue>1</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The index (in order of increasing base halo mass) of the tree at which to begin when building merger trees.
       !@   </description>
       !@   <type>integer</type>
       !@   <cardinality>1</cardinality>
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
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildTreesHaloMassDistribution',mergerTreeBuildTreesHaloMassDistribution,defaultValue="uniform")
       !@ <inputParameter>
       !@   <name>mergerTreeBuildTreesHaloMassExponent</name>
       !@   <defaultValue>1</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Halo masses will be (pseudo-)uniformly distributed in $[\log(M)]^{1/(1+\alpha)}$ where $\alpha=${\tt mergerTreeBuildTreesHaloMassExponent}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildTreesHaloMassExponent',mergerTreeBuildTreesHaloMassExponent,defaultValue=1.0d0)
       !@ <inputParameter>
       !@   <name>mergerTreeBuildTreesProcessDescending</name>
       !@   <defaultValue>false</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     If true, causes merger trees to be processed in order of decreasing mass.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildTreesProcessDescending',mergerTreeBuildTreesProcessDescending,defaultValue=.false.)
       !@ <inputParameter>
       !@   <name>mergerTreeBuildTreeMassesFile</name>
       !@   <defaultValue>null</defaultValue>       
       !@   <description>
       !@     Specifies the name of a file from which to read the masses of merger tree root halos when building merger trees.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildTreeMassesFile',mergerTreeBuildTreeMassesFile,defaultValue='null')

       ! Generate a randomly sampled set of halo masses.
       treeCount=max(2,int(dlog10(mergerTreeBuildHaloMassMaximum/mergerTreeBuildHaloMassMinimum)*dble(mergerTreeBuildTreesPerDecade)))
       call Alloc_Array(treeHaloMass,[treeCount])
       call Alloc_Array(treeWeight  ,[treeCount])

       ! Determine how to compute the tree root masses.
       select case (char(mergerTreeBuildTreesHaloMassDistribution))
       case ("quasi","uniform")
          ! Generate a randomly sampled set of halo masses.
          treeCount=max(2,int(dlog10(mergerTreeBuildHaloMassMaximum/mergerTreeBuildHaloMassMinimum)*dble(mergerTreeBuildTreesPerDecade)))
          call Alloc_Array(treeHaloMass,[treeCount])
          call Alloc_Array(treeWeight  ,[treeCount])
          
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
               &/mergerTreeBuildHaloMassMinimum)+dlog(mergerTreeBuildHaloMassMinimum))
       case ("read")
          ! Read masses from a file.
          !$omp critical (FoX_DOM_Access)
          doc => parseFile(char(mergerTreeBuildTreeMassesFile),iostat=ioErr)
          if (ioErr /= 0) call Galacticus_Error_Report('Merger_Tree_Build_Initialize','unable to read or parse merger tree root mass file')
          ! Get a list of all defined tree root masses.
          rootMassList => getElementsByTagname(doc,"treeRootMass")
          ! Get a count of the number of trees.
          treeCount=getLength(rootMassList)
          ! Allocate arrays for tree masses.
          call Alloc_Array(treeHaloMass,[treeCount])
          call Alloc_Array(treeWeight  ,[treeCount])
          ! Extract the tree masses from the XML data.
          do iTree=1,treeCount
             thisTree => item(rootMassList,iTree-1)
             call extractDataContent(thisTree,treeHaloMass(iTree))
          end do
          ! Finished - destroy the XML document.
          call destroy(doc)
          !$omp end critical (FoX_DOM_Access)
       case default
          call Galacticus_Error_Report('Merger_Tree_Build_Initialize','unknown halo mass distribution option')
       end select
       
       ! Find the cosmic time at which the trees are based.
       expansionFactor=Expansion_Factor_from_Redshift(mergerTreeBuildTreesBaseRedshift)
       mergerTreeBuildTreesBaseTime=Cosmology_Age(expansionFactor)       
       ! Compute the weight (number of trees per unit volume) for each tree.
       do iTree=1,treeCount
          ! Get the minimum mass of the interval occupied by this tree.
          if (iTree==1) then
             if (char(mergerTreeBuildTreesHaloMassDistribution) == "read") then
                MassMinimum=treeHaloMass(iTree)*dsqrt(treeHaloMass(iTree)/treeHaloMass(iTree+1))
             else
                MassMinimum=mergerTreeBuildHaloMassMinimum
             end if
          else
             MassMinimum=dsqrt(treeHaloMass(iTree)*treeHaloMass(iTree-1))
          end if
          ! Get the maximum mass of the interval occupied by this tree.
          if (iTree==treeCount) then
             if (char(mergerTreeBuildTreesHaloMassDistribution) == "read") then
                MassMaximum=treeHaloMass(iTree)*dsqrt(treeHaloMass(iTree)/treeHaloMass(iTree-1))
             else
                MassMaximum=mergerTreeBuildHaloMassMaximum
             end if
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
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
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

  subroutine Merger_Tree_Build_Do(thisTree,skipTree)
    !% Build a merger tree.
    use Tree_Nodes
    use Galacticus_State
    use Kind_Numbers
    use String_Handling
    use ISO_Varying_String
    implicit none
    type(mergerTree),        intent(inout) :: thisTree
    logical,                 intent(in)    :: skipTree
    integer(kind=kind_int8), parameter     :: baseNodeIndex=1
    integer(kind=kind_int8)                :: thisTreeIndex
    type(varying_string)                   :: message

    ! Get a base halo mass and initialize. Do this within an OpenMP critical section so that threads don't try to get the same
    ! tree.
    !$omp critical (Merger_Tree_Build_Do)
    if (nextTreeIndex<=treeCount) then
       ! Retrieve stored internal state if possible.
       call Galacticus_State_Retrieve
       ! Take a snapshot of the internal state and store it.
       call Galacticus_State_Snapshot
       message='Storing state for tree #'
       message=message//nextTreeIndex
       call Galacticus_State_Store(message)
       ! Determine the index of the tree to process.
       select case (mergerTreeBuildTreesProcessDescending)
       case(.false.)
          ! Processing trees in ascending order, to just use nextTreeIndex as the index of the tree to process.
          thisTreeIndex=nextTreeIndex
       case(.true. )
          ! Processing trees in descending order, so begin from the final index and work back.
          thisTreeIndex=treeCount+1-nextTreeIndex
       end select
       ! Give the tree an index.
       thisTree%index=thisTreeIndex
       ! Create the base node.
       call thisTree%createNode(thisTree%baseNode,baseNodeIndex)
       ! Assign a weight to the tree.
       thisTree%volumeWeight=treeWeight(thisTreeIndex)
       ! Assign a mass to it.
       call Tree_Node_Mass_Set(thisTree%baseNode,treeHaloMass(thisTreeIndex))
       ! Assign a time.
       call Tree_Node_Time_Set(thisTree%baseNode,mergerTreeBuildTreesBaseTime)
       ! Increment the tree index counter.
       nextTreeIndex=nextTreeIndex+1
    end if
    !$omp end critical (Merger_Tree_Build_Do)
    ! If we got a tree, we can now process it (in paralell if running under OpenMP).
    if (associated(thisTree%baseNode).and..not.skipTree) then
       ! Call routine to actually build the tree.
       call Merger_Tree_Builder(thisTree)
    end if
    return
  end subroutine Merger_Tree_Build_Do

end module Merger_Tree_Build
