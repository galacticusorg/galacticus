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


!% Contains a module which implements building of simple merger trees with smooth mass accretion histories and no branches using
!% the fitting function of \cite{wechsler_concentrations_2002}.

module Merger_Tree_Smooth_Accretion
  !% Implements building of simple merger trees with smooth mass accretion histories and no branches using the fitting function of
  !% \cite{wechsler_concentrations_2002}.
  use Merger_Trees
  use ISO_Varying_String
  implicit none
  private
  public :: Merger_Tree_Smooth_Accretion_Initialize

  ! Variables giving properties of the merger tree.
  double precision :: mergerTreeHaloMass,mergerTreeBuildTreesBaseTime,mergerTreeHaloMassResolution&
       &,mergerTreeHaloMassDeclineFactor,mergerTreeBaseRedshift

  ! Flag indicating whether or not we've already built the tree.
  logical          :: treeWasBuilt=.false.
  
contains

  !# <mergerTreeConstructMethod>
  !#  <unitName>Merger_Tree_Smooth_Accretion_Initialize</unitName>
  !# </mergerTreeConstructMethod>
  subroutine Merger_Tree_Smooth_Accretion_Initialize(mergerTreeConstructMethod,Merger_Tree_Construct)
    !% Initializes the smooth accretion merger tree module.
    use Input_Parameters
    use Cosmology_Functions
    implicit none
    type(varying_string),          intent(in)    :: mergerTreeConstructMethod
    procedure(),          pointer, intent(inout) :: Merger_Tree_Construct

    ! Check if our method is to be used.
    if (mergerTreeConstructMethod == 'smoothAccretion') then
       ! Assign pointer to our merger tree construction subroutine.
       Merger_Tree_Construct => Merger_Tree_Smooth_Accretion_Do
       ! Read the mass of the halo to construct.
       !@ <inputParameter>
       !@   <name>mergerTreeHaloMass</name>
       !@   <defaultValue>$10^{12}$</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The final mass of the merger tree base halo to consider when building a smoothly accreting merger tree, in units of $M_\odot$.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeHaloMass'                ,mergerTreeHaloMass                ,defaultValue=1.0d12)
       !@ <inputParameter>
       !@   <name>mergerTreeHaloMassResolution</name>
       !@   <defaultValue>$10^{12}$</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The final mass of the merger tree base halo to consider when building a smoothly accreting merger tree, in units of $M_\odot$.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeHaloMassResolution'      ,mergerTreeHaloMassResolution      ,defaultValue=1.0d9 )
       !@ <inputParameter>
       !@   <name>mergerTreeHaloMassDeclineFactor</name>
       !@   <defaultValue>$0.9$</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The factor by which halo mass should decrease in each step back in time building a smoothly accreting merger tree, in units of $M_\odot$.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeHaloMassDeclineFactor'   ,mergerTreeHaloMassDeclineFactor   ,defaultValue=0.9d0 )
       !@ <inputParameter>
       !@   <name>mergerTreeBaseRedshift</name>
       !@   <defaultValue>0</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The redshift at which to plant the base node when building the smoothly accreting merger tree.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBaseRedshift'            ,mergerTreeBaseRedshift            ,defaultValue=0.0d0 )
    end if
    return
  end subroutine Merger_Tree_Smooth_Accretion_Initialize

  subroutine Merger_Tree_Smooth_Accretion_Do(thisTree,skipTree)
    !% Build a merger tree with a smooth mass accretion history using the fitting function of \cite{wechsler_concentrations_2002}.
    use Tree_Nodes
    use Cosmology_Functions
    use Kind_Numbers
    use Dark_Matter_Halo_Mass_Accretion_Histories
    implicit none
    type(mergerTree),       intent(inout) :: thisTree
    logical,                intent(in)    :: skipTree
    type(treeNode),         pointer       :: currentNode,newNode
    integer(kind=kind_int8)               :: nodeIndex
    double precision                      :: mergerTreeBaseTime,expansionFactorBase,nodeMass,nodeTime

    ! Build the merger tree.
    !$omp critical (Merger_Tree_Build_Do)
    if (.not.treeWasBuilt) then
       ! Give the tree an index.
       thisTree%index=1
       ! Create the base node.
       nodeIndex=1
       call thisTree%createNode(thisTree%baseNode,nodeIndex)
       ! Assign an arbitrary weight to the tree.
       thisTree%volumeWeight=1.0
       ! Assign a mass to the base node.
       call Tree_Node_Mass_Set(thisTree%baseNode,mergerTreeHaloMass)
       ! Find the cosmic time at which the tree is based.
       expansionFactorBase=Expansion_Factor_from_Redshift(mergerTreeBaseRedshift)
       mergerTreeBaseTime =Cosmology_Age                 (expansionFactorBase   )
       ! Assign a time to the base node.
       call Tree_Node_Time_Set(thisTree%baseNode,mergerTreeBaseTime)
       ! Get a pointer to the current node (i.e. the base node).       
       currentNode => thisTree%baseNode
       ! Initialize current node mass.
       nodeMass=mergerTreeHaloMass
       ! Initialize node index counter.
       nodeIndex=1
       ! Step backwards, creating nodes until a sufficiently low mass has been reached.
       do while (nodeMass > mergerTreeHaloMassResolution)
          ! Increment node index.
          nodeIndex=nodeIndex+1
          ! Create a node.
          call thisTree%createNode(newNode,nodeIndex)
          ! Adjust the mass by the specified factor.
          nodeMass=nodeMass*mergerTreeHaloMassDeclineFactor
          ! Set the mass of the node.
          call Tree_Node_Mass_Set (newNode,nodeMass )
          ! Find the time corresponding to this expansion factor.
          nodeTime=Dark_Matter_Halo_Mass_Accretion_Time(thisTree%baseNode,nodeMass)
          ! Set the time for the new node.
          call Tree_Node_Time_Set (newNode,nodeTime )
          ! Create parent and child links.
          currentNode%childNode  => newNode
          newNode    %parentNode => currentNode
          ! Move the current node to the new node.
          currentNode => newNode
       end do
       ! Flag that the tree is now built.
       treeWasBuilt=.true.
    end if
    !$omp end critical (Merger_Tree_Build_Do)

    return
  end subroutine Merger_Tree_Smooth_Accretion_Do

end module Merger_Tree_Smooth_Accretion
