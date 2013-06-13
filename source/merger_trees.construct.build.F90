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
       &,mergerTreeBuildTreesBaseTime,mergerTreeBuildTreesPerDecade
  integer              :: mergerTreeBuildTreesBeginAtTree
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
    use, intrinsic :: ISO_C_Binding
    use Input_Parameters
    use Memory_Management
    use Cosmology_Functions
    use FGSL
    use Quasi_Random
    use Halo_Mass_Function
    use Sort
    use Galacticus_Error
    use Numerical_Ranges
    use Numerical_Integration
    use Numerical_Interpolation
    use FoX_dom
    !# <include directive="mergerTreeBuildMethod" type="moduleUse">
    include 'merger_trees.build.modules.inc'
    !# </include>
    implicit none
    type(varying_string),          intent(in)       :: mergerTreeConstructMethod
    procedure(),          pointer, intent(inout)    :: Merger_Tree_Construct
    type(Node),           pointer                   :: doc,thisTree
    type(NodeList),       pointer                   :: rootMassList
    integer,              parameter                 :: massFunctionSamplePerDecade=100
    double precision,     parameter                 :: toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3
    double precision,     allocatable, dimension(:) :: massFunctionSampleLogMass,massFunctionSampleProbability,massFunctionSampleLogMassMonotonic
    integer                                         :: iTree,ioErr,iSample,jSample,massFunctionSampleCount
    type(fgsl_qrng)                                 :: quasiSequenceObject
    logical                                         :: quasiSequenceReset=.true.
    double precision                                :: expansionFactor,massMinimum,massMaximum,massFunctionSampleLogPrevious,probability
    type(fgsl_function)                             :: integrandFunction
    type(fgsl_integration_workspace)                :: integrationWorkspace
    type(fgsl_interp)                               :: interpolationObject
    type(fgsl_interp_accel)                         :: interpolationAccelerator
    type(c_ptr)                                     :: parameterPointer
    logical                                         :: integrandReset=.true.,interpolationReset=.true.

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
       call Get_Input_Parameter('mergerTreeBuildTreesPerDecade'   ,mergerTreeBuildTreesPerDecade   ,defaultValue=10.0d0)
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
       !@   <name>mergerTreeBuildTreesProcessDescending</name>
       !@   <defaultValue>true</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     If true, causes merger trees to be processed in order of decreasing mass.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildTreesProcessDescending',mergerTreeBuildTreesProcessDescending,defaultValue=.true.)
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
       
       ! Find the cosmic time at which the trees are based.
       expansionFactor=Expansion_Factor_from_Redshift(mergerTreeBuildTreesBaseRedshift)
       mergerTreeBuildTreesBaseTime=Cosmology_Age(expansionFactor)       

       ! Generate a randomly sampled set of halo masses.
       treeCount=max(2,int(log10(mergerTreeBuildHaloMassMaximum/mergerTreeBuildHaloMassMinimum)*mergerTreeBuildTreesPerDecade))
       call Alloc_Array(treeHaloMass,[treeCount])
       call Alloc_Array(treeWeight  ,[treeCount])

       ! Determine how to compute the tree root masses.
       select case (char(mergerTreeBuildTreesHaloMassDistribution))
       case ("quasi","uniform")
          ! Generate a randomly sampled set of halo masses.
          treeCount=max(2,int(log10(mergerTreeBuildHaloMassMaximum/mergerTreeBuildHaloMassMinimum)*mergerTreeBuildTreesPerDecade))
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
          ! Create a cumulative probability for sampling halo masses.
          massFunctionSampleCount=max(2,int(log10(mergerTreeBuildHaloMassMaximum/mergerTreeBuildHaloMassMinimum)*massFunctionSamplePerDecade))
          call Alloc_Array(massFunctionSampleLogMass         ,[massFunctionSampleCount])
          call Alloc_Array(massFunctionSampleLogMassMonotonic,[massFunctionSampleCount])
          call Alloc_Array(massFunctionSampleProbability     ,[massFunctionSampleCount])
          massFunctionSampleLogMass=Make_Range(log10(mergerTreeBuildHaloMassMinimum),log10(mergerTreeBuildHaloMassMaximum),massFunctionSampleCount,rangeType=rangeTypeLogarithmic)
          massFunctionSampleLogPrevious=log10(mergerTreeBuildHaloMassMinimum)
          jSample=0
          do iSample=1,massFunctionSampleCount
             probability=Integrate(massFunctionSampleLogPrevious&
                  &,massFunctionSampleLogMass(iSample),Mass_Function_Sampling_Integrand,parameterPointer,integrandFunction &
                  &,integrationWorkspace,toleranceAbsolute=toleranceAbsolute,toleranceRelative=toleranceRelative,reset=integrandReset)
             if (iSample == 1 .or. probability > 0.0d0) then
                jSample=jSample+1
                massFunctionSampleProbability     (jSample)=probability
                massFunctionSampleLogMassMonotonic(jSample)=massFunctionSampleLogMass(iSample)
                if (jSample > 1) massFunctionSampleProbability(jSample)=massFunctionSampleProbability(jSample)+massFunctionSampleProbability(jSample-1)
             end if
             massFunctionSampleLogPrevious=massFunctionSampleLogMass(iSample)
          end do
          call Integrate_Done(integrandFunction,integrationWorkspace)
          massFunctionSampleCount=jSample
          if (massFunctionSampleCount < 2) call Galacticus_Error_Report('Merger_Tree_Build_Initialize','tabulated mass function sampling density has fewer than 2 non-zero points')
          ! Normalize the cumulative probability distribution.
          massFunctionSampleProbability=massFunctionSampleProbability/massFunctionSampleProbability(massFunctionSampleCount)
          ! Compute the corresponding halo masses by interpolation in the cumulative probability distribution function.
          do iTree=1,treeCount
             treeHaloMass(iTree)=Interpolate(massFunctionSampleCount,massFunctionSampleProbability(1:massFunctionSampleCount)&
                  &,massFunctionSampleLogMassMonotonic(1:massFunctionSampleCount),interpolationObject,interpolationAccelerator&
                  &,treeHaloMass(iTree) ,reset=interpolationReset,extrapolationType=extrapolationTypeFixed)
          end do
          treeHaloMass=10.0d0**treeHaloMass
          call Interpolate_Done(interpolationObject,interpolationAccelerator,interpolationReset)
          call Dealloc_Array(massFunctionSampleLogMass         )
          call Dealloc_Array(massFunctionSampleProbability     )
          call Dealloc_Array(massFunctionSampleLogMassMonotonic)
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

       ! Compute the weight (number of trees per unit volume) for each tree.
       do iTree=1,treeCount
          ! Get the minimum mass of the interval occupied by this tree.
          if (iTree==1) then
             if (char(mergerTreeBuildTreesHaloMassDistribution) == "read") then
                massMinimum=treeHaloMass(iTree)*sqrt(treeHaloMass(iTree)/treeHaloMass(iTree+1))
             else
                massMinimum=mergerTreeBuildHaloMassMinimum
             end if
          else
             massMinimum=sqrt(treeHaloMass(iTree)*treeHaloMass(iTree-1))
          end if
          ! Get the maximum mass of the interval occupied by this tree.
          if (iTree==treeCount) then
             if (char(mergerTreeBuildTreesHaloMassDistribution) == "read") then
                massMaximum=treeHaloMass(iTree)*sqrt(treeHaloMass(iTree)/treeHaloMass(iTree-1))
             else
                massMaximum=mergerTreeBuildHaloMassMaximum
             end if
          else
             massMaximum=sqrt(treeHaloMass(iTree)*treeHaloMass(iTree+1))
          end if
          ! For distributions of masses, adjust the masses at the end points so that they are at the
          ! geometric mean of their range.
          if     (                                                          &
               &   (iTree == 1 .or. iTree == treeCount)                     &
               &  .and.                                                     &
               &   char(mergerTreeBuildTreesHaloMassDistribution) /= "read" &
               & ) treeHaloMass(iTree)=sqrt(massMinimum*massMaximum)
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
       !# <include directive="mergerTreeBuildMethod" type="functionCall" functionType="void">
       !#  <functionArgs>mergerTreeBuildMethod,Merger_Tree_Builder</functionArgs>
       include 'merger_trees.build.inc'
       !# </include>
       if (.not.associated(Merger_Tree_Builder)) call Galacticus_Error_Report('Merger_Tree_Build','method '&
            &//char(mergerTreeBuildMethod)//' is unrecognized')
    end if
    return
  end subroutine Merger_Tree_Build_Initialize

  function Mass_Function_Sampling_Integrand(logMass,parameterPointer) bind(c)
    !% The integrand over the mass function sampling density function.
    use, intrinsic :: ISO_C_Binding
    use Merger_Trees_Mass_Function_Sampling
    implicit none
    real(c_double)        :: Mass_Function_Sampling_Integrand
    real(c_double), value :: logMass
    type(c_ptr),    value :: parameterPointer

    Mass_Function_Sampling_Integrand=Merger_Tree_Construct_Mass_Function_Sampling(10.0d0**logMass,mergerTreeBuildTreesBaseTime,mergerTreeBuildHaloMassMinimum,mergerTreeBuildHaloMassMaximum)
    return
  end function Mass_Function_Sampling_Integrand
  
  subroutine Merger_Tree_Build_Do(thisTree,skipTree)
    !% Build a merger tree.
    use Galacticus_Nodes
    use Galacticus_State
    use Kind_Numbers
    use String_Handling
    implicit none
    type(mergerTree),          intent(inout) :: thisTree
    logical,                   intent(in)    :: skipTree
    class(nodeComponentBasic), pointer       :: baseNodeBasicComponent
    integer(kind=kind_int8),   parameter     :: baseNodeIndex=1
    integer(kind=kind_int8)                  :: thisTreeIndex
    type(varying_string)                     :: message

    ! Get a base halo mass and initialize. Do this within an OpenMP critical section so that threads don't try to get the same
    ! tree.
    !$omp critical (Merger_Tree_Build_Do)
    if (nextTreeIndex <= treeCount) then
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
       ! Get the basic component of the base node.
       baseNodeBasicComponent => thisTree%baseNode%basic(autoCreate=.true.)
       ! Assign a mass to it.
       call baseNodeBasicComponent%massSet(treeHaloMass(thisTreeIndex) )
       ! Assign a time.
       call baseNodeBasicComponent%timeSet(mergerTreeBuildTreesBaseTime)
       ! Increment the tree index counter.
       nextTreeIndex=nextTreeIndex+1
    end if
    !$omp end critical (Merger_Tree_Build_Do)
    ! If we got a tree, we can now process it (in parallel if running under OpenMP).
    if (associated(thisTree%baseNode).and..not.skipTree) then
       ! Call routine to actually build the tree.
       call Merger_Tree_Builder(thisTree)
    end if
    return
  end subroutine Merger_Tree_Build_Do

end module Merger_Tree_Build
