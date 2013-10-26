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

!% Contains a module which computes conditional mass functions in merger trees.

module Merger_Trees_Conditional_Mass_Function
  !% Computes conditional mass functions in merger trees.
  implicit none
  private
  public :: Merger_Tree_Conditional_Mass_Function, Merger_Tree_Conditional_Mass_Function_Output

  ! Flag indicating if module is initialized.
  logical                                           :: moduleInitialized=.false.

  ! Arrays to store parent halo masses and times.
  double precision, allocatable, dimension(:      ) :: timeParents,timeProgenitors&
       &,mergerTreeComputeConditionalMassFunctionParentRedshifts,mergerTreeComputeConditionalMassFunctionProgenitorRedshifts&
       &,massParents,massRatios

  ! Array to store the conditional mass function.
  double precision, allocatable, dimension(:,:    ) :: normalization
  double precision, allocatable, dimension(:,:,:  ) :: conditionalMassFunction      ,conditionalMassFunctionError
  double precision, allocatable, dimension(:,:,:,:) :: primaryProgenitorMassFunction,primaryProgenitorMassFunctionError
  double precision, allocatable, dimension(:,:,:,:) :: formationRateFunction        ,formationRateFunctionError

  ! Extent and size of the conditional mass function arrays.
  integer                                           :: mergerTreeComputeConditionalMassFunctionParentMassCount&
       &,mergerTreeComputeConditionalMassFunctionMassRatioCount,massFunctionTimeCount&
       &,mergerTreeComputeConditionalMassFunctionPrimaryProgenitorDepth
  double precision                                  :: massParentLogarithmicMinimum,massRatioLogarithmicMinimum &
       &,massParentLogarithmicBinWidthInverse,massRatioLogarithmicBinWidthInverse&
       &,mergerTreeConditionalMassFunctionFormationRateTimeFraction

  ! Option indicating if conditional mass function calculation is required.
  logical                                           :: mergerTreeComputeConditionalMassFunction

contains

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Merger_Tree_Conditional_Mass_Function</unitName>
  !#   <after>Merger_Tree_Monotonic_Mass_Growth</after>
  !#   <after>Merger_Tree_Prune_Branches</after>
  !#   <after>Merger_Tree_Prune_Hierarchy</after>
  !#   <after>Merger_Tree_Regrid_Time</after>
  !# </mergerTreePreEvolveTask>
  subroutine Merger_Tree_Conditional_Mass_Function(thisTree)
    !% Compute conditional mass function on {\tt thisTree}.
    use Galacticus_Nodes
    use Input_Parameters
    use Memory_Management
    use Cosmology_Functions
    use Numerical_Ranges
    use Numerical_Comparison
    use Galacticus_Error
    implicit none
    type            (mergerTree             ), intent(in   ), target                 :: thisTree
    type            (treeNode               ), pointer                               :: thisNode,childNode,parentNode,parentsChildNode&
         &,descendentNode
    type            (mergerTree             ), pointer                               :: currentTree
    class           (nodeComponentBasic     ), pointer                               :: thisBasicComponent,childBasicComponent&
         &,parentBasicComponent,descendentBasic,parentsChildBasic
    class           (cosmologyFunctionsClass), pointer                               :: cosmologyFunctionsDefault
    double precision                         , allocatable  , dimension(:,:,:), save :: primaryProgenitorMass
    !$omp threadprivate(primaryProgenitorMass)
    integer                                                                          :: i,binMassParent,binMassRatio,iPrimary,jPrimary&
         &,binMassRatioCreation,binMassRatioDescendent
    double precision                                                                 :: branchBegin,branchEnd,parentBranchBegin &
         &,parentBranchEnd ,massProgenitor,massParent,branchMassInitial,branchMassFinal,parentBranchMassInitial &
         &,parentBranchMassFinal ,massRatioLogarithmic,massParentLogarithmic,massRatio &
         &,mergerTreeComputeConditionalMassFunctionParentMassMinimum,mergerTreeComputeConditionalMassFunctionParentMassMaximum &
         &,mergerTreeComputeConditionalMassFunctionMassRatioMinimum ,mergerTreeComputeConditionalMassFunctionMassRatioMaximum&
         &,massRatioCreation,massRatioCreationLogarithmic,massRatioDescendent,massRatioDescendentLogarithmic

    ! Check if module is initialized.
    if (.not.moduleInitialized) then
       !$omp critical (Merger_Tree_Conditional_Mass_Function_Initialize)
       if (.not.moduleInitialized) then
          ! Get parameter specifying if pruning is required.
          !@ <inputParameter>
          !@   <name>mergerTreeComputeConditionalMassFunction</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not the conditional mass function of merger trees should be computed.
          !@   </description>
          !@   <type>logical</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeComputeConditionalMassFunction',mergerTreeComputeConditionalMassFunction,defaultValue=.false.)
          if (mergerTreeComputeConditionalMassFunction) then
             !@ <inputParameter>
             !@   <name>mergerTreeComputeConditionalMassFunctionParentMassCount</name>
             !@   <defaultValue>10</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The number of bins in parent mass when constructing conditional halo mass functions.
             !@   </description>
             !@   <type>logical</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeComputeConditionalMassFunctionParentMassCount',mergerTreeComputeConditionalMassFunctionParentMassCount,defaultValue=10)
             !@ <inputParameter>
             !@   <name>mergerTreeComputeConditionalMassFunctionParentMassMinimum</name>
             !@   <defaultValue>$10^{10}M_\odot$</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The minimum parent halo mass to bin when constructing conditional halo mass functions.
             !@   </description>
             !@   <type>logical</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeComputeConditionalMassFunctionParentMassMinimum',mergerTreeComputeConditionalMassFunctionParentMassMinimum,defaultValue=1.0d10)
             !@ <inputParameter>
             !@   <name>mergerTreeComputeConditionalMassFunctionParentMassMaximum</name>
             !@   <defaultValue>$10^{15}M_\odot$</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The maximum parent halo mass to bin when constructing conditional halo mass functions.
             !@   </description>
             !@   <type>logical</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeComputeConditionalMassFunctionParentMassMaximum',mergerTreeComputeConditionalMassFunctionParentMassMaximum,defaultValue=1.0d15)
             !@ <inputParameter>
             !@   <name>mergerTreeComputeConditionalMassFunctionMassRatioCount</name>
             !@   <defaultValue>10</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The number of bins in mass ratio when constructing conditional halo mass functions.
             !@   </description>
             !@   <type>logical</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeComputeConditionalMassFunctionMassRatioCount',mergerTreeComputeConditionalMassFunctionMassRatioCount,defaultValue=10)
             !@ <inputParameter>
             !@   <name>mergerTreeComputeConditionalMassFunctionMassRatioMinimum</name>
             !@   <defaultValue>$10^{10}M_\odot$</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The minimum mass ratio to bin when constructing conditional halo mass functions.
             !@   </description>
             !@   <type>logical</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeComputeConditionalMassFunctionMassRatioMinimum',mergerTreeComputeConditionalMassFunctionMassRatioMinimum,defaultValue=1.0d-4)
             !@ <inputParameter>
             !@   <name>mergerTreeComputeConditionalMassFunctionMassRatioMaximum</name>
             !@   <defaultValue>$10^{15}M_\odot$</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The maximum mass ratio to bin when constructing conditional halo mass functions.
             !@   </description>
             !@   <type>logical</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeComputeConditionalMassFunctionMassRatioMaximum',mergerTreeComputeConditionalMassFunctionMassRatioMaximum,defaultValue=1.0d+1)
             !@ <inputParameter>
             !@   <name>mergerTreeComputeConditionalMassFunctionParentRedshifts</name>
             !@   <defaultValue>0.0</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The set of parent halo redshifts to use when constructing conditional halo mass functions.
             !@   </description>
             !@   <type>logical</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             massFunctionTimeCount=Get_Input_Parameter_Array_Size('mergerTreeComputeConditionalMassFunctionParentRedshifts')
             if (Get_Input_Parameter_Array_Size('mergerTreeComputeConditionalMassFunctionParentRedshifts') /= massFunctionTimeCount) &
                  & call Galacticus_Error_Report('Merger_Tree_Conditional_Mass_Function','mismatch in sizes of parent and progenitor redshift arrays')
             call Alloc_Array(mergerTreeComputeConditionalMassFunctionParentRedshifts    ,[massFunctionTimeCount])
             call Alloc_Array(mergerTreeComputeConditionalMassFunctionProgenitorRedshifts,[massFunctionTimeCount])
             call Alloc_Array(timeProgenitors                                            ,[massFunctionTimeCount])
             call Alloc_Array(timeParents                                                ,[massFunctionTimeCount])
             call Get_Input_Parameter('mergerTreeComputeConditionalMassFunctionParentRedshifts',mergerTreeComputeConditionalMassFunctionParentRedshifts,defaultValue=[0.0d0])
             !@ <inputParameter>
             !@   <name>mergerTreeComputeConditionalMassFunctionProgenitorRedshifts</name>
             !@   <defaultValue>1.0</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The set of progenitor halo redshifts to use when constructing conditional halo mass functions.
             !@   </description>
             !@   <type>logical</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeComputeConditionalMassFunctionProgenitorRedshifts',mergerTreeComputeConditionalMassFunctionProgenitorRedshifts,defaultValue=[1.0d0])
             !@ <inputParameter>
             !@   <name>mergerTreeComputeConditionalMassFunctionPrimaryProgenitorDepth</name>
             !@   <defaultValue>2</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The depth in progenitor ranking for which to store ranked progenitor mass functions. For example, a value of 2 means store mass functions for the most massive, and second most massive progenitor.
             !@   </description>
             !@   <type>logical</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeComputeConditionalMassFunctionPrimaryProgenitorDepth',mergerTreeComputeConditionalMassFunctionPrimaryProgenitorDepth,defaultValue=2)
             !@ <inputParameter>
             !@   <name>mergerTreeConditionalMassFunctionFormationRateTimeFraction</name>
             !@   <defaultValue>0.01</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The fraction of the current time over which to estimate the formation rate of halos when computing merger tree statistics.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeConditionalMassFunctionFormationRateTimeFraction',mergerTreeConditionalMassFunctionFormationRateTimeFraction,defaultValue=0.01d0)
             ! Construct bins for parent node mass.
             massParentLogarithmicMinimum        =  log( mergerTreeComputeConditionalMassFunctionParentMassMinimum)
             massParentLogarithmicBinWidthInverse= dble( mergerTreeComputeConditionalMassFunctionParentMassCount  ) &
                  &                               / log(                                                            &
                  &                                      mergerTreeComputeConditionalMassFunctionParentMassMaximum  &
                  &                                     /mergerTreeComputeConditionalMassFunctionParentMassMinimum  &
                  &                                    )
             call Alloc_Array(massParents,[mergerTreeComputeConditionalMassFunctionParentMassCount+1])
             massParents=Make_Range(                                                           &
                  &                 mergerTreeComputeConditionalMassFunctionParentMassMinimum, &
                  &                 mergerTreeComputeConditionalMassFunctionParentMassMaximum, &
                  &                 mergerTreeComputeConditionalMassFunctionParentMassCount  , &
                  &                 rangeType=rangeTypeLogarithmic                             &
                  &                )
             ! Construct bins for mass ratio.
             massRatioLogarithmicMinimum        =  log( mergerTreeComputeConditionalMassFunctionMassRatioMinimum)
             massRatioLogarithmicBinWidthInverse= dble( mergerTreeComputeConditionalMassFunctionMassRatioCount  ) &
                  &                              / log(                                                           &
                  &                                     mergerTreeComputeConditionalMassFunctionMassRatioMaximum  &
                  &                                    /mergerTreeComputeConditionalMassFunctionMassRatioMinimum  &
                  &                                   )
             call Alloc_Array(massRatios,[mergerTreeComputeConditionalMassFunctionMassRatioCount+1])
             massRatios=Make_Range(&
                  &                mergerTreeComputeConditionalMassFunctionMassRatioMinimum, &
                  &                mergerTreeComputeConditionalMassFunctionMassRatioMaximum, &
                  &                mergerTreeComputeConditionalMassFunctionMassRatioCount  , &
                  &                rangeType=rangeTypeLogarithmic                            &
                  &               )
             ! Get the default cosmology functions object.
             cosmologyFunctionsDefault => cosmologyFunctions()
             ! Construct arrays of times for progenitors.
             do i=1,massFunctionTimeCount
                timeProgenitors(i)=&
                     & cosmologyFunctionsDefault%cosmicTime(                            &
                     &  cosmologyFunctionsDefault%expansionFactorFromRedshift(          &
                     &   mergerTreeComputeConditionalMassFunctionProgenitorRedshifts(i) &
                     &  )                                                               &
                     & )
             end do
             ! Construct arrays of times for parents.
             do i=1,massFunctionTimeCount
                timeParents(i)=                                                     & 
                     & cosmologyFunctionsDefault%cosmicTime(                        &
                     &  cosmologyFunctionsDefault%expansionFactorFromRedshift(      &
                     &   mergerTreeComputeConditionalMassFunctionParentRedshifts(i) &
                     &  )                                                           &
                     & )
             end do
             ! Allocate and initialize array to store the condtional mass function.
             call Alloc_Array(                                                          &
                  &           normalization                                           , &
                  &           [                                                         &
                  &            massFunctionTimeCount                                  , &
                  &            mergerTreeComputeConditionalMassFunctionParentMassCount  &
                  &           ]                                                         &
                  &          )
             call Alloc_Array(                                                          &
                  &           conditionalMassFunction                                 , &
                  &           [                                                         &
                  &            massFunctionTimeCount                                  , &
                  &            mergerTreeComputeConditionalMassFunctionParentMassCount, &
                  &            mergerTreeComputeConditionalMassFunctionMassRatioCount   &
                  &           ]                                                         &
                  &          )
             call Alloc_Array(                                                          &
                  &           conditionalMassFunctionError                            , &
                  &           [                                                         &
                  &            massFunctionTimeCount                                  , &
                  &            mergerTreeComputeConditionalMassFunctionParentMassCount, &
                  &            mergerTreeComputeConditionalMassFunctionMassRatioCount   &
                  &           ]                                                         &
                  &          )
             call Alloc_Array(                                                                 &
                  &           primaryProgenitorMassFunction                                  , &
                  &           [                                                                &
                  &            massFunctionTimeCount                                         , &
                  &            mergerTreeComputeConditionalMassFunctionParentMassCount       , &
                  &            mergerTreeComputeConditionalMassFunctionMassRatioCount        , &
                  &            mergerTreeComputeConditionalMassFunctionPrimaryProgenitorDepth  &
                  &           ]                                                                &
                  &          )
             call Alloc_Array(                                                                 &
                  &           primaryProgenitorMassFunctionError                             , &
                  &           [                                                                &
                  &            massFunctionTimeCount                                         , &
                  &            mergerTreeComputeConditionalMassFunctionParentMassCount       , &
                  &            mergerTreeComputeConditionalMassFunctionMassRatioCount        , &
                  &            mergerTreeComputeConditionalMassFunctionPrimaryProgenitorDepth  &
                  &           ]                                                                &
                  &          )
             call Alloc_Array(                                                                 &
                  &           formationRateFunction                                          , &
                  &           [                                                                &
                  &            massFunctionTimeCount                                         , &
                  &            mergerTreeComputeConditionalMassFunctionParentMassCount       , &
                  &            mergerTreeComputeConditionalMassFunctionMassRatioCount        , &
                  &            2                                                               &
                  &           ]                                                                &
                  &          )
             call Alloc_Array(                                                                 &
                  &           formationRateFunctionError                                     , &
                  &           [                                                                &
                  &            massFunctionTimeCount                                         , &
                  &            mergerTreeComputeConditionalMassFunctionParentMassCount       , &
                  &            mergerTreeComputeConditionalMassFunctionMassRatioCount        , &
                  &            2                                                               &
                  &           ]                                                                &
                  &          )
             normalization                     =0.0d0
             conditionalMassFunction           =0.0d0
             conditionalMassFunctionError      =0.0d0
             primaryProgenitorMassFunction     =0.0d0
             primaryProgenitorMassFunctionError=0.0d0
             formationRateFunction             =0.0d0
             formationRateFunctionError        =0.0d0
          end if
          ! Flag that module is initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical (Merger_Tree_Conditional_Mass_Function_Initialize)
    end if

    ! Return immediately if the conditional mass function is not to be computed.
    if (.not.mergerTreeComputeConditionalMassFunction) return

    ! Allocate primary progenitor mass array.
    if (.not.allocated(primaryProgenitorMass))                                               &
         & call Alloc_Array(                                                                 &
         &                  primaryProgenitorMass                                          , &
         &                  [                                                                &
         &                   massFunctionTimeCount                                         , &
         &                   mergerTreeComputeConditionalMassFunctionParentMassCount       , &
         &                   mergerTreeComputeConditionalMassFunctionPrimaryProgenitorDepth  &
         &                  ]                                                                &
         &                 )

    ! Iterate over trees.
    currentTree => thisTree
    do while (associated(currentTree))
       ! Initialize primary progenitor masses to zero.
       primaryProgenitorMass=0.0d0
       ! Get root node of the tree.       
       thisNode => currentTree%baseNode
       ! Walk the tree, pruning hierarchy.
       do while (associated(thisNode))          
          ! Get the child node, and process if child exists.
          childNode => thisNode%firstChild
          do while (associated(childNode))
             ! Get the basic components.
             thisBasicComponent  => thisNode %basic()
             childBasicComponent => childNode%basic()
             ! Determine range of times spanned by this branch.
             branchBegin=childBasicComponent%time()
             branchEnd  = thisBasicComponent%time()
             ! Does the branch span a progenitor node time?
             do i=1,massFunctionTimeCount
                ! Does the branch span a parent node time?
                if     (                                                          &
                     &     branchBegin <= timeParents(i)                          &
                     &  .and.                                                     &
                     &   (                                                        &
                     &     branchEnd   >  timeParents(i)                          &
                     &    .or.                                                    &
                     &     (                                                      &
                     &       .not.associated(thisNode%parent)                     &
                     &      .and.                                                 &
                     &       Values_Agree(branchEnd,timeParents(i),relTol=1.0d-6) &
                     &     )                                                      &
                     &   )                                                        &
                     & ) then
                   ! Get the masses on the branch.
                   branchMassInitial=childBasicComponent%mass()
                   if (childNode%isPrimaryProgenitor()) then
                      branchMassFinal=thisBasicComponent%mass()
                   else
                      branchMassFinal=branchMassInitial
                   end if
                   ! Interpolate to get the mass at the required time.
                   massParent=branchMassInitial+(branchMassFinal-branchMassInitial)*(timeParents(i)-branchBegin)&
                        &/(branchEnd-branchBegin)
                   massParentLogarithmic=log(massParent)
                   binMassParent=int(                                                      &
                        &             (massParentLogarithmic-massParentLogarithmicMinimum) &
                        &            *                                                     &
                        &             massParentLogarithmicBinWidthInverse                 &
                        &           )                                                      &
                        &        +1
                   if    (binMassParent >= 1 .and. binMassParent <= mergerTreeComputeConditionalMassFunctionParentMassCount) then
                      !$omp atomic
                      normalization(i,binMassParent)=normalization(i,binMassParent)+currentTree%volumeWeight
                   end if
                end if
                ! Check if the branch spans the progenitor time.
                if     (                                   &
                     &   branchBegin <= timeProgenitors(i) &
                     &  .and.                              &
                     &   branchEnd   >  timeProgenitors(i) &
                     & ) then
                   ! Get the masses on the branch.
                   branchMassInitial=childBasicComponent%mass()
                   if (childNode%isPrimaryProgenitor()) then
                      branchMassFinal=thisBasicComponent%mass()
                   else
                      branchMassFinal=branchMassInitial
                   end if
                   ! Interpolate to get the mass at the required time.
                   massProgenitor=branchMassInitial+(branchMassFinal-branchMassInitial)*(timeProgenitors(i)-branchBegin)&
                        &/(branchEnd-branchBegin)
                   ! Walk up the tree to find parents.
                   parentNode => thisNode
                   parentWalk : do while (associated(parentNode))
                      ! Get the parent's child.
                      parentsChildNode => parentNode%firstChild
                      ! Get the basic components.
                      parentBasicComponent => parentNode      %basic()
                      parentsChildBasic    => parentsChildNode%basic()
                      ! Determine range of times spanned by this branch.
                      parentBranchBegin=parentsChildBasic   %time()
                      parentBranchEnd  =parentBasicComponent%time()
                      ! Does the branch span a parent node time?
                      if     (                                                                &
                           &   parentBranchBegin <= timeParents(i)                            &
                           &  .and.                                                           &
                           &   (                                                              &
                           &     parentBranchEnd   >  timeParents(i)                          &
                           &    .or.                                                          &
                           &     (                                                            &
                           &       .not.associated(parentNode%parent)                         &
                           &      .and.                                                       &
                           &       Values_Agree(parentBranchEnd,timeParents(i),relTol=1.0d-6) &
                           &     )                                                            &
                           &   )                                                              &
                           & ) then
                         ! Get the masses on the parent branch.
                         parentBranchMassInitial=parentsChildBasic%mass()
                         parentBranchMassFinal  =      parentBasicComponent%mass()
                         ! Find the parent mass at the required time.
                         massParent=parentBranchMassInitial+(parentBranchMassFinal-parentBranchMassInitial)*(timeParents(i)&
                              &-parentBranchBegin)/(parentBranchEnd-parentBranchBegin)
                         ! Accumulate to mass function array.
                         massParentLogarithmic=log(               massParent)
                         massRatio            =    massProgenitor/massParent
                         massRatioLogarithmic =log(               massRatio )
                         binMassParent        =int(                                                      &
                              &                     (massParentLogarithmic-massParentLogarithmicMinimum) &
                              &                    *                                                     &
                              &                     massParentLogarithmicBinWidthInverse                 &
                              &                   )                                                      &
                              &                +1
                         binMassRatio         =int(                                                      &
                              &                     (massRatioLogarithmic-massRatioLogarithmicMinimum)   &
                              &                    *                                                     &
                              &                     massRatioLogarithmicBinWidthInverse                  &
                              &                   )                                                      &
                              &                +1
                         ! Check if within binned ranges.
                         if    (binMassParent >= 1 .and. binMassParent <= mergerTreeComputeConditionalMassFunctionParentMassCount) then
                            if (binMassRatio  >= 1 .and. binMassRatio  <= mergerTreeComputeConditionalMassFunctionMassRatioCount) then
                               !$omp atomic
                               conditionalMassFunction             (i,binMassParent,binMassRatio) &
                                    & =conditionalMassFunction     (i,binMassParent,binMassRatio) &
                                    & + massRatio*currentTree%volumeWeight
                               !$omp atomic
                               conditionalMassFunctionError        (i,binMassParent,binMassRatio) &
                                    & =conditionalMassFunctionError(i,binMassParent,binMassRatio) &
                                    & +(massRatio*currentTree%volumeWeight)**2
                            end if
                            ! Check for formation.
                            if (branchBegin > timeProgenitors(i)*(1.0d0-mergerTreeConditionalMassFunctionFormationRateTimeFraction) .and. .not.associated(childNode%firstChild)) then
                               ! This is a newly formed halo, accumulate to formation rate arrays.
                               ! Find the mass at creation.
                               massRatioCreation           =branchMassInitial/massParent
                               massRatioCreationLogarithmic=log(massRatioCreation)
                               binMassRatioCreation=int(                                                              &
                                    &                      (massRatioCreationLogarithmic-massRatioLogarithmicMinimum) &
                                    &                   *                                                             &
                                    &                    massRatioLogarithmicBinWidthInverse                          &
                                    &                  )                                                              &
                                    &               +1
                               if (binMassRatioCreation >= 1 .and. binMassRatioCreation <= mergerTreeComputeConditionalMassFunctionMassRatioCount) then
                                  !$omp atomic
                                  formationRateFunction     (i,binMassParent,binMassRatioCreation,1)=formationRateFunction     (i,binMassParent,binMassRatioCreation,1)+ massRatioCreation*currentTree%volumeWeight
                                  !$omp atomic
                                  formationRateFunctionError(i,binMassParent,binMassRatioCreation,1)=formationRateFunctionError(i,binMassParent,binMassRatioCreation,1)+(massRatioCreation*currentTree%volumeWeight)**2
                               end if
                               ! Find the mass of this node just prior to it becoming a subhalo.
                               descendentNode => thisNode
                               do while (associated(descendentNode%parent).and.associated(descendentNode%parent%firstChild,descendentNode))
                                  descendentNode => descendentNode%parent
                               end do
                               descendentBasic => descendentNode%basic()
                               massRatioDescendent=descendentBasic%mass()/massParent
                               massRatioDescendentLogarithmic=log(massRatioDescendent)
                               binMassRatioDescendent=int(                                                              &
                                    &                      (massRatioDescendentLogarithmic-massRatioLogarithmicMinimum) &
                                    &                     *                                                             &
                                    &                      massRatioLogarithmicBinWidthInverse                          &
                                    &                    )                                                              &
                                    &                 +1
                               if (binMassRatioDescendent >= 1 .and. binMassRatioDescendent <= mergerTreeComputeConditionalMassFunctionMassRatioCount) then
                                  !$omp atomic
                                  formationRateFunction     (i,binMassParent,binMassRatioDescendent,2)=formationRateFunction     (i,binMassParent,binMassRatioDescendent,2)+ massRatioDescendent*currentTree%volumeWeight
                                  !$omp atomic
                                  formationRateFunctionError(i,binMassParent,binMassRatioDescendent,2)=formationRateFunctionError(i,binMassParent,binMassRatioDescendent,2)+(massRatioDescendent*currentTree%volumeWeight)**2
                               end if
                            end if
                            ! Accumulate to the primary progenitor mass array if necessary.
                            iPrimary=1
                            do while (massRatio < primaryProgenitorMass(i,binMassParent,iPrimary))
                               iPrimary=iPrimary+1
                               if (iPrimary > mergerTreeComputeConditionalMassFunctionPrimaryProgenitorDepth) exit
                            end do
                            if (iPrimary <= mergerTreeComputeConditionalMassFunctionPrimaryProgenitorDepth) then
                               if (iPrimary < mergerTreeComputeConditionalMassFunctionPrimaryProgenitorDepth) then
                                  do jPrimary=mergerTreeComputeConditionalMassFunctionPrimaryProgenitorDepth,iPrimary+1,-1
                                     primaryProgenitorMass        (i,binMassParent,jPrimary  ) &
                                          & =primaryProgenitorMass(i,binMassParent,jPrimary-1)
                                  end do
                               end if
                               primaryProgenitorMass(i,binMassParent,iPrimary)=massRatio
                            end if
                         end if
                      end if
                      parentNode => parentNode%parent
                   end do parentWalk
                end if
                ! Record the mass of the branch at the parent time.
             end do
             ! Move to the next child.
             childNode => childNode%sibling
          end do
          ! Move to the next node.
          call thisNode%walkTree(thisNode)
       end do
       ! Store the computed primary progenitor mass functions.
       do i=1,massFunctionTimeCount
          do binMassParent=1,mergerTreeComputeConditionalMassFunctionParentMassCount
             do iPrimary=1,mergerTreeComputeConditionalMassFunctionPrimaryProgenitorDepth
                if (primaryProgenitorMass(i,binMassParent,iPrimary) > 0.0d0) then
                   massRatioLogarithmic=log(primaryProgenitorMass(i,binMassParent,iPrimary))
                   binMassRatio         =int(                                                      &
                        &                     (massRatioLogarithmic-massRatioLogarithmicMinimum)   &
                        &                    *                                                     &
                        &                     massRatioLogarithmicBinWidthInverse                  &
                        &                   )                                                      &
                        &                +1
                   if (binMassRatio  >= 1 .and. binMassRatio  <= mergerTreeComputeConditionalMassFunctionMassRatioCount) then
                      !$omp atomic
                      primaryProgenitorMassFunction             (i,binMassParent,binMassRatio,iPrimary)= &
                           &  primaryProgenitorMassFunction     (i,binMassParent,binMassRatio,iPrimary)  &
                           & +primaryProgenitorMass             (i,binMassParent,             iPrimary)  &
                           & *currentTree%volumeWeight
                      !$omp atomic
                      primaryProgenitorMassFunctionError          (i,binMassParent,binMassRatio,iPrimary)= &
                           &    primaryProgenitorMassFunctionError(i,binMassParent,binMassRatio,iPrimary)  &
                           & +(                                                                            &
                           &    primaryProgenitorMass             (i,binMassParent,             iPrimary)  &
                           &   *currentTree%volumeWeight                                                   &
                           &  )**2
                   end if
                end if
             end do
          end do
       end do
       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do
    return
  end subroutine Merger_Tree_Conditional_Mass_Function

  !# <hdfPreCloseTask>
  !#  <unitName>Merger_Tree_Conditional_Mass_Function_Output</unitName>
  !# </hdfPreCloseTask>
  subroutine Merger_Tree_Conditional_Mass_Function_Output
    !% Outputs conditional mass function.
    use ISO_Varying_String
    use IO_HDF5
    use Galacticus_HDF5
    use Numerical_Constants_Astronomical
    implicit none
    type   (hdf5Object) :: conditionalMassFunctionGroup,massDataset
    integer             :: i,j,iPrimary

    ! Return immediately if the conditional mass function is not to be computed.
    if (.not.mergerTreeComputeConditionalMassFunction) return
    ! Normalize the conditional mass functions.
    normalization=normalization/massRatioLogarithmicBinWidthInverse/log(10.0d0)
    do i=1,massFunctionTimeCount
       do j=1,mergerTreeComputeConditionalMassFunctionParentMassCount
          if (normalization(i,j) > 0.0d0) then
             conditionalMassFunction     (i,j,:  )=     conditionalMassFunction     (i,j,:  ) /normalization(i,j)
             conditionalMassFunctionError(i,j,:  )=sqrt(conditionalMassFunctionError(i,j,:  ))/normalization(i,j)
             formationRateFunction       (i,j,:,:)=     formationRateFunction       (i,j,:,:) /normalization(i,j)
             formationRateFunctionError  (i,j,:,:)=sqrt(formationRateFunctionError  (i,j,:,:))/normalization(i,j)
             do iPrimary=1,mergerTreeComputeConditionalMassFunctionPrimaryProgenitorDepth
                primaryProgenitorMassFunction     (i,j,:,iPrimary)=     primaryProgenitorMassFunction     (i,j,:,iPrimary) /normalization(i,j)
                primaryProgenitorMassFunctionError(i,j,:,iPrimary)=sqrt(primaryProgenitorMassFunctionError(i,j,:,iPrimary))/normalization(i,j)
             end do
          end if
       end do
    end do
    ! Output the data.
    conditionalMassFunctionGroup=galacticusOutputFile%openGroup('conditionalMassFunction','Conditional mass functions of merger trees.')
    call conditionalMassFunctionGroup%writeDataset(massParents,"massParent","Mass of parent node [Msolar]",datasetReturned=massDataset)
    call massDataset%writeAttribute(massSolar,"unitsInSI")
    call massDataset%close()
    call conditionalMassFunctionGroup%writeDataset(massRatios ,"massRatio" ,"Mass of ratio node [Msolar]" ,datasetReturned=massDataset)
    call massDataset%writeAttribute(massSolar,"unitsInSI")
    call massDataset%close()
    call conditionalMassFunctionGroup%writeDataset(mergerTreeComputeConditionalMassFunctionParentRedshifts    ,"redshiftParent"                    ,"Redshift of parent node []"                )
    call conditionalMassFunctionGroup%writeDataset(mergerTreeComputeConditionalMassFunctionProgenitorRedshifts,"redshiftProgenitor"                ,"Redshift of progenitor node []"            )
    call conditionalMassFunctionGroup%writeDataset(conditionalMassFunction                                    ,"conditionalMassFunction"           ,"Conditional mass functions []"             )
    call conditionalMassFunctionGroup%writeDataset(conditionalMassFunctionError                               ,"conditionalMassFunctionError"      ,"Conditional mass function errors []"       )
    call conditionalMassFunctionGroup%writeDataset(primaryProgenitorMassFunction                              ,"primaryProgenitorMassFunction"     ,"Primary progenitor mass functions []"      )
    call conditionalMassFunctionGroup%writeDataset(primaryProgenitorMassFunctionError                         ,"primaryProgenitorMassFunctionError","Primary progenitor mass function errors []")
    call conditionalMassFunctionGroup%writeDataset(formationRateFunction                                      ,"formationRateFunction"             ,"Formation rate functions []"               )
    call conditionalMassFunctionGroup%writeDataset(formationRateFunctionError                                 ,"formationRateFunctionError"        ,"Formation rate function errors []"         )
    call conditionalMassFunctionGroup%close()    
    return
  end subroutine Merger_Tree_Conditional_Mass_Function_Output

end module Merger_Trees_Conditional_Mass_Function
