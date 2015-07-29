!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Contains a module which implements a merger tree operator which accumulates conditional mass functions for trees.

  !# <mergerTreeOperator name="mergerTreeOperatorConditionalMF" defaultThreadPrivate="no">
  !#  <description>Provides a merger tree operator which accumulates conditional mass functions for trees.</description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorConditionalMF
     !% A merger tree operator class which accumulates conditional mass functions for trees.
     private
     double precision                , allocatable, dimension(:      ) :: timeParents                         , timeProgenitors                    , &
          &                                                               parentRedshifts                     , progenitorRedshifts                , &
          &                                                               massParents                         , massRatios
     double precision                , allocatable, dimension(:,:    ) :: normalization
     double precision                , allocatable, dimension(:,:,:  ) :: conditionalMassFunction             , conditionalMassFunctionError
     double precision                , allocatable, dimension(:,:,:,:) :: primaryProgenitorMassFunction       , primaryProgenitorMassFunctionError
     double precision                , allocatable, dimension(:,:,:,:) :: formationRateFunction               , formationRateFunctionError
     integer                                                           :: parentMassCount                     , primaryProgenitorDepth             , &
          &                                                               massRatioCount                      , timeCount
     double precision                                                  :: massParentLogarithmicMinimum        , massRatioLogarithmicMinimum        , &
          &                                                               massParentLogarithmicBinWidthInverse, massRatioLogarithmicBinWidthInverse, &
          &                                                               formationRateTimeFraction
     type            (varying_string)                                  :: outputGroupName
   contains
     final     ::             conditionalMFDestructor
     procedure :: operate  => conditionalMFOperate
     procedure :: finalize => conditionalMFFinalize
  end type mergerTreeOperatorConditionalMF
  
  interface mergerTreeOperatorConditionalMF
     !% Constructors for the conditional mass function merger tree operator class.
     module procedure conditionalMFConstructorParameters
     module procedure conditionalMFConstructorInternal
  end interface mergerTreeOperatorConditionalMF

contains

  function conditionalMFConstructorParameters(parameters)
    !% Constructor for the conditional mass function merger tree operator class which takes a parameter set as input.
    use Memory_Management
    implicit none
    type            (mergerTreeOperatorConditionalMF)                            :: conditionalMFConstructorParameters
    type            (inputParameters                ), intent(in   )             :: parameters
    double precision                                 , allocatable, dimension(:) :: progenitorRedshifts               , parentRedshifts
    integer                                                                      :: parentMassCount                   , massRatioCount   , &
         &                                                                          primaryProgenitorDepth
    double precision                                                             :: parentMassMinimum                 , parentMassMaximum, &
         &                                                                          massRatioMinimum                  , massRatioMaximum , &
         &                                                                          formationRateTimeFraction
    type            (varying_string                 )                            :: outputGroupName
    !# <inputParameterList label="allowedParameterNames" />

    !# <inputParameter>
    !#   <name>parentMassCount</name>
    !#   <source>parameters</source>
    !#   <defaultValue>10</defaultValue>
    !#   <description>The number of bins in parent mass when constructing conditional halo mass functions.</description>
    !#   <type>integer</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>parentMassMinimum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d10</defaultValue>
    !#   <description>The minimum parent halo mass to bin when constructing conditional halo mass functions.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>parentMassMaximum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d15</defaultValue>
    !#   <description>The maximum parent halo mass to bin when constructing conditional halo mass functions.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massRatioCount</name>
    !#   <source>parameters</source>
    !#   <defaultValue>10</defaultValue>
    !#   <description>The number of bins in mass ratio when constructing conditional halo mass functions.</description>
    !#   <type>integer</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massRatioMinimum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d-4</defaultValue>
    !#   <description>The minimum mass ratio to bin when constructing conditional halo mass functions.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massRatioMaximum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d1</defaultValue>
    !#   <description>The maximum mass ratio to bin when constructing conditional halo mass functions.</description>
    !#   <type>logical</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    call Alloc_Array(parentRedshifts    ,[max(1,parameters%count('parentRedshifts'    ,zeroIfNotPresent=.true.))])
    !# <inputParameter>
    !#   <name>parentRedshifts</name>
    !#   <source>parameters</source>
    !#   <defaultValue>[0.0d0]</defaultValue>
    !#   <description>The set of parent halo redshifts to use when constructing conditional halo mass functions.</description>
    !#   <type>real</type>
    !#   <cardinality>1..</cardinality>
    !# </inputParameter>
    call Alloc_Array(progenitorRedshifts,[max(1,parameters%count('progenitorRedshifts',zeroIfNotPresent=.true.))])
    !# <inputParameter>
    !#   <name>progenitorRedshifts</name>
    !#   <source>parameters</source>
    !#   <defaultValue>[1.0d0]</defaultValue>
    !#   <description>The set of progenitor halo redshifts to use when constructing conditional halo mass functions.</description>
    !#   <type>real</type>
    !#   <cardinality>1..</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>primaryProgenitorDepth</name>
    !#   <source>parameters</source>
    !#   <defaultValue>2</defaultValue>
    !#   <description>The depth in progenitor ranking for which to store ranked progenitor mass functions. For example, a value of 2 means store mass functions for the most massive, and second most massive progenitor.</description>
    !#   <type>integer</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>formationRateTimeFraction</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.01d0</defaultValue>
    !#   <description>The fraction of the current time over which to estimate the formation rate of halos when computing merger tree statistics.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>outputGroupName</name>
    !#   <source>parameters</source>
    !#   <defaultValue>var_str('conditionalMassFunction')</defaultValue>
    !#   <description>The name of the HDF5 group to which the conditional mass function should be output.</description>
    !#   <type>string</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    ! Construct the instance.
    conditionalMFConstructorParameters=conditionalMFConstructorInternal(                           &
         &                                                              parentMassCount          , &
         &                                                              parentMassMinimum        , &
         &                                                              parentMassMaximum        , &
         &                                                              massRatioCount           , &
         &                                                              massRatioMinimum         , &
         &                                                              massRatioMaximum         , &
         &                                                              parentRedshifts          , &
         &                                                              progenitorRedshifts      , &
         &                                                              primaryProgenitorDepth   , &
         &                                                              formationRateTimeFraction, &
         &                                                              outputGroupName            &
         &                                                             )
    return
  end function conditionalMFConstructorParameters

  function conditionalMFConstructorInternal(parentMassCount,parentMassMinimum,parentMassMaximum,massRatioCount,massRatioMinimum,massRatioMaximum,parentRedshifts,progenitorRedshifts,primaryProgenitorDepth,formationRateTimeFraction,outputGroupName)
    !% Internal constructor for the conditional mass function merger tree operator class.
    use Cosmology_Functions
    use Numerical_Ranges
    use Memory_Management
    use Galacticus_Error
    implicit none
    type            (mergerTreeOperatorConditionalMF)                              :: conditionalMFConstructorInternal
    double precision                                 , intent(in   ), dimension(:) :: parentRedshifts                 , progenitorRedshifts
    integer                                          , intent(in   )               :: parentMassCount                 , massRatioCount     , &
         &                                                                            primaryProgenitorDepth
    double precision                                 , intent(in   )               :: parentMassMinimum               , parentMassMaximum  , &
         &                                                                            massRatioMinimum                , massRatioMaximum   , &
         &                                                                            formationRateTimeFraction
    type            (varying_string                 ), intent(in   )               :: outputGroupName
    class           (cosmologyFunctionsClass        ), pointer                     :: cosmologyFunctions_
    integer                                                                        :: i
    
    ! Store array sizes.
    conditionalMFConstructorInternal%timeCount             =size(parentRedshifts)
    conditionalMFConstructorInternal%parentMassCount       =parentMassCount
    conditionalMFConstructorInternal%massRatioCount        =massRatioCount
    conditionalMFConstructorInternal%primaryProgenitorDepth=primaryProgenitorDepth
    conditionalMFConstructorInternal%outputGroupName       =outputGroupName
    if (size(progenitorRedshifts) /= conditionalMFConstructorInternal%timeCount) &
         & call Galacticus_Error_Report('conditionalMFConstructorInternal','mismatch in sizes of parent and progenitor redshift arrays')
    ! Allocate arrays.
    call Alloc_Array(conditionalMFConstructorInternal%parentRedshifts    ,[conditionalMFConstructorInternal%timeCount        ])
    call Alloc_Array(conditionalMFConstructorInternal%progenitorRedshifts,[conditionalMFConstructorInternal%timeCount        ])
    call Alloc_Array(conditionalMFConstructorInternal%timeProgenitors    ,[conditionalMFConstructorInternal%timeCount        ])
    call Alloc_Array(conditionalMFConstructorInternal%timeParents        ,[conditionalMFConstructorInternal%timeCount        ])
    call Alloc_Array(conditionalMFConstructorInternal%massParents        ,[conditionalMFConstructorInternal%parentMassCount+1])
    call Alloc_Array(conditionalMFConstructorInternal%massRatios         ,[conditionalMFConstructorInternal%massRatioCount +1])
    call Alloc_Array(                                                                      &
         &            conditionalMFConstructorInternal%normalization                     , &
         &           [                                                                     &
         &            conditionalMFConstructorInternal%timeCount                         , &
         &            conditionalMFConstructorInternal%parentMassCount                     &
         &           ]                                                                     &
         &          )
    call Alloc_Array(                                                                      &
         &            conditionalMFConstructorInternal%conditionalMassFunction           , &
         &           [                                                                     &
         &            conditionalMFConstructorInternal%timeCount                         , &
         &            conditionalMFConstructorInternal%parentMassCount                   , &
         &            conditionalMFConstructorInternal%massRatioCount                      &
         &           ]                                                                     &
         &          )
    call Alloc_Array(                                                                      &
         &            conditionalMFConstructorInternal%conditionalMassFunctionError      , &
         &           [                                                                     &
         &            conditionalMFConstructorInternal%timeCount                         , &
         &            conditionalMFConstructorInternal%parentMassCount                   , &
         &            conditionalMFConstructorInternal%massRatioCount                      &
         &           ]                                                                     &
         &          )
    call Alloc_Array(                                                                      &
         &            conditionalMFConstructorInternal%primaryProgenitorMassFunction     , &
         &           [                                                                     &
         &            conditionalMFConstructorInternal%timeCount                         , &
         &            conditionalMFConstructorInternal%parentMassCount                   , &
         &            conditionalMFConstructorInternal%massRatioCount                    , &
         &            conditionalMFConstructorInternal%primaryProgenitorDepth              &
         &           ]                                                                     &
         &          )
    call Alloc_Array(                                                                      &
         &            conditionalMFConstructorInternal%primaryProgenitorMassFunctionError, &
         &           [                                                                     &
         &            conditionalMFConstructorInternal%timeCount                         , &
         &            conditionalMFConstructorInternal%parentMassCount                   , &
         &            conditionalMFConstructorInternal%massRatioCount                    , &
         &            conditionalMFConstructorInternal%primaryProgenitorDepth              &
         &           ]                                                                     &
         &          )
    call Alloc_Array(                                                                      &
         &            conditionalMFConstructorInternal%formationRateFunction             , &
         &           [                                                                     &
         &            conditionalMFConstructorInternal%timeCount                         , &
         &            conditionalMFConstructorInternal%parentMassCount                   , &
         &            conditionalMFConstructorInternal%massRatioCount                    , &
         &            2                                                                    &
         &           ]                                                                     &
         &          )
    call Alloc_Array(                                                                      &
         &            conditionalMFConstructorInternal%formationRateFunctionError        , &
         &           [                                                                     &
         &            conditionalMFConstructorInternal%timeCount                         , &
         &            conditionalMFConstructorInternal%parentMassCount                   , &
         &            conditionalMFConstructorInternal%massRatioCount                    , &
         &            2                                                                    &
         &           ]                                                                     &
         &          )
 
    ! Construct bins for parent node mass.
    conditionalMFConstructorInternal%massParentLogarithmicMinimum        =  log( parentMassMinimum)
    conditionalMFConstructorInternal%massParentLogarithmicBinWidthInverse= dble( parentMassCount  ) &
         &                                                                / log(                    &
         &                                                                      +parentMassMaximum  &
         &                                                                      /parentMassMinimum  &
         &                                                                     )
    conditionalMFConstructorInternal%massParents=Make_Range(                                        &
         &                                                  parentMassMinimum                     , &
         &                                                  parentMassMaximum                     , &
         &                                                  parentMassCount                       , &
         &                                                  rangeType        =rangeTypeLogarithmic  &
         &                                                 )
    ! Construct bins for mass ratio.
    conditionalMFConstructorInternal%massRatioLogarithmicMinimum        =  log( massRatioMinimum)
    conditionalMFConstructorInternal%massRatioLogarithmicBinWidthInverse= dble( massRatioCount  ) &
         &                                                               / log(                   &
         &                                                                     +massRatioMaximum  &
         &                                                                     /massRatioMinimum  &
         &                                                                    )
    conditionalMFConstructorInternal%massRatios=Make_Range(                                       &
         &                                                 massRatioMinimum                     , &
         &                                                 massRatioMaximum                     , &
         &                                                 massRatioCount                       , &
         &                                                 rangeType       =rangeTypeLogarithmic  &
         &                                                )
    ! Get the default cosmology functions object.
    cosmologyFunctions_ => cosmologyFunctions()
    ! Construct arrays of times for progenitors.
    do i=1,conditionalMFConstructorInternal%timeCount
       conditionalMFConstructorInternal%timeProgenitors(i)=              &
                     & cosmologyFunctions_%cosmicTime(                   &
                     &  cosmologyFunctions_%expansionFactorFromRedshift( &
                     &   progenitorRedshifts(i)                          &
                     &  )                                                &
                     & )
    end do
    ! Construct arrays of times for parents.
    do i=1,conditionalMFConstructorInternal%timeCount
       conditionalMFConstructorInternal%timeParents(i)=         & 
            & cosmologyFunctions_%cosmicTime(                   &
            &  cosmologyFunctions_%expansionFactorFromRedshift( &
            &   parentRedshifts(i)                              &
            &  )                                                &
            & )
    end do
    ! Initialize mass function arrays.
    conditionalMFConstructorInternal%normalization                     =0.0d0
    conditionalMFConstructorInternal%conditionalMassFunction           =0.0d0
    conditionalMFConstructorInternal%conditionalMassFunctionError      =0.0d0
    conditionalMFConstructorInternal%primaryProgenitorMassFunction     =0.0d0
    conditionalMFConstructorInternal%primaryProgenitorMassFunctionError=0.0d0
    conditionalMFConstructorInternal%formationRateFunction             =0.0d0
    conditionalMFConstructorInternal%formationRateFunctionError        =0.0d0
    return
  end function conditionalMFConstructorInternal

  elemental subroutine conditionalMFDestructor(self)
    !% Destructor for the merger tree operator function class.
    implicit none
    type(mergerTreeOperatorConditionalMF), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine conditionalMFDestructor

  subroutine conditionalMFOperate(self,tree)
    !% Compute conditional mass function on {\normalfont \ttfamily tree}.
    use Galacticus_Nodes
    use Input_Parameters
    use Memory_Management
    use Numerical_Comparison
    use Galacticus_Error
    implicit none
    class           (mergerTreeOperatorConditionalMF), intent(inout)                                  :: self
    type            (mergerTree                     ), intent(inout)                         , target :: tree
    type            (treeNode                       ), pointer                                        :: node                          , nodeChild            , &
         &                                                                                               nodeParent                    , nodeParentChild      , &
         &                                                                                               descendentNode
    type            (mergerTree                     ), pointer                                        :: treeCurrent
    class           (nodeComponentBasic             ), pointer                                        :: basic                         , basicChild           , &
         &                                                                                               basicParent                   , descendentBasic      , &
         &                                                                                               basicParentChild
    integer                                                                                           :: i                             , binMassParent        , &
         &                                                                                               binMassRatio                  , iPrimary             , &
         &                                                                                               jPrimary                      , binMassRatioCreation , &
         &                                                                                               binMassRatioDescendent
    double precision                                                                                  :: branchBegin                   , branchEnd            , &
         &                                                                                               parentBranchBegin             , parentBranchEnd      , &
         &                                                                                               massProgenitor                , massParent           , &
         &                                                                                               branchMassInitial             , branchMassFinal      , &
         &                                                                                               parentBranchMassInitial       , parentBranchMassFinal, &
         &                                                                                               massRatioLogarithmic          , massParentLogarithmic, &
         &                                                                                               massRatio                     , massRatioCreation    , &
         &                                                                                               massRatioCreationLogarithmic  , massRatioDescendent  , &
         &                                                                                               massRatioDescendentLogarithmic
    double precision                                  , dimension(                                                                                              &
         &                                                        self%timeCount             ,                                                                  &
         &                                                        self%parentMassCount       ,                                                                  &
         &                                                        self%primaryProgenitorDepth                                                                   &
         &                                                       )                                    :: primaryProgenitorMass
    
   ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Initialize primary progenitor masses to zero.
       primaryProgenitorMass=0.0d0
       ! Get root node of the tree.       
       node => treeCurrent%baseNode
       ! Walk the tree, accumulating statistics.
       do while (associated(node))          
          ! Get the child node, and process if child exists.
          nodeChild => node%firstChild
          do while (associated(nodeChild))
             ! Get the basic components.
             basic      => node     %basic()
             basicChild => nodeChild%basic()
             ! Determine range of times spanned by this branch.
             branchBegin=basicChild%time()
             branchEnd  =basic     %time()
             ! Does the branch span a progenitor node time?
             do i=1,self%timeCount
                ! Does the branch span a parent node time?
                if     (                                                               &
                     &     branchBegin <= self%timeParents(i)                          &
                     &  .and.                                                          &
                     &   (                                                             &
                     &     branchEnd   >  self%timeParents(i)                          &
                     &    .or.                                                         &
                     &     (                                                           &
                     &       .not.associated(node%parent)                              &
                     &      .and.                                                      &
                     &       Values_Agree(branchEnd,self%timeParents(i),relTol=1.0d-6) &
                     &     )                                                           &
                     &   )                                                             &
                     & ) then
                   ! Get the masses on the branch.
                   branchMassInitial=basicChild%mass()
                   if (nodeChild%isPrimaryProgenitor()) then
                      branchMassFinal=basic%mass()
                   else
                      branchMassFinal=branchMassInitial
                   end if
                   ! Interpolate to get the mass at the required time.
                   massParent=                     +branchMassInitial  &
                        &     +(branchMassFinal    -branchMassInitial) &
                        &     *(self%timeParents(i)-branchBegin      ) &
                        &     /(branchEnd          -branchBegin      )
                   massParentLogarithmic=log(massParent)
                   binMassParent=int(                                                            &
                        &            +(+massParentLogarithmic-self%massParentLogarithmicMinimum) &
                        &            *self%massParentLogarithmicBinWidthInverse                  &
                        &           )                                                            &
                        &        +1
                   if (binMassParent >= 1 .and. binMassParent <= self%parentMassCount) then
                      !$omp critical(conditionalMassFunctionAccumulate)
                      self%normalization(i,binMassParent)=self%normalization(i,binMassParent)+treeCurrent%volumeWeight
                      !$omp end critical(conditionalMassFunctionAccumulate)
                   end if
                end if
                ! Check if the branch spans the progenitor time.
                if     (                                        &
                     &   branchBegin <= self%timeProgenitors(i) &
                     &  .and.                                   &
                     &   branchEnd   >  self%timeProgenitors(i) &
                     & ) then
                   ! Get the masses on the branch.
                   branchMassInitial=basicChild%mass()
                   if (nodeChild%isPrimaryProgenitor()) then
                      branchMassFinal=basic%mass()
                   else
                      branchMassFinal=branchMassInitial
                   end if
                   ! Interpolate to get the mass at the required time.
                   massProgenitor=                         +branchMassInitial  &
                        &         +(branchMassFinal        -branchMassInitial) &
                        &         *(self%timeProgenitors(i)-branchBegin      ) &
                        &         /(branchEnd              -branchBegin      )
                   ! Walk up the tree to find parents.
                   nodeParent => node
                   parentWalk : do while (associated(nodeParent))
                      ! Get the parent's child.
                      nodeParentChild   => nodeParent      %firstChild
                      ! Get the basic components.
                      basicParent       => nodeParent      %basic     ()
                      basicParentChild  => nodeParentChild %basic     ()
                      ! Determine range of times spanned by this branch.
                      parentBranchBegin =  basicParentChild%time      ()
                      parentBranchEnd   =  basicParent     %time      ()
                      ! Does the branch span a parent node time?
                      if     (                                                                     &
                           &   parentBranchBegin <= self%timeParents(i)                            &
                           &  .and.                                                                &
                           &   (                                                                   &
                           &     parentBranchEnd >  self%timeParents(i)                            &
                           &    .or.                                                               &
                           &     (                                                                 &
                           &       .not.associated(nodeParent%parent)                              &
                           &      .and.                                                            &
                           &       Values_Agree(parentBranchEnd,self%timeParents(i),relTol=1.0d-6) &
                           &     )                                                                 &
                           &   )                                                                   &
                           & ) then
                         ! Get the masses on the parent branch.
                         parentBranchMassInitial=basicParentChild%mass()
                         parentBranchMassFinal  =     basicParent%mass()
                         ! Find the parent mass at the required time.
                         massParent=                       +parentBranchMassInitial  &
                              &     +(parentBranchMassFinal-parentBranchMassInitial) &
                              &     *(self%timeParents(i)  -parentBranchBegin      ) &
                              &     /(parentBranchEnd      -parentBranchBegin      )
                         ! Accumulate to mass function array.
                         massParentLogarithmic=log(               massParent)
                         massRatio            =    massProgenitor/massParent
                         massRatioLogarithmic =log(               massRatio )
                         binMassParent        =int(                                                           &
                              &                     (massParentLogarithmic-self%massParentLogarithmicMinimum) &
                              &                    *self%massParentLogarithmicBinWidthInverse                 &
                              &                   )                                                           &
                              &                +1
                         binMassRatio         =int(                                                           &
                              &                     (massRatioLogarithmic -self%massRatioLogarithmicMinimum ) &
                              &                    *self%massRatioLogarithmicBinWidthInverse                  &
                              &                   )                                                           &
                              &                +1
                         ! Check if within binned ranges.
                         if    (binMassParent >= 1 .and. binMassParent <= self%parentMassCount) then
                            if (binMassRatio  >= 1 .and. binMassRatio  <= self%massRatioCount) then
                               !$omp critical(conditionalMassFunctionAccumulate)
                               self%conditionalMassFunction             (i,binMassParent,binMassRatio) &
                                    & =self%conditionalMassFunction     (i,binMassParent,binMassRatio) &
                                    &   +massRatio                                                     &
                                    &   *treeCurrent%volumeWeight
                               self%conditionalMassFunctionError        (i,binMassParent,binMassRatio) &
                                    & =self%conditionalMassFunctionError(i,binMassParent,binMassRatio) &
                                    & +(                                                               &
                                    &   +massRatio                                                     &
                                    &   *treeCurrent%volumeWeight                                      &
                                    &  )**2
                               !$omp end critical(conditionalMassFunctionAccumulate)
                            end if
                            ! Check for formation.
                            if (branchBegin > self%timeProgenitors(i)*(1.0d0-self%formationRateTimeFraction) .and. .not.associated(nodeChild%firstChild)) then
                               ! This is a newly formed halo, accumulate to formation rate arrays.
                               ! Find the mass at creation.
                               massRatioCreation           =branchMassInitial/massParent
                               massRatioCreationLogarithmic=log(massRatioCreation)
                               binMassRatioCreation=int(                                                                   &
                                    &                      (massRatioCreationLogarithmic-self%massRatioLogarithmicMinimum) &
                                    &                   *self%massRatioLogarithmicBinWidthInverse                          &
                                    &                  )                                                                   &
                                    &               +1
                               if (binMassRatioCreation >= 1 .and. binMassRatioCreation <= self%massRatioCount) then
                                  !$omp critical(conditionalMassFunctionAccumulate)
                                  self%formationRateFunction     (i,binMassParent,binMassRatioCreation,1)=self%formationRateFunction     (i,binMassParent,binMassRatioCreation,1)+ massRatioCreation*treeCurrent%volumeWeight
                                  self%formationRateFunctionError(i,binMassParent,binMassRatioCreation,1)=self%formationRateFunctionError(i,binMassParent,binMassRatioCreation,1)+(massRatioCreation*treeCurrent%volumeWeight)**2
                                  !$omp end critical(conditionalMassFunctionAccumulate)
                               end if
                               ! Find the mass of this node just prior to it becoming a subhalo.
                               descendentNode => node
                               do while (associated(descendentNode%parent).and.associated(descendentNode%parent%firstChild,descendentNode))
                                  descendentNode => descendentNode%parent
                               end do
                               descendentBasic     => descendentNode %basic()
                               massRatioDescendent =  descendentBasic%mass ()/massParent
                               massRatioDescendentLogarithmic=log(massRatioDescendent)
                               binMassRatioDescendent=int(                                                                   &
                                    &                      (massRatioDescendentLogarithmic-self%massRatioLogarithmicMinimum) &
                                    &                     *self%massRatioLogarithmicBinWidthInverse                          &
                                    &                    )                                                                   &
                                    &                 +1
                               if (binMassRatioDescendent >= 1 .and. binMassRatioDescendent <= self%massRatioCount) then
                                  !$omp critical(conditionalMassFunctionAccumulate)
                                  self%formationRateFunction     (i,binMassParent,binMassRatioDescendent,2)=self%formationRateFunction     (i,binMassParent,binMassRatioDescendent,2)+ massRatioDescendent*treeCurrent%volumeWeight
                                  self%formationRateFunctionError(i,binMassParent,binMassRatioDescendent,2)=self%formationRateFunctionError(i,binMassParent,binMassRatioDescendent,2)+(massRatioDescendent*treeCurrent%volumeWeight)**2
                                  !$omp end critical(conditionalMassFunctionAccumulate)
                               end if
                            end if
                            ! Accumulate to the primary progenitor mass array if necessary.
                            iPrimary=1
                            do while (massRatio < primaryProgenitorMass(i,binMassParent,iPrimary))
                               iPrimary=iPrimary+1
                               if (iPrimary > self%primaryProgenitorDepth) exit
                            end do
                            if (iPrimary <= self%primaryProgenitorDepth) then
                               if (iPrimary < self%primaryProgenitorDepth) then
                                  do jPrimary=self%primaryProgenitorDepth,iPrimary+1,-1
                                     primaryProgenitorMass        (i,binMassParent,jPrimary  ) &
                                          & =primaryProgenitorMass(i,binMassParent,jPrimary-1)
                                  end do
                               end if
                               primaryProgenitorMass(i,binMassParent,iPrimary)=massRatio
                            end if
                         end if
                      end if
                      nodeParent => nodeParent%parent
                   end do parentWalk
                end if
                ! Record the mass of the branch at the parent time.
             end do
             ! Move to the next child.
             nodeChild => nodeChild%sibling
          end do
          ! Move to the next node.
          call node%walkTree(node)
       end do
       ! Store the computed primary progenitor mass functions.
       do i=1,self%timeCount
          do binMassParent=1,self%parentMassCount
             do iPrimary=1,self%primaryProgenitorDepth
                if (primaryProgenitorMass(i,binMassParent,iPrimary) > 0.0d0) then
                   massRatioLogarithmic=log(primaryProgenitorMass(i,binMassParent,iPrimary))
                   binMassRatio         =int(                                                           &
                        &                     (massRatioLogarithmic-self%massRatioLogarithmicMinimum)   &
                        &                    *                                                          &
                        &                     self%massRatioLogarithmicBinWidthInverse                  &
                        &                   )                                                           &
                        &                +1
                   if (binMassRatio  >= 1 .and. binMassRatio  <= self%massRatioCount) then
                      !$omp critical(conditionalMassFunctionAccumulate)
                      self          %primaryProgenitorMassFunction     (i,binMassParent,binMassRatio,iPrimary)= &
                           & +  self%primaryProgenitorMassFunction     (i,binMassParent,binMassRatio,iPrimary)  &
                           & +       primaryProgenitorMass             (i,binMassParent,             iPrimary)  &
                           & *treeCurrent%volumeWeight
                      self          %primaryProgenitorMassFunctionError(i,binMassParent,binMassRatio,iPrimary)= &
                           &    self%primaryProgenitorMassFunctionError(i,binMassParent,binMassRatio,iPrimary)  &
                           & +(                                                                                 &
                           &   +     primaryProgenitorMass             (i,binMassParent,             iPrimary)  &
                           &   *treeCurrent%volumeWeight                                                        &
                           &  )**2
                      !$omp end critical(conditionalMassFunctionAccumulate)
                   end if
                end if
             end do
          end do
       end do
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine conditionalMFOperate

  subroutine conditionalMFFinalize(self)
    !% Outputs conditional mass function.
    use ISO_Varying_String
    use IO_HDF5
    use Galacticus_HDF5
    use Numerical_Constants_Astronomical
    use Memory_Management
    implicit none
    class           (mergerTreeOperatorConditionalMF), intent(inout)                     :: self
    type            (hdf5Object                     )                                    :: conditionalMassFunctionGroup , massDataset
    double precision                                 , allocatable  , dimension(:,:,:  ) :: conditionalMassFunction      , conditionalMassFunctionError
    double precision                                 , allocatable  , dimension(:,:,:,:) :: primaryProgenitorMassFunction, primaryProgenitorMassFunctionError
    double precision                                 , allocatable  , dimension(:,:,:,:) :: formationRateFunction        , formationRateFunctionError
    integer                                                                              :: i                            , j                                 , &
         &                                                                                  iPrimary

    ! Normalize the conditional mass functions.
    self%normalization=self%normalization/self%massRatioLogarithmicBinWidthInverse/log(10.0d0)
    do i=1,self%timeCount
       do j=1,self%parentMassCount
          if (self%normalization(i,j) > 0.0d0) then
             self%conditionalMassFunction     (i,j,:  )=     self%conditionalMassFunction     (i,j,:  ) /self%normalization(i,j)
             self%conditionalMassFunctionError(i,j,:  )=sqrt(self%conditionalMassFunctionError(i,j,:  ))/self%normalization(i,j)
             self%formationRateFunction       (i,j,:,:)=     self%formationRateFunction       (i,j,:,:) /self%normalization(i,j)
             self%formationRateFunctionError  (i,j,:,:)=sqrt(self%formationRateFunctionError  (i,j,:,:))/self%normalization(i,j)
             do iPrimary=1,self%primaryProgenitorDepth
                self%primaryProgenitorMassFunction     (i,j,:,iPrimary)=     self%primaryProgenitorMassFunction     (i,j,:,iPrimary) /self%normalization(i,j)
                self%primaryProgenitorMassFunctionError(i,j,:,iPrimary)=sqrt(self%primaryProgenitorMassFunctionError(i,j,:,iPrimary))/self%normalization(i,j)
             end do
          end if
       end do
    end do
    ! Output the data.
    !$omp critical(HDF5_Access)
    ! Check if our output group already exists.
    if (galacticusOutputFile%hasGroup('char(self%outputGroupName)')) then
       ! Our group does exist. Read existing mass functions, add them to our own, then write back to file.
       conditionalMassFunctionGroup=galacticusOutputFile%openGroup(char(self%outputGroupName),'Conditional mass functions of merger trees.')
       call Alloc_Array(conditionalMassFunction           ,shape(self%conditionalMassFunction           ))
       call Alloc_Array(conditionalMassFunctionError      ,shape(self%conditionalMassFunctionError      ))
       call Alloc_Array(primaryProgenitorMassFunction     ,shape(self%primaryProgenitorMassFunction     ))
       call Alloc_Array(primaryProgenitorMassFunctionError,shape(self%primaryProgenitorMassFunctionError))
       call Alloc_Array(formationRateFunction             ,shape(self%formationRateFunction             ))
       call Alloc_Array(formationRateFunctionError        ,shape(self%formationRateFunctionError        ))
       call conditionalMassFunctionGroup%readDataset('conditionalMassFunction'           ,conditionalMassFunction           )
       call conditionalMassFunctionGroup%readDataset('conditionalMassFunctionError'      ,conditionalMassFunctionError      )
       call conditionalMassFunctionGroup%readDataset('primaryProgenitorMassFunction'     ,primaryProgenitorMassFunction     )
       call conditionalMassFunctionGroup%readDataset('primaryProgenitorMassFunctionError',primaryProgenitorMassFunctionError)
       call conditionalMassFunctionGroup%readDataset('formationRateFunction'             ,formationRateFunction             )
       call conditionalMassFunctionGroup%readDataset('formationRateFunctionError'        ,formationRateFunctionError        )
       call weightedAverage(self%conditionalMassFunction      ,conditionalMassFunction      ,self%conditionalMassFunctionError      ,conditionalMassFunctionError      ,self%conditionalMassFunction      ,conditionalMassFunction      )
       call weightedAverage(self%primaryProgenitorMassFunction,primaryProgenitorMassFunction,self%primaryProgenitorMassFunctionError,primaryProgenitorMassFunctionError,self%primaryProgenitorMassFunction,primaryProgenitorMassFunction)
       call weightedAverage(self%formationRateFunction        ,formationRateFunction        ,self%formationRateFunctionError        ,formationRateFunctionError        ,self%formationRateFunction        ,formationRateFunction        )
       call Dealloc_Array(conditionalMassFunction           )
       call Dealloc_Array(conditionalMassFunctionError      )
       call Dealloc_Array(primaryProgenitorMassFunction     )
       call Dealloc_Array(primaryProgenitorMassFunctionError)
       call Dealloc_Array(formationRateFunction             )
       call Dealloc_Array(formationRateFunctionError        )
    else
       ! Our group does not already exist. Simply write the data.
       conditionalMassFunctionGroup=galacticusOutputFile%openGroup(char(self%outputGroupName),'Conditional mass functions of merger trees.')
       call conditionalMassFunctionGroup%writeDataset  (self%massParents                       ,"massParent"                        ,"Mass of parent node [Msolar]"              ,datasetReturned=massDataset)
       call massDataset                 %writeAttribute(massSolar                              ,"unitsInSI"                                                                                                  )
       call massDataset                 %close         (                                                                                                                                                     )
       call conditionalMassFunctionGroup%writeDataset  (self%massRatios                        ,"massRatio"                         ,"Mass of ratio node [Msolar]"               ,datasetReturned=massDataset)
       call massDataset                 %writeAttribute(massSolar                              ,"unitsInSI"                                                                                                  )
       call massDataset                 %close         (                                                                                                                                                     )
       call conditionalMassFunctionGroup%writeDataset  (self%parentRedshifts                   ,"redshiftParent"                    ,"Redshift of parent node []"                                            )
       call conditionalMassFunctionGroup%writeDataset  (self%progenitorRedshifts               ,"redshiftProgenitor"                ,"Redshift of progenitor node []"                                        )
    end if
    call conditionalMassFunctionGroup%writeDataset  (self%conditionalMassFunction           ,"conditionalMassFunction"           ,"Conditional mass functions []"                                         )
    call conditionalMassFunctionGroup%writeDataset  (self%conditionalMassFunctionError      ,"conditionalMassFunctionError"      ,"Conditional mass function errors []"                                   )
    call conditionalMassFunctionGroup%writeDataset  (self%primaryProgenitorMassFunction     ,"primaryProgenitorMassFunction"     ,"Primary progenitor mass functions []"                                  )
    call conditionalMassFunctionGroup%writeDataset  (self%primaryProgenitorMassFunctionError,"primaryProgenitorMassFunctionError","Primary progenitor mass function errors []"                            )
    call conditionalMassFunctionGroup%writeDataset  (self%formationRateFunction             ,"formationRateFunction"             ,"Formation rate functions []"                                           )
    call conditionalMassFunctionGroup%writeDataset  (self%formationRateFunctionError        ,"formationRateFunctionError"        ,"Formation rate function errors []"                                     )
    call conditionalMassFunctionGroup%close         (                                                                                                                                                     )    
    !$omp end critical(HDF5_Access)
    return

  contains

    elemental subroutine weightedAverage(x,y,xError,yError,z,zError)
      !% Computed a minimum-variance weighted average of two values, {\normalfont \ttfamily x} and {\normalfont \ttfamily y},
      !% given errors on those quantities.
      implicit none
      double precision, intent(in   ) :: x     , y     , &
           &                             xError, yError
      double precision, intent(  out) :: z     , zError

      ! Compute the weighted mean.      
      z      =+    (                 &
           &        +    x/xError**2 &
           &        +    y/yError**2 &
           &       )                 &
           &  /    (                 &
           &        +1.0d0/xError**2 &
           &        +1.0d0/yError**2 &
           &       )
      ! Compute the variance on the weighted mean.
      zError =+1.0d0                 &
           &  /sqrt(                 &
           &        +1.0d0/xError**2 &
           &        +1.0d0/yError**2 &
           &       )
      return
    end subroutine weightedAverage
      
  end subroutine conditionalMFFinalize
