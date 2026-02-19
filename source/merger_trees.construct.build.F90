!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
!!    Andrew Benson <abenson@carnegiescience.edu>
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

  !!{
  Implements a merger tree constructor class which builds merger trees after drawing masses at random from a mass distribution.
  !!}

  use :: Cosmology_Functions      , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters     , only : cosmologyParametersClass
  use :: Halo_Mass_Functions      , only : haloMassFunctionClass
  use :: Merger_Trees_Build_Masses, only : mergerTreeBuildMassesClass
  use :: Merger_Trees_Builders    , only : mergerTreeBuilderClass
  use :: Merger_Tree_Seeds        , only : mergerTreeSeedsClass
  use :: Numerical_Random_Numbers , only : randomNumberGeneratorClass
  use :: Output_Times             , only : outputTimesClass

  !![
  <mergerTreeConstructor name="mergerTreeConstructorBuild">
   <description>
    A merger tree constructor class which builds merger trees. This class first creates a distribution of tree root halo masses
    and then builds a merger tree from each root halo.
   </description>
  </mergerTreeConstructor>
  !!]
  type, extends(mergerTreeConstructorClass) :: mergerTreeConstructorBuild
     !!{
     A class implementing merger tree construction by building trees.
     !!}
     private
     class           (cosmologyParametersClass  ), pointer                   :: cosmologyParameters_   => null()
     class           (cosmologyFunctionsClass   ), pointer                   :: cosmologyFunctions_    => null()
     class           (mergerTreeBuildMassesClass), pointer                   :: mergerTreeBuildMasses_ => null()
     class           (mergerTreeBuilderClass    ), pointer                   :: mergerTreeBuilder_     => null()
     class           (mergerTreeSeedsClass      ), pointer                   :: mergerTreeSeeds_       => null()
     class           (haloMassFunctionClass     ), pointer                   :: haloMassFunction_      => null()
     class           (outputTimesClass          ), pointer                   :: outputTimes_           => null()
     class           (randomNumberGeneratorClass), pointer                   :: randomNumberGenerator_ => null()
     ! Variables giving the mass range and sampling frequency for mass function sampling.
     double precision                                                        :: timeBase                        , timeSnapTolerance, &
          &                                                                     redshiftBase
     integer                                                                 :: treeBeginAt
     ! Direction in which to process trees.
     logical                                                                 :: processDescending
     ! Array of halo masses to use.
     integer         (c_size_t                  )                            :: treeCount                       , treeNumberOffset
     double precision                            , allocatable, dimension(:) :: treeMass                        , treeWeight       , &
          &                                                                     treeMassMinimum                 , treeMassMaximum
     integer         (c_size_t                  ), allocatable, dimension(:) :: treeMassCount                   , rankMass
     logical                                                                 :: computeTreeWeights
   contains
     !![
     <methods>
       <method description="Construct the set of tree masses to be built." method="constructMasses" />
     </methods>
     !!]
     final     ::                    buildDestructor
     procedure :: construct       => buildConstruct
     procedure :: constructMasses => buildConstructMasses
  end type mergerTreeConstructorBuild

  interface mergerTreeConstructorBuild
     !!{
     Constructors for the \refClass{mergerTreeConstructorBuild} merger tree constructor class.
     !!}
     module procedure buildConstructorParameters
     module procedure buildConstructorInternal
  end interface mergerTreeConstructorBuild

contains

  function buildConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeConstructorBuild} merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter         , inputParameters
    implicit none
    type            (mergerTreeConstructorBuild)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (cosmologyParametersClass  ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass   ), pointer       :: cosmologyFunctions_
    class           (mergerTreeBuilderClass    ), pointer       :: mergerTreeBuilder_
    class           (mergerTreeSeedsClass      ), pointer       :: mergerTreeSeeds_
    class           (haloMassFunctionClass     ), pointer       :: haloMassFunction_
    class           (mergerTreeBuildMassesClass), pointer       :: mergerTreeBuildMasses_
    class           (randomNumberGeneratorClass), pointer       :: randomNumberGenerator_
    class           (outputTimesClass          ), pointer       :: outputTimes_
    double precision                                            :: redshiftBase          , timeSnapTolerance
    integer                                                     :: treeBeginAt
    logical                                                     :: processDescending

    !![
    <inputParameter>
      <name>redshiftBase</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The redshift at which to plant the base node when building merger trees.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>timeSnapTolerance</name>
      <defaultValue>1.0d-6</defaultValue>
      <description>The fractional tolerance within which the tree base time will be snapped to a nearby output time.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>treeBeginAt</name>
      <defaultValue>0</defaultValue>
      <description>The index (in order of increasing base halo mass) of the tree at which to begin when building merger trees. A value of ``0'' means to begin with tree number 1 (if processing trees in ascending order), or equal to the number of trees (otherwise).</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>processDescending</name>
      <defaultValue>.true.</defaultValue>
      <description>If true, causes merger trees to be processed in order of decreasing mass.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="mergerTreeBuilder"     name="mergerTreeBuilder_"     source="parameters"/>
    <objectBuilder class="mergerTreeSeeds"       name="mergerTreeSeeds_"       source="parameters"/>
    <objectBuilder class="haloMassFunction"      name="haloMassFunction_"      source="parameters"/>
    <objectBuilder class="mergerTreeBuildMasses" name="mergerTreeBuildMasses_" source="parameters"/>
    <objectBuilder class="outputTimes"           name="outputTimes_"           source="parameters"/>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=mergerTreeConstructorBuild(                                                                                                        &
         &                          cosmologyFunctions_    %cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftBase     )), &
         &                                                                                                             timeSnapTolerance  , &
         &                                                                                                             treeBeginAt        , &
         &                                                                                                             processDescending  , &
         &                          cosmologyParameters_                                                                                  , &
         &                          cosmologyFunctions_                                                                                   , &
         &                          mergerTreeBuildMasses_                                                                                , &
         &                          mergerTreeBuilder_                                                                                    , &
         &                          mergerTreeSeeds_                                                                                      , &
         &                          haloMassFunction_                                                                                     , &
         &                          outputTimes_                                                                                          , &
         &                          randomNumberGenerator_                                                                                  &
         &                         )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="mergerTreeBuilder_"    />
    <objectDestructor name="mergerTreeSeeds_"      />
    <objectDestructor name="haloMassFunction_"     />
    <objectDestructor name="mergerTreeBuildMasses_"/>
    <objectDestructor name="outputTimes_"          />
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function buildConstructorParameters

  function buildConstructorInternal(timeBase,timeSnapTolerance,treeBeginAt,processDescending,cosmologyParameters_,cosmologyFunctions_,mergerTreeBuildMasses_,mergerTreeBuilder_,mergerTreeSeeds_,haloMassFunction_,outputTimes_,randomNumberGenerator_) result(self)
    !!{
    Initializes the merger tree building module.
    !!}
    use, intrinsic :: ISO_C_Binding       , only : c_size_t
    use            :: Numerical_Comparison, only : Values_Agree
    implicit none
    type            (mergerTreeConstructorBuild)                        :: self
    double precision                            , intent(in   )         :: timeBase            , timeSnapTolerance
    integer                                     , intent(in   )         :: treeBeginAt
    logical                                     , intent(in   )         :: processDescending
    class           (cosmologyParametersClass  ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass   ), intent(in   ), target :: cosmologyFunctions_
    class           (mergerTreeBuilderClass    ), intent(in   ), target :: mergerTreeBuilder_
    class           (mergerTreeSeedsClass      ), intent(in   ), target :: mergerTreeSeeds_
    class           (haloMassFunctionClass     ), intent(in   ), target :: haloMassFunction_
    class           (mergerTreeBuildMassesClass), intent(in   ), target :: mergerTreeBuildMasses_
    class           (outputTimesClass          ), intent(in   ), target :: outputTimes_
    class           (randomNumberGeneratorClass), intent(in   ), target :: randomNumberGenerator_
    integer         (c_size_t                  )                        :: i
    !![
    <constructorAssign variables="timeBase, timeSnapTolerance, treeBeginAt, processDescending, *cosmologyParameters_, *cosmologyFunctions_, *mergerTreeBuildMasses_, *mergerTreeBuilder_, *mergerTreeSeeds_, *haloMassFunction_, *outputTimes_, *randomNumberGenerator_"/>
    !!]

    ! Store base redshift.
    self%redshiftBase=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeBase))
    ! Set offset for tree numbers.
    if (self%treeBeginAt == 0) then
       self%treeNumberOffset=0_c_size_t
    else
       self%treeNumberOffset=self%treeBeginAt-1_c_size_t
    end if
    ! Snap tree base time if necessary.
    if (self%timeSnapTolerance > 0.0d0) then
       do i=1_c_size_t,self%outputTimes_%count()
          if (Values_Agree(self%timeBase,self%outputTimes_%time(i),relTol=self%timeSnapTolerance)) self%timeBase=self%outputTimes_%time(i)
       end do
    end if
    return
  end function buildConstructorInternal

  subroutine buildDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeConstructorBuild} merger tree constructor class.
    !!}
    implicit none
    type(mergerTreeConstructorBuild), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%mergerTreeBuilder_"    />
    <objectDestructor name="self%mergerTreeSeeds_"      />
    <objectDestructor name="self%haloMassFunction_"     />
    <objectDestructor name="self%mergerTreeBuildMasses_"/>
    <objectDestructor name="self%outputTimes_"          />
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine buildDestructor

  function buildConstruct(self,treeNumber,finished) result(tree)
    !!{
    Build a merger tree.
    !!}
    use :: Functions_Global       , only : State_Retrieve_       , State_Store_
    use :: Galacticus_Nodes       , only : mergerTree            , nodeComponentBasic, treeNode
    use :: Kind_Numbers           , only : kind_int8
    use :: Merger_Tree_State_Store, only : treeStateStoreSequence
    use :: String_Handling        , only : operator(//)
    implicit none
    type   (mergerTree                ), pointer       :: tree
    class  (mergerTreeConstructorBuild), intent(inout) :: self
    integer(c_size_t                  ), intent(in   ) :: treeNumber
    logical                            , intent(  out) :: finished
    class  (nodeComponentBasic        ), pointer       :: basicBase
    integer(kind_int8                 ), parameter     :: indexNodeBase=1
    integer(kind_int8                 )                :: treeIndex
    type   (varying_string            )                :: message

    ! Ensure masses are constructed.
    call self%constructMasses()
    ! Prepare to store/restore internal state.
    treeStateStoreSequence=-1_c_size_t
    ! Retrieve stored internal state if possible.
    call State_Retrieve_()
    if (treeStateStoreSequence > 0_c_size_t) then
       if (self%processDescending) then
          self%treeNumberOffset=self%treeCount-treeStateStoreSequence
       else
          self%treeNumberOffset=              +treeStateStoreSequence-1_c_size_t
       end if
    end if
    ! Determine the index of the tree to process.
    if (self%processDescending) then
       ! Processing trees in descending order, so begin from the final index and work back.
       treeIndex=self%treeCount+1-(treeNumber+self%treeNumberOffset)
    else
       ! Processing trees in ascending order, to just use treeNumber as the index of the tree to process.
       treeIndex=                +(treeNumber+self%treeNumberOffset)
    end if
    if     (                             &
         &   treeIndex >  0_kind_int8    &
         &  .and.                        &
         &   treeIndex <= self%treeCount &
         & ) then
       ! Allocate the tree.
       allocate(tree)
       ! Give the tree an index.
       tree%index=treeIndex
       ! Restart the random number sequence.
       allocate(tree%randomNumberGenerator_,mold=self%randomNumberGenerator_)
       !$omp critical(mergerTreeConstructBuildDeepCopyReset)
       !![
       <deepCopyReset variables="self%randomNumberGenerator_"/>
       <deepCopy source="self%randomNumberGenerator_" destination="tree%randomNumberGenerator_"/>
       <deepCopyFinalize variables="tree%randomNumberGenerator_"/>
       !!]
       !$omp end critical(mergerTreeConstructBuildDeepCopyReset)
       call self                 %randomSequenceNonDeterministicWarn(tree)
       call self%mergerTreeSeeds_%set                               (tree)
       ! Store the internal state.
       if (treeStateStoreSequence == -1_c_size_t) treeStateStoreSequence=treeNumber
       message=var_str('Storing state for tree #')//treeNumber
       call State_Store_(message)
       ! Initialize.
       tree%event            => null()
       tree%initializedUntil =  0.0d0
       call tree%properties%initialize()
       ! Create the base node.
       tree%firstTree => tree
       tree%nodeBase  => treeNode(indexNodeBase,tree)
       ! Get the basic component of the base node.
       basicBase     => tree%nodeBase%basic(autoCreate=.true.)
       ! Assign a mass to it.
       call basicBase%massSet(self%treeMass(self%rankMass(treeIndex)))
       ! Assign a time.
       call basicBase%timeSet(self%timeBase                          )
       ! Assign a weight to the tree, computing it if necessary.
       if (self%computeTreeWeights) then
          ! The weight is computed by finding the total mass in halos per unit volume within the mass range represented by this
          ! tree, and dividing that by the tree mass. This ensures that we account for all mass that should be in halos. (The
          ! alternative, just computing the number density of halos within the mass range by integrating the mass function, gives
          ! a slightly different answer as the contribution from each mass in the range is then not mass-weighted.) We finally
          ! divide through by the number of trees with this mass.
          self%treeWeight(self%rankMass(treeIndex))=+self%haloMassFunction_   %massFraction   (                                                &
               &                                                                               self%timeBase                                 , &
               &                                                                               self%treeMassMinimum(self%rankMass(treeIndex)), &
               &                                                                               self%treeMassMaximum(self%rankMass(treeIndex)), &
               &                                                                               tree%nodeBase                                   &
               &                                                                              )                                                &
               &                                    *self%cosmologyParameters_%densityCritical(                                                &
               &                                                                              )                                                &
               &                                    *self%cosmologyParameters_%OmegaMatter    (                                                &
               &                                                                              )                                                &
               &                                    /                                          self%treeMass       (self%rankMass(treeIndex))  &
               &                                    /dble                                     (                                                &
               &                                                                               self%treeMassCount  (              treeIndex)   &
               &                                                                              )
       end if
       tree%volumeWeight=self%treeWeight(self%rankMass(treeIndex))
       ! Build the tree.
       call self%mergerTreeBuilder_%build(tree)
    else
       nullify(tree)
    end if
    finished=.not.associated(tree)
    return
  end function buildConstruct

  subroutine buildConstructMasses(self)
    !!{
    Construct the set of tree masses to be built.
    !!}
    use :: Error  , only : Error_Report
    use :: Sorting, only : sortIndex
    implicit none
    class  (mergerTreeConstructorBuild), intent(inout) :: self
    integer(c_size_t                  )                :: iTreeFirst, iTreeLast

    if (.not.allocated(self%treeMass)) then
       ! Generate set of merger tree masses.
       call self%mergerTreeBuildMasses_%construct(self%timeBase,self%treeMass,self%treeMassMinimum,self%treeMassMaximum,self%treeWeight)
       self%treeCount=size(self%treeMass,kind=c_size_t)
       ! Sort halos by mass.
       allocate(self%rankMass(self%treeCount))
       self%rankMass=sortIndex(self%treeMass)
       ! Compute the weight (number of trees per unit volume) for each tree if weights were not supplied.
       self%computeTreeWeights=.not.allocated(self%treeWeight)
       if (.not.self%computeTreeWeights) then
          if     (                                 &
               &   allocated(self%treeMassMinimum) &
               &  .or.                             &
               &   allocated(self%treeMassMaximum) &
               & ) call Error_Report('mass interval should not be set if tree weights are provided'//{introspection:location})
       else
          allocate(self%treeWeight   (self%treeCount))
          allocate(self%treeMassCount(self%treeCount))
          iTreeFirst=0
          do while (iTreeFirst < self%treeCount)
             iTreeFirst=iTreeFirst+1
             ! Find the last tree with the same mass.
             iTreeLast=iTreeFirst
             if (iTreeLast < self%treeCount) then
                do while (self%treeMass(self%rankMass(iTreeLast+1)) == self%treeMass(self%rankMass(iTreeFirst)))
                   iTreeLast=iTreeLast+1
                   if (iTreeLast == self%treeCount) exit
                end do
             end if
             ! Store the number of trees at this mass.
             self%treeMassCount(iTreeFirst:iTreeLast)=iTreeLast-iTreeFirst+1
             ! Update to the last tree processed.
             iTreeFirst=iTreeLast
          end do
       end if
    end if
    return
  end subroutine buildConstructMasses

