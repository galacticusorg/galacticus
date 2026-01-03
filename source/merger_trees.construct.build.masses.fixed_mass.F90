!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  Implementation of a merger tree masses class which uses a fixed mass for trees.
  !!}
  use :: Cosmology_Parameters     , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales  , only : darkMatterHaloScaleClass
  use :: Nodes_Operators          , only : nodeOperatorClass
  use :: Numerical_Random_Numbers , only : randomNumberGeneratorClass

  !![
  <mergerTreeBuildMasses name="mergerTreeBuildMassesFixedMass">
   <description>A merger tree masses class which uses a fixed mass for trees.</description>
  </mergerTreeBuildMasses>
  !!]
  type, extends(mergerTreeBuildMassesClass) :: mergerTreeBuildMassesFixedMass
     !!{
     Implementation of a merger tree masses class which samples masses from a distribution.
     !!}
     private
     class           (cosmologyParametersClass  ), pointer                   :: cosmologyParameters_   => null()
     class           (darkMatterHaloScaleClass  ), pointer                   :: darkMatterHaloScale_   => null()
     class           (nodeOperatorClass         ), pointer                   :: nodeOperator_          => null()
     class           (randomNumberGeneratorClass), pointer                   :: randomNumberGenerator_ => null()
     double precision                            , allocatable, dimension(:) :: massTree                        , radiusTree
     integer                                     , allocatable, dimension(:) :: treeCount
     double precision                                                        :: massIntervalFractional
   contains
     final     ::              fixedMassDestructor
     procedure :: construct => fixedMassConstruct
  end type mergerTreeBuildMassesFixedMass

  interface mergerTreeBuildMassesFixedMass
     module procedure fixedMassConstructorParameters
     module procedure fixedMassConstructorInternal
  end interface mergerTreeBuildMassesFixedMass

contains

  function fixedMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeBuildMassesFixedMass} merger tree masses class which takes a parameter set as
    input.
    !!}
    use :: Error            , only : Error_Report
    use :: Input_Parameters , only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeBuildMassesFixedMass)                            :: self
    type            (inputParameters               ), intent(inout)             :: parameters
    double precision                                , allocatable, dimension(:) :: massTree              , radiusTree
    integer                                         , allocatable, dimension(:) :: treeCount
    class           (cosmologyParametersClass      ), pointer                   :: cosmologyParameters_
    class           (darkMatterHaloScaleClass      ), pointer                   :: darkMatterHaloScale_
    class           (nodeOperatorClass             ), pointer                   :: nodeOperator_
    class           (randomNumberGeneratorClass    ), pointer                   :: randomNumberGenerator_
    double precision                                                            :: massIntervalFractional
    integer                                                                     :: fixedHalosCount

    fixedHalosCount=-1
    if (parameters%isPresent('massTree' )) then
       if     (                                                       &
            &   fixedHalosCount /= -1                                 &
            &  .and.                                                  &
            &   fixedHalosCount /= parameters%count('massTree'  )     &
            & )                                                       &
            & call Error_Report(                                      &
            &                   'parameter cardinality mismatch'//    &
            &                   {introspection:location}              &
            &                  )
       if (fixedHalosCount == -1) fixedHalosCount=parameters%count('massTree' )
    end if
    if (parameters%isPresent('treeCount')) then
       if     (                                                       &
            &   fixedHalosCount /= -1                                 &
            &  .and.                                                  &
            &   fixedHalosCount /= parameters%count('treeCount' )     &
            & )                                                       &
            & call Error_Report(                                      &
            &                   'parameter cardinality mismatch'//    &
            &                   {introspection:location}              &
            &                  )
       if (fixedHalosCount == -1) fixedHalosCount=parameters%count('treeCount')
    end if
    if (parameters%isPresent('radiusTree')) then
       if     (                                                       &
            &   fixedHalosCount /= -1                                 &
            &  .and.                                                  &
            &   fixedHalosCount /= parameters%count('radiusTree')     &
            & )                                                       &
            & call Error_Report(                                      &
            &                   'parameter cardinality mismatch'//    &
            &                   {introspection:location}              &
            &                  )
       if (fixedHalosCount == -1) fixedHalosCount=parameters%count('radiusTree')
    end if
    if (fixedHalosCount == -1) fixedHalosCount=1
    allocate(massTree  (fixedHalosCount))
    allocate(treeCount (fixedHalosCount))
    allocate(radiusTree(fixedHalosCount))
    !![
    <inputParameter>
      <name>massTree</name>
      <defaultValue>spread(1.0d12,1,fixedHalosCount)</defaultValue>
      <description>Specifies the masses of halos to use when building halos.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>treeCount</name>
      <defaultValue>spread(1,1,fixedHalosCount)</defaultValue>
      <description>Specifies the number of halos to use when building halos.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusTree</name>
      <defaultValue>spread(-1.0d0,1,fixedHalosCount)</defaultValue>
      <description>Specifies the radii within which halo masses are specified when building halos.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massIntervalFractional</name>
      <defaultValue>0.1d0</defaultValue>
      <description>The fractional mass interval occupied by the trees. Where the intervals of trees of different mass would overlap this interval will be truncated.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"   name="darkMatterHaloScale_"   source="parameters"/>
    <objectBuilder class="nodeOperator"          name="nodeOperator_"          source="parameters"/>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=mergerTreeBuildMassesFixedMass(massTree,radiusTree,treeCount,massIntervalFractional,cosmologyParameters_,darkMatterHaloScale_,nodeOperator_,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="darkMatterHaloScale_"  />
    <objectDestructor name="nodeOperator_"         />
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function fixedMassConstructorParameters

  function fixedMassConstructorInternal(massTree,radiusTree,treeCount,massIntervalFractional,cosmologyParameters_,darkMatterHaloScale_,nodeOperator_,randomNumberGenerator_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeBuildMassesFixedMass} merger tree masses class.
    !!}
    implicit none
    type            (mergerTreeBuildMassesFixedMass)                              :: self
    double precision                                , intent(in   ), dimension(:) :: massTree              , radiusTree
    integer                                         , intent(in   ), dimension(:) :: treeCount
    double precision                                , intent(in   )               :: massIntervalFractional
    class           (cosmologyParametersClass      ), intent(in   ), target       :: cosmologyParameters_
    class           (darkMatterHaloScaleClass      ), intent(in   ), target       :: darkMatterHaloScale_
    class           (nodeOperatorClass             ), intent(in   ), target       :: nodeOperator_
    class           (randomNumberGeneratorClass    ), intent(in   ), target       :: randomNumberGenerator_
    !![
    <constructorAssign variables="massTree, radiusTree, treeCount, massIntervalFractional, *cosmologyParameters_, *darkMatterHaloScale_, *nodeOperator_, *randomNumberGenerator_"/>
    !!]

    return
  end function fixedMassConstructorInternal

  subroutine fixedMassDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeBuildMassesFixedMass} merger tree masses class.
    !!}
    implicit none
    type(mergerTreeBuildMassesFixedMass), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%darkMatterHaloScale_"  />
    <objectDestructor name="self%nodeOperator_"         />
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine fixedMassDestructor

  subroutine fixedMassConstruct(self,time,mass,massMinimum,massMaximum,weight)
    !!{
    Construct a set of merger tree masses by sampling from a distribution.
    !!}
    use :: Calculations_Resets, only : Calculations_Reset
    use :: Galacticus_Nodes   , only : nodeComponentBasic           , treeNode  , mergerTree
    use :: Root_Finder        , only : rangeExpandMultiplicative    , rootFinder
    use :: Sorting            , only : sort
    implicit none
    class           (mergerTreeBuildMassesFixedMass), intent(inout)                            :: self
    double precision                                , intent(in   )                            :: time
    double precision                                , intent(  out), allocatable, dimension(:) :: mass       , weight     , &
         &                                                                                        massMinimum, massMaximum
    type            (rootFinder                    )                                           :: finder
    type            (mergerTree                    )                                           :: tree
    type            (treeNode                      ), pointer                                  :: node
    class           (nodeComponentBasic            ), pointer                                  :: basic
    integer                                                                                    :: indexStart , i
    !$GLC attributes unused :: weight

    finder=rootFinder(                                               &
         &            rootFunction       =massEnclosed             , &
         &            toleranceAbsolute  =1.0d-6                   , &
         &            toleranceRelative  =1.0d-6                   , &
         &            rangeExpandUpward  =2.0d+0                   , &
         &            rangeExpandDownward=0.5d+0                   , &
         &            rangeExpandType    =rangeExpandMultiplicative  &
         &           )
    ! Restart the random number sequence.
    tree%index=1
    allocate(tree%randomNumberGenerator_,mold=self%randomNumberGenerator_)
    !$omp critical(mergerTreeConstructBuildMassesFixedDeepCopyReset)
    !![
    <deepCopyReset variables="self%randomNumberGenerator_"/>
    <deepCopy source="self%randomNumberGenerator_" destination="tree%randomNumberGenerator_"/>
    <deepCopyFinalize variables="tree%randomNumberGenerator_"/>
    !!]
    !$omp end critical(mergerTreeConstructBuildMassesFixedDeepCopyReset)
    call tree%randomNumberGenerator_%seedSet(seed=tree%index,offset=.true.)
    node          => treeNode      (hostTree  =tree  )
    tree%nodeBase => node
    basic         => node    %basic(autoCreate=.true.)
    call basic%timeSet            (time)
    call basic%timeLastIsolatedSet(time)
    do i=1,size(self%treeCount)
       if (self%radiusTree(i) > 0.0d0) then
          ! Set the halo mass.
          call Calculations_Reset(node)
          ! Convert masses to virial masses.        
          self%massTree(i)=finder%find(rootGuess=self%massTree(i))
       end if
    end do
    call tree%destroy()
    allocate(mass       (sum(self%treeCount)))
    allocate(massMinimum(sum(self%treeCount)))
    allocate(massMaximum(sum(self%treeCount)))
    indexStart=1
    do i=1,size(self%treeCount)
       mass      (indexStart:indexStart+self%treeCount(i)-1)=+self%massTree (i)
       indexStart                                           =+self%treeCount(i) &
            &                                                +indexStart
    end do
    call sort(mass)
    indexStart=1
    do i=1,size(self%treeCount)
       massMinimum(indexStart:indexStart+self%treeCount(i)-1)=+self%massTree (i)/sqrt(1.0d0+self%massIntervalFractional)
       massMaximum(indexStart:indexStart+self%treeCount(i)-1)=+self%massTree (i)*sqrt(1.0d0+self%massIntervalFractional)
       if (i > 1 .and. massMinimum(indexStart) < massMaximum(indexStart-1)) then
          massMinimum(indexStart                    :indexStart+self%treeCount(i)-1)=sqrt(                    &
               &                                                                          +self%massTree(i-1) &
               &                                                                          *self%massTree(i  ) &
               &                                                                         )
          massMaximum(indexStart-self%treeCount(i-1):indexStart                  -1)=sqrt(                    &
               &                                                                          +self%massTree(i-1) &
               &                                                                          *self%massTree(i  ) &
               &                                                                         )
       end if
       indexStart=+self%treeCount(i)+indexStart
    end do
    return

  contains

    double precision function massEnclosed(massTree)
      !!{
      Root finding function used to set the halo mass given the halo radius.
      !!}
      use :: Galactic_Structure_Options, only : massTypeDark
      use :: Mass_Distributions        , only : massDistributionClass
      implicit none
      double precision                       , intent(in   ) :: massTree
      class           (massDistributionClass), pointer       :: massDistribution_

      call basic              %massSet           (massTree)
      call self %nodeOperator_%nodeTreeInitialize(node    )
      call self %nodeOperator_%nodeInitialize    (node    )
      call Calculations_Reset(node)
      massDistribution_ =>    node                                  %massDistribution    (massType=     massTypeDark   )
      massEnclosed      =  +  massDistribution_                     %massEnclosedBySphere(         self%radiusTree  (i)) &
           &               *  self             %cosmologyParameters_%OmegaMatter()                                       &
           &               /(                                                                                            &
           &                 +self             %cosmologyParameters_%OmegaMatter()                                       &
           &                 -self             %cosmologyParameters_%OmegaBaryon()                                       &
           &                )                                                                                            &
           &               -                                                                       self%massTree    (i)
      !![
      <objectDestructor name="massDistribution_"/>
      !!]
      return
    end function massEnclosed

  end subroutine fixedMassConstruct
