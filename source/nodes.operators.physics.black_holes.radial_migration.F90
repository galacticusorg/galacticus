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
  Implements a node operator class that implements radial migration of non-central black holes.
  !!}

  use :: Black_Hole_Binary_Mergers          , only : blackHoleBinaryMergerClass
  use :: Black_Hole_Binary_Recoil_Velocities, only : blackHoleBinaryRecoilClass
  use :: Black_Hole_Binary_Separations      , only : blackHoleBinarySeparationGrowthRateClass
  use :: Black_Hole_Seeds                   , only : blackHoleSeedsClass
  use :: Dark_Matter_Halo_Scales            , only : darkMatterHaloScaleClass

  !![
  <nodeOperator name="nodeOperatorBlackHolesRadialMigration">
   <description>A node operator class that implements radial migration of non-central black holes.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorBlackHolesRadialMigration
     !!{
     A node operator class that implements radial migration of non-central black holes.
     !!}
     private
     class(darkMatterHaloScaleClass                ), pointer :: darkMatterHaloScale_                 => null()
     class(blackHoleBinaryRecoilClass              ), pointer :: blackHoleBinaryRecoil_               => null()
     class(blackHoleBinaryMergerClass              ), pointer :: blackHoleBinaryMerger_               => null()
     class(blackHoleBinarySeparationGrowthRateClass), pointer :: blackHoleBinarySeparationGrowthRate_ => null()
     class(blackHoleSeedsClass                     ), pointer :: blackHoleSeeds_                      => null()
   contains
     final     ::                                blackHolesRadialMigrationDestructor
     procedure :: differentialEvolutionScales => blackHolesRadialMigrationDifferentialEvolutionScales
     procedure :: differentialEvolution       => blackHolesRadialMigrationDifferentialEvolution
  end type nodeOperatorBlackHolesRadialMigration
  
  interface nodeOperatorBlackHolesRadialMigration
     !!{
     Constructors for the \refClass{nodeOperatorBlackHolesRadialMigration} node operator class.
     !!}
     module procedure blackHolesRadialMigrationConstructorParameters
     module procedure blackHolesRadialMigrationConstructorInternal
  end interface nodeOperatorBlackHolesRadialMigration

  ! Submodule-scope objects used in callback functions.
  class  (nodeOperatorBlackHolesRadialMigration), pointer :: self_
  integer                                                 :: mergingInstance
  !$omp threadprivate(self_,mergingInstance)

contains

  function blackHolesRadialMigrationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorBlackHolesRadialMigration} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorBlackHolesRadialMigration   )                :: self
    type (inputParameters                         ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass                ), pointer       :: darkMatterHaloScale_
    class(blackHoleBinaryRecoilClass              ), pointer       :: blackHoleBinaryRecoil_
    class(blackHoleBinaryMergerClass              ), pointer       :: blackHoleBinaryMerger_
    class(blackHoleBinarySeparationGrowthRateClass), pointer       :: blackHoleBinarySeparationGrowthRate_
    class(blackHoleSeedsClass                     ), pointer       :: blackHoleSeeds_

    !![
    <objectBuilder class="darkMatterHaloScale"                 name="darkMatterHaloScale_"                 source="parameters"/>
    <objectBuilder class="blackHoleBinaryRecoil"               name="blackHoleBinaryRecoil_"               source="parameters"/>
    <objectBuilder class="blackHoleBinaryMerger"               name="blackHoleBinaryMerger_"               source="parameters"/>
    <objectBuilder class="blackHoleBinarySeparationGrowthRate" name="blackHoleBinarySeparationGrowthRate_" source="parameters"/>
    <objectBuilder class="blackHoleSeeds"                      name="blackHoleSeeds_"                      source="parameters"/>
    !!]
    self=nodeOperatorBlackHolesRadialMigration(darkMatterHaloScale_,blackHoleBinaryRecoil_,blackHoleBinaryMerger_,blackHoleBinarySeparationGrowthRate_,blackHoleSeeds_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"                />
    <objectDestructor name="blackHoleBinaryRecoil_"              />
    <objectDestructor name="blackHoleBinaryMerger_"              />
    <objectDestructor name="blackHoleBinarySeparationGrowthRate_"/>
    <objectDestructor name="blackHoleSeeds_"                     />
    !!]
    return
  end function blackHolesRadialMigrationConstructorParameters

  function blackHolesRadialMigrationConstructorInternal(darkMatterHaloScale_,blackHoleBinaryRecoil_,blackHoleBinaryMerger_,blackHoleBinarySeparationGrowthRate_,blackHoleSeeds_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorBlackHolesRadialMigration} node operator class.
    !!}
    implicit none
    type (nodeOperatorBlackHolesRadialMigration   )                        :: self
    class(darkMatterHaloScaleClass                ), intent(in   ), target :: darkMatterHaloScale_
    class(blackHoleBinaryRecoilClass              ), intent(in   ), target :: blackHoleBinaryRecoil_
    class(blackHoleBinaryMergerClass              ), intent(in   ), target :: blackHoleBinaryMerger_
    class(blackHoleBinarySeparationGrowthRateClass), intent(in   ), target :: blackHoleBinarySeparationGrowthRate_
    class(blackHoleSeedsClass                     ), intent(in   ), target :: blackHoleSeeds_
    !![
    <constructorAssign variables="*darkMatterHaloScale_, *blackHoleBinaryRecoil_, *blackHoleBinaryMerger_, *blackHoleBinarySeparationGrowthRate_, *blackHoleSeeds_"/>
    !!]
    
    return
  end function blackHolesRadialMigrationConstructorInternal

  subroutine blackHolesRadialMigrationDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorBlackHolesRadialMigration} node operator class.
    !!}
    implicit none
    type(nodeOperatorBlackHolesRadialMigration), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"                />
    <objectDestructor name="self%blackHoleBinaryRecoil_"              />
    <objectDestructor name="self%blackHoleBinaryMerger_"              />
    <objectDestructor name="self%blackHoleBinarySeparationGrowthRate_"/>
    <objectDestructor name="self%blackHoleSeeds_"                     />
    !!]
    return
  end subroutine blackHolesRadialMigrationDestructor

  subroutine blackHolesRadialMigrationDifferentialEvolutionScales(self,node)
    !!{
    Set absolute ODE solver scale for black hole radial migration.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, nodeComponentSpheroid
    implicit none
    class           (nodeOperatorBlackHolesRadialMigration), intent(inout) :: self
    type            (treeNode                             ), intent(inout) :: node
    double precision                                       , parameter     :: scaleSizeRelative=1.0d-4
    double precision                                       , parameter     :: scaleSizeAbsolute=1.0d-6
    class           (nodeComponentSpheroid                ), pointer       :: spheroid
    class           (nodeComponentBlackHole               ), pointer       :: blackHole
    integer                                                                :: instance

    ! Determine if at least one black hole exists.
    if (node%blackHoleCount() <= 0) return
    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Loop over instances.
    do instance=1,node%blackHoleCount()
       ! Get the black hole.
       blackHole => node%blackHole(instance=instance)
       ! Set scale for radius.
       call blackHole%radialPositionScale(                                                      &
            &                             maxval(                                               &
            &                                    [                                              &
            &                                     scaleSizeAbsolute,                            &
            &                                     scaleSizeRelative*spheroid %halfMassRadius(), &
            &                                                       blackHole%radialPosition()  &
            &                                    ]                                              &
            &                                   )                                               &
            &                            )
    end do
    return
  end subroutine blackHolesRadialMigrationDifferentialEvolutionScales
  
  subroutine blackHolesRadialMigrationDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Account for radial migration of black holes.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBlackHole        , propertyInactive
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (nodeOperatorBlackHolesRadialMigration), intent(inout), target  :: self
    type            (treeNode                             ), intent(inout), target  :: node
    logical                                                , intent(inout)          :: interrupt
    procedure       (interruptTask                        ), intent(inout), pointer :: functionInterrupt
    integer                                                , intent(in   )          :: propertyType
    class           (nodeComponentBlackHole               )               , pointer :: blackHole         , blackHoleCentral
    integer                                                                         :: iInstance         , instanceCount
    double precision                                                                :: radiusHardBinary  , radialMigrationRate
    
    if (propertyInactive(propertyType)) return
    ! Get a count of the number of black holes associated with this node.
    instanceCount=node%blackHoleCount()
    ! Get the central black hole.
    blackHoleCentral => node%blackHole(instance=1)
    ! Do radial migration for non-central black holes.
    do iInstance=2,instanceCount
       ! Get the black hole.
       blackHole => node%blackHole(instance=iInstance)
       ! Compute the hard binary radius.
       radiusHardBinary= (                                                     &
            &              gravitationalConstant_internal                      &
            &             *(                                                   &
            &               +blackHoleCentral%mass()                           &
            &               +blackHole       %mass()                           &
            &              )                                                   &
            &            )                                                     &
            &           /(                                                     &
            &              4.0d0                                               &
            &             *  self%darkMatterHaloScale_%velocityVirial(node)**2 &
            &            )
       ! Places a new black hole in the center of the galaxy in case there is no central one.
       if     (                                                       &
            &   blackHoleCentral%mass          () == 0.0d0            &
            &  .and.                                                  &
            &   blackHole       %radialPosition() <= radiusHardBinary &
            & ) then
          self_             => self
          mergingInstance   =  iInstance
          interrupt         =  .true.
          functionInterrupt => mergeBlackHoles
          return
       end if
       ! Check for a black hole that is about to merge.
       if (blackHole%radialPosition() <= 0.0d0) then
          ! Record which instance is merging, then trigger an interrupt.
          self_             => self
          mergingInstance   =  iInstance
          interrupt         =  .true.
          functionInterrupt => mergeBlackHoles
          return
       end if
       ! Set the rate of radial migration.
       radialMigrationRate=self%blackHoleBinarySeparationGrowthRate_%growthRate(blackHole)
       call blackHole%radialPositionRate(radialMigrationRate)
    end do
    return
  end subroutine blackHolesRadialMigrationDifferentialEvolution
  
  subroutine mergeBlackHoles(node,timeEnd)
    !!{
    Merge two black holes.
    !!}
    use :: Galacticus_Nodes                     , only : nodeComponentBlackHole , treeNode
    use :: Events_Black_Hole_Merger             , only : Event_Black_Hole_Merger
    use :: Nodes_Operators_Black_Holes_Utilities, only : blackHolesRecoilEscapes
    implicit none
    type            (treeNode              ), intent(inout), target   :: node
    double precision                        , intent(in   ), optional :: timeEnd
    class           (nodeComponentBlackHole)               , pointer  :: blackHole1      , blackHole2        , &
         &                                                               blackHolePrimary, blackHoleSecondary
    double precision                                                  :: massBlackHoleNew, spinBlackHoleNew  , &
         &                                                               massBlackHole1  , massBlackHole2    , &
         &                                                               velocityRecoil  , spinBlackHole1    , &
         &                                                               spinBlackHole2
    !$GLC attributes unused :: timeEnd
    
    ! Get the black holes.
    blackHole1 => node%blackHole(instance=              1)
    blackHole2 => node%blackHole(instance=mergingInstance)
    ! Process the merger to get the mass and spin of the merged black hole.
    call self_%blackHoleBinaryMerger_%merge(                   &
         &                                  blackHole2%mass(), &
         &                                  blackHole1%mass(), &
         &                                  blackHole2%spin(), &
         &                                  blackHole1%spin(), &
         &                                  massBlackHoleNew , &
         &                                  spinBlackHoleNew   &
         &                                 )
    ! Check which black hole is more massive in order to compute an appropriate recoil velocity.
    if (blackHole1%mass() >= blackHole2%mass()) then
       blackHolePrimary   => blackHole1
       blackHoleSecondary => blackHole2
    else
       blackHolePrimary   => blackHole2
       blackHoleSecondary => blackHole1
    end if
    massBlackHole1=blackHolePrimary  %mass()
    massBlackHole2=blackHoleSecondary%mass()
    spinBlackHole1=blackHolePrimary  %spin()
    spinBlackHole2=blackHoleSecondary%spin()
    ! Process the black hole merger.
    call Event_Black_Hole_Merger(blackHolePrimary,blackHoleSecondary,blackHole1)
    ! Calculate the recoil velocity of the binary black hole and check whether it escapes the galaxy
    velocityRecoil=self_%blackHoleBinaryRecoil_%velocity(blackHolePrimary,blackHoleSecondary)
    ! Compare the recoil velocity to the potential and determine whether the binary is ejected or stays in the galaxy.
    if     (                                                                                             &
         &  blackHolesRecoilEscapes(                                                                     &
         &                          node                 =node                                         , &
         &                          radius               =0.0d0                                        , &
         &                          radiusEscape         =self_%darkMatterHaloScale_%radiusVirial(node), &
         &                          velocityRecoil       =velocityRecoil                               , &
         &                          ignoreCentralBlackHole=.true.                                        &
         &                         )                                                                     &
         & ) then
       massBlackHoleNew=self_%blackHoleSeeds_%mass(node)
       spinBlackHoleNew=self_%blackHoleSeeds_%spin(node)
    end if
    ! Set the mass and spin of the central black hole.
    call blackHole1%massSet(massBlackHoleNew)
    call blackHole1%spinSet(spinBlackHoleNew)
    ! Remove the merging black hole from the list.
    call node%blackHoleRemove(mergingInstance)
    return
  end subroutine mergeBlackHoles
