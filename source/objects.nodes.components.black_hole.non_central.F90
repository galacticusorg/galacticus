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
Contains a module which implements the standard black hole node component.
!!}

module Node_Component_Black_Hole_Noncentral
  !!{
  Implement black hole tree node methods.
  !!}
  use :: Black_Hole_Binary_Mergers          , only : blackHoleBinaryMergerClass
  use :: Black_Hole_Binary_Recoil_Velocities, only : blackHoleBinaryRecoilClass
  use :: Black_Hole_Binary_Separations      , only : blackHoleBinarySeparationGrowthRateClass
  use :: Black_Hole_Seeds                   , only : blackHoleSeedsClass
  use :: Dark_Matter_Halo_Scales            , only : darkMatterHaloScaleClass
  implicit none
  private
  public :: Node_Component_Black_Hole_Noncentral_Rate_Compute       , Node_Component_Black_Hole_Noncentral_Scale_Set        , &
       &    Node_Component_Black_Hole_Noncentral_Initialize         , Node_Component_Black_Hole_Noncentral_Thread_Initialize, &
       &    Node_Component_Black_Hole_Noncentral_Thread_Uninitialize, Node_Component_Black_Hole_Noncentral_State_Store      , &
       &    Node_Component_Black_Hole_Noncentral_State_Restore

  !![
  <component>
   <class>blackHole</class>
   <name>nonCentral</name>
   <extends>
    <class>blackHole</class>
    <name>standard</name>
   </extends>
   <isDefault>false</isDefault>
   <output instances="first"/>
   <properties>
    <property>
      <name>radialPosition</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
    </property>
   </properties>
  </component>
  !!]

  ! Objects used by this component.
  class(darkMatterHaloScaleClass                ), pointer :: darkMatterHaloScale_
  class(blackHoleBinaryRecoilClass              ), pointer :: blackHoleBinaryRecoil_
  class(blackHoleBinaryMergerClass              ), pointer :: blackHoleBinaryMerger_
  class(blackHoleBinarySeparationGrowthRateClass), pointer :: blackHoleBinarySeparationGrowthRate_
  class(blackHoleSeedsClass                     ), pointer :: blackHoleSeeds_
  !$omp threadprivate(darkMatterHaloScale_,blackHoleBinaryRecoil_,blackHoleBinaryMerger_,blackHoleBinarySeparationGrowthRate_,blackHoleSeeds_)

  ! Option specifying whether the triple black hole interaction should be used.
  logical :: tripleInteraction

  ! Index of black hole instance about to merge.
  integer :: mergingInstance
  !$omp threadprivate(mergingInstance)

  ! Index of black hole involved in three-body interactions
  integer :: binaryInstance           , tripleInstance
  !$omp threadprivate(binaryInstance,tripleInstance)

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Black_Hole_Noncentral_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Black_Hole_Noncentral_Initialize(parameters)
    !!{
    Initializes the noncentral black hole component module.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters
    type(inputParameters)                :: subParameters

    ! Find our parameters.
    subParameters=parameters%subParameters('componentBlackHole')
    !![
    <inputParameter>
      <name>tripleInteraction</name>
      <defaultValue>.false.</defaultValue>
      <description>Determines whether or not triple black hole interactions will be accounted for.</description>
      <source>subParameters</source>
    </inputParameter>
    !!]
    return
  end subroutine Node_Component_Black_Hole_Noncentral_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Black_Hole_Noncentral_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Black_Hole_Noncentral_Thread_Initialize(parameters)
    !!{
    Initializes the tree node random spin module.
    !!}
    use :: Galacticus_Nodes, only : defaultBlackHoleComponent
    use :: Input_Parameters, only : inputParameter           , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters
    type(inputParameters)                :: subParameters

    if (defaultBlackHoleComponent%noncentralIsActive()) then
       ! Find our parameters.
       subParameters=parameters%subParameters('componentBlackHole')
       !![
       <objectBuilder class="darkMatterHaloScale"                 name="darkMatterHaloScale_"                 source="subParameters"/>
       <objectBuilder class="blackHoleBinaryRecoil"               name="blackHoleBinaryRecoil_"               source="subParameters"/>
       <objectBuilder class="blackHoleBinaryMerger"               name="blackHoleBinaryMerger_"               source="subParameters"/>
       <objectBuilder class="blackHoleBinarySeparationGrowthRate" name="blackHoleBinarySeparationGrowthRate_" source="subParameters"/>
       <objectBuilder class="blackHoleSeeds"                      name="blackHoleSeeds_"                      source="subParameters"/>
       !!]
    end if
    return
  end subroutine Node_Component_Black_Hole_Noncentral_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Black_Hole_Noncentral_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Black_Hole_Noncentral_Thread_Uninitialize()
    !!{
    Uninitializes the tree node random spin module.
    !!}
    use :: Galacticus_Nodes, only : defaultBlackHoleComponent
    implicit none

    if (defaultBlackHoleComponent%noncentralIsActive()) then
       !![
       <objectDestructor name="darkMatterHaloScale_"                />
       <objectDestructor name="blackHoleBinaryRecoil_"              />
       <objectDestructor name="blackHoleBinaryMerger_"              />
       <objectDestructor name="blackHoleBinarySeparationGrowthRate_"/>
       <objectDestructor name="blackHoleSeeds_"                     />
       !!]
    end if
    return
  end subroutine Node_Component_Black_Hole_Noncentral_Thread_Uninitialize

  !![
  <rateComputeTask>
   <unitName>Node_Component_Black_Hole_Noncentral_Rate_Compute</unitName>
  </rateComputeTask>
  !!]
  subroutine Node_Component_Black_Hole_Noncentral_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !!{
    Compute the black hole node mass rate of change.
    !!}
    use :: Galacticus_Nodes                , only : defaultBlackHoleComponent     , interruptTask, nodeComponentBlackHole, propertyInactive, &
          &                                         treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    type            (treeNode              ), intent(inout)          :: node
    logical                                 , intent(inout)          :: interrupt
    procedure       (interruptTask         ), intent(inout), pointer :: interruptProcedure
    integer                                 , intent(in   )          :: propertyType
    class           (nodeComponentBlackHole)               , pointer :: blackHoleBinary   , blackHoleCentral   , &
         &                                                              blackHole
    integer                                                          :: iInstance         , instanceCount
    double precision                                                 :: binaryRadius      , radialMigrationRate, &
         &                                                              radiusHardBinary
    logical                                                          :: binaryRadiusFound

    ! Return immediately if inactive variables are requested.
    if (propertyInactive(propertyType)) return
    if (defaultBlackHoleComponent%noncentralIsActive()) then
       ! Get a count of the number of black holes associated with this node.
       instanceCount=node%blackHoleCount()
       ! Get the central black hole.
       blackHoleCentral => node%blackHole(instance=1)
       ! Do radial migration for non-central black holes.
       do iInstance=2,instanceCount
          ! Get the black hole.
          blackHole => node%blackHole(instance=iInstance)
          ! Compute the hard binary radius.
          radiusHardBinary= (                                                &
               &              gravitationalConstant_internal                 &
               &             *(                                              &
               &               +blackHoleCentral%mass()                      &
               &               +blackHole       %mass()                      &
               &              )                                              &
               &            )                                                &
               &           /(                                                &
               &              4.0d0                                          &
               &             *  darkMatterHaloScale_%velocityVirial(node)**2 &
               &            )
          ! Places a new black hole in the center of the galaxy in case there is no central one.
          if     (                                                       &
               &   blackHoleCentral%mass          () == 0.0d0            &
               &  .and.                                                  &
               &   blackHole       %radialPosition() <= radiusHardBinary &
               & ) then
             mergingInstance=iInstance
             interrupt=.true.
             interruptProcedure => Node_Component_Black_Hole_Noncentral_Merge_Black_Holes
             return
          end if
          ! Check for a black hole that is about to merge.
          if (blackHole%radialPosition() <= 0.0d0) then
             ! Record which instance is merging, then trigger an interrupt.
             mergingInstance=iInstance
             interrupt=.true.
             interruptProcedure => Node_Component_Black_Hole_Noncentral_Merge_Black_Holes
             return
          end if
          ! Set the rate of radial migration.
          radialMigrationRate=blackHoleBinarySeparationGrowthRate_%growthRate(blackHole)
          call blackHole%radialPositionRate(radialMigrationRate)
       end do
       ! Loop over black holes, testing for triple black hole interactions. Find the three closest black holes then check if a
       ! three body interaction occurs using the radial condition derived in Hoffman and Loeb (2007).
       binaryRadiusFound=.false.
       binaryRadius     =huge(1.0d0)
       if (tripleInteraction) then
          if (instanceCount >= 3 .and. blackHoleCentral%mass() > 0.0d0) then
             do iInstance=2,instanceCount
                ! Get the black hole.
                blackHole => node%blackHole(instance=iInstance)
                if     (                                                      &
                     &  (          blackHole%radialPosition() <= binaryRadius &
                     &   .or. .not.binaryRadiusFound                          &
                     &  )                                                     &
                     &  .and. .not.blackHole%mass          () <= 0.0d0        &
                     &  .and. .not.blackHole%radialPosition() <= 0.0d0        &
                     & ) then
                   binaryRadius     =blackHole%radialPosition()
                   binaryInstance   =iInstance
                   binaryRadiusFound=.true.
                end if
             end do
             if (binaryRadiusFound) then
                ! Get the binary black hole.
                blackHoleBinary => node%blackHole(instance=binaryInstance)
                ! Compute the hard binary radius.
                radiusHardBinary= (                                              &
                     &              gravitationalConstant_internal               &
                     &             *(                                            &
                     &                blackHoleCentral%mass()                    &
                     &               + blackHoleBinary%mass()                    &
                     &              )                                            &
                     &            )                                              &
                     &           /(                                              &
                     &              4.0d0                                        &
                     &             *darkMatterHaloScale_%velocityVirial(node)**2 &
                     &            )
                ! Search for a third black hole.
                do iInstance=2,instanceCount
                   ! Get the black hole.
                   blackHole => node%blackHole(instance=iInstance)
                   if     (      .not. iInstance                         == binaryInstance   &
                        &  .and. .not. blackHole%mass                 () <= 0.0d0            &
                        &  .and. .not. blackHole%radialPosition       () <= 0.0d0            &
                        &  .and.       blackHole%radialPosition       () <= radiusHardBinary &
                        &  .and.       blackHole%tripleInteractionTime() == 0.0d0            &
                        & ) then
                      tripleInstance=iInstance
                      interrupt=.true.
                      interruptProcedure => Node_Component_Black_Hole_Noncentral_Triple_Interaction
                      return
                   end if
                end do
             end if
          end if
       end if
    end if
    return
  end subroutine Node_Component_Black_Hole_Noncentral_Rate_Compute

  !![
  <scaleSetTask>
   <unitName>Node_Component_Black_Hole_Noncentral_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Black_Hole_Noncentral_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : defaultBlackHoleComponent, nodeComponentBlackHole, nodeComponentSpheroid, treeNode
    implicit none
    type            (treeNode              ), intent(inout), pointer :: node
    double precision                        , parameter              :: scaleSizeRelative=1.0d-4
    double precision                        , parameter              :: scaleSizeAbsolute=1.0d-6
    class           (nodeComponentSpheroid )               , pointer :: spheroid
    class           (nodeComponentBlackHole)               , pointer :: blackHole
    integer                                                          :: instance

    ! Determine if the noncentral implementation is active and at least one black hole exists.
    if (defaultBlackHoleComponent%noncentralIsActive().and.node%blackHoleCount() > 0) then
       ! Get the spheroid component.
       spheroid => node%spheroid()
       ! Loop over instances.
       do instance=1,node%blackHoleCount()
          ! Get the black hole.
          blackHole => node%blackHole(instance=instance)
          ! Set scale for radius.
          call blackHole%radialPositionScale(                                                    &
               &                             maxval(                                             &
               &                                  [                                              &
               &                                   scaleSizeAbsolute,                            &
               &                                   scaleSizeRelative*spheroid %halfMassRadius(), &
               &                                                     blackHole%radialPosition()  &
               &                                  ]                                              &
               &                                   )                                             &
               &                            )
       end do
    end if
    return
  end subroutine Node_Component_Black_Hole_Noncentral_Scale_Set

  subroutine Node_Component_Black_Hole_Noncentral_Merge_Black_Holes(node,timeEnd)
    !!{
    Merge two black holes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, treeNode
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
    call blackHoleBinaryMerger_%merge(blackHole2%mass(), &
         &                            blackHole1%mass(), &
         &                            blackHole2%spin(), &
         &                            blackHole1%spin(), &
         &                            massBlackHoleNew , &
         &                            spinBlackHoleNew   &
         &                           )
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
    ! Calculate the recoil velocity of the binary black hole and check whether it escapes the galaxy
    velocityRecoil=blackHoleBinaryRecoil_%velocity(blackHolePrimary,blackHoleSecondary)
    ! Compare the recoil velocity to the potential and determine whether the binary is ejected or stays in the galaxy.
    if (Node_Component_Black_Hole_Noncentral_Recoil_Escapes(node,velocityRecoil,radius=0.0d0,ignoreCentralBlackHole=.true.)) then
       massBlackHoleNew=blackHoleSeeds_%mass(node)
       spinBlackHoleNew=blackHoleSeeds_%spin(node)
    end if
    ! Set the mass and spin of the central black hole.
    call blackHole1%massSet(massBlackHoleNew)
    call blackHole1%spinSet(spinBlackHoleNew)
    ! Remove the merging black hole from the list.
    call node%blackHoleRemove(mergingInstance)
    return
  end subroutine Node_Component_Black_Hole_Noncentral_Merge_Black_Holes

  subroutine Node_Component_Black_Hole_Noncentral_Triple_Interaction(node,timeEnd)
    !!{
    Handles triple black holes interactions, using conditions similar to those of \cite{volonteri_assembly_2003}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic            , nodeComponentBlackHole, treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    type            (treeNode              ), intent(inout), target   :: node
    double precision                        , intent(in   ), optional :: timeEnd
    class           (nodeComponentBasic    )               , pointer  :: basic
    class           (nodeComponentBlackHole)               , pointer  :: blackHoleBinary          , blackHoleCentral           , &
         &                                                               ejectedBlackHoleComponent, newBinaryBlackHoleComponent, &
         &                                                               tripleBlackHoleComponent
    integer                                                           :: ejectedInstance          , newBinaryInstance
    double precision                                                  :: bindingEnergy            , kineticEnergyChange        , &
         &                                                               massBinary               , massEjected                , &
         &                                                               massRatioIntruder        , newRadius                  , &
         &                                                               velocityBinary           , velocityEjected
    logical                                                           :: removeBinary             , removeEjected
    !$GLC attributes unused :: timeEnd

    ! Get the basic component.
    basic            => node%basic    (                       )
    ! Get the black holes.
    blackHoleCentral => node%blackHole(instance=             1)
    blackHoleBinary  => node%blackHole(instance=binaryInstance)
    tripleBlackHoleComponent  => node%blackHole(instance=tripleInstance)
    ! We have to distinguish two cases, where a different black hole is ejected, the one with the lowest mass.
    massRatioIntruder=+  tripleBlackHoleComponent%mass() &
         &            /(                                 &
         &              +blackHoleCentral        %mass() &
         &              +blackHoleBinary         %mass() &
         &             )
    ! Set the triple interaction time for the triple black hole.
    call tripleBlackHoleComponent%tripleInteractionTimeSet(basic%time())
    ! Branch on intruder mass ratio.
    if (massRatioIntruder <= 2.0d0) then
       if (tripleBlackHoleComponent%mass() <= blackHoleBinary%mass()) then
          newRadius           = blackHoleBinary%radialPosition()/(1.0d0+0.4d0*massRatioIntruder)
          call blackHoleBinary%radialPositionSet(newRadius)
          bindingEnergy       =+gravitationalConstant_internal              &
               &               *(                                           &
               &                 +tripleBlackHoleComponent%mass          () &
               &                 *blackHoleCentral        %mass          () &
               &                )                                           &
               &               /  tripleBlackHoleComponent%radialPosition()
          kineticEnergyChange=0.4d0*massRatioIntruder*bindingEnergy
          ejectedInstance    =tripleInstance
          newBinaryInstance  =binaryInstance
          ejectedBlackHoleComponent   => tripleBlackHoleComponent
          newBinaryBlackHoleComponent => blackHoleBinary
       else
          newRadius          = tripleBlackHoleComponent%radialPosition()/(1.0d0+0.4d0*massRatioIntruder )
          call tripleBlackHoleComponent%radialPositionSet(newRadius)
          bindingEnergy      =+gravitationalConstant_internal      &
               &              *(                                   &
               &                +blackHoleBinary %mass          () &
               &                *blackHoleCentral%mass          () &
               &               )                                   &
               &              /  blackHoleBinary %radialPosition()
          kineticEnergyChange=0.4d0*massRatioIntruder*bindingEnergy
          ejectedInstance    =binaryInstance
          newBinaryInstance  =tripleInstance
          ejectedBlackHoleComponent   => blackHoleBinary
          newBinaryBlackHoleComponent => tripleBlackHoleComponent
       end if
    else
       ! This latter case can be referred to as head-on collision.
       newRadius             =0.53d0*tripleBlackHoleComponent%radialPosition()
       call tripleBlackHoleComponent%radialPositionSet(newRadius)
       bindingEnergy         =+gravitationalConstant_internal      &
            &                 *(                                   &
            &                   +blackHoleBinary %mass          () &
            &                   *blackHoleCentral%mass          () &
            &                  )                                   &
            &                 /  blackHoleBinary %radialPosition()
       kineticEnergyChange         =  0.9d0*massRatioIntruder*bindingEnergy
       ejectedInstance             =  binaryInstance
       newBinaryInstance           =  tripleInstance
       ejectedBlackHoleComponent   => blackHoleBinary
       newBinaryBlackHoleComponent => tripleBlackHoleComponent
    end if
    ! First we find the lightest black hole and tag it as being ejected.
    massEjected= ejectedBlackHoleComponent  %mass()
    massBinary = newBinaryBlackHoleComponent%mass() &
         &      +blackHoleCentral           %mass()
    velocityEjected=sqrt(kineticEnergyChange/(1.0d0+massEjected/massBinary )/massEjected*2.0d0)
    velocityBinary =sqrt(kineticEnergyChange/(1.0d0+massBinary /massEjected)/massBinary *2.0d0)
    ! Determine whether the ejected black hole is actually ejected.
    removeEjected=Node_Component_Black_Hole_Noncentral_Recoil_Escapes(node,velocityEjected,ejectedBlackHoleComponent%radialPosition(),ignoreCentralBlackHole=.false.)
    ! Determine whether the binary black hole is ejected.
    removeBinary=Node_Component_Black_Hole_Noncentral_Recoil_Escapes(node,velocityBinary,newBinaryBlackHoleComponent%radialPosition(),ignoreCentralBlackHole=.true. )
    ! Remove the binary black hole from the list if required.
    if (removeBinary) then
       ! Set the central black hole as a zero mass component.
       call blackHoleCentral%massSet(0.0d0)
       call blackHoleCentral%spinSet(0.0d0)
       ! Remove the binary black hole.
       call node%blackHoleRemove(newBinaryInstance)
       ! If this removal has changed the position of the ejected black hole in the list then update its index.
       if (ejectedInstance > newBinaryInstance) ejectedInstance=ejectedInstance-1
    end if
    ! Remove the ejected black hole from the list if required.
    if (removeEjected) call node%blackHoleRemove(ejectedInstance)
    return
  end subroutine Node_Component_Black_Hole_Noncentral_Triple_Interaction

  logical function Node_Component_Black_Hole_Noncentral_Recoil_Escapes(node,velocityRecoil,radius,ignoreCentralBlackHole)
    !!{
    Return true if the given recoil velocity is sufficient to eject a black hole from the halo.
    !!}
    use :: Coordinates               , only : coordinateSpherical   , assignment(=)
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Galactic_Structure_Options, only : componentTypeBlackHole
    use :: Galacticus_Nodes          , only : treeNode
    implicit none
    type            (treeNode             ), intent(inout) :: node
    double precision                       , intent(in   ) :: velocityRecoil        , radius
    logical                                , intent(in   ) :: ignoreCentralBlackHole
    class           (massDistributionClass), pointer       :: massDistribution_
    double precision                                       :: potential             , potentialSelf
    type            (coordinateSpherical  )                :: coordinates           , coordinatesVirial
    
    ! Compute relevant potentials.
    coordinates       =  [                     radius            ,0.0d0,0.0d0]
    coordinatesVirial =  [darkMatterHaloScale_%radiusVirial(node),0.0d0,0.0d0]
    massDistribution_ => node             %massDistribution   (                             )
    potential         =  massDistribution_%potentialDifference(coordinates,coordinatesVirial)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    if (ignoreCentralBlackHole) then
       ! Compute potential of central black hole to be subtracted off of total value.
       massDistribution_ => node             %massDistribution   (componentType=componentTypeBlackHole                  )
       potentialSelf     =  massDistribution_%potentialDifference(              coordinates           ,coordinatesVirial)
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
    else
       ! No correction for central black hole as it is to be included.
       potentialSelf=0.0d0
    end if
    ! Evaluate the escape condition.
    Node_Component_Black_Hole_Noncentral_Recoil_Escapes= &
         &  +0.5d0*velocityRecoil**2                     &
         &  +      potential                             &
         & >                                             &
         &  +      potentialSelf
    return
  end function Node_Component_Black_Hole_Noncentral_Recoil_Escapes

  !![
  <stateStoreTask>
   <unitName>Node_Component_Black_Hole_NonCentral_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_Black_Hole_NonCentral_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentBlackHole -> nonCentral',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="darkMatterHaloScale_ blackHoleBinaryRecoil_ blackHoleBinaryMerger_ blackHoleBinarySeparationGrowthRate_ blackHoleSeeds_"/>
    !!]
    return
  end subroutine Node_Component_Black_Hole_NonCentral_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Node_Component_Black_Hole_NonCentral_State_Restore</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_Black_Hole_NonCentral_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentBlackHole -> nonCentral',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="darkMatterHaloScale_ blackHoleBinaryRecoil_ blackHoleBinaryMerger_ blackHoleBinarySeparationGrowthRate_ blackHoleSeeds_"/>
    !!]
    return
  end subroutine Node_Component_Black_Hole_NonCentral_State_Restore

end module Node_Component_Black_Hole_Noncentral
