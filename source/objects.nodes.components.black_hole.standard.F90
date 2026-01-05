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
Contains a module which implements the standard black hole node component.
!!}

module Node_Component_Black_Hole_Standard
  !!{
  Implement black hole tree node methods.
  !!}
  use :: Black_Hole_Binary_Initial_Separation, only : blackHoleBinaryInitialSeparationClass
  use :: Black_Hole_Binary_Mergers           , only : blackHoleBinaryMergerClass
  use :: Black_Hole_Binary_Recoil_Velocities , only : blackHoleBinaryRecoilClass
  use :: Black_Hole_Seeds                    , only : blackHoleSeedsClass
  use :: Cosmology_Parameters                , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScaleClass
  implicit none
  private
  public :: Node_Component_Black_Hole_Standard_Thread_Uninitialize, Node_Component_Black_Hole_Standard_Post_Evolve  , &
       &    Node_Component_Black_Hole_Standard_Thread_Initialize  , Node_Component_Black_Hole_Standard_Initialize   , &
       &    Node_Component_Black_Hole_Standard_State_Store        , Node_Component_Black_Hole_Standard_State_Restore, &
       &    Node_Component_Black_Hole_Standard_Scale_Set  

  !![
  <component>
   <class>blackHole</class>
   <name>standard</name>
   <isDefault>true</isDefault>
   <output instances="first"/>
   <properties>
    <property>
      <name>mass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
      <output unitsInSI="massSolar" comment="Mass of the black hole."/>
    </property>
    <property>
      <name>spin</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
      <getFunction>Node_Component_Black_Hole_Standard_Spin</getFunction>
      <output unitsInSI="0.0d0" comment="Spin of the black hole."/>
    </property>
    <property>
      <name>radialPosition</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
    <property>
      <name>tripleInteractionTime</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
   </properties>
   <bindings>
    <binding method="massDistribution" function="Node_Component_Black_Hole_Standard_Mass_Distribution"/>
    <binding method="massBaryonic"     function="Node_Component_Black_Hole_Standard_Mass_Baryonic"    />
   </bindings>
   <functions>objects.nodes.components.black_hole.standard.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(blackHoleBinaryRecoilClass           ), pointer :: blackHoleBinaryRecoil_
  class(blackHoleBinaryInitialSeparationClass), pointer :: blackHoleBinaryInitialSeparation_
  class(blackHoleBinaryMergerClass           ), pointer :: blackHoleBinaryMerger_
  class(blackHoleSeedsClass                  ), pointer :: blackHoleSeeds_
  class(darkMatterHaloScaleClass             ), pointer :: darkMatterHaloScale_
  !$omp threadprivate(blackHoleBinaryRecoil_,blackHoleBinaryInitialSeparation_,blackHoleBinaryMerger_,blackHoleSeeds_,darkMatterHaloScale_)

  ! Output options.
  logical          :: outputMergers

  ! Record of whether cold mode is explicitly tracked.
  logical          :: coldModeTracked
  ! A threadprivate object used to track to which thread events are attached.
  integer :: thread
  !$omp threadprivate(thread)

contains
  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Black_Hole_Standard_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Black_Hole_Standard_Initialize(parameters)
    !!{
    Initializes the standard black hole component module.
    !!}
    use :: Galacticus_Nodes, only : defaultHotHaloComponent, nodeComponentBlackHoleStandard
    use :: Input_Parameters, only : inputParameter         , inputParameters
    implicit none
    type(inputParameters               ), intent(inout) :: parameters
    type(inputParameters               )                :: subParameters

    ! Find our parameters.
    subParameters=parameters%subParameters('componentBlackHole')

    ! Get options controlling output.
    !![
    <inputParameter>
      <name>outputMergers</name>
      <defaultValue>.false.</defaultValue>
      <description>Determines whether or not properties of black hole mergers will be output.</description>
      <source>subParameters</source>
    </inputParameter>
    !!]

    ! Check if cold mode is explicitly tracked.
    coldModeTracked=defaultHotHaloComponent%massColdIsGettable()
    return
  end subroutine Node_Component_Black_Hole_Standard_Initialize
  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Black_Hole_Standard_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Black_Hole_Standard_Thread_Initialize(parameters)
    !!{
    Initializes the tree node standard black hole module.
    !!}
    use :: Events_Hooks              , only : satelliteMergerEvent     , openMPThreadBindingAtLevel, dependencyRegEx, dependencyDirectionBefore
    use :: Galacticus_Nodes          , only : defaultBlackHoleComponent
    use :: Input_Parameters          , only : inputParameter           , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters
    type(dependencyRegEx), dimension(1)  :: dependencies
    type(inputParameters)                :: subParameters
    if (defaultBlackHoleComponent%standardIsActive()) then
       dependencies(1)=dependencyRegEx(dependencyDirectionBefore,'^remnantStructure:')
       call satelliteMergerEvent%attach(thread,satelliteMerger,openMPThreadBindingAtLevel,label='nodeComponentBlackHoleStandard',dependencies=dependencies)
       ! Find our parameters.
       subParameters=parameters%subParameters('componentBlackHole')
       !![
       <objectBuilder class="blackHoleSeeds"                   name="blackHoleSeeds_"                   source="subParameters"/>
       <objectBuilder class="blackHoleBinaryRecoil"            name="blackHoleBinaryRecoil_"            source="subParameters"/>
       <objectBuilder class="blackHoleBinaryInitialSeparation" name="blackHoleBinaryInitialSeparation_" source="subParameters"/>
       <objectBuilder class="blackHoleBinaryMerger"            name="blackHoleBinaryMerger_"            source="subParameters"/>
       <objectBuilder class="darkMatterHaloScale"              name="darkMatterHaloScale_"              source="subParameters"/>
       !!]     
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Thread_Initialize
  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Black_Hole_Standard_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Black_Hole_Standard_Thread_Uninitialize()
    !!{
    Uninitializes the tree node standard black hole module.
    !!}
    use :: Events_Hooks    , only : satelliteMergerEvent
    use :: Galacticus_Nodes, only : defaultBlackHoleComponent
    implicit none
    if (defaultBlackHoleComponent%standardIsActive()) then
       if (satelliteMergerEvent%isAttached(thread,satelliteMerger)) call satelliteMergerEvent%detach(thread,satelliteMerger)
       !![
       <objectDestructor name="blackHoleSeeds_"                  />
       <objectDestructor name="blackHoleBinaryRecoil_"           />
       <objectDestructor name="blackHoleBinaryInitialSeparation_"/>
       <objectDestructor name="blackHoleBinaryMerger_"           />
       <objectDestructor name="darkMatterHaloScale_"             />
       !!]
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Thread_Uninitialize
  !![

  <scaleSetTask>
   <unitName>Node_Component_Black_Hole_Standard_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Black_Hole_Standard_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : defaultBlackHoleComponent, nodeComponentBlackHole, nodeComponentSpheroid, nodeComponentNSC, &
      &                             treeNode
    implicit none
    type            (treeNode              ), intent(inout), pointer :: node
    double precision                        , parameter              :: scaleMassRelative=1.0d-4
    double precision                        , parameter              :: scaleSizeRelative=1.0d-4
    double precision                        , parameter              :: scaleSizeAbsolute=1.0d-6
    double precision                        , parameter              :: scaleMassAbsolute=1.0d+0
    class           (nodeComponentNSC      )               , pointer :: nuclearStarCluster
    class           (nodeComponentSpheroid )               , pointer :: spheroid
    class           (nodeComponentBlackHole)               , pointer :: blackHole
    integer                                                          :: instance
    double precision                                                 :: massStellar

    ! Determine if the standard implementation is active and at least one black hole exists.
    if (defaultBlackHoleComponent%standardIsActive().and.node%blackHoleCount() > 0) then
       ! Get the spheroid component.
       spheroid => node%spheroid()
       ! Get the nuclear star cluster component.
       nuclearStarCluster => node%NSC()
       ! Loop over instances.
       do instance=1,node%blackHoleCount()
          ! Get the black hole.
          blackHole => node%blackHole(instance=instance)
          ! Set scale for mass. We set an absolute mass scale that depends on the stellar masses of spheroid and nuclear star
          ! cluster components, if present. We use the smaller of the two masses (if both are present) to ensure accurate tracking
          ! of the black hole as it couples to each of these components.
          if      (spheroid%massStellar() > 0.0d0 .and. nuclearStarCluster%massStellar() > 0.0d0) then
             massStellar=min(spheroid%massStellar(),nuclearStarCluster%massStellar())
          else if (spheroid%massStellar() > 0.0d0                                               ) then
             massStellar=    spheroid%massStellar()
          else
             massStellar=                           nuclearStarCluster%massStellar()
          end if
          call blackHole%massScale(                                            &
               &                   maxval(                                     &
               &                          [                                    &
               &                           scaleMassAbsolute                 , &
               &                           scaleMassRelative*massStellar     , &
               &                                             blackHole%mass()  &
               &                           ]                                   &
               &                      )                                        &
               &                  )
          ! Set scale for spin.
          call blackHole%spinScale(1.0d0)
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
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Scale_Set

  subroutine satelliteMerger(self,node)
    !!{
    Merge any black hole associated with {\normalfont \ttfamily node} before it merges with its host halo.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBlackHole , treeNode
    use :: Events_Black_Hole_Merger, only : Event_Black_Hole_Merger
    implicit none
    class           (*                     ), intent(inout) :: self
    type            (treeNode              ), intent(inout) :: node
    type            (treeNode              ), pointer       :: hostNode
    class           (nodeComponentBlackHole), pointer       :: blackHoleHostCentral, blackHole         , &
         &                                                     blackHolePrimary    , blackHoleSecondary
    integer                                                 :: instance
    double precision                                        :: massBlackHoleNew    , spinBlackHoleNew  , &
         &                                                     massBlackHole1      , massBlackHole2    , &
         &                                                     radiusInitial       , velocityRecoil    , &
         &                                                     spinBlackHole1      , spinBlackHole2
    !$GLC attributes unused :: self
    
    ! Find the node to merge with.
    hostNode => node%mergesWith()
    ! Find the initial radius of the satellite black hole in the remnant.
    radiusInitial=blackHoleBinaryInitialSeparation_%separationInitial(node,hostNode)
    ! If the separation is non-positive, assume that the black holes merge instantaneously.
    if (radiusInitial <= 0.0d0) then
       ! Get the central black hole of the host galaxy.
       blackHoleHostCentral => hostNode%blackHole(instance=1,autoCreate=.true.)
       ! Loop over all black holes in the satellite galaxy.
       do instance=1,node%blackHoleCount()
          ! Get the black hole.
          blackHole => node%blackHole(instance=instance)
          ! Compute the outcome of the merger.
          call blackHoleBinaryMerger_%merge(                             &
               &                            blackHole           %mass(), &
               &                            blackHoleHostCentral%mass(), &
               &                            blackHole           %spin(), &
               &                            blackHoleHostCentral%spin(), &
               &                            massBlackHoleNew           , &
               &                            spinBlackHoleNew             &
               &                           )
          ! Merge the black holes instantaneously.
          ! Check which black hole is more massive in order to compute an appropriate recoil velocity
          if (blackHoleHostCentral%mass() >= blackHole%mass()) then
             blackHolePrimary   => blackHoleHostCentral
             blackHoleSecondary => blackHole
          else
             blackHolePrimary   => blackHole
             blackHoleSecondary => blackHoleHostCentral
          end if
          massBlackHole1=blackHolePrimary  %mass()
          massBlackHole2=blackHoleSecondary%mass()
          spinBlackHole1=blackHolePrimary  %spin()
          spinBlackHole2=blackHoleSecondary%spin()
          ! Process the black hole merger.
          call Event_Black_Hole_Merger(blackHolePrimary,blackHoleSecondary,blackHoleHostCentral)
          ! Now calculate the recoil velocity of the binary black hole and check whether it escapes the galaxy.
          velocityRecoil=blackHoleBinaryRecoil_%velocity(blackHolePrimary,blackHoleSecondary)
          if (Node_Component_Black_Hole_Standard_Recoil_Escapes(node,velocityRecoil,radius=0.0d0,ignoreCentralBlackHole=.true.)) then
             massBlackHoleNew=0.0d0
             spinBlackHoleNew=0.0d0
          end if
          ! Set the new black hole mass in the host.
          call Node_Component_Black_Hole_Standard_Output_Merger(node,massBlackHole1,massBlackHole2)
          call blackHoleHostCentral%massSet(massBlackHoleNew)
          call blackHoleHostCentral%spinSet(spinBlackHoleNew)
       end do
    else
       ! Adjust the radii of the black holes in the satellite galaxy.
       do instance=node%blackHoleCount(),1,-1
          blackHole => node%blackHole(instance=instance)
          call blackHole%radialPositionSet(radiusInitial)
          ! Declares them as not having interacted in a triple black hole interaction.
          call blackHole%tripleInteractionTimeSet(0.0d0)
          ! Remove this black hole if it has no mass.
          if (blackHole%mass() <= 0.0d0) call node%blackHoleRemove(instance)
       end do
       ! Move black holes from the satellite to the host.
       call node%blackHoleMove(hostNode)
    end if
    return
  end subroutine satelliteMerger
  
  logical function Node_Component_Black_Hole_Standard_Recoil_Escapes(node,velocityRecoil,radius,ignoreCentralBlackHole)
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

    ! Return false immediately if the recoil velocity is zero.
    if (velocityRecoil <= 0.0d0) then
       Node_Component_Black_Hole_Standard_Recoil_Escapes=.false.
       return
    end if
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
    Node_Component_Black_Hole_Standard_Recoil_Escapes= &
         &  +0.5d0*velocityRecoil**2                   &
         &  +      potential                           &
         & >                                           &
         &  +      potentialSelf
    return
  end function Node_Component_Black_Hole_Standard_Recoil_Escapes

  subroutine Node_Component_Black_Hole_Standard_Output_Merger(node,massBlackHole1,massBlackHole2)
    !!{
    Outputs properties of merging black holes.
    !!}
    use :: Output_HDF5     , only : outputFile
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    use :: HDF5_Access     , only : hdf5Access
    use :: IO_HDF5         , only : hdf5Object
    implicit none
    type            (treeNode          ), intent(inout) :: node
    double precision                    , intent(in   ) :: massBlackHole1    , massBlackHole2
    class           (nodeComponentBasic), pointer       :: basic
    type            (hdf5Object        )                :: mergersGroup
    ! Exit if merger data is not to be output.
    if (.not.outputMergers) return
    ! Ignore mergers with zero mass black holes.
    if (massBlackHole2 <= 0.0d0    ) return
    ! Get the basic component.
    basic => node%basic()
    ! Open the group to which black hole mergers should be written.
    !$ call hdf5Access%set()
    mergersGroup=outputFile%openGroup("blackHoleMergers","Black hole mergers data.")
    ! Append to the datasets.
    call    mergersGroup%writeDataset([massBlackHole1            ],"massBlackHole1","Mass of the first merging black hole." ,appendTo=.true.)
    call    mergersGroup%writeDataset([massBlackHole2            ],"massBlackHole2","Mass of the second merging black hole.",appendTo=.true.)
    call    mergersGroup%writeDataset([basic%time()              ],"timeOfMerger"  ,"The time of the black hole merger."    ,appendTo=.true.)
    call    mergersGroup%writeDataset([node%hostTree%volumeWeight],"volumeWeight"  ,"The weight for the black hole merger." ,appendTo=.true.)
    !$ call hdf5Access  %unset       (                                                                                                      )
    return
  end subroutine Node_Component_Black_Hole_Standard_Output_Merger

  !![
  <postStepTask>
    <unitName>Node_Component_Black_Hole_Standard_Post_Evolve</unitName>
  </postStepTask>
  !!]
  subroutine Node_Component_Black_Hole_Standard_Post_Evolve(node,status)
    !!{
    Keep black hole spin in physical range.
    !!}
    use :: Galacticus_Nodes, only : defaultBlackHoleComponent, nodeComponentBlackHole, treeNode
    use :: Interface_GSL   , only : GSL_Success              , GSL_Continue
    implicit none
    type            (treeNode              ), intent(inout), pointer :: node
    integer                                 , intent(inout)          :: status
    class           (nodeComponentBlackHole)               , pointer :: blackHole
    double precision                        , parameter              :: spinMaximum=0.9999d0
    integer                                                          :: i                   , instanceCount
    double precision                                                 :: spin
    ! Check if the standard component is active.
    if (defaultBlackHoleComponent%standardIsActive()) then
       ! Get a count of the number of black holes associated with this node.
       instanceCount=node%blackHoleCount()
       if (instanceCount > 0) then
          ! Iterate over instances.
          do i=1,instanceCount
             ! Get the black hole component.
             blackHole => node%blackHole(instance=i)
             if (blackHole%spin() > spinMaximum .or. blackHole%spin() < 0.0d0) then
                ! Note that "status" is not set to failure as this change in state of the black hole should not change any
                ! calculation of differential evolution rates as the spin was in an unphysical regime anyway
                spin=max(min(blackHole%spin(),spinMaximum),0.0d0)
                call blackHole%spinSet(spin)
                ! Indicate that ODE evolution should continue after this state change.
                if (status == GSL_Success) status=GSL_Continue
             end if
             if (blackHole%mass() < 0.0d0) then
                ! Note that "status" is not set to failure as this change in state of the black hole should not change any
                ! calculation of differential evolution rates as the mass was in an unphysical regime anyway.
                call blackHole%massSet(blackHoleSeeds_%mass(node))
                ! Indicate that ODE evolution should continue after this state change.
                if (status == GSL_Success) status=GSL_Continue
             end if
          end do
       end if
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Post_Evolve

  !![
  <stateStoreTask>
   <unitName>Node_Component_Black_Hole_Standard_State_Store</unitName>
  </stateStoreTask>
  !!]

  subroutine Node_Component_Black_Hole_Standard_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentBlackHole -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="blackHoleBinaryRecoil_ blackHoleBinaryInitialSeparation_ blackHoleBinaryMerger_ blackHoleSeeds_ darkMatterHaloScale_"/>
    !!]
    return
  end subroutine Node_Component_Black_Hole_Standard_State_Store
  !![
  <stateRetrieveTask>
   <unitName>Node_Component_Black_Hole_Standard_State_Restore</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_Black_Hole_Standard_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentBlackHole -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="blackHoleBinaryRecoil_ blackHoleBinaryInitialSeparation_ blackHoleBinaryMerger_ blackHoleSeeds_ darkMatterHaloScale_"/>
    !!]
    return
  end subroutine Node_Component_Black_Hole_Standard_State_Restore

end module Node_Component_Black_Hole_Standard
