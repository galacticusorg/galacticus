!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!+ Contributions to this file made by: Andrew Benson, Christoph Behrens.

!!{
Contains a module implementing a halo angular momentum component using the approach of
\cite{vitvitska_origin_2002}.
!!}

module Node_Component_Halo_Angular_Momentum_Vitvitska
  !!{
  Implements a node spin component using the approach of \cite{vitvitska_origin_2002}.
  !!}
  use :: Dark_Matter_Halo_Scales           , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profile_Scales        , only : darkMatterProfileScaleRadiusClass
  use :: Dark_Matter_Profiles_DMO          , only : darkMatterProfileDMOClass
  use :: Halo_Spin_Distributions           , only : haloSpinDistributionClass
  use :: Virial_Orbits                     , only : virialOrbitClass
  use :: Merger_Trees_Build_Mass_Resolution, only : mergerTreeMassResolutionClass
  implicit none
  private
  public :: Node_Component_Halo_Angular_Momentum_Vitvitska_Bindings     , Node_Component_Halo_Angular_Momentum_Vitvitska_Scale_Set  , &
       &    Node_Component_Halo_Angular_Momentum_Vitvitska_Rate_Compute , Node_Component_Halo_Angular_Momentum_Vitvitska_Thread_Init, &
       &    Node_Component_Halo_Angular_Momentum_Vitvitska_Thread_Uninit, Node_Component_Halo_Angular_Momentum_Vitvitska_Initialize

  !![
  <component>
   <class>spin</class>
   <name>vitvitska</name>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>angularMomentum</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="get" />
      <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum magnitude of the DMO halo."/>
    </property>
    <property>
      <name>angularMomentumVector</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="get" />
      <output labels="[X,Y,Z]" unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum vector of the DMO halo."/>
      <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
    </property>
    <property>
      <name>angularMomentumAccretionRate</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
    </property>
    <property>
      <name>angularMomentumGrowthRate</name>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" isDeferred="get" />
      <type>double</type>
      <rank>0</rank>
    </property>
   </properties>
  </component>
  !!]

  ! Objects used by this component.
  class(haloSpinDistributionClass        ), pointer :: haloSpinDistribution_
  class(darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_
  class(virialOrbitClass                 ), pointer :: virialOrbit_
  class(darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_
  class(darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_
  class(mergerTreeMassResolutionClass    ), pointer :: mergerTreeMassResolution_
  !$omp threadprivate(haloSpinDistribution_,darkMatterProfileDMO_,virialOrbit_,darkMatterHaloScale_,darkMatterProfileScaleRadius_,mergerTreeMassResolution_)

  ! Parameter controlling scaling of orbital angular momentum with mass ratio.
  double precision :: angularMomentumVitvitskaMassExponent

  ! Parameter controlling differential evolution of angular momentum.
  logical          :: angularMomentumVitvitskaEvolveDifferentially

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Vitvitska_Bindings</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Vitvitska_Bindings(parameters_)
    !!{
    Initializes the ``Vitvitska'' implementation of the spin component.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpinVitvitska
    use :: Input_Parameters, only : inputParameter            , inputParameters
    implicit none
    type(inputParameters           ), intent(inout) :: parameters_
    type(nodeComponentSpinVitvitska)                :: spin

    !![
    <inputParameter>
      <name>angularMomentumVitvitskaMassExponent</name>
      <defaultValue>1.0d0</defaultValue>
      <source>parameters_</source>
      <description>The exponent of mass ratio appearing in the orbital angular momentum term in the Vitvitska model.</description>
    </inputParameter>
    <inputParameter>
      <name>angularMomentumVitvitskaEvolveDifferentially</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters_</source>
      <description>If true, the angular momentum vector of the halo evolves under differential evolution due to the accretion of subresolution mass. If false, the angular momentum vector is held fixed during differential evolution and only updates during node mergers. While the former should be more realistic, the latter is currently preferred as the standard hot halo component treatment of angular momentum relies on this assumption.</description>
    </inputParameter>
    !!]
    ! Bind deferred functions.
    call spin%angularMomentumFunction          (Vitvitska_Halo_Angular_Momentum            )
    call spin%angularMomentumVectorFunction    (Vitvitska_Halo_Angular_Momentum_Vector     )
    call spin%angularMomentumGrowthRateFunction(Vitvitska_Halo_Angular_Momentum_Growth_Rate)
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Vitvitska_Bindings

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Vitvitska_Thread_Init</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Vitvitska_Thread_Init(parameters_)
    !!{
    Initializes the tree node Vitvitsake spin module.
    !!}
    use :: Dark_Matter_Profile_Scales, only : darkMatterProfileScaleRadius
    use :: Events_Hooks              , only : nodePromotionEvent                 , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes          , only : defaultSpinComponent
    use :: Input_Parameters          , only : inputParameter                     , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    if (defaultSpinComponent%vitvitskaIsActive()) then
       !![
       <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters_"/>
       <objectBuilder class="haloSpinDistribution"         name="haloSpinDistribution_"         source="parameters_"/>
       <objectBuilder class="darkMatterProfileDMO"         name="darkMatterProfileDMO_"         source="parameters_"/>
       <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters_"/>
       <objectBuilder class="virialOrbit"                  name="virialOrbit_"                  source="parameters_"/>
       <objectBuilder class="mergerTreeMassResolution"     name="mergerTreeMassResolution_"     source="parameters_"/>
       !!]
       call nodePromotionEvent%attach(defaultSpinComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentSpinVitvitska')
    end if
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Vitvitska_Thread_Init

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Vitvitska_Thread_Uninit</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Vitvitska_Thread_Uninit()
    !!{
    Uninitializes the Vitvitska angular momentum module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultSpinComponent
    implicit none

    if (defaultSpinComponent%vitvitskaIsActive()) then
       !![
       <objectDestructor name="darkMatterProfileScaleRadius_"/>
       <objectDestructor name="haloSpinDistribution_"        />
       <objectDestructor name="darkMatterProfileDMO_"        />
       <objectDestructor name="darkMatterHaloScale_"         />
       <objectDestructor name="virialOrbit_"                 />
       <objectDestructor name="mergerTreeMassResolution_"    />
       !!]
       if (nodePromotionEvent%isAttached(defaultSpinComponent,nodePromotion)) call nodePromotionEvent%detach(defaultSpinComponent,nodePromotion)
    end if
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Vitvitska_Thread_Uninit

  !![
  <mergerTreeInitializeTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Vitvitska_Initialize</unitName>
   <sortName>spin</sortName>
  </mergerTreeInitializeTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Vitvitska_Initialize(node)
    !!{
    Initialize the spin of {\normalfont \ttfamily node}.
    !!}
    use :: Dark_Matter_Halo_Spins  , only : Dark_Matter_Halo_Angular_Momentum_Scale
    use :: Galacticus_Nodes        , only : defaultSpinComponent                   , nodeComponentBasic                 , nodeComponentDarkMatterProfile, nodeComponentSatellite, &
          &                                 nodeComponentSpin                      , nodeComponentSpinVitvitska         , treeNode
    use :: Beta_Functions          , only : Beta_Function                          , Beta_Function_Incomplete_Normalized
    use :: Numerical_Constants_Math, only : Pi
    use :: Vectors                 , only : Vector_Magnitude
    implicit none
    type            (treeNode                      ), intent(inout), pointer   :: node
    type            (treeNode                      )               , pointer   :: nodeChild                         , nodeSibling                       , &
         &                                                                        nodeUnresolved
    class           (nodeComponentBasic            )               , pointer   :: basicChild                        , basicSibling                      , &
         &                                                                        basic                             , basicUnresolved
    class           (nodeComponentSpin             )               , pointer   :: spin                              , spinSibling                       , &
         &                                                                        spinChild
    class           (nodeComponentSatellite        )               , pointer   :: satelliteSibling
    class           (nodeComponentDarkMatterProfile)               , pointer   :: darkMatterProfileUnresolved
    double precision                                , dimension(3)             :: angularMomentumOrbital            , angularMomentumTotal              , &
         &                                                                        angularMomentumUnresolved
    double precision                                               , parameter :: massFunctionSlopeLogarithmic=1.9d0
    double precision                                                           :: angularMomentumValue              , massRatio                         , &
         &                                                                        theta                             , phi                               , &
         &                                                                        massUnresolved                    , radiusScaleUnresolved             , &
         &                                                                        massResolution                    , angularMomentumSubresolutionFactor, &
         &                                                                        a                                 , b

    ! Check if we are the default method.
    if (defaultSpinComponent%vitvitskaIsActive()) then
       ! Get the spin component.
       spin => node%spin(autoCreate=.true.)
       ! Ensure that the spin has not yet been assigned for this node.
       select type (spin)
       class is (nodeComponentSpinVitvitska)
          if (spin%angularMomentumValue() == 0.0d0) then
             basic => node%basic()
             ! If this node has no children, draw its spin from a distribution, and assign a direction which is isotropically
             ! distributed.
             if (.not.associated(node%firstChild)) then
                theta               =acos(2.0d0   *node%hostTree%randomNumberGenerator_%uniformSample()-1.0d0)
                phi                 =     2.0d0*Pi*node%hostTree%randomNumberGenerator_%uniformSample()
                angularMomentumValue=haloSpinDistribution_%sample(node)*Dark_Matter_Halo_Angular_Momentum_Scale(node,darkMatterProfileDMO_)
                angularMomentumTotal=angularMomentumValue*[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
             else
                nodeChild                =>  node      %firstChild
                spinChild =>  nodeChild %spin()
                basicChild               =>  nodeChild %basic              ()
                massUnresolved           =  +basic     %mass               () &
                     &                      -basicChild%mass               ()
                ! Get the resolution of the tree.
                massResolution=mergerTreeMassResolution_%resolution(node%hostTree)
                ! Node has multiple progenitors - iterate over them and sum their angular momenta.
                angularMomentumTotal =   0.0d0
                nodeSibling          =>  nodeChild
                do while(associated(nodeSibling%sibling))
                   nodeSibling                =>  nodeSibling   %sibling
                   basicSibling               =>  nodeSibling   %basic              (                      )
                   spinSibling                =>  nodeSibling   %spin               (                      )
                   satelliteSibling           =>  nodeSibling   %satellite          (autoCreate=.true.     )
                   massRatio                  =  +basicSibling  %mass               (                      ) &
                        &                        /basicChild    %mass               (                      )
                   angularMomentumOrbital     =  +Orbital_Angular_Momentum          (           nodeSibling)
                   massUnresolved             =  +massUnresolved                                             &
                        &                        -basicSibling  %mass               (                      )
                   ! Add orbital angular momentum of this sibling scaled by the reduced mass to correct to the center of mass of the
                   ! sibling-child binary system.
                   angularMomentumTotal=+angularMomentumTotal                    &
                        &               +angularMomentumOrbital                  &
                        &               /(                                       &
                        &                 +1.0d0                                 &
                        &                 +massRatio                             &
                        &                )**angularMomentumVitvitskaMassExponent
                   ! Add the spin angular momentum of the sibling.
                   angularMomentumTotal=+angularMomentumTotal                    &
                        &               +spinSibling%angularMomentumVector()
                end do
                ! Add in the spin angular momentum of the primary child.
                angularMomentumTotal=+angularMomentumTotal                       &
                     &               +spinChild%angularMomentumVector()
                ! Account for unresolved accretion. The assumption is that unresolved accretion has the mean specific angular momentum averaged over the distribution of virial orbits.
                nodeUnresolved              => treeNode                                       (                         )
                basicUnresolved             => nodeUnresolved               %basic            (autoCreate=.true.        )
                darkMatterProfileUnresolved => nodeUnresolved               %darkMatterProfile(autoCreate=.true.        )
                call basicUnresolved            %massSet            (massResolution       )
                call basicUnresolved            %timeSet            (basicChild%time()    )
                call basicUnresolved            %timeLastIsolatedSet(basicChild%time()    )
                radiusScaleUnresolved       =  darkMatterProfileScaleRadius_%radius           (           nodeUnresolved)
                call darkMatterProfileUnresolved%scaleSet           (radiusScaleUnresolved)
                ! Compute a correction factor to the orbital angular momentum which takes into account the mass dependence of the
                ! 1/(1+m/M)ᵅ term that is applied to the angular momentum, and the reduced mass factor that appears in the orbital
                ! angular momentum. Averaging this over a power-law mass function gives the result below. In the case that α=0 the
                ! result is identically 1 - in this case we avoid computing beta functions.
                massRatio=+basicUnresolved%mass()       &
                     &    /basicChild     %mass()
                a        =+2.0d0                        &
                     &    -massFunctionSlopeLogarithmic
                if (angularMomentumVitvitskaMassExponent == 0.0d0) then
                   angularMomentumSubresolutionFactor=+1.0d0
                else
                   b                                 =+massFunctionSlopeLogarithmic                                         &
                        &                             +angularMomentumVitvitskaMassExponent                                 &
                        &                             -2.0d0
                   angularMomentumSubresolutionFactor=+Beta_Function_Incomplete_Normalized(a,b,massRatio/(1.0d0+massRatio)) &
                        &                             *Beta_Function                      (a,b                            ) &
                        &                             *           (2.0d0-massFunctionSlopeLogarithmic)                      &
                        &                             /massRatio**(2.0d0-massFunctionSlopeLogarithmic)
                end if
                ! Accumulate the angular momentum of the unresolved mass. Set the angular momentum accretion rate of the child
                ! node such that it accretes the correct amount over its lifetime.
                angularMomentumUnresolved=+massUnresolved                                                                         &
                     &                    *virialOrbit_                      %angularMomentumVectorMean(nodeUnresolved,nodeChild) &
                     &                    *angularMomentumSubresolutionFactor      
                angularMomentumTotal     =+angularMomentumTotal                                                                   &
                     &                    +angularMomentumUnresolved
                call spinChild%angularMomentumAccretionRateSet(                           &
                     &                                         +angularMomentumUnresolved &
                     &                                         /(                         &
                     &                                           +basic     %time()       &
                     &                                           -basicChild%time()       &
                     &                                          )                         &
                     &                                        )
                call nodeUnresolved%destroy()
                deallocate(nodeUnresolved)
                ! Compute the magnitude of the angular momentum.
                angularMomentumValue=Vector_Magnitude(angularMomentumTotal)
             end if
             call spin%angularMomentumSet      (angularMomentumValue)
             call spin%angularMomentumVectorSet(angularMomentumTotal)
          end if
       end select
    end if
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Vitvitska_Initialize

  subroutine Node_Component_Halo_Angular_Momentum_Vitvitska_Branch_Init(self)
    !!{
    Ensure a branch of a merger tree is initialized with spins in the Vitvitska model.
    !!}
    use :: Galacticus_Nodes   , only : nodeComponentSpinVitvitska, treeNode
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodesBranch
    implicit none
    class(nodeComponentSpinVitvitska          ), intent(inout) :: self
    type (treeNode                            ), pointer       :: node
    type (mergerTreeWalkerIsolatedNodesBranch )                :: treeWalker

    treeWalker=mergerTreeWalkerIsolatedNodesBranch(self%hostNode)
    do while (treeWalker%next(node))
       call Node_Component_Halo_Angular_Momentum_Vitvitska_Initialize(node)
    end do
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Vitvitska_Branch_Init

  double precision function Vitvitska_Halo_Angular_Momentum(self)
    !!{
    Return the spin parameter.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpinVitvitska
    implicit none
    class(nodeComponentSpinVitvitska), intent(inout) :: self

    if (self%angularMomentumValue() == 0.0d0) call Node_Component_Halo_Angular_Momentum_Vitvitska_Branch_Init(self)
    Vitvitska_Halo_Angular_Momentum=self%angularMomentumValue()
    return
  end function Vitvitska_Halo_Angular_Momentum

  function Vitvitska_Halo_Angular_Momentum_Vector(self)
    !!{
    Return the spin parameter vector.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpinVitvitska
    implicit none
    double precision                            , dimension(:) , allocatable :: Vitvitska_Halo_Angular_Momentum_Vector
    class           (nodeComponentSpinVitvitska), intent(inout)              :: self

    if (all(self%angularMomentumVectorValue() == 0.0d0)) call Node_Component_Halo_Angular_Momentum_Vitvitska_Branch_Init(self)
    Vitvitska_Halo_Angular_Momentum_Vector=self%angularMomentumVectorValue()
    return
  end function Vitvitska_Halo_Angular_Momentum_Vector

  double precision function Vitvitska_Halo_Angular_Momentum_Growth_Rate(self)
    !!{
    Return the growth rate of the spin parameter.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSpinVitvitska
    implicit none
    class(nodeComponentSpinVitvitska), intent(inout) :: self
    class(nodeComponentBasic        ), pointer       :: basic

    if (angularMomentumVitvitskaEvolveDifferentially) then
       ! The vector angular momentum accretion rate is available for each node. For the rate of change of the spin
       ! magnitude what we care about if the projection of this vector rate of change onto the instantaneous spin
       ! unit vector.
       basic                                       =>             self%hostNode%basic                ()
       Vitvitska_Halo_Angular_Momentum_Growth_Rate = +Dot_Product(self         %angularMomentumVector(),self%angularMomentumAccretionRate()) &
            &                                        /            self         %angularMomentum      () 
    else
       ! No differential evolution is applied.
       Vitvitska_Halo_Angular_Momentum_Growth_Rate = +0.0d0
    end if
    return
  end function Vitvitska_Halo_Angular_Momentum_Growth_Rate

  subroutine nodePromotion(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this
    case, we simply update the spin of {\normalfont \ttfamily node} to be consistent with the
    merging event.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic     , nodeComponentSpin, nodeComponentSpinVitvitska, treeNode
    implicit none
    class(*                 ), intent(inout)          :: self
    type (treeNode          ), intent(inout), target  :: node
    type (treeNode          )               , pointer :: nodeParent
    class(nodeComponentSpin )               , pointer :: spinParent , spin
    class(nodeComponentBasic)               , pointer :: basicParent, basic
    !$GLC attributes unused :: self

    spin        => node      %spin  ()
    nodeParent  => node      %parent
    basic       => node      %basic ()
    basicParent => nodeParent%basic ()
    spinParent  => nodeParent%spin  ()
    if (basic%time() /= basicParent%time()) &
         & call Galacticus_Error_Report('node has not been evolved to its parent'//{introspection:location})
    ! Adjust the spin to that of the parent node.
    call spin%angularMomentumSet             (spinParent%angularMomentum             ())
    call spin%angularMomentumVectorSet       (spinParent%angularMomentumVector       ())
    call spin%angularMomentumAccretionRateSet(spinParent%angularMomentumAccretionRate())
    return
  end subroutine nodePromotion

  !![
  <rateComputeTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Vitvitska_Rate_Compute</unitName>
  </rateComputeTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Vitvitska_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !!{
    Compute rates of change of properties in the Vitvitska implementation of the spin component.
    !!}
    use :: Galacticus_Nodes, only : defaultSpinComponent, nodeComponentSpin, nodeComponentSpinVitvitska, propertyTypeInactive, &
          &                         treeNode
    implicit none
    type            (treeNode         ), intent(inout)          :: node
    logical                            , intent(inout)          :: interrupt
    procedure       (                 ), intent(inout), pointer :: interruptProcedure
    integer                            , intent(in   )          :: propertyType
    class           (nodeComponentSpin)               , pointer :: spin
    double precision                                            :: angularMomentumMagnitude, angularMomentumGrowthRate
    double precision                   , dimension(3)           :: angularMomentumVector
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Return immediately if this class is not in use.
    if (.not.defaultSpinComponent%vitvitskaIsActive()) return
    ! Get the spin component.
    spin => node%spin()
    ! Ensure that it is of the Vitvitska class.
    select type (spin)
    class is (nodeComponentSpinVitvitska)
       ! Rate of change of the spin magnitude is determined from the virtual property. For the rate of change of the spin vector
       ! we split the results into two. The first part is the change in the magnitude of the spin vector, which is simply
       ! proportional to the rate of change of the scalar spin. The second is the perpendicular component of the rate of change
       ! (which therefore does not change the magnitude of the spin, and so does not depend on changes in the mass or energy of
       ! the halo). This second term is computed by taking the vector rate of angular momentum accretion and subtracting off that
       ! component along the instantaneous spin unit vector.
       angularMomentumGrowthRate=spin%angularMomentumGrowthRate()
       if (angularMomentumGrowthRate /= 0.0d0) then
          angularMomentumMagnitude =spin%angularMomentum      ()
          angularMomentumVector    =spin%angularMomentumVector()
          call spin%angularMomentumRate      (                                                                                          &
               &                              +angularMomentumGrowthRate                                                                &
               &                             )
          call spin%angularMomentumVectorRate(                                                                                          &
               &                              +  angularMomentumGrowthRate                                                              &
               &                              *  angularMomentumVector                                                                  &
               &                              /  angularMomentumMagnitude                                                               &
               &                              + (                                                                                       &
               &                                +                            spin%angularMomentumAccretionRate()                        &
               &                                -Dot_Product                (spin%angularMomentumAccretionRate(),angularMomentumVector) &
               &                                *angularMomentumVector                                                                  &
               &                                /angularMomentumMagnitude**2                                                            &
               &                              )                                                                                         &
               &                             )
       end if
    end select
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Vitvitska_Rate_Compute

  !![
  <scaleSetTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Vitvitska_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Vitvitska_Scale_Set(node)
    !!{
    Set scales for properties in the Vitvitska implementation of the spin component.
    !!}
    use :: Galacticus_Nodes      , only : nodeComponentSpin                      , nodeComponentSpinVitvitska, treeNode, defaultSpinComponent
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum_Scale
    implicit none
    type            (treeNode         ), intent(inout), pointer :: node
    double precision                   , parameter              :: spinScaleAbsolute=1.0d-4
    class           (nodeComponentSpin)               , pointer :: spin

    ! Return immediately if this class is not in use.
    if (.not.defaultSpinComponent%vitvitskaIsActive()) return
    ! Get the spin component.
    spin => node%spin()
    ! Ensure that it is of the Vitvitska class.
    select type (spin)
    class is (nodeComponentSpinVitvitska)
       ! Set scale for spin.
       call spin%angularMomentumScale(spinScaleAbsolute*Dark_Matter_Halo_Angular_Momentum_Scale(node,darkMatterProfileDMO_))
    end select
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Vitvitska_Scale_Set

  function Orbital_Angular_Momentum(node)
    !!{
    Returns the orbital angular momentum vector associated with a satellite by drawing a
    random position towards the host at virial radius distance and a random velocity vector
    consistent with the orbital parameters of the satellite.
    !!}
    use :: Coordinates     , only : assignment(=)     , coordinateCartesian
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSatellite, treeNode
    use :: Kepler_Orbits   , only : keplerOrbit
    use :: Vectors         , only : Vector_Product
    implicit none
    double precision                        , dimension(3)                :: Orbital_Angular_Momentum
    type            (treeNode              ), pointer     , intent(inout) :: node
    class           (nodeComponentSatellite), pointer                     :: satellite
    class           (nodeComponentBasic    ), pointer                     :: basic
    double precision                        , dimension(3)                :: haloVelocity            , haloPosition
    type            (keplerOrbit           )                              :: orbit
    type            (coordinateCartesian   )                              :: coordinates

    ! Get the orbital properties.
    basic        => node       %basic      (                 )
    satellite    => node       %satellite  (autoCreate=.true.)
    orbit        =  satellite  %virialOrbit(                 )
    coordinates  =  orbit      %position   (                 )
    haloPosition =  coordinates
    coordinates  =  orbit      %velocity   (                 )
    haloVelocity =  coordinates
    ! Calculate the orbital angular momentum vector.
    Orbital_Angular_Momentum=+               basic%mass()  &
         &                   *Vector_Product(              &
         &                                   haloPosition, &
         &                                   haloVelocity  &
         &                                  )
    return
  end function Orbital_Angular_Momentum

end module Node_Component_Halo_Angular_Momentum_Vitvitska
