!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!+ Contributions to this file made by: Chrsitoph Behrens.

!% Contains a module implementing a node spin component using the approach of
!% \cite{vitvitska_origin_2002}.

module Node_Component_Spin_Vitvitska
  !% Implements a node spin component using the approach of \cite{vitvitska_origin_2002}.
  use Galacticus_Nodes
  use Kepler_Orbits
  use Vectors
  use Dark_Matter_Halo_Spins
  use Numerical_Constants_Physical
  use Dark_Matter_Profiles
  use ISO_Varying_String
  implicit none
  private
  public :: Node_Component_Spin_Vitvitska_Promote     , Node_Component_Spin_Vitvitska_Initialize_Spins, &
       &    Node_Component_Spin_Vitvitska_Bindings    , Node_Component_Spin_Vitvitska_Scale_Set       , &
       &    Node_Component_Spin_Vitvitska_Rate_Compute

  !# <component>
  !#  <class>spin</class>
  !#  <name>vitvitska</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>spin</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="get" />
  !#     <output unitsInSI="0.0d0" comment="Spin parameter of the node."/>
  !#   </property>
  !#   <property>
  !#     <name>spinVector</name>
  !#     <type>double</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="get" />
  !#     <output labels="[X,Y,Z]" unitsInSI="0.0d0" comment="Spin vector of the node."/>
  !#     <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>spinGrowthRate</name>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" isDeferred="get" />
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#   </property>
  !#  </properties>
  !# </component>

  ! Module initialization state.
  logical :: moduleInitialized=.false.
  
contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Spin_Vitvitska_Bindings</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Spin_Vitvitska_Bindings()
    !% Initializes the ``Vitvitskae'' implementation of the spin component.
    implicit none
    type(nodeComponentspinVitvitska) :: spin

    ! Initialize the bindings.
    if (.not.moduleInitialized) then
       !$omp critical (Node_Component_Spin_Vitvitska_Bindings)
       if (.not.moduleInitialized) then
          ! Bind deferred functions.
          call spin%spinFunction          (Node_Component_Spin_Vitvitska_Spin            )
          call spin%spinVectorFunction    (Node_Component_Spin_Vitvitska_Spin_Vector     )
          call spin%spinGrowthRateFunction(Node_Component_Spin_Vitvitska_Spin_Growth_Rate)
          ! Record that the module is now initialize.
          moduleInitialized=.true.
       end if
       !$omp end critical (Node_Component_Spin_Vitvitska_Bindings)
    end if
    return
  end subroutine Node_Component_Spin_Vitvitska_Bindings
  
  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Spin_Vitvitska_Initialize_Spins</unitName>
  !#  <sortName>spin</sortName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Spin_Vitvitska_Initialize_Spins(node)
    !% Initialize the spin of {\normalfont \ttfamily node}.
    use Halo_Spin_Distributions
    use Dark_Matter_Halo_Spins
    use Dark_Matter_Profiles
    use ISO_Varying_String
    implicit none
    type            (treeNode                 ), intent(inout), pointer :: node
    type            (treeNode                 )               , pointer :: nodeChild             , nodeSibling
    class           (nodeComponentBasic       )               , pointer :: basicChild            , basicSibling        , &
         &                                                                 basic
    class           (nodeComponentSpin        )               , pointer :: spin                  , spinSibling         , &
         &                                                                 spinChild
    class           (nodeComponentSatellite   )               , pointer :: satelliteSibling
    class           (haloSpinDistributionClass)               , pointer :: haloSpinDistribution_
    class           (darkMatterProfileClass   )               , pointer :: darkMatterProfile_
    double precision                           , dimension(3)           :: angularMomentumOrbital, angularMomentumTotal, &
         &                                                                 spinVector
    double precision                                                    :: spinValue             , massRatio           , &
         &                                                                 theta                 , phi
    
    ! Check if we are the default method.
    if (defaultSpinComponent%vitvitskaIsActive()) then
       ! Get required objects.
       darkMatterProfile_    => darkMatterProfile   ()
       haloSpinDistribution_ => haloSpinDistribution()
       ! Get the spin component.
       spin => node%spin(autoCreate=.true.)
       ! Ensure that the spin has not yet been assigned for this node.
       select type (spin)
          class is (nodeComponentSpinVitvitska)
          if (spin%spinValue() == 0.0d0) then
             basic => node%basic()
             ! If this node has no children, draw its spin from a distribution, and assign a direction which is isotropically
             ! distributed.
             if (.not.associated(node%firstChild)) then
                theta     =acos(2.0d0   *node%hostTree%randomNumberGenerator%sample()-1.0d0)
                phi       =     2.0d0*Pi*node%hostTree%randomNumberGenerator%sample()
                spinValue =haloSpinDistribution_%sample(node)             
                spinVector=spinValue*[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
             else
                nodeChild  => node     %firstChild
                spinChild  => nodeChild%spin      ()
                basicChild => nodeChild%basic     ()
                ! Node has multiple progenitors - iterate over them and sum their angular momenta.
                nodeSibling          =>                                    nodeChild
                angularMomentumTotal =  +Dark_Matter_Halo_Angular_Momentum(nodeChild) &
                     &                  *spinChild%spinVector()                       &
                     &                  /spinChild%spin      ()
                do while(associated(nodeSibling%sibling))
                   nodeSibling            =>  nodeSibling %sibling
                   basicSibling           =>  nodeSibling %basic    (                 )
                   spinSibling            =>  nodeSibling %spin     (                 )
                   satelliteSibling       =>  nodeSibling %satellite(autoCreate=.true.)
                   massRatio              =  +basicSibling%mass     (                 ) &
                        &                    /basicChild  %mass     (                 )
                   angularMomentumOrbital =  +Orbital_Angular_Momentum(nodeSibling)
                   ! Add orbital angular momentum of this sibling scaled by the reduced mass to correct to the center of mass of the
                   ! sibling-child binary system.
                   angularMomentumTotal=+angularMomentumTotal   &
                        &               +angularMomentumOrbital &
                        &               /(                      &
                        &                 +1.0d0                &
                        &                 +massRatio            &
                        &               )
                   ! Add the spin angular momentum of the sibling.
                   angularMomentumTotal=+angularMomentumTotal                           &
                        &               +Dark_Matter_Halo_Angular_Momentum(nodeSibling) &
                        &               *spinSibling%spinVector()                       &
                        &               /spinSibling%spin      ()
                end do
                ! Convert angular momentum back to spin.
                spinValue =+Dark_Matter_Halo_Spin(node,Vector_Magnitude(angularMomentumTotal))
                spinVector=+spinValue                              &
                     &     *                 angularMomentumTotal  &
                     &     /Vector_Magnitude(angularMomentumTotal)
             end if
             call spin%spinSet      (spinValue )
             call spin%spinVectorSet(spinVector)
          end if
       end select
    end if
    return
  end subroutine Node_Component_Spin_Vitvitska_Initialize_Spins

  subroutine Node_Component_Spin_Vitvitska_Branch_Initialize(self)
    !% Ensure a branch of a merger tree is initialized with spins in the Vitvitska model.
    implicit none
    class  (nodeComponentSpinVitvitska), intent(inout) :: self
    type   (treeNode                  ), pointer       :: node
    logical                                            :: finished

    node     => self%hostNode
    finished =  .false.
    do while (.not.finished)
       node => node%walkBranch(self%hostNode)
       if (.not.associated(node)) then
          ! When a null pointer is returned, the full branch has been walked. We then have to
          ! handle the base node of the branch as a special case.
          node     => self%hostNode
          finished =  .true.
       end if
       call Node_Component_Spin_Vitvitska_Initialize_Spins(node)
    end do
    return
  end subroutine Node_Component_Spin_Vitvitska_Branch_Initialize
  
  double precision function Node_Component_Spin_Vitvitska_Spin(self)
    !% Return the spin parameter.
    implicit none
    class(nodeComponentSpinVitvitska), intent(inout) :: self
    
    if (self%spinValue() == 0.0d0) call Node_Component_Spin_Vitvitska_Branch_Initialize(self)
    Node_Component_Spin_Vitvitska_Spin=self%spinValue()
    return
  end function Node_Component_Spin_Vitvitska_Spin

  function Node_Component_Spin_Vitvitska_Spin_Vector(self)
    !% Return the spin parameter vector.
    implicit none
    double precision                            , dimension(:) , allocatable :: Node_Component_Spin_Vitvitska_Spin_Vector
    class           (nodeComponentSpinVitvitska), intent(inout)              :: self
    
    if (all(self%spinVectorValue() == 0.0d0)) call Node_Component_Spin_Vitvitska_Branch_Initialize(self)
    Node_Component_Spin_Vitvitska_Spin_Vector=self%spinVectorValue()
    return
  end function Node_Component_Spin_Vitvitska_Spin_Vector

  double precision function Node_Component_Spin_Vitvitska_Spin_Growth_Rate(self)
    !% Return the growth rate of the spin parameter.
    use Dark_Matter_Profiles
    implicit none
    class(nodeComponentSpinVitvitska), intent(inout) :: self
    class(nodeComponentBasic        ), pointer       :: basic
    class(darkMatterProfileClass    ), pointer       :: darkMatterProfile_

    ! Assumes that the angular momentum of the halo remains unchanged during differential evolution, so changes in spin arise only
    ! from changes in the mass and energy of the halo.
    basic                                          =>  self             %hostNode%basic     (             )
    darkMatterProfile_                             =>  darkMatterProfile                    (             )
    Node_Component_Spin_Vitvitska_Spin_Growth_Rate =  +self             %spin               (             ) &
         &                                            *(                                                    &
         &                                              +0.5d0                                              &
         &                                              *darkMatterProfile_%energyGrowthRate(self%hostNode) &
         &                                              /darkMatterProfile_%energy          (self%hostNode) &
         &                                              -2.5d0                                              &
         &                                              *basic             %accretionRate   (             ) &
         &                                              /basic             %mass            (             ) &
         &                                             )
    return
  end function Node_Component_Spin_Vitvitska_Spin_Growth_Rate

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Spin_Vitvitska_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Spin_Vitvitska_Promote(node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this
    !% case, we simply update the spin of {\normalfont \ttfamily node} to be consistent with the
    !% merging event.
    use Galacticus_Error
    implicit none
    type (treeNode          ), intent(inout), pointer :: node
    type (treeNode          )               , pointer :: nodeParent
    class(nodeComponentSpin )               , pointer :: spinParent , spin
    class(nodeComponentBasic)               , pointer :: basicParent, basic

    ! Ensure that the spin component is of the vitvitska class.
    spin => node%spin()
    select type (spin)
    class is (nodeComponentSpinVitvitska)
       nodeParent  => node      %parent
       basic       => node      %basic ()
       basicParent => nodeParent%basic ()
       spinParent  => nodeParent%spin  ()
       if (basic%time() /= basicParent%time()) &
            & call Galacticus_Error_Report('node has not been evolved to its parent'//{introspection:location})   
       ! Adjust the spin to that of the parent node.
       call spin%spinSet      (spinParent%spin      ())
       call spin%spinVectorSet(spinParent%spinVector())
    end select
    return
  end subroutine Node_Component_Spin_Vitvitska_Promote
   
  !# <rateComputeTask>
  !#  <unitName>Node_Component_Spin_Vitvitska_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Spin_Vitvitska_Rate_Compute(node,odeConverged,interrupt,interruptProcedure,propertyType)
    !% Compute rates of change of properties in the Vitvitska implementation of the spin component.
    implicit none
    type            (treeNode         ), intent(inout), pointer :: node
    logical                            , intent(in   )          :: odeConverged
    logical                            , intent(inout)          :: interrupt
    procedure       (                 ), intent(inout), pointer :: interruptProcedure
    integer                            , intent(in   )          :: propertyType
    class           (nodeComponentSpin)               , pointer :: spin
    double precision                                            :: spinMagnitude     , spinGrowthRate
    double precision                   , dimension(3)           :: spinVector
    !GCC$ attributes unused :: interrupt, interruptProcedure, odeConverged
    
    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Return immediately if this class is not in use.
    if (.not.defaultSpinComponent%vitvitskaIsActive()) return
    ! Get the spin component.
    spin => node%spin()
    ! Ensure that it is of the Vitvitska class.
    select type (spin)
    class is (nodeComponentSpinVitvitska)
       ! Rate of change is set equal to the precomputed growth rate.
       spinGrowthRate=spin%spinGrowthRate()
       if (spinGrowthRate /= 0.0d0) then
          spinMagnitude =spin%spin          ()
          spinVector    =spin%spinVector    ()
          call spin%spinRate      (spinGrowthRate                         )
          call spin%spinVectorRate(spinGrowthRate*spinVector/spinMagnitude)
       end if
    end select
    return
  end subroutine Node_Component_Spin_Vitvitska_Rate_Compute

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Spin_Vitvitska_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Spin_Vitvitska_Scale_Set(node)
    !% Set scales for properties in the Vitvitska implementation of the spin component.
    implicit none
    type            (treeNode         ), intent(inout), pointer :: node
    double precision                   , parameter              :: spinScaleAbsolute=1.0d-4
    class           (nodeComponentSpin)               , pointer :: spin

    ! Get the spin component.
    spin => node%spin()
    ! Ensure that it is of the Vitvitska class.
    select type (spin)
    class is (nodeComponentSpinVitvitska)
       ! Set scale for spin.
       call spin%spinScale(spinScaleAbsolute)
    end select
    return
  end subroutine Node_Component_Spin_Vitvitska_Scale_Set

  function Orbital_Angular_Momentum(node)
    !% Returns the orbital angular momentum vector associated with a satellite by drawing a
    !% random position towards the host at virial radius distance and a random velocity vector
    !% consistent with the orbital parameters of the satellite.
    use Virial_Orbits
    use Satellite_Merging_Timescales
    use Vectors
    implicit none
    double precision                        , dimension(3)                :: Orbital_Angular_Momentum
    type            (treeNode              ), pointer     , intent(inout) :: node
    class           (nodeComponentSatellite), pointer                     :: satellite
    class           (nodeComponentBasic    ), pointer                     :: basic
    double precision                        , dimension(3)                :: haloVelocity            , haloPosition     , &
         &                                                                   vectorRadial            , vectorTangential1, &
         &                                                                   vectorTangential2       , p
    double precision                        , dimension(2)                :: q
    type            (keplerOrbit           )                              :: orbit
    double precision                                                      :: velocityRadial          , radius           , &
         &                                                                   velocityTangential    
    integer                                                               :: i

    ! Get the orbital properties.
    basic              => node     %basic             (                 )
    satellite          => node     %satellite         (autoCreate=.true.)
    orbit              =  satellite%virialOrbit       (                 )
    velocityRadial     =  orbit    %velocityRadial    (                 )
    velocityTangential =  orbit    %velocityTangential(                 )
    radius             =  orbit    %radius            (                 )
    ! Find the position vector of the subhalo on the virial sphere.
    p=1.0d0
    do while(Vector_Magnitude(p) > 1.0d0)
       do i=1,3
          p(i)=node%hostTree%randomNumberGenerator%sample()*2.0d0-1.0d0
       end do
    end do
    vectorRadial=+                 p  &
         &       /Vector_Magnitude(p)
    haloPosition=+vectorRadial        &
         &       *radius
    ! Find two tangential normal vectors. The first we construct by hand.
    if (vectorRadial(1) == 0.0d0 .and. vectorRadial(2) == 0.0d0) then
       vectorTangential1   =[+0.0d0                                ,+1.0d0                                ,+0.0d0]
    else 
       if (.not.vectorRadial(1) == 0.0d0) then
          vectorTangential1=[-1.0d0*vectorRadial(2)/vectorRadial(1),+1.0d0                                ,+0.0d0]
       else
          vectorTangential1=[+1.0d0                                ,-1.0d0*vectorRadial(1)/vectorRadial(2),+0.0d0]
       end if
    end if
    vectorTangential1=+                 vectorTangential1  &
         &            /Vector_Magnitude(vectorTangential1)
   ! The second vector is obtained from the cross product.
   vectorTangential2=Vector_Product(vectorRadial,vectorTangential1)
   ! Choose a random direction for the tangential velocity.
   q=1.0d0
   do while(sqrt(sum(q**2)) > 1.0d0)
      do i=1,2
         q(i)=node%hostTree%randomNumberGenerator%sample()*2.0d0-1.0d0
      end do
   end do
   q=q/sqrt(sum(q**2))
   ! Construct the full velocity vector.
   haloVelocity=+       velocityRadial     &
        &       *       vectorRadial       &
        &       +       velocityTangential &
        &       *(                         &
        &         +q(1)*vectorTangential1  &
        &         +q(2)*vectorTangential2  &
        &       )
   ! Calculate the orbital angular momentum vector.
   Orbital_Angular_Momentum=+               basic%mass()  &
        &                   *Vector_Product(              &
        &                                   haloPosition, &
        &                                   haloVelocity  &
        &                                  )
   return   
 end function Orbital_Angular_Momentum

end module Node_Component_Spin_Vitvitska
