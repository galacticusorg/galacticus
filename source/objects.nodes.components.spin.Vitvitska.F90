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

!+ Contributions to this file made by: Chrsitoph Behrens.

!% Contains a module implementing a node spin component using the approach of
!% \cite{vitvitska_origin_2002}.

module Node_Component_Spin_Vitvitska
  !% Implements a node spin component using the approach of \cite{vitvitska_origin_2002}.
  use Galacticus_Nodes
  use FGSL
  use Kepler_Orbits
  use Vectors
  use Dark_Matter_Halo_Spins
  use Numerical_Constants_Physical
  use Dark_Matter_Profiles
  use ISO_Varying_String
  implicit none
  private
  public :: Node_Component_Spin_Vitvitska_Initialize, Node_Component_Spin_Vitvitska_Initialize_Spins, &
       &    Node_Component_Spin_Vitvitska_Promote

  !# <component>
  !#  <class>spin</class>
  !#  <name>vitvitska</name>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>spin</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="0.0d0" comment="Spin parameter of the node."/>
  !#   </property>
  !#   <property>
  !#     <name>spinGrowthRate</name>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" />
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#     <getFunction bindsTo="component">spinVitvitskaSpinGrowthRate</getFunction>
  !#   </property>
  !#  </properties>
  !#  <functions>objects.nodes.components.spin.Vitvitska.bound_functions.inc</functions>
  !# </component>

  ! Record of whether the module has been initialized.
  logical                    :: moduleInitialized              =.false.

  ! Random number sequence used for choosing satellite orbits.
  type            (fgsl_rng) :: pseudoSequenceObject
  logical                    :: reset                          =.true.
  !$omp threadprivate(pseudoSequenceObject,reset)

  ! Parameter controlling weight given to mergers.
  double precision           :: spinVitvitskaMergerRatioExponent

contains

  subroutine Node_Component_Spin_Vitvitska_Initialize()
    !% Initializes the ``vitvitska'' spin component module.
    use Input_Parameters
    implicit none
    
    ! Test whether module is already initialized.
    !$omp critical (Node_Component_Spin_Vitvitska_Initialize)
    if (.not.moduleInitialized) then    
       !@ <inputParameter>
       !@   <name>spinVitvitskaMergerRatioExponent</name>
       !@   <defaultValue>0.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The angular momentum of a merging halo is weighted by $(1+M_{\rm satellite}/M_{\rm host})^\alpha$ where $\alpha=${\normalfont \ttfamily [spinVitvitskaMaximumMergerRatio]}.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('spinVitvitskaMergerRatioExponent',spinVitvitskaMergerRatioExponent,defaultValue=2.0d0)
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Spin_Vitvitska_Initialize)
    return
  end subroutine Node_Component_Spin_Vitvitska_Initialize

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
    class           (nodeComponentBasic       )               , pointer :: basicChild            , basicSibling             , &
         &                                                                 basic
    class           (nodeComponentSpin        )               , pointer :: spin                  , spinNew                  , &
         &                                                                 spinChild
    class           (nodeComponentSatellite   )               , pointer :: satelliteSibling
    class           (haloSpinDistributionClass)               , pointer :: haloSpinDistribution_
    class           (darkMatterProfileClass   )               , pointer :: darkMatterProfile_
    double precision                           , dimension(3)           :: angularMomentumOrbital, angularMomentumTotal
    double precision                                                    :: spinValue             , massRatio
    
    ! Check if we are the default method.
    if (defaultSpinComponent%vitvitskaIsActive()) then
       ! Ensure the module is initialized.
       call Node_Component_Spin_Vitvitska_Initialize()
       ! Get required objects.
       darkMatterProfile_    => darkMatterProfile   ()
       haloSpinDistribution_ => haloSpinDistribution()
       ! Get the spin component.
       spin => node%spin()
       ! Ensure that the spin has not yet been assigned for this node.
       select type (spin)
       type is (nodeComponentSpin)
          basic => node%basic()
          ! If this node has no children, draw its spin from a distribution.
          spinNew => node%spin(autoCreate=.true.)
          if (.not.associated(node%firstChild)) then           
             spinValue=haloSpinDistribution_%sample(node)
          else
             nodeChild  => node     %firstChild
             spinChild  => nodeChild%spin      ()
             basicChild => nodeChild%basic     ()
             ! Node has multiple progenitors - iterate over them and sum their angular momenta.
             angularMomentumTotal =  [0.0d0,0.0d0,Dark_Matter_Halo_Angular_Momentum(nodeChild)]
             nodeSibling          =>                                                nodeChild
             do while(associated(nodeSibling%sibling))
                nodeSibling            =>  nodeSibling %sibling
                basicSibling           =>  nodeSibling %basic    (                 )
                satelliteSibling       =>  nodeSibling %satellite(autoCreate=.true.)
                massRatio              =  +basicSibling%mass     (                 ) &
                     &                    /basicChild  %mass     (                 )
                angularMomentumOrbital =  +Orbital_Angular_Momentum(nodeSibling)
                ! Sum angular momenta with a weight dependent on the mass ratio.
                angularMomentumTotal=+angularMomentumTotal               &
                     &               +angularMomentumOrbital             &
                     &               /(                                  &
                     &                 +1.0d0                            &
                     &                 +massRatio                        &
                     &               )**spinVitvitskaMergerRatioExponent
             end do
             ! Convert angular momentum back to spin.
             spinValue=Dark_Matter_Halo_Spin(node,Vector_Magnitude(angularMomentumTotal))
          end if
          call spinNew%spinSet(spinValue)
       end select
    end if
    return
  end subroutine Node_Component_Spin_Vitvitska_Initialize_Spins
  
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
            & call Galacticus_Error_Report('Node_Component_Spin_Vitvitska_Promote','node has not been evolved to its parent')   
       ! Adjust the spin to that of the parent node.
       call spin%spinSet(spinParent%spin())
    end select
    return
  end subroutine Node_Component_Spin_Vitvitska_Promote
   
  function Orbital_Angular_Momentum(node)
    !% Returns the orbital angular momentum vector associated with a satellite by drawing a
    !% random position towards the host at virial radius distance and a random velocity vector
    !% consistent with the orbital parameters of the satellite.
    use Virial_Orbits
    use Satellite_Merging_Timescales
    use Pseudo_Random
    use Vectors
    use FGSL
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
    basic              => node     %basic             ()
    satellite          => node     %satellite         ()
    orbit              =  satellite%virialOrbit       ()
    velocityRadial     =  orbit    %velocityRadial    ()
    velocityTangential =  orbit    %velocityTangential()
    radius             =  orbit    %radius            ()
    ! Find the position vector of the subhalo on the virial sphere.
    p=1.0d0
    do while(Vector_Magnitude(p) > 1.0d0)
       do i=1,3
          p(i)=Pseudo_Random_Get(pseudoSequenceObject,reset)*2.0d0-1.0d0
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
         q(i)=Pseudo_Random_Get(pseudoSequenceObject,reset)*2.0d0-1.0d0
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
