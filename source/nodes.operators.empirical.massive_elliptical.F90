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
  Implements a node operator class that inserts an empirical model of the formation history of a massive elliptical galaxy.
  !!}

  !![
  <nodeOperator name="nodeOperatorEmpiricalMassiveElliptical">
   <description>A node operator class that inserts an empirical model of the formation history of a massive elliptical galaxy.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorEmpiricalMassiveElliptical
     !!{     
     A node operator class that inserts an empirical model of the formation history of a massive elliptical galaxy. The galaxy is
     assumed to grow in the main branch of the tree with a constant specific star formation rate, such that it mass is given by:
     \begin{equation}
       M_\star(t) = M_{\star,0} \exp(-\phi_\star [t-t_0]),
     \end{equation}
     where $M_{\star,0}=${\normalfont \ttfamily [massStellarFinal]} is the stellar mass in the root node of the tree,
     $\phi_\star=${\normalfont \ttfamily [rateStarFormationSpecific]}, and $t_0$ is the cosmic time at the root node of the tree.
     !!}
     private
     double precision :: massStellarFinal                  , rateStarFormationSpecific                , &
          &              angularMomentumPseudoSpecificFinal, rateAngularMomentumPseudoSpecificSpecific, &
          &              radiusFinal                       , rateRadiusSpecific
     logical          :: useAngularMomentum
   contains
     procedure :: nodeInitialize                      => empiricalMassiveEllipticalNodeInitialize
     procedure :: differentialEvolution               => empiricalMassiveEllipticalDifferentialEvolution
     procedure :: differentialEvolutionSolveAnalytics => empiricalMassiveEllipticalSolveAnalytics
  end type nodeOperatorEmpiricalMassiveElliptical
  
  interface nodeOperatorEmpiricalMassiveElliptical
     !!{
     Constructors for the \refClass{nodeOperatorEmpiricalMassiveElliptical} node operator class.
     !!}
     module procedure empiricalMassiveEllipticalConstructorParameters
     module procedure empiricalMassiveEllipticalConstructorInternal
  end interface nodeOperatorEmpiricalMassiveElliptical
  
contains

  function empiricalMassiveEllipticalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorEmpiricalMassiveElliptical} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorEmpiricalMassiveElliptical)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    double precision                                                        :: massStellarFinal                  , rateStarFormationSpecific                , &
         &                                                                     angularMomentumPseudoSpecificFinal, rateAngularMomentumPseudoSpecificSpecific, &
         &                                                                     radiusFinal                       , rateRadiusSpecific
    
    !![
    <inputParameter>
      <name>massStellarFinal</name>
      <source>parameters</source>
      <description>The final stellar mass of the elliptical galaxy.</description>
    </inputParameter>
    <inputParameter>
      <name>rateStarFormationSpecific</name>
      <source>parameters</source>
      <description>The specific star formation rate of the elliptical galaxy.</description>
    </inputParameter>
    !!]
    if (parameters%isPresent('angularMomentumPseudoSpecificFinal')) then
       !![
       <inputParameter>
	 <name>angularMomentumPseudoSpecificFinal</name>
	 <source>parameters</source>
	 <description>The final specific pseudo-angular momentum of the elliptical galaxy.</description>
       </inputParameter>
       <inputParameter>
	 <name>rateAngularMomentumPseudoSpecificSpecific</name>
	 <source>parameters</source>
	 <description>The specific growth rate of the specific pseudo-angular momentum of the elliptical galaxy.</description>
       </inputParameter>
       !!]
    end if
    if (parameters%isPresent(                       'radiusFinal')) then
       !![
       <inputParameter>
	 <name>radiusFinal</name>
	 <source>parameters</source>
	 <description>The final radius of the elliptical galaxy.</description>
       </inputParameter>
       <inputParameter>
	 <name>rateRadiusSpecific</name>
	 <source>parameters</source>
	 <description>The specific growth rate of the radius of the elliptical galaxy.</description>
       </inputParameter>
       !!]
    end if
    !![
    <conditionalCall>
      <call>self=nodeOperatorEmpiricalMassiveElliptical(massStellarFinal,rateStarFormationSpecific{conditions})</call>
      <argument name="angularMomentumPseudoSpecificFinal"        value="angularMomentumPseudoSpecificFinal"        parameterPresent="parameters"/>
      <argument name="rateAngularMomentumPseudoSpecificSpecific" value="rateAngularMomentumPseudoSpecificSpecific" parameterPresent="parameters"/>
      <argument name="radiusFinal"                               value="radiusFinal"                               parameterPresent="parameters"/>
      <argument name="rateRadiusSpecific"                        value="rateRadiusSpecific"                        parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function empiricalMassiveEllipticalConstructorParameters

  function empiricalMassiveEllipticalConstructorInternal(massStellarFinal,rateStarFormationSpecific,angularMomentumPseudoSpecificFinal,rateAngularMomentumPseudoSpecificSpecific,radiusFinal,rateRadiusSpecific) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorEmpiricalMassiveElliptical} node operator class.
    !!}
    implicit none
    type            (nodeOperatorEmpiricalMassiveElliptical)                          :: self
    double precision                                        , intent(in   )           :: massStellarFinal                  , rateStarFormationSpecific
    double precision                                        , intent(in   ), optional :: angularMomentumPseudoSpecificFinal, rateAngularMomentumPseudoSpecificSpecific, &
         &                                                                               radiusFinal                       , rateRadiusSpecific
    !![
    <constructorAssign variables="massStellarFinal, rateStarFormationSpecific, angularMomentumPseudoSpecificFinal, rateAngularMomentumPseudoSpecificSpecific, radiusFinal, rateRadiusSpecific"/>
    !!]

    self%useAngularMomentum=.true.
    if      (present(angularMomentumPseudoSpecificFinal)) then
       if     (                                                                                                             &
            &   .not.present(rateAngularMomentumPseudoSpecificSpecific)                                                     &
            & ) call Error_Report('must provide angular momentum specific growth rate'          //{introspection:location})
       if     (                                                                                                             &
            &        present(                              radiusFinal)                                                     &
            &  .or.&
            &        present(                       rateRadiusSpecific)                                                     &
            & ) call Error_Report('either angular momentum or radius must be provided, not both'//{introspection:location})
       self%useAngularMomentum=.true.
    else if (present(                       radiusFinal)) then
       if     (                                                                                                             &
            &   .not.present(                       rateRadiusSpecific)                                                     &
            & ) call Error_Report('must provide radius specific growth rate'                    //{introspection:location})
       if     (                                                                                                             &
            &        present(       angularMomentumPseudoSpecificFinal)                                                     &
            &  .or.&
            &        present(rateAngularMomentumPseudoSpecificSpecific)                                                     &
            & ) call Error_Report('either angular momentum or radius must be provided, not both'//{introspection:location})
       self%useAngularMomentum=.false.
    else
       self%useAngularMomentum=.false.
       call Error_Report('either angular momentum or radius must be provided'//{introspection:location})
    end if
    return
  end function empiricalMassiveEllipticalConstructorInternal

  subroutine empiricalMassiveEllipticalNodeInitialize(self,node)
    !!{
    Initialize nodes for the massive elliptical.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSpheroid
    implicit none
    class           (nodeOperatorEmpiricalMassiveElliptical), intent(inout), target  :: self
    type            (treeNode                              ), intent(inout), target  :: node
    type            (treeNode                              )               , pointer :: nodeRoot
    class           (nodeComponentBasic                    )               , pointer :: basicLeaf  , basicRoot
    class           (nodeComponentSpheroid                 )               , pointer :: spheroid
    double precision                                                                 :: timeLeaf   , timeRoot       , &
         &                                                                              massStellar, angularMomentum

    ! Initialize only the leaf node on the main branch.
    if (associated(node%firstChild).or..not.node%isOnMainBranch()) return
    ! Get times in the leaf and root nodes.
    nodeRoot  =>  node     %hostTree   %nodeBase
    basicLeaf =>  node                 %basic   (                 )
    basicRoot =>  nodeRoot             %basic   (                 )
    spheroid  =>  node                 %spheroid(autoCreate=.true.)
    timeLeaf  =   basicLeaf            %time    (                 )
    timeRoot  =   basicRoot            %time    (                 )
    ! Compute the initial mass of the elliptical galaxy.
    massStellar    =+self%massStellarFinal                               &
         &          *exp(                                                &
         &               +self%rateStarFormationSpecific                 &
         &               *(                                              &
         &                 +timeLeaf                                     &
         &                 -timeRoot                                     &
         &                )                                              &
         &              )
    ! Compute the initial angular momentum the elliptical galaxy.
    angularMomentum=+self%angularMomentumPseudoSpecificFinal             &
         &          *exp(                                                &
         &               +self%rateAngularMomentumPseudoSpecificSpecific &
         &               *(                                              &
         &                 +timeLeaf                                     &
         &                 -timeRoot                                     &
         &                )                                              &
         &              )                                                &
         &          *massStellar
    ! Set the properties
    call spheroid%    massStellarSet(massStellar    )
    call spheroid%angularMomentumSet(angularMomentum)
    return
  end subroutine empiricalMassiveEllipticalNodeInitialize

  subroutine empiricalMassiveEllipticalDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Impose the star formation rate for the massive elliptical.
    !!}
    use :: Galacticus_Nodes, only : propertyInactive, nodeComponentSpheroid
    implicit none
    class           (nodeOperatorEmpiricalMassiveElliptical), intent(inout), target  :: self
    type            (treeNode                              ), intent(inout), target  :: node
    logical                                                 , intent(inout)          :: interrupt
    procedure       (interruptTask                         ), intent(inout), pointer :: functionInterrupt
    integer                                                 , intent(in   )          :: propertyType
    class           (nodeComponentSpheroid                 )               , pointer :: spheroid
    double precision                                                                 :: rateStarFormation, rateAngularMomentum

    ! Return immediately for inactive property evolution.
    if (propertyInactive(propertyType).or..not.node%isOnMainBranch()) return
    ! Compute and set the rates of star formation and angular momentum growth.
    spheroid               =>    node    %spheroid                                 ()
    rateStarFormation      =  +  spheroid%massStellar                              () &
         &                    *  self    %rateStarFormationSpecific
    call    spheroid%    massStellarRate(rateStarFormation  )
    if (self%useAngularMomentum) then
       rateAngularMomentum =  +  spheroid%angularMomentum                          () &
            &                 *(                                                      &
            &                   +self    %rateStarFormationSpecific                   &
            &                   +self    %rateAngularMomentumPseudoSpecificSpecific   &
            &                   )
       call spheroid%angularMomentumRate(rateAngularMomentum)
    end if
    return
  end subroutine empiricalMassiveEllipticalDifferentialEvolution

  subroutine empiricalMassiveEllipticalSolveAnalytics(self,node,time)
    !!{
    Set radii of empirical elliptical galaxies.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSpheroid
    implicit none
    class           (nodeOperatorEmpiricalMassiveElliptical), intent(inout) :: self
    type            (treeNode                              ), intent(inout) :: node
    double precision                                        , intent(in   ) :: time
    type            (treeNode                              ), pointer       :: nodeRoot
    class           (nodeComponentBasic                    ), pointer       :: basicLeaf, basicRoot
    class           (nodeComponentSpheroid                 ), pointer       :: spheroid
    double precision                                                        :: timeLeaf , timeRoot , &
         &                                                                     radius

    if (self%useAngularMomentum.or..not.node%isOnMainBranch()) return
    ! Get times in the leaf and root nodes.
    nodeRoot  =>  node     %hostTree   %nodeBase
    basicLeaf =>  node                 %basic   ()
    basicRoot =>  nodeRoot             %basic   ()
    spheroid  =>  node                 %spheroid()
    timeLeaf  =   basicLeaf            %time    ()
    timeRoot  =   basicRoot            %time    ()
    ! Compute the initial mass of the elliptical galaxy.
    radius    =  +self     %radiusFinal        &
         &       *exp(                         &
         &            +self%rateRadiusSpecific &
         &            *(                       &
         &              +timeLeaf              &
         &              -timeRoot              &
         &             )                       &
         &           )
    call spheroid%radiusSet(radius)
    return
  end subroutine empiricalMassiveEllipticalSolveAnalytics
