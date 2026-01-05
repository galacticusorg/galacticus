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
  Implements a node operator class that inserts an empirical model of the formation history of a central disk galaxy.
  !!}

  !![
  <nodeOperator name="nodeOperatorEmpiricalCentralDisk">
   <description>A node operator class that inserts an empirical model of the formation history of a central disk galaxy.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorEmpiricalCentralDisk
     !!{     
     A node operator class that inserts an empirical model of the formation history of a central disk galaxy. The galaxy is
     assumed to grow in the main branch of the tree with a constant specific star formation rate, such that it mass is given by:
     \begin{equation}
       M_\star(t) = M_{\star,0} \exp(-\phi_\star [t-t_0]),
     \end{equation}
     where $M_{\star,0}=${\normalfont \ttfamily [massStellarFinal]} is the stellar mass in the root node of the tree,
     $\phi_\star=${\normalfont \ttfamily [rateStarFormationSpecific]}, and $t_0$ is the cosmic time at the root node of the tree.
     !!}
     private
     double precision :: massStellarFinal            , rateStarFormationSpecific          , &
          &              angularMomentumSpecificFinal, rateAngularMomentumSpecificSpecific
   contains
     procedure :: nodeInitialize        => empiricalCentralDiskNodeInitialize
     procedure :: differentialEvolution => empiricalCentralDiskDifferentialEvolution
  end type nodeOperatorEmpiricalCentralDisk
  
  interface nodeOperatorEmpiricalCentralDisk
     !!{
     Constructors for the \refClass{nodeOperatorEmpiricalCentralDisk} node operator class.
     !!}
     module procedure empiricalCentralDiskConstructorParameters
     module procedure empiricalCentralDiskConstructorInternal
  end interface nodeOperatorEmpiricalCentralDisk
  
contains

  function empiricalCentralDiskConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorEmpiricalCentralDisk} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorEmpiricalCentralDisk)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    double precision                                                  :: massStellarFinal            , rateStarFormationSpecific          , &
         &                                                               angularMomentumSpecificFinal, rateAngularMomentumSpecificSpecific
    
    !![
    <inputParameter>
      <name>massStellarFinal</name>
      <source>parameters</source>
      <description>The final stellar mass of the disk galaxy.</description>
    </inputParameter>
    <inputParameter>
      <name>rateStarFormationSpecific</name>
      <source>parameters</source>
      <description>The specific star formation rate of the disk galaxy.</description>
    </inputParameter>
    <inputParameter>
      <name>angularMomentumSpecificFinal</name>
      <source>parameters</source>
      <description>The final specific angular momentum of the disk galaxy.</description>
    </inputParameter>
    <inputParameter>
      <name>rateAngularMomentumSpecificSpecific</name>
      <source>parameters</source>
      <description>The specific growth rate of the specific angular momentum of the disk galaxy.</description>
    </inputParameter>
    !!]
    self=nodeOperatorEmpiricalCentralDisk(massStellarFinal,rateStarFormationSpecific,angularMomentumSpecificFinal,rateAngularMomentumSpecificSpecific)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function empiricalCentralDiskConstructorParameters

  function empiricalCentralDiskConstructorInternal(massStellarFinal,rateStarFormationSpecific,angularMomentumSpecificFinal,rateAngularMomentumSpecificSpecific) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorEmpiricalCentralDisk} node operator class.
    !!}
    implicit none
    type            (nodeOperatorEmpiricalCentralDisk)                :: self
    double precision                                  , intent(in   ) :: massStellarFinal            , rateStarFormationSpecific          , &
         &                                                               angularMomentumSpecificFinal, rateAngularMomentumSpecificSpecific
    !![
    <constructorAssign variables="massStellarFinal, rateStarFormationSpecific, angularMomentumSpecificFinal, rateAngularMomentumSpecificSpecific"/>
    !!]

    return
  end function empiricalCentralDiskConstructorInternal

  subroutine empiricalCentralDiskNodeInitialize(self,node)
    !!{
    Initialize nodes for the massive elliptical.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDisk
    implicit none
    class           (nodeOperatorEmpiricalCentralDisk), intent(inout), target  :: self
    type            (treeNode                        ), intent(inout), target  :: node
    type            (treeNode                        )               , pointer :: nodeRoot
    class           (nodeComponentBasic              )               , pointer :: basicLeaf  , basicRoot
    class           (nodeComponentDisk               )               , pointer :: disk
    double precision                                                           :: timeLeaf   , timeRoot       , &
         &                                                                        massStellar, angularMomentum

    ! Initialize only the leaf node on the main branch.
    if (associated(node%firstChild).or..not.node%isOnMainBranch()) return
    ! Get times in the leaf and root nodes.
    nodeRoot => node
    do while (associated(nodeRoot%parent))
       nodeRoot => nodeRoot%parent
    end do
    basicLeaf => node     %basic(                 )
    basicRoot => nodeRoot %basic(                 )
    disk      => node     %disk (autoCreate=.true.)
    timeLeaf  =  basicLeaf%time (                 )
    timeRoot  =  basicRoot%time (                 )
    ! Compute the initial mass of the elliptical galaxy.
    massStellar    =+self%massStellarFinal                         &
         &          *exp(                                          &
         &               +self%rateStarFormationSpecific           &
         &               *(                                        &
         &                 +timeLeaf                               &
         &                 -timeRoot                               &
         &                )                                        &
         &              )
    ! Compute the initial angular momentum the elliptical galaxy.
    angularMomentum=+self%angularMomentumSpecificFinal             &
         &          *exp(                                          &
         &               +self%rateAngularMomentumSpecificSpecific &
         &               *(                                        &
         &                 +timeLeaf                               &
         &                 -timeRoot                               &
         &                )                                        &
         &              )                                          &
         &          *massStellar
    ! Set the properties
    call disk%    massStellarSet(massStellar    )
    call disk%angularMomentumSet(angularMomentum)
    return
  end subroutine empiricalCentralDiskNodeInitialize

  subroutine empiricalCentralDiskDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Impose the star formation rate for the massive elliptical.
    !!}
    use :: Galacticus_Nodes, only : propertyInactive, nodeComponentDisk
    implicit none
    class           (nodeOperatorEmpiricalCentralDisk), intent(inout), target  :: self
    type            (treeNode                        ), intent(inout), target  :: node
    logical                                           , intent(inout)          :: interrupt
    procedure       (interruptTask                   ), intent(inout), pointer :: functionInterrupt
    integer                                           , intent(in   )          :: propertyType
    class           (nodeComponentDisk               )               , pointer :: disk
    double precision                                                           :: rateStarFormation, rateAngularMomentum

    ! Return immediately for inactive property evolution.
    if (propertyInactive(propertyType)) return
    ! Compute and set the rates of star formation and angular momentum growth.
    disk                =>    node%disk                               ()
    rateStarFormation   =  +  disk%massStellar                        () &
         &                 *  self%rateStarFormationSpecific
    rateAngularMomentum =  +  disk%angularMomentum                    () &
         &                 *(                                            &
         &                   +self%rateStarFormationSpecific             &
         &                   +self%rateAngularMomentumSpecificSpecific   &
         &                   )
    call disk%    massStellarRate(rateStarFormation  )
    call disk%angularMomentumRate(rateAngularMomentum)
    return
  end subroutine empiricalCentralDiskDifferentialEvolution

