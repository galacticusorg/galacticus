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
  Implements a node operator class that implements triple interactions between black holes.
  !!}
  
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <nodeOperator name="nodeOperatorBlackHolesTripleInteraction">
   <description>A node operator class that implements triple interactions between black holes.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorBlackHolesTripleInteraction
     !!{
     A node operator class that implements triple interactions between black holes.
     !!}
     private
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::                          blackHolesTripleInteractionDestructor
     procedure :: differentialEvolution => blackHolesTripleInteractionDifferentialEvolution
  end type nodeOperatorBlackHolesTripleInteraction
  
  interface nodeOperatorBlackHolesTripleInteraction
     !!{
     Constructors for the \refClass{nodeOperatorBlackHolesTripleInteraction} node operator class.
     !!}
     module procedure blackHolesTripleInteractionConstructorParameters
     module procedure blackHolesTripleInteractionConstructorInternal
  end interface nodeOperatorBlackHolesTripleInteraction
  
  ! Submodule-scope objects used in callback functions.
  integer                                                   :: instanceBinary, instanceTriple
  class  (nodeOperatorBlackHolesTripleInteraction), pointer :: self_
  !$omp threadprivate(self_,instanceBinary,instanceTriple)

contains

  function blackHolesTripleInteractionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorBlackHolesTripleInteraction} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorBlackHolesTripleInteraction)                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass               ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodeOperatorBlackHolesTripleInteraction(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function blackHolesTripleInteractionConstructorParameters

  function blackHolesTripleInteractionConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorBlackHolesTripleInteraction} node operator class.
    !!}
    implicit none
    type (nodeOperatorBlackHolesTripleInteraction)                        :: self
    class(darkMatterHaloScaleClass               ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]
    
    return
  end function blackHolesTripleInteractionConstructorInternal

  subroutine blackHolesTripleInteractionDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorBlackHolesTripleInteraction} node operator class.
    !!}
    implicit none
    type(nodeOperatorBlackHolesTripleInteraction), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine blackHolesTripleInteractionDestructor

  subroutine blackHolesTripleInteractionDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Account for accretion onto black holes.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBlackHole        , propertyInactive
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (nodeOperatorBlackHolesTripleInteraction), intent(inout), target  :: self
    type            (treeNode                               ), intent(inout), target  :: node
    logical                                                  , intent(inout)          :: interrupt
    procedure       (interruptTask                          ), intent(inout), pointer :: functionInterrupt
    integer                                                  , intent(in   )          :: propertyType
    class           (nodeComponentBlackHole                 )               , pointer :: blackHoleBinary  , blackHoleCentral, &
         &                                                                               blackHole
    integer                                                                           :: iInstance        , instanceCount
    double precision                                                                  :: binaryRadius     , radiusHardBinary
    logical                                                                           :: binaryRadiusFound

    if (propertyInactive(propertyType)) return
    ! Get a count of the number of black holes associated with this node.
    instanceCount=node%blackHoleCount()
    ! Get the central black hole.
    blackHoleCentral => node%blackHole(instance=1)
    ! Iterate over black holes, testing for triple black hole interactions. Find the three closest black holes then check if a
    ! three body interaction occurs using the radial condition derived in Hoffman and Loeb (2007;
    ! http://adsabs.harvard.edu/abs/2007MNRAS.377..957H).
    binaryRadiusFound=.false.
    binaryRadius     =huge(1.0d0)
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
             instanceBinary   =iInstance
             binaryRadiusFound=.true.
          end if
       end do
       if (binaryRadiusFound) then
          ! Get the binary black hole.
          blackHoleBinary => node%blackHole(instance=instanceBinary)
          ! Compute the hard binary radius.
          radiusHardBinary= (                                                   &
               &              gravitationalConstant_internal                    &
               &             *(                                                 &
               &                blackHoleCentral%mass()                         &
               &               + blackHoleBinary%mass()                         &
               &              )                                                 &
               &            )                                                   &
               &           /(                                                   &
               &              4.0d0                                             &
               &             *self%darkMatterHaloScale_%velocityVirial(node)**2 &
               &            )
          ! Search for a third black hole.
          do iInstance=2,instanceCount
             ! Get the black hole.
             blackHole => node%blackHole(instance=iInstance)
             if     (      .not. iInstance                         == instanceBinary   &
                  &  .and. .not. blackHole%mass                 () <= 0.0d0            &
                  &  .and. .not. blackHole%radialPosition       () <= 0.0d0            &
                  &  .and.       blackHole%radialPosition       () <= radiusHardBinary &
                  &  .and.       blackHole%tripleInteractionTime() == 0.0d0            &
                  & ) then
                self_             => self
                instanceTriple    =  iInstance
                interrupt         =  .true.
                functionInterrupt => interactionTriple
                return
             end if
          end do
       end if
    end if
    return
  end subroutine blackHolesTripleInteractionDifferentialEvolution
  
  subroutine interactionTriple(node,timeEnd)
    !!{
    Handles triple black holes interactions, using conditions similar to those of \cite{volonteri_assembly_2003}.
    !!}
    use :: Galacticus_Nodes                     , only : nodeComponentBasic            , nodeComponentBlackHole, treeNode
    use :: Numerical_Constants_Astronomical     , only : gravitationalConstant_internal
    use :: Nodes_Operators_Black_Holes_Utilities, only : blackHolesRecoilEscapes
    implicit none
    type            (treeNode              ), intent(inout), target   :: node
    double precision                        , intent(in   ), optional :: timeEnd
    class           (nodeComponentBasic    )               , pointer  :: basic
    class           (nodeComponentBlackHole)               , pointer  :: blackHoleBinary          , blackHoleCentral           , &
         &                                                               ejectedBlackHoleComponent, newBinaryBlackHoleComponent, &
         &                                                               blackHoleTriple
    integer                                                           :: ejectedInstance          , newBinaryInstance
    double precision                                                  :: bindingEnergy            , kineticEnergyChange        , &
         &                                                               massBinary               , massEjected                , &
         &                                                               massRatioIntruder        , newRadius                  , &
         &                                                               velocityBinary           , velocityEjected            , &
         &                                                               radiusVirial
    logical                                                           :: removeBinary             , removeEjected
    !$GLC attributes unused :: timeEnd

    ! Get the basic component.
    basic            => node%basic    (                       )
    ! Get the black holes.
    blackHoleCentral => node%blackHole(instance=             1)
    blackHoleBinary  => node%blackHole(instance=instanceBinary)
    blackHoleTriple  => node%blackHole(instance=instanceTriple)
    ! We have to distinguish two cases, where a different black hole is ejected, the one with the lowest mass.
    massRatioIntruder=+  blackHoleTriple %mass() &
         &            /(                         &
         &              +blackHoleCentral%mass() &
         &              +blackHoleBinary %mass() &
         &             )
    ! Set the triple interaction time for the triple black hole.
    call blackHoleTriple%tripleInteractionTimeSet(basic%time())
    ! Branch on intruder mass ratio.
    if (massRatioIntruder <= 2.0d0) then
       if (blackHoleTriple%mass() <= blackHoleBinary%mass()) then
          newRadius           = blackHoleBinary%radialPosition()/(1.0d0+0.4d0*massRatioIntruder)
          call blackHoleBinary%radialPositionSet(newRadius)
          bindingEnergy       =+gravitationalConstant_internal     &
               &               *(                                  &
               &                 +blackHoleTriple %mass         () &
               &                 *blackHoleCentral%mass         () &
               &                )                                  &
               &               /  blackHoleTriple%radialPosition()
          kineticEnergyChange=0.4d0*massRatioIntruder*bindingEnergy
          ejectedInstance    =instanceTriple
          newBinaryInstance  =instanceBinary
          ejectedBlackHoleComponent   => blackHoleTriple
          newBinaryBlackHoleComponent => blackHoleBinary
       else
          newRadius          = blackHoleTriple%radialPosition()/(1.0d0+0.4d0*massRatioIntruder )
          call blackHoleTriple%radialPositionSet(newRadius)
          bindingEnergy      =+gravitationalConstant_internal      &
               &              *(                                   &
               &                +blackHoleBinary %mass          () &
               &                *blackHoleCentral%mass          () &
               &               )                                   &
               &              /  blackHoleBinary %radialPosition()
          kineticEnergyChange=0.4d0*massRatioIntruder*bindingEnergy
          ejectedInstance    =instanceBinary
          newBinaryInstance  =instanceTriple
          ejectedBlackHoleComponent   => blackHoleBinary
          newBinaryBlackHoleComponent => blackHoleTriple
       end if
    else
       ! This latter case can be referred to as head-on collision.
       newRadius             =0.53d0*blackHoleTriple%radialPosition()
       call blackHoleTriple%radialPositionSet(newRadius)
       bindingEnergy         =+gravitationalConstant_internal      &
            &                 *(                                   &
            &                   +blackHoleBinary %mass          () &
            &                   *blackHoleCentral%mass          () &
            &                  )                                   &
            &                 /  blackHoleBinary %radialPosition()
       kineticEnergyChange         =  0.9d0*massRatioIntruder*bindingEnergy
       ejectedInstance             =  instanceBinary
       newBinaryInstance           =  instanceTriple
       ejectedBlackHoleComponent   => blackHoleBinary
       newBinaryBlackHoleComponent => blackHoleTriple
    end if
    ! First we find the lightest black hole and tag it as being ejected.
    massEjected= ejectedBlackHoleComponent  %mass()
    massBinary = newBinaryBlackHoleComponent%mass() &
         &      +blackHoleCentral           %mass()
    velocityEjected=sqrt(kineticEnergyChange/(1.0d0+massEjected/massBinary )/massEjected*2.0d0)
    velocityBinary =sqrt(kineticEnergyChange/(1.0d0+massBinary /massEjected)/massBinary *2.0d0)
    ! Determine whether the ejected black hole is actually ejected.
    radiusVirial =self_%darkMatterHaloScale_%radiusVirial(node)
    removeEjected=blackHolesRecoilEscapes(node,ejectedBlackHoleComponent  %radialPosition(),radiusVirial,velocityEjected,ignoreCentralBlackHole=.false.)
    ! Determine whether the binary black hole is ejected.
    removeBinary =blackHolesRecoilEscapes(node,newBinaryBlackHoleComponent%radialPosition(),radiusVirial,velocityBinary ,ignoreCentralBlackHole=.true. )
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
  end subroutine interactionTriple
