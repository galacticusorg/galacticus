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
  Implements a model of the ram pressure stripping force from hot halos based on orbital position within the host halo.
  !!}

  use :: Hot_Halo_Mass_Distributions, only : hotHaloMassDistributionClass

  !![
  <hotHaloRamPressureForce name="hotHaloRamPressureForceRelativePosition">
   <description>
    A hot halo ram pressure force class which computes the force based on the current position relative to the host
    halo. Specifically, the ram pressure force is    
    \begin{equation}
    \mathcal{F}_\mathrm{ram, hot, host} = \rho_\mathrm{hot, host}(r) v^2(r),
    \end{equation}
    where $\rho_\mathrm{hot, host}(r)$ is the hot halo density profile of the nodes host halo, $v(r)$ is the orbital velocity
    of the node in that host, and $r$ is the instantaneous distance to the host halo.
   </description>
  </hotHaloRamPressureForce>
  !!]
  type, extends(hotHaloRamPressureForceClass) :: hotHaloRamPressureForceRelativePosition
     !!{
     Implementation of a hot halo ram pressure force class based on orbital position within the host halo.
     !!}
     private
     class(hotHaloMassDistributionClass), pointer :: hotHaloMassDistribution_ => null()
   contains
     final     ::          relativePositionDestructor
     procedure :: force => relativePositionForce
  end type hotHaloRamPressureForceRelativePosition

  interface hotHaloRamPressureForceRelativePosition
     !!{
     Constructors for the \refClass{hotHaloRamPressureForceRelativePosition} hot halo ram pressure force class.
     !!}
     module procedure relativePositionConstructorParameters
     module procedure relativePositionConstructorInternal
  end interface hotHaloRamPressureForceRelativePosition

contains

  function relativePositionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{hotHaloRamPressureForceRelativePosition} hot halo ram pressure force class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (hotHaloRamPressureForceRelativePosition)                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(hotHaloMassDistributionClass           ), pointer       :: hotHaloMassDistribution_

    !![
    <objectBuilder class="hotHaloMassDistribution" name="hotHaloMassDistribution_" source="parameters"/>
    !!]
    self=hotHaloRamPressureForceRelativePosition(hotHaloMassDistribution_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="hotHaloMassDistribution_"/>
    !!]
    return
  end function relativePositionConstructorParameters

  function relativePositionConstructorInternal(hotHaloMassDistribution_) result(self)
    !!{
    Internal constructor for the \refClass{hotHaloRamPressureForceRelativePosition} hot halo ram pressure force class.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Error_Report             , Component_List
    use :: Galacticus_Nodes, only : defaultPositionComponent
    implicit none
    type (hotHaloRamPressureForceRelativePosition)                        :: self
    class(hotHaloMassDistributionClass           ), intent(in   ), target :: hotHaloMassDistribution_
    !![
    <constructorAssign variables="*hotHaloMassDistribution_"/>
    !!]

    ! Ensure that required methods are supported.
    if     (                                                                                                                          &
         &  .not.                                                                                                                     &
         &       (                                                                                                                    &
         &        defaultPositionComponent%positionIsGettable().and.                                                                  &
         &        defaultPositionComponent%velocityIsGettable()                                                                       &
         &  )                                                                                                                         &
         & ) call Error_Report                                                                                                        &
         &        (                                                                                                                   &
         &         'this method requires that position, and velocity properties must all be gettable for the `position` component.'// &
         &         Component_List(                                                                                                    &
         &                        'position'                                                                                       ,  &
         &                        defaultPositionComponent%positionAttributeMatch(requireGettable=.true.).intersection.               &
         &                        defaultPositionComponent%velocityAttributeMatch(requireGettable=.true.)                             &
         &                       )                                                                                                 // &
         &         {introspection:location}                                                                                           &
         &        )
    
    return
  end function relativePositionConstructorInternal

  subroutine relativePositionDestructor(self)
    !!{
    Destructor for the \refClass{hotHaloRamPressureForceRelativePosition} hot halo ram pressure force class.
    !!}
    implicit none
    type(hotHaloRamPressureForceRelativePosition), intent(inout) :: self

    !![
    <objectDestructor name="self%hotHaloMassDistribution_"/>
    !!]
    return
  end subroutine relativePositionDestructor

  double precision function relativePositionForce(self,node) result(force)
    !!{
    Return a ram pressure force due to the hot halo based on orbital position within the host halo.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentPosition, nodeComponentBasic
    use :: Vectors                   , only : Vector_Magnitude
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Coordinates               , only : coordinateSpherical  , assignment(=)
    use :: Galactic_Structure_Options, only : componentTypeHotHalo , massTypeGaseous
    implicit none
    class           (hotHaloRamPressureForceRelativePosition), intent(inout) :: self
    type            (treeNode                               ), intent(inout) :: node
    class           (nodeComponentPosition                  ), pointer       :: position         , positionHost
    class           (nodeComponentBasic                     ), pointer       :: basic            , basicPrevious   , &
         &                                                                      basicCurrent     , basicHost
    type            (treeNode                               ), pointer       :: nodeHost         , nodeHostPrevious, &
         &                                                                      nodeHostCurrent
    class           (massDistributionClass                  ), pointer       :: massDistribution_
    type            (coordinateSpherical                    )                :: coordinates
    double precision                                                         :: radiusRelative   , velocityRelative

    ! Find the host node. Seek the descendant of the node closest in time to our satellite node. This is necessary as satellites
    ! can evolve ahead of their hosts.
    nodeHostPrevious => node%parent
    nodeHostCurrent  => node%parent
    basic            => node            %basic()
    basicCurrent     => nodeHostCurrent %basic()
    basicPrevious    => nodeHostPrevious%basic()
    do while (associated(nodeHostCurrent))
       if (basicCurrent%time() > basic%time()) exit
       nodeHostPrevious    => nodeHostCurrent
       nodeHostCurrent     => nodeHostCurrent%parent
       basicPrevious       => basicCurrent
       if (associated(nodeHostCurrent))                &
            & basicCurrent => nodeHostCurrent%basic ()
    end do
    if     (                                         &
         &    abs(basicPrevious%time()-basic%time()) &
         &   <                                       &
         &    abs(basicCurrent %time()-basic%time()) &
         &  .or.                                     &
         &   .not.associated(nodeHostCurrent)        &
         & ) then
       nodeHost  => nodeHostPrevious
       basicHost => basicPrevious
    else
       nodeHost  => nodeHostCurrent
       basicHost => basicCurrent
    end if
    ! Get the position components.
    position         =>  node    %position()
    positionHost     =>  nodeHost%position()
    ! Compute orbital position and velocity.
    radiusRelative   =  +Vector_Magnitude(position%position()-positionHost%position())
    velocityRelative =  +Vector_Magnitude(position%velocity()-positionHost%velocity())
    ! Find the ram pressure force this orbital radius.
    coordinates       =  [radiusRelative,0.0d0,0.0d0]
    massDistribution_ =>  nodeHost         %massDistribution(componentTypeHotHalo,massTypeGaseous)
    force             =  +massDistribution_%density         (coordinates                         )    &
         &               *velocityRelative                                                        **2
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function relativePositionForce
