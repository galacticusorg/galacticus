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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a galactic inclination class giving the inclination of a component whose axis is misaligned from that of
  its galaxy.
  !!}

  !![
  <galacticInclination name="galacticInclinationMisaligned" docformat="rst">
   <description>
   A galactic inclination class giving the inclination to the line of sight not of a galaxy but of a component of it
   whose symmetry axis is misaligned from that of the galaxy---the dusty torus around an AGN, for example, whose axis
   need not be that of the disk of its host.

   The orientation of the galaxy is given by the ``[galacticInclination]`` object which this class wraps. The axis of
   the component lies at an angle :math:`\psi` from that of the galaxy, at an azimuth :math:`\phi` about it measured
   from the plane containing the axis of the galaxy and the line of sight. If the galaxy is seen at inclination
   :math:`i` the component is seen at the inclination :math:`\theta` given by

   .. math::

      \cos\theta = \left| \cos i \cos\psi + \sin i \sin\psi \cos\phi \right|,

   the absolute value folding :math:`\theta` into :math:`[0,\pi/2]`, since only the angle between the axis and the
   line of sight matters.

   **Distribution of misalignments.** The axes of components are distributed isotropically within a cone of
   half-angle ``[angleMisalignmentMaximum]`` (in degrees) about the axis of the galaxy, so that :math:`\cos\psi` is
   uniform on :math:`[\cos\psi_\mathrm{max},1]` and :math:`\phi` uniform on :math:`[0,2\pi)`. A half-angle of
   :math:`90^\circ` makes the axis of the component independent of that of the galaxy---isotropically distributed, and
   uncorrelated with it---while a half-angle of zero aligns the two.

   **For each galaxy.** When the wrapped class supplies an inclination for each galaxy, so does this class, and the
   misalignment of each galaxy is then drawn once, by :galacticus-class:`nodeOperatorMisalignmentRandom`, which must be
   used with it. That operator stores a uniform deviate from which :math:`\psi` is found here, rather than
   :math:`\psi` itself, so that the distribution of misalignments is set in one place, by this class.

   **Averaging over orientation.** Consumers which average over the orientation of the galaxy supply its inclination
   to ``inclinationsMisaligned``, which returns the inclinations of the component to average over at that inclination
   of the galaxy. The azimuth :math:`\phi`, like :math:`i`, depends on the direction from which the galaxy is seen,
   so it is always averaged over, by the midpoint rule with ``[orderAzimuth]`` points on :math:`[0,\pi]` (the
   integrand being even in :math:`\phi`). If no misalignment is available for each galaxy---because the wrapped class
   supplies no inclination---the average is also taken over the distribution of :math:`\psi`, by a Gauss-Legendre rule
   of order ``[orderMisalignment]`` in :math:`\cos\psi`. Otherwise the galaxy's own :math:`\psi` is used, so that,
   for example, the luminosity absorbed by a torus and by the dust of its host galaxy is that for its own
   misalignment.
   </description>
  </galacticInclination>
  !!]
  type, extends(galacticInclinationClass) :: galacticInclinationMisaligned
     !!{RST
     A galactic inclination class giving the inclination of a component whose axis is misaligned from that of its
     galaxy.
     !!}
     private
     class           (galacticInclinationClass), pointer                   :: galacticInclination_         => null()
     ! The maximum misalignment is stored in degrees, as the parameter is given, so that a descriptor round-trips in
     ! the same units.
     double precision                                                      :: angleMisalignmentMaximum              , cosineMisalignmentMaximum
     integer                                                               :: orderMisalignment                     , orderAzimuth                , &
          &                                                                   misalignmentDeviateID                 , misalignmentAzimuthID
     logical                                                               :: perGalaxy
     ! The Gauss-Legendre rule in the cosine of the misalignment, used when averaging over its distribution.
     double precision                          , allocatable, dimension(:) :: cosineMisalignmentQuadrature          , weightMisalignmentQuadrature
   contains
     final     ::                           misalignedDestructor
     procedure :: inclination            => misalignedInclination
     procedure :: isAvailable            => misalignedIsAvailable
     procedure :: inclinationsMisaligned => misalignedInclinationsMisaligned
  end type galacticInclinationMisaligned

  interface galacticInclinationMisaligned
     !!{RST
     Constructors for the :galacticus-class:`galacticInclinationMisaligned` galactic inclination class.
     !!}
     module procedure misalignedConstructorParameters
     module procedure misalignedConstructorInternal
  end interface galacticInclinationMisaligned

contains

  function misalignedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`galacticInclinationMisaligned` galactic inclination class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticInclinationMisaligned)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (galacticInclinationClass     ), pointer       :: galacticInclination_
    double precision                                               :: angleMisalignmentMaximum
    integer                                                        :: orderMisalignment       , orderAzimuth

    !![
    <inputParameter docformat="rst">
      <name>angleMisalignmentMaximum</name>
      <defaultValue>90.0d0</defaultValue>
      <description>
      The half-angle, in degrees, of the cone about the axis of the galaxy within which the axes of components are
      isotropically distributed. 90 makes the axis of the component independent of that of the galaxy, and zero aligns
      them.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>orderMisalignment</name>
      <defaultValue>8</defaultValue>
      <description>
      The order of the Gauss-Legendre rule in the cosine of the misalignment used to average over its distribution.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>orderAzimuth</name>
      <defaultValue>16</defaultValue>
      <description>
      The number of points of the midpoint rule in azimuth used to average over the azimuth of the axis of the
      component about that of the galaxy.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="galacticInclination" name="galacticInclination_" source="parameters"/>
    !!]
    self=galacticInclinationMisaligned(angleMisalignmentMaximum,orderMisalignment,orderAzimuth,galacticInclination_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticInclination_"/>
    !!]
    return
  end function misalignedConstructorParameters

  function misalignedConstructorInternal(angleMisalignmentMaximum,orderMisalignment,orderAzimuth,galacticInclination_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`galacticInclinationMisaligned` galactic inclination class. The
    maximum misalignment is given in degrees.
    !!}
    use :: Error                              , only : Error_Report
    use :: Numerical_Constants_Astronomical   , only : degreesToRadians
    use :: Numerical_Quadrature_Gauss_Legendre, only : gaussLegendreRule
    implicit none
    type            (galacticInclinationMisaligned)                        :: self
    double precision                               , intent(in   )         :: angleMisalignmentMaximum
    integer                                        , intent(in   )         :: orderMisalignment       , orderAzimuth
    class           (galacticInclinationClass     ), intent(in   ), target :: galacticInclination_
    !![
    <constructorAssign variables="angleMisalignmentMaximum, orderMisalignment, orderAzimuth, *galacticInclination_"/>
    !!]

    if (angleMisalignmentMaximum < 0.0d0 .or. angleMisalignmentMaximum > 90.0d0) call Error_Report('`angleMisalignmentMaximum` must be in the range [0,90] degrees'//{introspection:location})
    if (orderMisalignment        < 1                                           ) call Error_Report('`orderMisalignment` must be positive'                          //{introspection:location})
    if (orderAzimuth             < 1                                           ) call Error_Report('`orderAzimuth` must be positive'                               //{introspection:location})
    self%cosineMisalignmentMaximum=cos(angleMisalignmentMaximum*degreesToRadians)
    ! The rule is on [0,1], with weights summing to one; map it onto [cos(ψ_max),1].
    call gaussLegendreRule(orderMisalignment,self%cosineMisalignmentQuadrature,self%weightMisalignmentQuadrature)
    self%cosineMisalignmentQuadrature=+self%cosineMisalignmentMaximum                                           &
         &                            +(1.0d0-self%cosineMisalignmentMaximum)*self%cosineMisalignmentQuadrature
    ! A misalignment is needed for each galaxy only if the galaxy itself has an orientation; `nodeOperatorMisalignmentRandom`
    ! creates and sets it.
    self%perGalaxy=self%galacticInclination_%isAvailable()
    if (self%perGalaxy) then
       !![
       <addMetaProperty component="basic" name="misalignmentDeviate" id="self%misalignmentDeviateID" isEvolvable="no"/>
       <addMetaProperty component="basic" name="misalignmentAzimuth" id="self%misalignmentAzimuthID" isEvolvable="no"/>
       !!]
    end if
    return
  end function misalignedConstructorInternal

  subroutine misalignedDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`galacticInclinationMisaligned` galactic inclination class.
    !!}
    implicit none
    type(galacticInclinationMisaligned), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticInclination_"/>
    !!]
    return
  end subroutine misalignedDestructor

  logical function misalignedIsAvailable(self) result(isAvailable)
    !!{RST
    Return true if an inclination of the component is available for each galaxy, which it is if the wrapped class
    supplies one for the galaxy.
    !!}
    implicit none
    class(galacticInclinationMisaligned), intent(inout) :: self

    isAvailable=self%perGalaxy
    return
  end function misalignedIsAvailable

  double precision function misalignedInclination(self,node) result(inclination)
    !!{RST
    Return the inclination of the component of the galaxy in ``node``, from the inclination of the galaxy and the
    misalignment of the component stored for it.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (galacticInclinationMisaligned), intent(inout)         :: self
    type            (treeNode                     ), intent(inout), target :: node
    class           (nodeComponentBasic           ), pointer               :: basic
    double precision                                                       :: cosineMisalignment

    basic              =>  node%basic()
    cosineMisalignment =  +1.0d0                                                       &
         &                -(1.0d0-self%cosineMisalignmentMaximum)                      &
         &                *basic%floatRank0MetaPropertyGet(self%misalignmentDeviateID)
    inclination        =  misalignedInclinationComponent(                                                             &
         &                                               self %galacticInclination_%inclination(node)               , &
         &                                               cosineMisalignment                                         , &
         &                                               basic%floatRank0MetaPropertyGet(self%misalignmentAzimuthID)  &
         &                                              )
    return
  end function misalignedInclination

  subroutine misalignedInclinationsMisaligned(self,node,inclinationHost,inclinations,weights)
    !!{RST
    Return the inclinations of the component, and their weights, over which to average at inclination
    ``inclinationHost`` of the galaxy: over the azimuth of the component's axis, and over the distribution of its
    misalignment unless the galaxy has its own.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBasic
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (galacticInclinationMisaligned), intent(inout)                            :: self
    type            (treeNode                     ), intent(inout), target                    :: node
    double precision                               , intent(in   )                            :: inclinationHost
    double precision                               , intent(inout), allocatable, dimension(:) :: inclinations      , weights
    class           (nodeComponentBasic           ), pointer                                  :: basic
    double precision                                                                          :: cosineMisalignment, weightMisalignment
    integer                                                                                   :: countMisalignment , countAzimuth      , &
         &                                                                                       j                 , l                 , &
         &                                                                                       k

    if (self%cosineMisalignmentMaximum >= 1.0d0) then
       ! The component is aligned with the galaxy.
       countMisalignment=1
       countAzimuth     =1
    else if (self%perGalaxy) then
       countMisalignment=1
       countAzimuth     =self%orderAzimuth
    else
       countMisalignment=self%orderMisalignment
       countAzimuth     =self%orderAzimuth
    end if
    if (allocated(inclinations)) deallocate(inclinations)
    if (allocated(weights     )) deallocate(weights     )
    allocate(inclinations(countMisalignment*countAzimuth))
    allocate(weights     (countMisalignment*countAzimuth))
    if (self%perGalaxy) basic => node%basic()
    k=0
    do j=1,countMisalignment
       if (self%cosineMisalignmentMaximum >= 1.0d0) then
          cosineMisalignment=1.0d0
          weightMisalignment=1.0d0
       else if (self%perGalaxy) then
          cosineMisalignment=+1.0d0                                                       &
               &             -(1.0d0-self%cosineMisalignmentMaximum)                      &
               &             *basic%floatRank0MetaPropertyGet(self%misalignmentDeviateID)
          weightMisalignment=1.0d0
       else
          cosineMisalignment=self%cosineMisalignmentQuadrature(j)
          weightMisalignment=self%weightMisalignmentQuadrature(j)
       end if
       do l=1,countAzimuth
          k              =+k &
               &          +1
          inclinations(k)=misalignedInclinationComponent(inclinationHost,cosineMisalignment,Pi*(dble(l)-0.5d0)/dble(countAzimuth))
          weights     (k)=+weightMisalignment &
               &          /dble(countAzimuth)
       end do
    end do
    return
  end subroutine misalignedInclinationsMisaligned

  double precision function misalignedInclinationComponent(inclinationHost,cosineMisalignment,azimuth) result(inclination)
    !!{RST
    Return the inclination of a component whose axis is at an angle :math:`\psi` (given by its cosine,
    ``cosineMisalignment``) from that of its galaxy, and at an ``azimuth`` :math:`\phi` about it, when the galaxy is
    seen at inclination ``inclinationHost``, folded into :math:`[0,\pi/2]`.
    !!}
    implicit none
    double precision, intent(in   ) :: inclinationHost, cosineMisalignment, &
         &                             azimuth

    inclination=acos(                                                                                         &
         &           min(                                                                                     &
         &               1.0d0                                                                              , &
         &               abs(                                                                                 &
         &                   +cos(inclinationHost)*                     cosineMisalignment                    &
         &                   +sin(inclinationHost)*sqrt(max(0.0d0,1.0d0-cosineMisalignment**2))*cos(azimuth)  &
         &                  )                                                                                 &
         &              )                                                                                     &
         &          )
    return
  end function misalignedInclinationComponent
