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
  Implements a dust attenuation class which averages another over orientation.
  !!}

  !![
  <dustAttenuation name="dustAttenuationInclinationAveraged" docformat="rst">
   <description>
   A dust attenuation class which averages the transmission of another attenuator over orientation, for models in
   which no inclination is available for each galaxy---which is the default, see
   :galacticus-class:`galacticInclinationNull`.

   Orientations are isotropic, so uniform in :math:`\cos i`, and the average is

   .. math::

      \bar{T}(\lambda) = \int_0^1 T(\lambda,i) \, \mathrm{d}\cos i.

   It is the *transmission* which is averaged, not the optical depth. Those are not the same thing: attenuation is
   exponential in optical depth, so a galaxy seen at a range of orientations transmits more light than one seen at
   the mean optical depth, and averaging the depth instead would systematically over-attenuate. What an observer
   averaging over an orientation-blind sample measures is the mean flux, which is the mean transmission.

   The wrapped attenuator is evaluated at each abscissa through the optional ``inclination`` argument of its
   ``transmission`` method, so nothing is mutated and the arrangement is safe under threading. Wrapping a
   :galacticus-class:`dustAttenuationSequence` averages the sequence as a whole---the product of its members
   evaluated at a common orientation---rather than the product of separately averaged members, which is a different
   and incorrect quantity.

   The integral is evaluated by Gauss-Legendre quadrature of order ``order``. Fixed-order quadrature is used in
   preference to an adaptive rule because the integrand is smooth and bounded in :math:`[0,1]`, so it converges
   quickly, and because the cost of the average multiplies the cost of every luminosity: a fixed order makes that
   cost predictable and puts it under the user's control. Order 8 is accurate to better than one part in
   :math:`10^{6}` for the attenuation curves of interest.
   </description>
  </dustAttenuation>
  !!]
  type, extends(dustAttenuationClass) :: dustAttenuationInclinationAveraged
     !!{RST
     A dust attenuation class which averages another over orientation.
     !!}
     private
     class           (dustAttenuationClass), pointer                   :: dustAttenuation_ => null()
     integer                                                           :: order
     ! Abscissae in cos(i) and their weights, computed once at construction.
     double precision                      , allocatable, dimension(:) :: cosineInclination         , weight
   contains
     final     ::                      inclinationAveragedDestructor
     procedure :: transmission      => inclinationAveragedTransmission
     procedure :: request           => inclinationAveragedRequest
     procedure :: supportsComponent => inclinationAveragedSupportsComponent
  end type dustAttenuationInclinationAveraged

  interface dustAttenuationInclinationAveraged
     !!{RST
     Constructors for the :galacticus-class:`dustAttenuationInclinationAveraged` dust attenuation class.
     !!}
     module procedure inclinationAveragedConstructorParameters
     module procedure inclinationAveragedConstructorInternal
  end interface dustAttenuationInclinationAveraged

contains

  function inclinationAveragedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustAttenuationInclinationAveraged` dust attenuation class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (dustAttenuationInclinationAveraged)                :: self
    type   (inputParameters                   ), intent(inout) :: parameters
    class  (dustAttenuationClass              ), pointer       :: dustAttenuation_
    integer                                                    :: order

    !![
    <inputParameter docformat="rst">
      <name>order</name>
      <defaultValue>8</defaultValue>
      <description>
      The order of the Gauss-Legendre quadrature used to average over orientation. The cost of the average is
      proportional to this.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="dustAttenuation" name="dustAttenuation_" source="parameters"/>
    !!]
    self=dustAttenuationInclinationAveraged(order,dustAttenuation_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="dustAttenuation_"/>
    !!]
    return
  end function inclinationAveragedConstructorParameters

  function inclinationAveragedConstructorInternal(order,dustAttenuation_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustAttenuationInclinationAveraged` dust attenuation class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (dustAttenuationInclinationAveraged)                        :: self
    integer                                    , intent(in   )         :: order
    class  (dustAttenuationClass              ), intent(in   ), target :: dustAttenuation_
    !![
    <constructorAssign variables="order, *dustAttenuation_"/>
    !!]

    if (self%order < 1) call Error_Report('`order` must be positive'//{introspection:location})
    call gaussLegendreRule(self%order,self%cosineInclination,self%weight)
    return
  end function inclinationAveragedConstructorInternal

  subroutine inclinationAveragedDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`dustAttenuationInclinationAveraged` dust attenuation class.
    !!}
    implicit none
    type(dustAttenuationInclinationAveraged), intent(inout) :: self

    !![
    <objectDestructor name="self%dustAttenuation_"/>
    !!]
    return
  end subroutine inclinationAveragedDestructor

  subroutine gaussLegendreRule(order,abscissae,weights)
    !!{RST
    Return the abscissae and weights of the Gauss-Legendre rule of the given ``order`` on the interval
    :math:`[0,1]`.

    The nodes are the roots of the Legendre polynomial of that order, found by Newton iteration from the standard
    Chebyshev-like starting guess, with the polynomial and its derivative evaluated by the usual recurrence. Both are
    then mapped from :math:`[-1,1]` onto :math:`[0,1]`, and the weights scaled by the half-width of the interval so
    that they sum to unity---which is what makes the result an average rather than an integral.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    integer                       , intent(in   )               :: order
    double precision, allocatable , intent(inout), dimension(:) :: abscissae                     , weights
    double precision              , parameter                   :: toleranceRelative     =1.0d-15
    integer                       , parameter                   :: countIterationsMaximum=100
    double precision                                            :: root                          , rootPrevious , &
         &                                                         legendre                      , legendreLower, &
         &                                                         legendreLowerLower            , derivative
    integer                                                     :: i                             , j            , &
         &                                                         countIterations

    if (allocated(abscissae)) deallocate(abscissae)
    if (allocated(weights  )) deallocate(weights  )
    allocate(abscissae(order))
    allocate(weights  (order))
    do i=1,order
       ! Initial guess for the i'th root of the Legendre polynomial of this order.
       root=cos(Pi*(dble(i)-0.25d0)/(dble(order)+0.5d0))
       countIterations=0
       do
          ! Evaluate the Legendre polynomial and its derivative at the current estimate by recurrence.
          legendre     =1.0d0
          legendreLower=0.0d0
          do j=1,order
             legendreLowerLower=legendreLower
             legendreLower     =legendre
             legendre          =(dble(2*j-1)*root*legendreLower-dble(j-1)*legendreLowerLower)/dble(j)
          end do
          derivative  =dble(order)*(root*legendre-legendreLower)/(root**2-1.0d0)
          rootPrevious=root
          root        =rootPrevious-legendre/derivative
          countIterations=countIterations+1
          if (abs(root-rootPrevious) <= toleranceRelative*abs(root) .or. countIterations >= countIterationsMaximum) exit
       end do
       ! Map from [-1,1] onto [0,1]. The weights are halved along with the interval, so that they sum to unity and
       ! the quadrature returns a mean.
       abscissae(i)=0.5d0*(1.0d0-root)
       weights  (i)=1.0d0/((1.0d0-root**2)*derivative**2)
    end do
    return
  end subroutine gaussLegendreRule

  function inclinationAveragedTransmission(self,node,descriptors,inclination) result(transmission)
    !!{RST
    Return the transmission of the wrapped attenuator, averaged over orientation.

    Any ``inclination`` supplied by a caller is ignored: this class exists precisely to remove the dependence on
    orientation, so honoring an imposed angle would defeat it. Nesting one of these inside another is therefore
    harmless but pointless.
    !!}
    implicit none
    class           (dustAttenuationInclinationAveraged), intent(inout)                               :: self
    type            (treeNode                          ), intent(inout), target                       :: node
    type            (emissionDescriptor                ), intent(in   ), dimension(:                ) :: descriptors
    double precision                                    , intent(in   ), optional                     :: inclination
    double precision                                                   , dimension(size(descriptors)) :: transmission
    integer                                                                                           :: i
    !$GLC attributes unused :: inclination

    transmission=0.0d0
    do i=1,self%order
       transmission=+transmission                                                        &
            &       +self%weight(i)                                                      &
            &       *self%dustAttenuation_%transmission(                                 &
            &                                           node                           , &
            &                                           descriptors                    , &
            &                                           acos(self%cosineInclination(i))  &
            &                                          )
    end do
    return
  end function inclinationAveragedTransmission

  logical function inclinationAveragedSupportsComponent(self,componentType) result(supportsComponent)
    !!{RST
    Return whether the wrapped attenuator supports the given component: averaging over orientation changes nothing
    about which components can be attenuated.

    Forwarding this matters. The base class accepts any individual component, so a decorator which did not forward
    would make the construction-time check pass for a component its wrapped attenuator rejects, and that component
    would then be attenuated anyway.
    !!}
    implicit none
    class(dustAttenuationInclinationAveraged), intent(inout) :: self
    type (enumerationComponentTypeType      ), intent(in   ) :: componentType

    supportsComponent=self%dustAttenuation_%supportsComponent(componentType)
    return
  end function inclinationAveragedSupportsComponent

  function inclinationAveragedRequest(self) result(request)
    !!{RST
    Return the decomposition request of the wrapped attenuator: averaging over orientation changes nothing about
    which parcels must be distinguished.
    !!}
    implicit none
    type (decompositionRequest              )                :: request
    class(dustAttenuationInclinationAveraged), intent(inout) :: self

    request=self%dustAttenuation_%request()
    return
  end function inclinationAveragedRequest
