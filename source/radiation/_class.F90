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

!!{RST
Contains a module which implements a class to describe radiation fields.
!!}

module Radiation_Fields
  !!{RST
  Implements a class to describe radiation fields.
  !!}
  use :: Galacticus_Nodes       , only : treeNode
  use :: Numerical_Interpolation, only : interpolator
  use :: Numerical_Ranges       , only : rangeLattice
  implicit none
  private
  public :: crossSectionFunctionTemplate
  
  !![
  <deepCopyActions class="rateCoefficient">
   <rateCoefficient>
    <methodCall method="interpolatorDeepCopy"/>
   </rateCoefficient>
  </deepCopyActions>
  !!]

  type :: rateCoefficient
     !!{RST
     Type used to store tables of rate coefficients.
     !!}
     private
     procedure       (crossSectionFunctionTemplate), pointer     , nopass      :: crossSectionFunction => null()
     double precision                              , dimension(2)              :: wavelengthRange
     ! Lattice to which the tabulated times are pinned, and the rate coefficients tabulated on it. The lattice is the source of
     ! truth for the extent of the tabulation, and is undefined until the first tabulation is made.
     type            (rangeLattice                )                            :: latticeTime
     double precision                              , dimension(:), allocatable :: rateCoefficient_
     type            (interpolator                )                            :: interpolator_
   contains
     !![
     <methods docformat="rst">
       <method description="Perform deep copy actions on interpolators." method="interpolatorDeepCopy" />
     </methods>
     !!]
     procedure :: interpolatorDeepCopy => rateCoefficientInterpolatorDeepCopy
  end type rateCoefficient
  
  !![
  <functionClass docformat="rst">
   <name>radiationField</name>
   <descriptiveName>Radiation Fields</descriptiveName>
   <description>
   Class providing radiation fields---the specific intensity (flux per unit frequency per steradian, in units of ergs cm\ :math:`^{-2}` s\ :math:`^{-1}` Hz\ :math:`^{-1}` ster\ :math:`^{-1}`) of a radiation background as a function of wavelength and cosmic time. Radiation fields are used to compute photoionization and photodissociation rates in the :term:`IGM` and :term:`CGM` by integrating the flux weighted by relevant cross-sections over wavelength. Implementations include the cosmic microwave background, ultraviolet and X-ray ionizing backgrounds, and stellar radiation fields.
   </description>
   <default>null</default>
   <method name="flux">
    <description>
    Return the flux (in units of ergs cm\ :math:`^{-2}` s\ :math:`^{-1}` Hz\ :math:`^{-1}` ster\ :math:`^{-1}`) of the given radiation field.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision          , intent(in   ) :: wavelength</argument>
    <argument>type            (treeNode), intent(inout) :: node</argument>
   </method>
   <method name="integrateOverCrossSection">
    <description>
    Integrates the flux (in units of ergs cm\ :math:`^{-2}` s\ :math:`^{-1}` Hz\ :math:`^{-1}` ster\ :math:`^{-1}`) of the given radiation structure between the wavelengths given in ``wavelengthRange`` over a cross section specified by the function ``crossSectionFunction``.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                              , dimension(2), intent(in   ) :: wavelengthRange</argument>
    <argument>procedure       (crossSectionFunctionTemplate), pointer                     :: crossSectionFunction</argument>
    <argument>type            (treeNode                    )              , intent(inout) :: node</argument>
    <code>
     radiationFieldIntegrateOverCrossSection=radiationFieldIntegrateOverCrossSection_(self,wavelengthRange,crossSectionFunction,node)
    </code>
   </method>
   <method name="time">
     <description>
     Return the time for which the radiation field is currently set.
     </description>
     <type>double precision</type>
     <pass>yes</pass>
   </method>
   <method name="timeSet">
     <description>
     Set the cosmic time (in Gyr) at which the radiation field properties---such as the CMB temperature or the UV background intensity---should be evaluated for subsequent flux queries.
     </description>
     <type>void</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="timeDependentOnly">
     <description>
     Return true if the radiation field depends upon time, but upon no other variables.
     </description>
     <type>logical</type>
     <pass>yes</pass>
   </method>
   <data>type(rateCoefficient), allocatable, dimension(:) :: rateCoefficients</data>
  </functionClass>
  !!]

  ! Module global variables for use in integrand routines.
  class    (radiationFieldClass         ), pointer :: self_
  type     (treeNode                    ), pointer :: node_
  procedure(crossSectionFunctionTemplate), pointer :: crossSectionFunction_
  !$omp threadprivate(self_,node_,crossSectionFunction_)

  abstract interface
     double precision function crossSectionFunctionTemplate(wavelength)
       double precision, intent(in   ) :: wavelength
     end function crossSectionFunctionTemplate
  end interface   
  
contains

  double precision function radiationFieldIntegrateOverCrossSection_(self,wavelengthRange,crossSectionFunction,node)
    !!{RST
    Integrate the photon number of the radiation field over a given cross-section function (which should return the cross section in units of cm\ :math:`^2`), i.e.:

    .. math::

       {4 \pi \over \mathrm{h}} \int_{\lambda_1}^{\lambda_2} \sigma(\lambda) j_{\nu}(\lambda) {\mathrm{d}\lambda \over \lambda},

    where :math:`j_{\nu}` is the flux of energy per unit area per unit solid angle and per unit frequency.
    !!}
    use :: Numerical_Constants_Math    , only : Pi
    use :: Numerical_Constants_Physical, only : plancksConstant
    use :: Numerical_Constants_Units   , only : ergs
    use :: Numerical_Integration       , only : integrator          , GSL_Integ_Gauss15
    use :: Numerical_Ranges            , only : Range_Lattice_Extend, Range_Pinned     , gridSchemePerDecade
    implicit none
    class           (radiationFieldClass         ), target      , intent(inout) :: self
    double precision                              , dimension(2), intent(in   ) :: wavelengthRange
    procedure       (crossSectionFunctionTemplate), pointer                     :: crossSectionFunction
    type            (treeNode                    ), target      , intent(inout) :: node
    type            (integrator                  ), save                        :: integrator_
    type            (rateCoefficient             ), dimension(:), allocatable   :: rateCoefficients
    double precision                              , dimension(:), allocatable   :: time
    logical                                       , dimension(:), allocatable   :: isComputed
    type            (rangeLattice                )                              :: latticeTime
    integer                                       , parameter                   :: countTimesPerDecade  =30
    logical                                                                     :: integratorInitialized=.false., matched             , &
         &                                                                         recompute
    double precision                                                            :: timeCurrent
    integer                                                                     :: i                            , indexRateCoefficient
    !$omp threadprivate(integrator_,integratorInitialized)

    ! Construct the integrator if necessary.
    if (.not.integratorInitialized) then
       ! A low order integrator is used here - this is efficient since typically these integrands can involve cross-sections with
       ! sharp edges, and the radiation field can also have sharp edges (and is often approximated using a tabulated function).
       integrator_          =integrator(                                         &
            &                           integrand        =crossSectionIntegrand, &
            &                           toleranceRelative=1.0d-3               , &
            &                           integrationRule  =GSL_Integ_Gauss15      &
            &                          )
       integratorInitialized=.true.
    end if
    ! Set module-scope pointers to self and the cross-section function for use in the integrand routine.
    self_                 => self
    node_                 => node
    crossSectionFunction_ => crossSectionFunction
    ! If this radiation field is time-dependent only we may be able to use a tabulated solution.
    if (self%timeDependentOnly()) then
       ! Look for an existing solution.
       matched=.false.
       if (allocated(self%rateCoefficients)) then
          do i=1,size(self%rateCoefficients)
             if     (                                                                                  &
                  &   associated(self%rateCoefficients(i)%crossSectionFunction,  crossSectionFunction) &
                  &  .and.                                                                             &
                  &   all       (self%rateCoefficients(i)%wavelengthRange     == wavelengthRange     ) &
                  & ) then
                matched             =.true.
                indexRateCoefficient=i
                exit
             end if
          end do
       end if
       ! Store the current time for the radiation field.
       timeCurrent=self%time()
       ! Determine if we need to (re)compute the rate coefficient.
       if (.not.matched) then
          if (allocated(self%rateCoefficients)) then
             call move_alloc(self%rateCoefficients,rateCoefficients)
             allocate(self%rateCoefficients(size(rateCoefficients)+1))
             do i=1,size(rateCoefficients)
                self%rateCoefficients(i)=rateCoefficients(i)
                call self%rateCoefficients(i)%interpolator_%GSLReallocate()
             end do
             deallocate(rateCoefficients)
          else
             allocate(self%rateCoefficients(                       1))
          end if
          indexRateCoefficient                                             =  size(self%rateCoefficients)
          self%rateCoefficients(indexRateCoefficient)%crossSectionFunction => crossSectionFunction
          self%rateCoefficients(indexRateCoefficient)%wavelengthRange      =  wavelengthRange
       end if
       ! Find the range of times to tabulate, pinning it to an absolute lattice so that the times at which the integral is
       ! evaluated - and therefore every rate coefficient interpolated between them - depend only on which lattice points are
       ! spanned, and not on the sequence of times at which this rate coefficient happened to be asked for. The request is passed
       ! as the target and the range already tabulated is unioned in through `latticeCurrent`; folding the latter into the target
       ! instead - as the `min`/`max` against the current bounds formerly did - would apply the safety margin to an already
       ! margined bound and so ratchet the range outward on every retabulation. The margin is the factor of two either side of
       ! the request that was applied before this tabulation was pinned. The bounds are anchored to half decades rather than to
       ! whole decades because this is a cosmic time axis, on which a whole decade is a large fraction of the history of the
       ! universe, and each additional point costs an integration of the radiation field over the cross section.
       latticeTime=Range_Pinned(                                                                         &
            &                                  [timeCurrent]                                           , &
            &                                   countTimesPerDecade                                    , &
            &                                   gridSchemePerDecade                                    , &
            &                   marginFactor  = 2.0d0                                                  , &
            &                   anchorEvery   = countTimesPerDecade/2                                  , &
            &                   latticeCurrent= self%rateCoefficients(indexRateCoefficient)%latticeTime  &
            &                  )
       ! Decide whether the tabulation must be built or extended. The decision is taken from the pinned lattice rather than from
       ! the bounds directly: the safety margin is applied to the request, so testing the request against the bounds would apply
       ! that margin only when the request happened to arrive outside them. The range reached would then depend on the order in
       ! which times were asked for - which is precisely the dependence that pinning the lattice exists to remove.
       recompute=.not.self%rateCoefficients(indexRateCoefficient)%latticeTime%isDefined()
       if (.not.recompute) recompute=.not.self%rateCoefficients(indexRateCoefficient)%latticeTime%covers(latticeTime)
       ! (Re)compute the table if necessary.
       if (recompute) then
          ! Extend the tabulation onto the new lattice, preserving the rate coefficients already computed - each is an integral
          ! over the radiation field at its own time, which does not change when the range of times tabulated does.
          call Range_Lattice_Extend(                                                              &
               &                    self%rateCoefficients(indexRateCoefficient)%latticeTime     , &
               &                    latticeTime                                                 , &
               &                    self%rateCoefficients(indexRateCoefficient)%rateCoefficient_, &
               &                    isComputed                                                    &
               &                   )
          time=latticeTime%values()
          do i=1,latticeTime%count
             if (isComputed(i)) cycle
             call self%timeSet(time(i))
             self%rateCoefficients(indexRateCoefficient)%rateCoefficient_(i)=+integrator_%integrate(wavelengthRange(1),wavelengthRange(2)) &
                  &                                                          *4.0d0                                                        &
                  &                                                          *Pi                                                           &
                  &                                                          *ergs                                                         &
                  &                                                          /plancksConstant
          end do
          self%rateCoefficients(indexRateCoefficient)%latticeTime  =latticeTime
          self%rateCoefficients(indexRateCoefficient)%interpolator_=interpolator(time,self%rateCoefficients(indexRateCoefficient)%rateCoefficient_)
          ! Restore the time in the radiation field.
          call self%timeSet(timeCurrent)
       end if
       ! Compute the rate coefficient by interpolating in the tables.
       radiationFieldIntegrateOverCrossSection_=self%rateCoefficients(indexRateCoefficient)%interpolator_%interpolate(timeCurrent)
       return
    end if
    ! Perform the integration.
    radiationFieldIntegrateOverCrossSection_=integrator_%integrate(wavelengthRange(1),wavelengthRange(2))    
    ! Scale result by multiplicative prefactors to give answer in units of inverse seconds.
    radiationFieldIntegrateOverCrossSection_=+radiationFieldIntegrateOverCrossSection_ &
         &                                   *4.0d0                                    &
         &                                   *Pi                                       &
         &                                   *ergs                                     &
         &                                   /plancksConstant
    return
  end function radiationFieldIntegrateOverCrossSection_

  double precision function crossSectionIntegrand(wavelength)
    !!{RST
    Integrand function use in integrating a radiation field over a cross section function.
    !!}
    double precision, intent(in   ) :: wavelength

    if (wavelength > 0.0d0) then
       crossSectionIntegrand=+crossSectionFunction_     (wavelength      ) &
            &                *self_                %flux(wavelength,node_) &
            &                /                           wavelength
    else
       crossSectionIntegrand=0.0d0
    end if
    return
  end function crossSectionIntegrand

  subroutine rateCoefficientInterpolatorDeepCopy(self)
    !!{RST
    Perform deep copy actions on the interpolator.
    !!}
    implicit none
    class(rateCoefficient), intent(inout) :: self

    call self%interpolator_%GSLReallocate()
    return
  end subroutine rateCoefficientInterpolatorDeepCopy

end module Radiation_Fields
