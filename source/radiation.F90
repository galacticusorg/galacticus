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
Contains a module which implements a class to describe radiation fields.
!!}

module Radiation_Fields
  !!{
  Implements a class to describe radiation fields.
  !!}
  use :: Galacticus_Nodes       , only : treeNode
  use :: Numerical_Interpolation, only : interpolator
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
     !!{
     Type used to store tables of rate coefficients.
     !!}
     private
     procedure       (crossSectionFunctionTemplate), pointer     , nopass :: crossSectionFunction => null(     )
     double precision                  , dimension(2)         :: wavelengthRange
     double precision                                         :: timeMinimum          =  huge(0.0d0), timeMaximum=-huge(0.0d0)
     type            (interpolator    ), allocatable          :: interpolator_
   contains
     !![
     <methods>
       <method description="Perform deep copy actions on interpolators." method="interpolatorDeepCopy" />
     </methods>
     !!]
     procedure :: interpolatorDeepCopy => rateCoefficientInterpolatorDeepCopy
  end type rateCoefficient
  
  !![
  <functionClass>
   <name>radiationField</name>
   <descriptiveName>Radiation Fields</descriptiveName>
   <description>Class providing radiation fields.</description>
   <default>null</default>
   <data>type(rateCoefficient), allocatable, dimension(:) :: rateCoefficients</data>
   <method name="flux">
    <description>Return the flux (in units of ergs cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$) of the given radiation field.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision          , intent(in   ) :: wavelength</argument>
    <argument>type            (treeNode), intent(inout) :: node</argument>
   </method>
   <method name="integrateOverCrossSection">
    <description>Integrates the flux (in units of ergs cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$) of the given radiation structure between the wavelengths given in {\normalfont \ttfamily wavelengthRange} over a cross section specified by the function {\normalfont \ttfamily crossSectionFunction}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision                  , dimension(2), intent(in   ) :: wavelengthRange</argument>
    <argument>procedure       (crossSectionFunctionTemplate), pointer                     :: crossSectionFunction</argument>
    <argument>type            (treeNode        )              , intent(inout) :: node</argument>
    <code>
     radiationFieldIntegrateOverCrossSection=radiationFieldIntegrateOverCrossSection_(self,wavelengthRange,crossSectionFunction,node)
    </code>
   </method>
   <method name="time">
     <description>Return the time for which the radiation field is currently set.</description>
     <type>double precision</type>
     <pass>yes</pass>
   </method>
   <method name="timeSet">
     <description>Set the time of the radiation field.</description>
     <type>void</type>
     <pass>yes</pass>
     <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="timeDependentOnly">
     <description>Return true if the radiation field depends upon time, but upon no other variables.</description>
     <type>logical</type>
     <pass>yes</pass>
   </method>
  </functionClass>
  !!]

  ! Module global variables for use in integrand routines.
  class    (radiationFieldClass), pointer :: self_
  type     (treeNode           ), pointer :: node_
  procedure(crossSectionFunctionTemplate   ), pointer :: crossSectionFunction_
  !$omp threadprivate(self_,node_,crossSectionFunction_)

  abstract interface
     double precision function crossSectionFunctionTemplate(wavelength)
       double precision, intent(in   ) :: wavelength
     end function crossSectionFunctionTemplate
  end interface   
  
contains

  double precision function radiationFieldIntegrateOverCrossSection_(self,wavelengthRange,crossSectionFunction,node)
    !!{
    Integrate the photon number of the radiation field over a given cross-section function (which should return the cross
    section in units of cm$^2$), i.e.:
    \begin{equation}
    {4 \pi \over \mathrm{h}} \int_{\lambda_1}^{\lambda_2} \sigma(\lambda) j_{\nu}(\lambda) {\mathrm{d}\lambda \over \lambda},
    \end{equation}
    where $j_{\nu}$ is the flux of energy per unit area per unit solid angle and per unit frequency.
    !!}
    use :: Numerical_Constants_Math    , only : Pi
    use :: Numerical_Constants_Physical, only : plancksConstant
    use :: Numerical_Constants_Units   , only : ergs
    use :: Numerical_Integration       , only : integrator     , GSL_Integ_Gauss15
    use :: Numerical_Ranges            , only : Make_Range     , rangeTypeLogarithmic
    implicit none
    class           (radiationFieldClass), target      , intent(inout) :: self
    double precision                     , dimension(2), intent(in   ) :: wavelengthRange
    procedure       (crossSectionFunctionTemplate   ), pointer                     :: crossSectionFunction
    type            (treeNode           ), target      , intent(inout) :: node
    type            (integrator         ), save                        :: integrator_
    type            (rateCoefficient    ), dimension(:), allocatable   :: rateCoefficients
    double precision                     , dimension(:), allocatable   :: time                         , rateCoefficient_
    double precision                     , parameter                   :: countTimesPerDecade  =30.0d0
    logical                                                            :: integratorInitialized=.false., matched         , &
         &                                                                recompute
    double precision                                                   :: timeCurrent
    integer                                                            :: countTimes                   , i               , &
         &                                                                indexRateCoefficient
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
       if (matched) then
          recompute= timeCurrent < self%rateCoefficients(indexRateCoefficient)%timeMinimum &
               &    .or.                                                                   &
               &     timeCurrent > self%rateCoefficients(indexRateCoefficient)%timeMaximum
       else
          recompute=.true.
          if (allocated(self%rateCoefficients)) then
             call move_alloc(self%rateCoefficients,rateCoefficients)
             allocate(self%rateCoefficients(size(rateCoefficients)+1))
             self%rateCoefficients(1:size(rateCoefficients))=rateCoefficients
             do i=1,size(rateCoefficients)
                call self%rateCoefficients(i)%interpolator_%GSLReallocate(gslFree=.false.)
             end do
             deallocate(rateCoefficients)
          else
             allocate(self%rateCoefficients(                       1))
          end if
          indexRateCoefficient                                             =  size(self%rateCoefficients)
          self%rateCoefficients(indexRateCoefficient)%crossSectionFunction => crossSectionFunction
          self%rateCoefficients(indexRateCoefficient)%wavelengthRange      =  wavelengthRange
       end if
       ! (Re)compute the table if necessary.
       if (recompute) then
          if (allocated(self%rateCoefficients(indexRateCoefficient)%interpolator_)) deallocate(self%rateCoefficients(indexRateCoefficient)%interpolator_)
          allocate(self%rateCoefficients(indexRateCoefficient)%interpolator_)
          self%rateCoefficients(indexRateCoefficient)%timeMinimum=min(self%rateCoefficients(indexRateCoefficient)%timeMinimum,timeCurrent/2.0d0)
          self%rateCoefficients(indexRateCoefficient)%timeMaximum=max(self%rateCoefficients(indexRateCoefficient)%timeMaximum,timeCurrent*2.0d0)
          countTimes=int(log10(self%rateCoefficients(indexRateCoefficient)%timeMaximum/self%rateCoefficients(indexRateCoefficient)%timeMinimum)*countTimesPerDecade+1.0d0)
          allocate(time            (countTimes))
          allocate(rateCoefficient_(countTimes))
          time=Make_Range(self%rateCoefficients(indexRateCoefficient)%timeMinimum,self%rateCoefficients(indexRateCoefficient)%timeMaximum,countTimes,rangeType=rangeTypeLogarithmic)
          do i=1,countTimes
             call self%timeSet(time(i))
             rateCoefficient_(i)=+integrator_%integrate(wavelengthRange(1),wavelengthRange(2)) &
                  &              *4.0d0                                                        &
                  &              *Pi                                                           &
                  &              *ergs                                                         &
                  &              /plancksConstant
          end do
          self%rateCoefficients(indexRateCoefficient)%interpolator_=interpolator(time,rateCoefficient_)
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
    !!{
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
    !!{
    Perform deep copy actions on the interpolator.
    !!}
    implicit none
    class(rateCoefficient), intent(inout) :: self

    call self%interpolator_%GSLReallocate(gslFree=.false.)
    return
  end subroutine rateCoefficientInterpolatorDeepCopy

end module Radiation_Fields
