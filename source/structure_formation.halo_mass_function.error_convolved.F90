!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements a dark matter halo mass function class which modifies another mass function by convolving
with a mass-dependent error.
!!}

  use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassError, nbodyHaloMassErrorClass

  !![
  <haloMassFunction name="haloMassFunctionErrorConvolved">
   <description>
    The halo mass function is computed by convolving another halo mass function with a mass dependent error. Specifically, the
    mass function is convolved with a Gaussian random error distribution with width computed using the given {\normalfont
    \ttfamily nbodyHaloMassError} object.
   </description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionErrorConvolved
     !!{
     A halo mass function class convolves another halo mass function with a mass dependent error.
     !!}
     private
     double precision                                   :: errorFractionalMaximum
     class           (haloMassFunctionClass  ), pointer :: massFunctionIntrinsic  => null()
     class           (nbodyHaloMassErrorClass), pointer :: nbodyHaloMassError_    => null()
   contains
     final     ::                 errorConvolvedDestructor
     procedure :: differential => errorConvolvedDifferential
  end type haloMassFunctionErrorConvolved

  interface haloMassFunctionErrorConvolved
     !!{
     Constructors for the {\normalfont \ttfamily errorConvolved} halo mass function class.
     !!}
     module procedure errorConvolvedConstructorParameters
     module procedure errorConvolvedConstructorInternal
  end interface haloMassFunctionErrorConvolved

  ! Mass error scale used in convolution integral.
  class           (haloMassFunctionClass), pointer :: errorConvolvedIntrinsicMassFunction
  double precision                                 :: errorConvolvedErrorMass            , errorConvolvedTime, &
       &                                              errorConvolvedMass
  !$omp threadprivate(errorConvolvedErrorMass,errorConvolvedTime,errorConvolvedMass, errorConvolvedIntrinsicMassFunction)

contains

  function errorConvolvedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily errorConvolved} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloMassFunctionErrorConvolved)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (haloMassFunctionClass         ), pointer       :: massFunctionIntrinsic
    class           (cosmologyParametersClass      ), pointer       :: cosmologyParameters_
    class           (nbodyHaloMassErrorClass       ), pointer       :: nbodyHaloMassError_
    double precision                                                :: errorFractionalMaximum

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>errorFractionalMaximum</name>
      <source>parameters</source>
      <description>Maximum allowed fractional error in halo mass.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_"  source="parameters"/>
    <objectBuilder class="nbodyHaloMassError"  name="nbodyHaloMassError_"   source="parameters"/>
    <objectBuilder class="haloMassFunction"    name="massFunctionIntrinsic" source="parameters"/>
    !!]
    self=haloMassFunctionErrorConvolved(massFunctionIntrinsic,cosmologyParameters_,nbodyHaloMassError_,errorFractionalMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_" />
    <objectDestructor name="nbodyHaloMassError_"  />
    <objectDestructor name="massFunctionIntrinsic"/>
    !!]
    return
  end function errorConvolvedConstructorParameters

  function errorConvolvedConstructorInternal(massFunctionIntrinsic,cosmologyParameters_,nbodyHaloMassError_,errorFractionalMaximum) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily errorConvolved} halo mass function class.
    !!}
    implicit none
    type            (haloMassFunctionErrorConvolved)                        :: self
    class           (haloMassFunctionClass         ), target, intent(in   ) :: massFunctionIntrinsic
    class           (cosmologyParametersClass      ), target, intent(in   ) :: cosmologyParameters_
    class           (nbodyHaloMassErrorClass       ), target, intent(in   ) :: nbodyHaloMassError_
    double precision                                        , intent(in   ) :: errorFractionalMaximum
    !![
    <constructorAssign variables="*massFunctionIntrinsic, *cosmologyParameters_, *nbodyHaloMassError_, errorFractionalMaximum"/>
    !!]

    return
  end function errorConvolvedConstructorInternal

  subroutine errorConvolvedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily errorConvolved} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionErrorConvolved), intent(inout) :: self

    !![
    <objectDestructor name="self%massFunctionIntrinsic" />
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%nbodyHaloMassError_"   />
    !!]
    return
  end subroutine errorConvolvedDestructor

  double precision function errorConvolvedDifferential(self,time,mass,node)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    use :: Galacticus_Nodes     , only : nodeComponentBasic, treeNode
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (haloMassFunctionErrorConvolved), intent(inout)           :: self
    double precision                                , intent(in   )           :: time                    , mass
    type            (treeNode                      ), intent(inout), optional :: node
    double precision                                , parameter               :: rangeIntegralSigma=5.0d0
    type            (treeNode                      ), pointer                 :: nodeWork
    class           (nodeComponentBasic            ), pointer                 :: basic
    double precision                                                          :: massLow                 , massHigh
    type            (integrator                    )                          :: integrator_

     ! Create a work node.
    nodeWork => treeNode      (                 )
    basic    => nodeWork%basic(autoCreate=.true.)
    ! Set the properties of the work node.
    call basic%massSet(mass)
    call basic%timeSet(time)
    ! Get the error at this mass scale.
    errorConvolvedErrorMass=+mass                                                       &
         &                  *min(                                                       &
         &                       self%nbodyHaloMassError_   %errorFractional(nodeWork), &
         &                       self%errorFractionalMaximum                            &
         &                      )
    ! Perform the convolution integral.
    massLow                             =  max(1.0d-3*mass,mass-rangeIntegralSigma*errorConvolvedErrorMass)
    massHigh                            =                  mass+rangeIntegralSigma*errorConvolvedErrorMass
    errorConvolvedIntrinsicMassFunction => self%massFunctionIntrinsic
    errorConvolvedTime                  =  time
    errorConvolvedMass                  =  mass
    integrator_                         =  integrator           (                                             &
         &                                                                         errorConvolvedConvolution, &
         &                                                       toleranceAbsolute=1.0d-100                 , &
         &                                                       toleranceRelative=1.0d-006                   &
         &                                                      )
    errorConvolvedDifferential          =  integrator_%integrate(                                             &
         &                                                       massLow                                    , &
         &                                                       massHigh                                     &
         &                                                      )
    ! Clean up our work node.
    call nodeWork%destroy()
    deallocate(nodeWork)
    return

  contains

    double precision function errorConvolvedConvolution(massPrime)
      !!{
      Integrand function used in convolving the dark matter halo mass function.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: massPrime

      ! Get the error at this mass scale.
      call basic%massSet(massPrime)
      errorConvolvedErrorMass=+massPrime                                                  &
           &                  *min(                                                       &
           &                       self%nbodyHaloMassError_   %errorFractional(nodeWork), &
           &                       self%errorFractionalMaximum                            &
           &                      )
      ! Return the convolution integrand.
      errorConvolvedConvolution=+errorConvolvedIntrinsicMassFunction%differential(errorConvolvedTime,massPrime,node=node) &
           &                    *exp(                                                                                     &
           &                         -0.5d0                                                                               &
           &                         *(                                                                                   &
           &                           +(                                                                                 &
           &                             +massPrime                                                                       &
           &                             -errorConvolvedMass                                                              &
           &                            )                                                                                 &
           &                           /errorConvolvedErrorMass                                                           &
           &                          )**2                                                                                &
           &                        )                                                                                     &
           &                    /sqrt(                                                                                    &
           &                          +2.0d0                                                                              &
           &                          *Pi                                                                                 &
           &                         )                                                                                    &
           &                    /errorConvolvedErrorMass
      return
    end function errorConvolvedConvolution

  end function errorConvolvedDifferential
