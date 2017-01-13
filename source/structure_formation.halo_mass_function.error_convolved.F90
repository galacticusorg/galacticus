!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements a dark matter halo mass function class which modifies another mass function by convolving
!% with a mass-dependent error.

  use Statistics_NBody_Halo_Mass_Errors
  
  !# <haloMassFunction name="haloMassFunctionErrorConvolved">
  !#  <description>
  !#   The halo mass function is computed by convolving another halo mass function with a mass dependent error. Specifically, the
  !#   mass function is convolved with a Gaussian random error distribution with width computed using the given {\normalfont
  !#   \ttfamily nbodyHaloMassError} object.  
  !#  </description>
  !# </haloMassFunction>
  type, extends(haloMassFunctionClass) :: haloMassFunctionErrorConvolved
     !% A halo mass function class convolves another halo mass function with a mass dependent error. 
     private
     double precision                                   :: errorFractionalMaximum
     class           (haloMassFunctionClass  ), pointer :: massFunctionIntrinsic
     class           (nbodyHaloMassErrorClass), pointer :: nbodyHaloMassError_
   contains
     final     ::                 errorConvolvedDestructor
     procedure :: differential => errorConvolvedDifferential
  end type haloMassFunctionErrorConvolved

  interface haloMassFunctionErrorConvolved
     !% Constructors for the {\normalfont \ttfamily errorConvolved} halo mass function class.
     module procedure errorConvolvedConstructorParameters
     module procedure errorConvolvedConstructorInternal
  end interface haloMassFunctionErrorConvolved

  ! Mass error scale used in convolution integral.
  class           (haloMassFunctionClass), pointer :: errorConvolvedIntrinsicMassFunction
  double precision                                 :: errorConvolvedErrorMass            , errorConvolvedTime, &
       &                                              errorConvolvedMass
  !$omp threadprivate(errorConvolvedErrorMass,errorConvolvedTime,errorConvolvedMass, errorConvolvedIntrinsicMassFunction)
  
contains

  function errorConvolvedConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily errorConvolved} halo mass function class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(haloMassFunctionErrorConvolved)                :: errorConvolvedConstructorParameters
    type(inputParameters               ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>errorFractionalMaximum</name>
    !#   <source>parameters</source>
    !#   <variable>errorConvolvedConstructorParameters%errorFractionalMaximum</variable>
    !#   <description>Maximum allowed fractional error in halo mass.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="errorConvolvedConstructorParameters%cosmologyParameters_"  source="parameters"/>
    !# <objectBuilder class="nbodyHaloMassError"  name="errorConvolvedConstructorParameters%nBodyHaloMassError_"   source="parameters"/>
    !# <objectBuilder class="haloMassFunction"    name="errorConvolvedConstructorParameters%massFunctionIntrinsic" source="parameters"/>
    return
  end function errorConvolvedConstructorParameters

  function errorConvolvedConstructorInternal(massFunctionIntrinsic,cosmologyParameters_,nBodyHaloMassError_,errorFractionalMaximum)
    !% Internal constructor for the {\normalfont \ttfamily errorConvolved} halo mass function class.
    implicit none
    type            (haloMassFunctionErrorConvolved)                        :: errorConvolvedConstructorInternal
    class           (haloMassFunctionClass         ), target, intent(in   ) :: massFunctionIntrinsic
    class           (cosmologyParametersClass      ), target, intent(in   ) :: cosmologyParameters_
    class           (nBodyHaloMassErrorClass       ), target, intent(in   ) :: nBodyHaloMassError_
    double precision                                        , intent(in   ) :: errorFractionalMaximum
    
    errorConvolvedConstructorInternal%massFunctionIntrinsic  => massFunctionIntrinsic
    errorConvolvedConstructorInternal%cosmologyParameters_   => cosmologyParameters_
    errorConvolvedConstructorInternal%nBodyHaloMassError_    => nBodyHaloMassError_
    errorConvolvedConstructorInternal%errorFractionalMaximum =  errorFractionalMaximum
    return
  end function errorConvolvedConstructorInternal
  
  subroutine errorConvolvedDestructor(self)
    !% Destructor for the {\normalfont \ttfamily errorConvolved} halo mass function class.
    implicit none
    type(haloMassFunctionErrorConvolved), intent(inout) :: self

    !# <objectDestructor name="self%massFunctionIntrinsic" />
    !# <objectDestructor name="self%cosmologyParameters_"  />
    !# <objectDestructor name="self%nBodyHaloMassError_"   />
    return
  end subroutine errorConvolvedDestructor

  double precision function errorConvolvedDifferential(self,time,mass)
    !% Return the differential halo mass function at the given time and mass.
    use, intrinsic :: ISO_C_Binding
    use               Numerical_Integration
    use               Galacticus_Nodes
    implicit none
    class           (haloMassFunctionErrorConvolved), intent(inout) :: self
    double precision                                , intent(in   ) :: time                    , mass
    double precision                                , parameter     :: rangeIntegralSigma=5.0d0
    type            (treeNode                      ), pointer       :: node
    class           (nodeComponentBasic            ), pointer       :: basic
    double precision                                                :: massLow                 , massHigh
    type            (fgsl_function                 )                :: integrandFunction
    type            (fgsl_integration_workspace    )                :: integrationWorkspace

     ! Create a work node.
    node  => treeNode      (                 )
    basic => node    %basic(autoCreate=.true.)
    ! Set the properties of the work node.
    call basic%massSet(mass)
    call basic%timeSet(time)
    ! Get the error at this mass scale.
    errorConvolvedErrorMass=+mass                                                   &
         &                  *min(                                                   &
         &                       self%nBodyHaloMassError_   %errorFractional(node), &
         &                       self%errorFractionalMaximum                        &
         &                      )
    ! Clean up our work node.
    call node%destroy()
    deallocate(node)
    ! Perform the convolution integral.
    massLow                             =  max(0.01d0*mass,mass-rangeIntegralSigma*errorConvolvedErrorMass)
    massHigh                            =                  mass+rangeIntegralSigma*errorConvolvedErrorMass
    errorConvolvedIntrinsicMassFunction => self%massFunctionIntrinsic
    errorConvolvedTime                  =  time
    errorConvolvedMass                  =  mass
    errorConvolvedDifferential          =  Integrate(                               &
         &                                           massLow                      , &
         &                                           massHigh                     , &
         &                                           errorConvolvedConvolution    , &
         &                                           integrandFunction            , &
         &                                           integrationWorkspace         , &
         &                                           toleranceAbsolute   =1.0d-100, &
         &                                           toleranceRelative   =1.0d-006  &
         &                                          )
    call Integrate_Done(integrandFunction,integrationWorkspace)    
    return
  end function errorConvolvedDifferential
  
  double precision function errorConvolvedConvolution(massPrime)
    !% Integrand function used in convolving the dark matter halo mass function.
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in   ) :: massPrime

    ! Return the convolution integrand.
    errorConvolvedConvolution=+errorConvolvedIntrinsicMassFunction%differential(errorConvolvedTime,massPrime) &
         &                    *exp(                                                                           &
         &                         -0.5d0                                                                     &
         &                         *(                                                                         &
         &                           +(                                                                       &
         &                             +massPrime                                                             &
         &                             -errorConvolvedMass                                                    &
         &                            )                                                                       &
         &                           /errorConvolvedErrorMass                                                 &
         &                          )**2                                                                      &
         &                        )                                                                           &
         &                    /sqrt(                                                                          &
         &                          +2.0d0                                                                    &
         &                          *Pi                                                                       &
         &                         )                                                                          &
         &                    /errorConvolvedErrorMass
    return
  end function errorConvolvedConvolution
