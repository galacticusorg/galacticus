!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements a transfer function class based on the thermal \gls{wdm} modifier of \cite{bode_halo_2001}.

  use :: Cosmology_Functions  , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters , only : cosmologyParametersClass
  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !# <transferFunction name="transferFunctionBode2001">
  !#  <description>Provides a transfer function based on the thermal \gls{wdm} modifier of \cite{bode_halo_2001}.</description>
  !# </transferFunction>
  type, extends(transferFunctionClass) :: transferFunctionBode2001
     !% A transfer function class which modifies another transfer function using the thermal \gls{wdm} modifier of \cite{bode_halo_2001}.
     private
     double precision                                    :: epsilon                       , eta        , &
          &                                                 nu                            , scaleCutOff, &
          &                                                 time                          , redshift
     class           (transferFunctionClass   ), pointer :: transferFunctionCDM  => null()
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     class           (darkMatterParticleClass ), pointer :: darkMatterParticle_  => null()
   contains
     final     ::                          bode2001Destructor
     procedure :: value                 => bode2001Value
     procedure :: logarithmicDerivative => bode2001LogarithmicDerivative
     procedure :: halfModeMass          => bode2001HalfModeMass
     procedure :: epochTime             => bode2001EpochTime
  end type transferFunctionBode2001

  interface transferFunctionBode2001
     !% Constructors for the {\normalfont \ttfamily bode2001} transfer function class.
     module procedure bode2001ConstructorParameters
     module procedure bode2001ConstructorInternal
  end interface transferFunctionBode2001

contains

  function bode2001ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily bode2001} transfer function class which takes a parameter set as input.
    use :: Cosmology_Functions           , only : cosmologyFunctions        , cosmologyFunctionsClass
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    use :: Galacticus_Error              , only : Galacticus_Error_Report
    use :: Input_Parameters              , only : inputParameter            , inputParameters
    implicit none
    type            (transferFunctionBode2001)                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (transferFunctionClass   ), pointer       :: transferFunctionCDM
    class           (cosmologyParametersClass), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    class           (darkMatterParticleClass ), pointer       :: darkMatterParticle_
    double precision                                          :: epsilon             , eta     , &
         &                                                       nu                  , redshift

    ! Validate parameters.
    if (.not.parameters%isPresent('transferFunctionMethod')) call Galacticus_Error_Report("an explicit 'transferFunctionMethod' must be given"//{introspection:location})
    ! Read parameters.
    !# <inputParameter>
    !#   <name>epsilon</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.359d0</defaultValue>
    !#   <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
    !#   <description>The parameter $\epsilon$ appearing in the warm dark matter transfer function \citep{barkana_constraints_2001}.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>eta</name>
    !#   <source>parameters</source>
    !#   <defaultValue>3.81d0</defaultValue>
    !#   <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
    !#   <description>The parameter $\epsilon$ appearing in the warm dark matter transfer function \citep{barkana_constraints_2001}.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>nu</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.1d0</defaultValue>
    !#   <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
    !#   <description>The parameter $\epsilon$ appearing in the warm dark matter transfer function \citep{barkana_constraints_2001}.</description>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !# <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    !# <objectBuilder class="transferFunction"    name="transferFunctionCDM"  source="parameters"/>
    !# <inputParameter>
    !#   <name>redshift</name>
    !#   <source>parameters</source>
    !#   <defaultValue>cosmologyFunctions_%redshiftFromExpansionFactor(cosmologyFunctions_%equalityEpochMatterRadiation(requestTypeExpansionFactor))</defaultValue>
    !#   <description>The redshift of the epoch at which the transfer function is defined.</description>
    !# </inputParameter>
     self=transferFunctionBode2001(transferFunctionCDM,epsilon,eta,nu,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),cosmologyParameters_,darkMatterParticle_,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    !# <objectDestructor name="cosmologyFunctions_" />
    !# <objectDestructor name="darkMatterParticle_" />
    !# <objectDestructor name="transferFunctionCDM" />
    return
  end function bode2001ConstructorParameters

  function bode2001ConstructorInternal(transferFunctionCDM,epsilon,eta,nu,time,cosmologyParameters_,darkMatterParticle_,cosmologyFunctions_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily bode2001} transfer function class.
    use :: Cosmology_Parameters , only : hubbleUnitsLittleH
    use :: Dark_Matter_Particles, only : darkMatterParticleWDMThermal
    use :: Galacticus_Error     , only : Galacticus_Error_Report
    implicit none
    type            (transferFunctionBode2001)                        :: self
    class           (transferFunctionClass   ), target, intent(in   ) :: transferFunctionCDM
    double precision                                  , intent(in   ) :: epsilon                   , eta                            , &
         &                                                               nu                        , time
    class           (cosmologyParametersClass), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), target, intent(in   ) :: cosmologyFunctions_
    class           (darkMatterParticleClass ), target, intent(in   ) :: darkMatterParticle_
    double precision                          , parameter             :: massReference       =1.0d0, degreesOfFreedomReference=1.5d0
    !# <constructorAssign variables="*transferFunctionCDM, epsilon, eta, nu, time, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterParticle_"/>

    self%redshift=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    ! Compute the comoving cut-off scale. This uses equation (4) from Barkana et al. (2001;
    ! http://adsabs.harvard.edu/abs/2001ApJ...558..482B), with the prefactor of 0.932 to give the cut-off scale at the epoch of
    ! matter-radiation equality as discussed in the paragraph following their equation (4).
    select type (particle => self%darkMatterParticle_)
    class is (darkMatterParticleWDMThermal)
       self%scaleCutOff=+0.932d0                                                                  &
            &           *0.201d0                                                                  &
            &           *(                                                                        &
            &             +(                                                                      &
            &               +self%cosmologyParameters_%OmegaMatter   (                  )         &
            &               -self%cosmologyParameters_%OmegaBaryon   (                  )         &
            &              )                                                                      &
            &             *  self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2      &
            &             /0.15d0                                                                 &
            &            )                                                               **0.15d0 &
            &           /(particle%degreesOfFreedomEffective()/degreesOfFreedomReference)**0.29d0 &
            &           /(particle%mass                     ()/            massReference)**1.15d0
    class default
       call Galacticus_Error_Report('transfer function expects a thermal warm dark matter particle'//{introspection:location})
    end select
    return
  end function bode2001ConstructorInternal

  subroutine bode2001Destructor(self)
    !% Destructor for the {\normalfont \ttfamily bode2001} transfer function class.
    implicit none
    type(transferFunctionBode2001), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%darkMatterParticle_" />
    !# <objectDestructor name="self%transferFunctionCDM" />
    return
  end subroutine bode2001Destructor

  double precision function bode2001Value(self,wavenumber)
    !% Return the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionBode2001), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber

    bode2001Value=+self%transferFunctionCDM%value(wavenumber)
    if (self%scaleCutOff > 0.0d0)                    &
         & bode2001Value=+bode2001Value              &
         &               /(                          &
         &                  1.0d0                    &
         &                 +                         &
         &                  (                        &
         &                   +self%epsilon           &
         &                   *wavenumber             &
         &                   *self%scaleCutOff       &
         &                  )**(2.0d0   *self%nu)    &
         &                )  **(self%eta/self%nu)
    return
  end function bode2001Value

  double precision function bode2001LogarithmicDerivative(self,wavenumber)
    !% Return the logarithmic derivative of the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionBode2001), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber

    bode2001LogarithmicDerivative=+self%transferFunctionCDM%logarithmicDerivative(wavenumber)
    if (self%scaleCutOff > 0.0d0)                                       &
         & bode2001LogarithmicDerivative=+bode2001LogarithmicDerivative &
         &                               +2.0d0                         &
         &                               *self%eta                      &
         &                               *(                             &
         &                                 +1.0d0                       &
         &                                 /(                           &
         &                                   +1.0d0                     &
         &                                   +(                         &
         &                                     +self%epsilon            &
         &                                     *wavenumber              &
         &                                     *self%scaleCutOff        &
         &                                    )**(2.0d0*self%nu)        &
         &                                  )                           &
         &                                 -1.0d0                       &
         &                                )
    return
  end function bode2001LogarithmicDerivative

  double precision function bode2001HalfModeMass(self,status)
    !% Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    !% to a \gls{cdm} transfer function.
    use :: Galacticus_Error        , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionBode2001), intent(inout)           :: self
    integer                                   , intent(  out), optional :: status
    double precision                          , parameter               :: wavenumberHalfModeScaleFree=sqrt(0.25d0+2.0d0*log(2.0d0))-0.5d0
    double precision                                                    :: matterDensity                                                  , wavenumberHalfMode

    matterDensity       =+self%cosmologyParameters_%OmegaMatter    () &
         &               *self%cosmologyParameters_%densityCritical()
    wavenumberHalfMode  =+(                            &
         &                 +2.0d0**(+self%nu/self%eta) &
         &                 -1.0d0                      &
         &                )      **(+0.5d0  /self%nu ) &
         &                /self%epsilon                &
         &                /self%scaleCutOff
    bode2001HalfModeMass=+4.0d0                &
         &               *Pi                   &
         &               /3.0d0                &
         &               *matterDensity        &
         &               *(                    &
         &                 +Pi                 &
         &                 /wavenumberHalfMode &
         &               )**3
    if (present(status)) status=errorStatusSuccess
    return
  end function bode2001HalfModeMass

  double precision function bode2001EpochTime(self)
    !% Return the cosmic time at the epoch at which this transfer function is defined.
    implicit none
    class(transferFunctionBode2001), intent(inout) :: self

    bode2001EpochTime=self%time
    return
  end function bode2001EpochTime
