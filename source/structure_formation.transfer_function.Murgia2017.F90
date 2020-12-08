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

  !# <transferFunction name="transferFunctionMurgia2017">
  !#  <description>Provides a transfer function based on the thermal \gls{wdm} modifier of \cite{bode_halo_2001}.</description>
  !# </transferFunction>
  type, extends(transferFunctionClass) :: transferFunctionMurgia2017
     !% A transfer function class which modifies another transfer function using the thermal \gls{wdm} modifier of \cite{bode_halo_2001}.
     private
     double precision                                    :: n_alpha                         , n_beta        , &
          &                                                 n_gamma                         ,  &
          &                                                 time                          , redshift
     class           (transferFunctionClass   ), pointer :: transferFunctionCDM  => null()
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     class           (darkMatterParticleClass ), pointer :: darkMatterParticle_  => null()
   contains
     final     ::                          Murgia2017Destructor
     procedure :: value                 => Murgia2017Value
     procedure :: logarithmicDerivative => Murgia2017LogarithmicDerivative
     procedure :: halfModeMass          => Murgia2017HalfModeMass
     procedure :: epochTime             => Murgia2017EpochTime
  end type transferFunctionMurgia2017

  interface transferFunctionMurgia2017
     !% Constructors for the {\normalfont \ttfamily bode2001} transfer function class.
     module procedure Murgia2017ConstructorParameters
     module procedure Murgia2017ConstructorInternal
  end interface transferFunctionMurgia2017

contains

  function Murgia2017ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily bode2001} transfer function class which takes a parameter set as input.
    use :: Cosmology_Functions           , only : cosmologyFunctions        , cosmologyFunctionsClass
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    use :: Galacticus_Error              , only : Galacticus_Error_Report
    use :: Input_Parameters              , only : inputParameter            , inputParameters
    implicit none
    type            (transferFunctionMurgia2017)                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (transferFunctionClass   ), pointer       :: transferFunctionCDM
    class           (cosmologyParametersClass), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    class           (darkMatterParticleClass ), pointer       :: darkMatterParticle_
    double precision                                          :: n_alpha             , n_beta     , &
         &                                                       n_gamma                  , redshift

    ! Validate parameters.
    if (.not.parameters%isPresent('transferFunctionMethod')) call Galacticus_Error_Report("an explicit 'transferFunctionMethod' must be given"//{introspection:location})
    ! Read parameters.
    !# <inputParameter>
    !#   <name>n_alpha</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.0075d0</defaultValue>
    !#   <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
    !#   <description>The parameter $\epsilon$ appearing in the warm dark matter transfer function \citep{barkana_constraints_2001}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>n_beta</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.5d0</defaultValue>
    !#   <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
    !#   <description>The parameter $\epsilon$ appearing in the warm dark matter transfer function \citep{barkana_constraints_2001}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>n_gamma</name>
    !#   <source>parameters</source>
    !#   <defaultValue>-10d0</defaultValue>
    !#   <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
    !#   <description>The parameter $\epsilon$ appearing in the warm dark matter transfer function \citep{barkana_constraints_2001}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
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
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
     self=transferFunctionMurgia2017(transferFunctionCDM,n_alpha,n_beta,n_gamma,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),cosmologyParameters_,darkMatterParticle_,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    !# <objectDestructor name="cosmologyFunctions_" />
    !# <objectDestructor name="darkMatterParticle_" />
    !# <objectDestructor name="transferFunctionCDM" />
    return
  end function Murgia2017ConstructorParameters

  function Murgia2017ConstructorInternal(transferFunctionCDM,n_alpha,n_beta,n_gamma,time,cosmologyParameters_,darkMatterParticle_,cosmologyFunctions_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily bode2001} transfer function class.
    use :: Cosmology_Parameters , only : hubbleUnitsLittleH
    use :: Dark_Matter_Particles, only : darkMatterParticleWDMThermal
    use :: Galacticus_Error     , only : Galacticus_Error_Report
    implicit none
    type            (transferFunctionMurgia2017)                        :: self
    class           (transferFunctionClass   ), target, intent(in   ) :: transferFunctionCDM
    double precision                                  , intent(in   ) :: n_alpha                   , n_beta                            , &
         &                                                               n_gamma                        , time
    class           (cosmologyParametersClass), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), target, intent(in   ) :: cosmologyFunctions_
    class           (darkMatterParticleClass ), target, intent(in   ) :: darkMatterParticle_
    double precision                          , parameter             :: massReference       =1.0d0, degreesOfFreedomReference=1.5d0
    !# <constructorAssign variables="*transferFunctionCDM, n_alpha, n_beta, n_gamma, time, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterParticle_"/>

    self%redshift=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    
    return
  end function Murgia2017ConstructorInternal

  subroutine Murgia2017Destructor(self)
    !% Destructor for the {\normalfont \ttfamily bode2001} transfer function class.
    implicit none
    type(transferFunctionMurgia2017), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%darkMatterParticle_" />
    !# <objectDestructor name="self%transferFunctionCDM" />
    return
  end subroutine Murgia2017Destructor

  double precision function Murgia2017Value(self,wavenumber)
    !% Return the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionMurgia2017), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber

    Murgia2017Value=+self%transferFunctionCDM%value(wavenumber)
    if (self%n_alpha > 0.0d0)                    &
         & Murgia2017Value=+Murgia2017Value              &
         &               *(                          &
         &                  1.0d0                    &
         &                 +                         &
         &                  (                        &
         &                   +self%n_alpha             &
         &                   *wavenumber             &
         &                  )**(self%n_beta)           &
         &                )  **(self%n_gamma)
    return
  end function Murgia2017Value

  double precision function Murgia2017LogarithmicDerivative(self,wavenumber)
    !% Return the logarithmic derivative of the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionMurgia2017), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber

    Murgia2017LogarithmicDerivative=+self%transferFunctionCDM%logarithmicDerivative(wavenumber)
    if (self%n_alpha > 0.0d0)                                       &
         & Murgia2017LogarithmicDerivative=+Murgia2017LogarithmicDerivative &
         &                               +self%n_beta                     &
         &                               *self%n_gamma                    &
         &                               *(                             &
         &                                   self%n_alpha               &
         &                                    *wavenumber)**self%n_beta      &
         &                                 /(                           &
         &                                  (                           &
         &                                   +1.0d0                     &
         &                                   +(                         &
         &                                     +self%n_alpha            &
         &                                     *wavenumber)**(self%n_beta) &
         &                                  )*wavenumber                &
         &                                )                             
    return
  end function Murgia2017LogarithmicDerivative

  double precision function Murgia2017HalfModeMass(self,status)
    !% Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    !% to a \gls{cdm} transfer function.
    use :: Galacticus_Error        , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionMurgia2017), intent(inout)           :: self
    integer                                   , intent(  out), optional :: status
    double precision                          , parameter               :: wavenumberHalfModeScaleFree=sqrt(0.25d0+2.0d0*log(2.0d0))-0.5d0
    double precision                                                    :: matterDensity                                                  , wavenumberHalfMode

    matterDensity       =+self%cosmologyParameters_%OmegaMatter    () &
         &               *self%cosmologyParameters_%densityCritical()
    wavenumberHalfMode  =+(                            &
         &                 +(1.0d0/self%n_alpha)        &
         &                 *((1.0d0/sqrt(2.0d0))**(1.0d0/self%n_gamma)  &
         &                 -1.0d0)**(1/self%n_beta))
    Murgia2017HalfModeMass=+4.0d0              &
         &               *Pi                   &
         &               /3.0d0                &
         &               *matterDensity        &
         &               *(                    &
         &                 +Pi                 &
         &                 /wavenumberHalfMode &
         &               )**3
    if (present(status)) status=errorStatusSuccess
    return
  end function Murgia2017HalfModeMass

  double precision function Murgia2017EpochTime(self)
    !% Return the cosmic time at the epoch at which this transfer function is defined.
    implicit none
    class(transferFunctionMurgia2017), intent(inout) :: self

    Murgia2017EpochTime=self%time
    return
  end function Murgia2017EpochTime
