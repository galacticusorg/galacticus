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

  !# <transferFunction name="transferFunctionfuzzyDM">
  !#  <description>Provides a transfer function based on the thermal \gls{wdm} modifier of \cite{bode_halo_2001}.</description>
  !# </transferFunction>
  type, extends(transferFunctionClass) :: transferFunctionfuzzyDM
     !% A transfer function class which modifies another transfer function using the thermal \gls{wdm} modifier of \cite{bode_halo_2001}.
     private
     double precision                                    :: m_22                         , n_beta        , &
          &                                                 n_gamma                         ,  &
          &                                                 time                          , redshift
     class           (transferFunctionClass   ), pointer :: transferFunctionCDM  => null()
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     class           (darkMatterParticleClass ), pointer :: darkMatterParticle_  => null()
   contains
     final     ::                          fuzzyDMDestructor
     procedure :: value                 => fuzzyDMValue
     procedure :: logarithmicDerivative => fuzzyDMLogarithmicDerivative
     procedure :: halfModeMass          => fuzzyDMHalfModeMass
     procedure :: epochTime             => fuzzyDMEpochTime
  end type transferFunctionfuzzyDM

  interface transferFunctionfuzzyDM
     !% Constructors for the {\normalfont \ttfamily bode2001} transfer function class.
     module procedure fuzzyDMConstructorParameters
     module procedure fuzzyDMConstructorInternal
  end interface transferFunctionfuzzyDM

contains

  function fuzzyDMConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily bode2001} transfer function class which takes a parameter set as input.
    use :: Cosmology_Functions           , only : cosmologyFunctions        , cosmologyFunctionsClass
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    use :: Galacticus_Error              , only : Galacticus_Error_Report
    use :: Input_Parameters              , only : inputParameter            , inputParameters
    implicit none
    type            (transferFunctionfuzzyDM)                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (transferFunctionClass   ), pointer       :: transferFunctionCDM
    class           (cosmologyParametersClass), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    class           (darkMatterParticleClass ), pointer       :: darkMatterParticle_
    double precision                                          :: m_22             , n_beta     , &
         &                                                       n_gamma                  , redshift

    ! Validate parameters.
    if (.not.parameters%isPresent('transferFunctionMethod')) call Galacticus_Error_Report("an explicit 'transferFunctionMethod' must be given"//{introspection:location})
    ! Read parameters.
    !# <inputParameter>
    !#   <name>m_22</name>
    !#   <source>parameters</source>
    !#   <defaultValue>40.0d0</defaultValue>
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
     self=transferFunctionfuzzyDM(transferFunctionCDM,m_22,n_beta,n_gamma,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),cosmologyParameters_,darkMatterParticle_,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    !# <objectDestructor name="cosmologyFunctions_" />
    !# <objectDestructor name="darkMatterParticle_" />
    !# <objectDestructor name="transferFunctionCDM" />
    return
  end function fuzzyDMConstructorParameters

  function fuzzyDMConstructorInternal(transferFunctionCDM,m_22,n_beta,n_gamma,time,cosmologyParameters_,darkMatterParticle_,cosmologyFunctions_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily bode2001} transfer function class.
    use :: Cosmology_Parameters , only : hubbleUnitsLittleH
    use :: Dark_Matter_Particles, only : darkMatterParticleWDMThermal
    use :: Galacticus_Error     , only : Galacticus_Error_Report
    implicit none
    type            (transferFunctionfuzzyDM)                         :: self
    class           (transferFunctionClass   ), target, intent(in   ) :: transferFunctionCDM
    double precision                                  , intent(in   ) :: m_22                   , n_beta                            , &
         &                                                               n_gamma                        , time
    class           (cosmologyParametersClass), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), target, intent(in   ) :: cosmologyFunctions_
    class           (darkMatterParticleClass ), target, intent(in   ) :: darkMatterParticle_
    double precision                          , parameter             :: massReference       =1.0d0, degreesOfFreedomReference=1.5d0
    !# <constructorAssign variables="*transferFunctionCDM, m_22, n_beta, n_gamma, time, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterParticle_"/>

    self%redshift=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    
    return
  end function fuzzyDMConstructorInternal

  subroutine fuzzyDMDestructor(self)
    !% Destructor for the {\normalfont \ttfamily bode2001} transfer function class.
    implicit none
    type(transferFunctionfuzzyDM), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%darkMatterParticle_" />
    !# <objectDestructor name="self%transferFunctionCDM" />
    return
  end subroutine fuzzyDMDestructor

  double precision function fuzzyDMValue(self,wavenumber)
    !% Return the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionfuzzyDM), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber
    double precision                                          :: n_alpha			      

    fuzzyDMValue=+self%transferFunctionCDM%value(wavenumber)
    n_alpha = 0.11d0*(self%m_22)**(-4.0d0/9.0d0) 
    if (n_alpha > 0.0d0)                    &
         & fuzzyDMValue=+fuzzyDMValue              &
         &               *(                          &
         &                  1.0d0                    &
         &                 +                         &
         &                  (                        &
         &                   +n_alpha             &
         &                   *wavenumber             &
         &                  )**(self%n_beta)           &
         &                )  **(self%n_gamma)
    return
  end function fuzzyDMValue

  double precision function fuzzyDMLogarithmicDerivative(self,wavenumber)
    !% Return the logarithmic derivative of the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionfuzzyDM), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber
    double precision                                          :: n_alpha

    n_alpha = 0.11d0*(self%m_22)**(-4.0d0/9.0d0)
    fuzzyDMLogarithmicDerivative=+self%transferFunctionCDM%logarithmicDerivative(wavenumber)
    if (n_alpha > 0.0d0)                                       &
         & fuzzyDMLogarithmicDerivative=+fuzzyDMLogarithmicDerivative &
         &                               +self%n_beta                     &
         &                               *self%n_gamma                    &
         &                               *(                             &
         &                                   n_alpha                    &
         &                                    *wavenumber)**self%n_beta      &
         &                                 /(                           &
         &                                  (                           &
         &                                   +1.0d0                     &
         &                                   +(                         &
         &                                     +n_alpha            &
         &                                     *wavenumber)**(self%n_beta) &
         &                                  )*wavenumber                &
         &                                )                             
    return
  end function fuzzyDMLogarithmicDerivative

  double precision function fuzzyDMHalfModeMass(self,status)
    !% Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    !% to a \gls{cdm} transfer function.
    use :: Galacticus_Error        , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionfuzzyDM), intent(inout)            :: self
    integer                                   , intent(  out), optional :: status
    double precision                          , parameter               :: wavenumberHalfModeScaleFree=sqrt(0.25d0+2.0d0*log(2.0d0))-0.5d0
    double precision                                                    :: matterDensity                                                  , wavenumberHalfMode
    double precision                                          :: n_alpha

    n_alpha = 0.11d0*(self%m_22)**(-4.0d0/9.0d0)
    matterDensity       =+self%cosmologyParameters_%OmegaMatter    () &
         &               *self%cosmologyParameters_%densityCritical()
    wavenumberHalfMode  =+(                            &
         &                 +(1.0d0/n_alpha)        &
         &                 *((1.0d0/sqrt(2.0d0))**(1.0d0/self%n_gamma)  &
         &                 -1.0d0)**(1/self%n_beta))
    fuzzyDMHalfModeMass=+4.0d0              &
         &               *Pi                   &
         &               /3.0d0                &
         &               *matterDensity        &
         &               *(                    &
         &                 +Pi                 &
         &                 /wavenumberHalfMode &
         &               )**3
    if (present(status)) status=errorStatusSuccess
    return
  end function fuzzyDMHalfModeMass

  double precision function fuzzyDMEpochTime(self)
    !% Return the cosmic time at the epoch at which this transfer function is defined.
    implicit none
    class(transferFunctionfuzzyDM), intent(inout) :: self

    fuzzyDMEpochTime=self%time
    return
  end function fuzzyDMEpochTime
