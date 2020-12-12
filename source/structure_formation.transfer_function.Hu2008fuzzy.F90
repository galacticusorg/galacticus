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

!% Contains a module which implements a transfer function class for fuzzy dark matter based on Hu et al 2008. More precise than fuzzyDM transfer fuction using Murgia nCDM framework. 

  use :: Cosmology_Functions  , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters , only : cosmologyParametersClass
  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !# <transferFunction name="transferFunctionHu2008fuzzy">
  !#  <description>Provides a transfer function based on the fuzzy dark matter paper Hu et al 2008.</description>
  !# </transferFunction>
  type, extends(transferFunctionClass) :: transferFunctionHu2008fuzzy
     !% A transfer function class for fuzzy dark matter
     private
     double precision                                    :: m_22                         , n_beta        , &
          &                                                 n_gamma                         ,  &
          &                                                 time                          , redshift
     class           (transferFunctionClass   ), pointer :: transferFunctionCDM  => null()
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     class           (darkMatterParticleClass ), pointer :: darkMatterParticle_  => null()
   contains
     final     ::                          Hu2008fuzzyDestructor
     procedure :: value                 => Hu2008fuzzyValue
     procedure :: logarithmicDerivative => Hu2008fuzzyLogarithmicDerivative
     procedure :: halfModeMass          => Hu2008fuzzyHalfModeMass
     procedure :: epochTime             => Hu2008fuzzyEpochTime
  end type transferFunctionHu2008fuzzy

  interface transferFunctionHu2008fuzzy
     !% Constructors for the Hu 2008 fuzzy dark matter transfer function class.
     module procedure Hu2008fuzzyConstructorParameters
     module procedure Hu2008fuzzyConstructorInternal
  end interface transferFunctionHu2008fuzzy

contains

  function Hu2008fuzzyConstructorParameters(parameters) result(self)
    !% Constructor for fuzzy DM  transfer function class which takes a parameter set as input.
    use :: Cosmology_Functions           , only : cosmologyFunctions        , cosmologyFunctionsClass
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    use :: Galacticus_Error              , only : Galacticus_Error_Report
    use :: Input_Parameters              , only : inputParameter            , inputParameters
    implicit none
    type            (transferFunctionHu2008fuzzy)                :: self
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
    !#   <description>Normalized dark matter particle mass. m_22 is the actual particle mass divided by 10^-22 eV.</description>
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
     self=transferFunctionHu2008fuzzy(transferFunctionCDM,m_22,n_beta,n_gamma,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),cosmologyParameters_,darkMatterParticle_,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    !# <objectDestructor name="cosmologyFunctions_" />
    !# <objectDestructor name="darkMatterParticle_" />
    !# <objectDestructor name="transferFunctionCDM" />
    return
  end function Hu2008fuzzyConstructorParameters

  function Hu2008fuzzyConstructorInternal(transferFunctionCDM,m_22,n_beta,n_gamma,time,cosmologyParameters_,darkMatterParticle_,cosmologyFunctions_) result(self)
    !% Internal constructor for the Hu fuzzy dark matter transfer function class.
    use :: Cosmology_Parameters , only : hubbleUnitsLittleH
    use :: Dark_Matter_Particles, only : darkMatterParticleWDMThermal
    use :: Galacticus_Error     , only : Galacticus_Error_Report
    implicit none
    type            (transferFunctionHu2008fuzzy)                         :: self
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
  end function Hu2008fuzzyConstructorInternal

  subroutine Hu2008fuzzyDestructor(self)
    !% Destructor for the Hu fuzzy DM transfer function class.
    implicit none
    type(transferFunctionHu2008fuzzy), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%darkMatterParticle_" />
    !# <objectDestructor name="self%transferFunctionCDM" />
    return
  end subroutine Hu2008fuzzyDestructor

  double precision function Hu2008fuzzyValue(self,wavenumber)
    !% Return the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionHu2008fuzzy), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber
    double precision                                          :: x
			      
    Hu2008fuzzyValue=+self%transferFunctionCDM%value(wavenumber)
    x                = (1.61d0/9.00d0)*(self%m_22**(-4.0d0/9.0d0))*wavenumber

    Hu2008fuzzyValue       =+Hu2008fuzzyValue    &
         &                         *(            &
         &             cos(x**3)/(1+x**8))   
    return
  end function Hu2008fuzzyValue

  double precision function Hu2008fuzzyLogarithmicDerivative(self,wavenumber)
    !% Return the logarithmic derivative of the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionHu2008fuzzy), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber
    double precision                                          :: x

    Hu2008fuzzyLogarithmicDerivative=+self%transferFunctionCDM%logarithmicDerivative(wavenumber)
    x                = (1.61d0/9.00d0)*(self%m_22**(-4.0d0/9.0d0))*wavenumber

    Hu2008fuzzyLogarithmicDerivative = Hu2008fuzzyLogarithmicDerivative + (-8.0d0*x**7 - 3*(x**10 + x**2)*tan(x**3))  / (1+x**8)                    
  end function Hu2008fuzzyLogarithmicDerivative

  double precision function Hu2008fuzzyHalfModeMass(self,status)
    !% Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    !% to a \gls{cdm} transfer function.
    use :: Galacticus_Error        , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionHu2008fuzzy), intent(inout)            :: self
    integer                                   , intent(  out), optional :: status
    double precision                          , parameter               :: wavenumberHalfModeScaleFree=sqrt(0.25d0+2.0d0*log(2.0d0))-0.5d0
    double precision                                                    :: matterDensity                                                  , wavenumberHalfMode

    matterDensity       =+self%cosmologyParameters_%OmegaMatter    () &
         &               *self%cosmologyParameters_%densityCritical()

    wavenumberHalfMode  = 4.5*self%m_22**(4.0d0/9.0d0)
                        
    Hu2008fuzzyHalfModeMass=+4.0d0              &
         &               *Pi                   &
         &               /3.0d0                &
         &               *matterDensity        &
         &               *(                    &
         &                 +Pi                 &
         &                 /wavenumberHalfMode &
         &               )**3
    if (present(status)) status=errorStatusSuccess
    return
  end function Hu2008fuzzyHalfModeMass

  double precision function Hu2008fuzzyEpochTime(self)
    !% Return the cosmic time at the epoch at which this transfer function is defined.
    implicit none
    class(transferFunctionHu2008fuzzy), intent(inout) :: self

    Hu2008fuzzyEpochTime=self%time
    return
  end function Hu2008fuzzyEpochTime
