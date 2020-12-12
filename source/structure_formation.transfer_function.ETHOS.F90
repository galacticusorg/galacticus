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

!% Contains a module which implements a transfer function class based on the ETHOS model.

  use :: Cosmology_Functions  , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters , only : cosmologyParametersClass
  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !# <transferFunction name="transferFunctionETHOSDM">
  !#  <description>Provides a transfer function based on the ETHOS model </description>
  !# </transferFunction>
  type, extends(transferFunctionClass) :: transferFunctionETHOSDM
     !% A transfer function class which modifies the CDM transfer function to fit the ETHOS model. 
     private
     double precision                                    :: n_alpha                         , n_beta        , &
          &                                                 n_gamma                         , n_sigma, n_tau, &
	  &						    k_peak, h_peak                  , h_2        , &
          &                                                 time                          , redshift
     class           (transferFunctionClass   ), pointer :: transferFunctionCDM  => null()
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     class           (darkMatterParticleClass ), pointer :: darkMatterParticle_  => null()
   contains
     final     ::                          ETHOSDMDestructor
     procedure :: value                 => ETHOSDMValue
     procedure :: logarithmicDerivative => ETHOSDMLogarithmicDerivative
     procedure :: halfModeMass          => ETHOSDMHalfModeMass
     procedure :: epochTime             => ETHOSDMEpochTime
  end type transferFunctionETHOSDM

  interface transferFunctionETHOSDM
     !% Constructors for the ETHOS transfer function class.
     module procedure ETHOSDMConstructorParameters
     module procedure ETHOSDMConstructorInternal
  end interface transferFunctionETHOSDM

contains

  function ETHOSDMConstructorParameters(parameters) result(self)
    !% Constructor for the ETHOS transfer function class which takes a parameter set as input.
    use :: Cosmology_Functions           , only : cosmologyFunctions        , cosmologyFunctionsClass
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    use :: Galacticus_Error              , only : Galacticus_Error_Report
    use :: Input_Parameters              , only : inputParameter            , inputParameters
    implicit none
    type            (transferFunctionETHOSDM)                 :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (transferFunctionClass   ), pointer       :: transferFunctionCDM
    class           (cosmologyParametersClass), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    class           (darkMatterParticleClass ), pointer       :: darkMatterParticle_
    double precision                                          :: n_alpha             , n_beta     , &
         &                                                       n_gamma             , n_sigma    , & 
	 &							 n_tau,  k_peak, h_peak, h_2	  , &
         &			                                 redshift

    ! Validate parameters.
    if (.not.parameters%isPresent('transferFunctionMethod')) call Galacticus_Error_Report("an explicit 'transferFunctionMethod' must be given"//{introspection:location})
    ! Read parameters.
    !# <inputParameter>
    !#   <name>n_alpha</name>
    !#   <source>parameters</source>
    !#   <defaultValue>40.0d0</defaultValue>
    !#   <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
    !#   <description>The parameter $\alpha$ from nCDM transfer function [Murgia et al (2017)], sets the cutoff scale length.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>n_beta</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.5d0</defaultValue>
    !#   <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
    !#   <description>The parameter $\beta$ from nCDM transfer function [Murgia et al (2017)], controls shape of cutoff .</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>n_gamma</name>
    !#   <source>parameters</source>
    !#   <defaultValue>-10d0</defaultValue>
    !#   <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
    !#   <description>The parameter $\gamma$ from nCDM transfer function [Murgia et al (2017)], controls shape of cutoff.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>n_sigma</name>
    !#   <source>parameters</source>
    !#   <defaultValue>-10d0</defaultValue>
    !#   <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
    !#   <description>From ETHOS paper [Bohr et al (2020)], determines width of first peak in transfer function.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>n_tau</name>
    !#   <source>parameters</source>
    !#   <defaultValue>-10d0</defaultValue>
    !#   <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
    !#   <description>From ETHOS paper [Bohr et al (2020)], determines damping of DAO.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>k_peak</name>
    !#   <source>parameters</source>
    !#   <defaultValue>-10d0</defaultValue>
    !#   <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
    !#   <description>Wavenumber of first peak in ETHOS transfer function.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>h_peak</name>
    !#   <source>parameters</source>
    !#   <defaultValue>-10d0</defaultValue>
    !#   <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
    !#   <description>Amplitude of first peak in ETHOS transfer function.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>h_2</name>
    !#   <source>parameters</source>
    !#   <defaultValue>-10d0</defaultValue>
    !#   <defaultSource>\citep[][for the transfer function at $z=z_\mathrm{eq}$]{barkana_constraints_2001}</defaultSource>
    !#   <description>Amplitude of second peak in ETHOS transfer function.</description>
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
     self=transferFunctionETHOSDM(transferFunctionCDM,n_alpha,n_beta,n_gamma,n_sigma,n_tau,k_peak,h_peak,h_2,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),cosmologyParameters_,darkMatterParticle_,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    !# <objectDestructor name="cosmologyFunctions_" />
    !# <objectDestructor name="darkMatterParticle_" />
    !# <objectDestructor name="transferFunctionCDM" />
    return
  end function ETHOSDMConstructorParameters

  function ETHOSDMConstructorInternal(transferFunctionCDM,n_alpha,n_beta,n_gamma,n_sigma,n_tau,k_peak,h_peak,h_2,time,cosmologyParameters_,darkMatterParticle_,cosmologyFunctions_) result(self)
    !% Internal constructor for the ETHOS transfer function class.
    use :: Cosmology_Parameters , only : hubbleUnitsLittleH
    use :: Dark_Matter_Particles, only : darkMatterParticleWDMThermal
    use :: Galacticus_Error     , only : Galacticus_Error_Report
    implicit none
    type            (transferFunctionETHOSDM)                         :: self
    class           (transferFunctionClass   ), target, intent(in   ) :: transferFunctionCDM
    double precision                                  , intent(in   ) :: n_alpha                   , n_beta                            , &
         &                                                               n_gamma, n_sigma, n_tau, k_peak, h_peak,h_2 , time
    class           (cosmologyParametersClass), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), target, intent(in   ) :: cosmologyFunctions_
    class           (darkMatterParticleClass ), target, intent(in   ) :: darkMatterParticle_
    double precision                          , parameter             :: massReference       =1.0d0, degreesOfFreedomReference=1.5d0
    !# <constructorAssign variables="*transferFunctionCDM, n_alpha, n_beta, n_gamma,n_sigma, n_tau, k_peak, h_peak, h_2, time, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterParticle_"/>

    self%redshift=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    
    return
  end function ETHOSDMConstructorInternal

  subroutine ETHOSDMDestructor(self)
    !% Destructor for the ETHOS transfer function class.
    implicit none
    type(transferFunctionETHOSDM), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%darkMatterParticle_" />
    !# <objectDestructor name="self%transferFunctionCDM" />
    return
  end subroutine ETHOSDMDestructor

  double precision function ETHOSDMValue(self,wavenumber)
    !% Return the transfer function at the given wavenumber.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionETHOSDM), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber			      

    ETHOSDMValue=+self%transferFunctionCDM%value(wavenumber)
    if (self%n_alpha > 0.0d0)                    &
         & ETHOSDMValue=+ETHOSDMValue              &
         &              *((                          &
         &                  1.0d0                    &
         &                 +                         &
         &                  (                        &
         &                   +self%n_alpha             &
         &                   *wavenumber             &
         &                  )**(self%n_beta)           &
         &                )  **(self%n_gamma)          &
	 &              - sqrt(self%h_peak)*exp(-.5d0*((wavenumber-self%k_peak)/(self%n_sigma*self%k_peak))**2.0d0) &
         &              + .25d0*sqrt(self%h_2)*erfc((wavenumber - 1.805d0*self%k_peak)/(self%n_tau*self%k_peak)-2.0d0)  &
	 &              * erf(-(wavenumber - 1.805d0*self%k_peak)/(self%n_sigma*self%k_peak)-2.0d0)               &   
         &              * cos(1.1083d0*Pi*wavenumber/self%k_peak))
    return
  end function ETHOSDMValue

  double precision function ETHOSDMLogarithmicDerivative(self,wavenumber)
    !% Return the logarithmic derivative of the transfer function at the given wavenumber.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionETHOSDM), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber

    
    ETHOSDMLogarithmicDerivative=+self%transferFunctionCDM%logarithmicDerivative(wavenumber)
    if (self%n_alpha > 0.0d0)                                       &
         & ETHOSDMLogarithmicDerivative=+ETHOSDMLogarithmicDerivative                         &
         &                 +    (self%n_alpha*(wavenumber*self%n_alpha)**(-1 + self%n_beta)   &
         &                 *    (1 + (wavenumber*self%n_alpha)**self%n_beta)**(-1 + self%n_gamma)*self%n_beta     &
         &                 *    self%n_gamma + (Sqrt(self%h_peak)*(wavenumber - self%k_peak)) &
         &                 /    (Exp((wavenumber - self%k_peak)**2                            &
         &                 /    (2.*self%k_peak**2*self%n_sigma**2))*self%k_peak**2*self%n_sigma**2)             & 
         &                 -    (Sqrt(self%h_2)*Cos((3.4818271379735677*wavenumber)/self%k_peak)                 &
         &                 *    Erfc(-2 + (-wavenumber + 1.805*self%k_peak)/(self%k_peak*self%n_sigma)))         &
         &                 /    (2.*Exp(-2 + (wavenumber - 1.805*self%k_peak)/(self%k_peak*self%n_tau))**2       &
         &    *    self%k_peak *self%n_tau*Sqrt(Pi)) + (Sqrt(self%h_2)*Cos((3.4818271379735677*wavenumber)/self%k_peak) &    
         &                 *    Erfc(-2 + (wavenumber - 1.805*self%k_peak)/(self%k_peak*self%n_tau)))            &
         &    /    (2.*Exp(-2 + (-wavenumber + 1.805*self%k_peak)/(self%k_peak*self%n_sigma))**2*self%k_peak     &
         &                 *    self%n_sigma*Sqrt(Pi))- (0.8704567844933919*Sqrt(self%h_2)                &
         &                 *    Erfc(-2 + (-wavenumber + 1.805*self%k_peak)/(self%k_peak*self%n_sigma))          &
         &                 *    Erfc(-2 + (wavenumber - 1.805*self%k_peak)/(self%k_peak*self%n_tau))             &
         &                 *    Sin((3.4818271379735677*wavenumber)/self%k_peak))/self%k_peak)/(-(Sqrt(self%h_peak)  &
         &                 /    Exp((wavenumber - self%k_peak)**2/(2.*self%k_peak**2*self%n_sigma**2))) + (1         &
         &                 +    (wavenumber*self%n_alpha)**self%n_beta)**self%n_gamma + (Sqrt(self%h_2)              &
         &                 *    Cos((3.4818271379735677*wavenumber)/self%k_peak)                             &
         &                 *    Erfc(-2 + (-wavenumber + 1.805*self%k_peak)/(self%k_peak*self%n_sigma))            &
         &                 *    Erfc(-2 + (wavenumber - 1.805*self%k_peak)/(self%k_peak*self%n_tau)))/4.)
    return
  end function ETHOSDMLogarithmicDerivative

  double precision function ETHOSDMHalfModeMass(self,status)
    !% Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    !% to a \gls{cdm} transfer function.
    use :: Galacticus_Error        , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionETHOSDM), intent(inout)            :: self
    integer                                   , intent(  out), optional :: status
    double precision                          , parameter               :: wavenumberHalfModeScaleFree=sqrt(0.25d0+2.0d0*log(2.0d0))-0.5d0
    double precision                                                    :: matterDensity                                                  , wavenumberHalfMode

    matterDensity       =+self%cosmologyParameters_%OmegaMatter    () &
         &               *self%cosmologyParameters_%densityCritical()
    wavenumberHalfMode  =+self%k_peak/2.5d0
    ETHOSDMHalfModeMass=+4.0d0              &
         &               *Pi                   &
         &               /3.0d0                &
         &               *matterDensity        &
         &               *(                    &
         &                 +Pi                 &
         &                 /wavenumberHalfMode &
         &               )**3
    if (present(status)) status=errorStatusSuccess
    return
  end function ETHOSDMHalfModeMass

  double precision function ETHOSDMEpochTime(self)
    !% Return the cosmic time at the epoch at which this transfer function is defined.
    implicit none
    class(transferFunctionETHOSDM), intent(inout) :: self

    ETHOSDMEpochTime=self%time
    return
  end function ETHOSDMEpochTime
