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

!+    Contributions to this file made by: Niusha Ahvazi
  
  !!{
  An  intergalactic medium state decorator class which provides a fit to the metallicity evolution of the \gls{igm} reported by \cite{ucci_astraeus_2023}.
  !!}

  !![
  <intergalacticMediumState name="intergalacticMediumStateMetallicityAstraeusV">
   <description>
    An intergalactic medium state class which provides a fit to the metallicity evolution of the \gls{igm} reported by \cite{ucci_astraeus_2023}.
   </description>
  </intergalacticMediumState>
  !!]
  type, extends(intergalacticMediumStateClass) :: intergalacticMediumStateMetallicityAstraeusV
     !!{
     An intergalactic medium state class which provides a fit to the metallicity evolution of the \gls{igm} reported by \cite{ucci_astraeus_2023}.
     !!}
     private
     class(intergalacticMediumStateClass), pointer :: intergalacticMediumState_ => null()
   contains
     final     ::                                metallicityAstraeusVDestructor
     procedure :: electronFraction            => metallicityAstraeusVElectronFraction
     procedure :: temperature                 => metallicityAstraeusVTemperature
     procedure :: neutralHydrogenFraction     => metallicityAstraeusVNeutralHydrogenFraction
     procedure :: neutralHeliumFraction       => metallicityAstraeusVNeutralHeliumFraction
     procedure :: singlyIonizedHeliumFraction => metallicityAstraeusVSinglyIonizedHeliumFraction
     procedure :: metallicity                 => metallicityAstraeusVMetallicity
  end type intergalacticMediumStateMetallicityAstraeusV

  interface intergalacticMediumStateMetallicityAstraeusV
     !!{
     Constructors for the \refClass{intergalacticMediumStateMetallicityAstraeusV} intergalactic medium state class.
     !!}
     module procedure metallicityAstraeusVIGMConstructorParameters
     module procedure metallicityAstraeusVIGMConstructorInternal
  end interface intergalacticMediumStateMetallicityAstraeusV

contains

  function metallicityAstraeusVIGMConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{intergalacticMediumStateMetallicityAstraeusV} \gls{igm} state class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (intergalacticMediumStateMetallicityAstraeusV)                :: self
    type (inputParameters                             ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                     ), pointer       :: cosmologyFunctions_
    class(intergalacticMediumStateClass               ), pointer       :: intergalacticMediumState_

    !![
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="intergalacticMediumState" name="intergalacticMediumState_" source="parameters"/>
    !!]
    ! Construct the object.
    self=intergalacticMediumStateMetallicityAstraeusV(cosmologyFunctions_,intergalacticMediumState_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="intergalacticMediumState_"/>
    !!]
    return
  end function metallicityAstraeusVIGMConstructorParameters

  function metallicityAstraeusVIGMConstructorInternal(cosmologyFunctions_,intergalacticMediumState_) result(self)
    !!{
    Constructor for the \refClass{intergalacticMediumStateMetallicityAstraeusV} \gls{igm} state class.
    !!}
    implicit none
    type (intergalacticMediumStateMetallicityAstraeusV)                        :: self
    class(cosmologyFunctionsClass                     ), intent(inout), target :: cosmologyFunctions_
    class(intergalacticMediumStateClass               ), intent(inout), target :: intergalacticMediumState_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *intergalacticMediumState_"/>
    !!]

    return
  end function metallicityAstraeusVIGMConstructorInternal

  subroutine metallicityAstraeusVDestructor(self)
    !!{
    Destructor for the metallicityAstraeusV \gls{igm} state class.
    !!}
    implicit none
    type(intergalacticMediumStateMetallicityAstraeusV), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%intergalacticMediumState_"/>
    !!]
    return
  end subroutine metallicityAstraeusVDestructor

  double precision function metallicityAstraeusVElectronFraction(self,time)
    !!{
    Return the electron fraction of the \gls{igm}.
    !!}
    implicit none
    class           (intergalacticMediumStateMetallicityAstraeusV), intent(inout) :: self
    double precision                                              , intent(in   ) :: time

    metallicityAstraeusVElectronFraction=self%intergalacticMediumState_%electronFraction(time)
    return
  end function metallicityAstraeusVElectronFraction

  double precision function metallicityAstraeusVNeutralHydrogenFraction(self,time)
    !!{
    Return the neutral hydrogen fraction of the \gls{igm}.
    !!}
    implicit none
    class           (intergalacticMediumStateMetallicityAstraeusV), intent(inout) :: self
    double precision                                              , intent(in   ) :: time

    metallicityAstraeusVNeutralHydrogenFraction=self%intergalacticMediumState_%neutralHydrogenFraction(time)  
    return
  end function metallicityAstraeusVNeutralHydrogenFraction

  double precision function metallicityAstraeusVNeutralHeliumFraction(self,time)
    !!{
    Return the neutral helium fraction of the \gls{igm}.
    !!}
    implicit none
    class           (intergalacticMediumStateMetallicityAstraeusV), intent(inout) :: self
    double precision                                              , intent(in   ) :: time

    metallicityAstraeusVNeutralHeliumFraction=self%intergalacticMediumState_%neutralHeliumFraction(time)  
    return
  end function metallicityAstraeusVNeutralHeliumFraction

  double precision function metallicityAstraeusVSinglyIonizedHeliumFraction(self,time)
    !!{
    Return the singly-ionized helium fraction of the \gls{igm}.
    !!}
    implicit none
    class           (intergalacticMediumStateMetallicityAstraeusV), intent(inout) :: self
    double precision                                              , intent(in   ) :: time

    metallicityAstraeusVSinglyIonizedHeliumFraction=self%intergalacticMediumState_%singlyIonizedHeliumFraction(time)  
    return
  end function metallicityAstraeusVSinglyIonizedHeliumFraction

  double precision function metallicityAstraeusVTemperature(self,time)
    !!{
    Return the temperature of the \gls{igm}.
    !!}
    implicit none
    class           (intergalacticMediumStateMetallicityAstraeusV), intent(inout) :: self
    double precision                                              , intent(in   ) :: time

    metallicityAstraeusVTemperature=self%intergalacticMediumState_%temperature(time)  
    return
  end function metallicityAstraeusVTemperature

  double precision function metallicityAstraeusVMetallicity(self,time)
    !!{
    Return the metallicity of the \gls{igm}.
    !!}
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    implicit none
    class           (intergalacticMediumStateMetallicityAstraeusV), intent(inout) :: self
    double precision                                              , intent(in   ) :: time
    double precision                                                              :: redshift

    ! Evaluate the IGM metallicity using a fit to the Astraeus V results (Ucci et al.; 2021, arXiv:2112.02115;
    ! https://ui.adsabs.harvard.edu/abs/2021arXiv211202115U; figure 9; solid blue curve "Photoionization" for all galaxies)
    ! computed by Niusha Ahvazi.
    redshift                       =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    metallicityAstraeusVMetallicity=+10.0d0**(                     &
         &                                    -4.30d-4*redshift**3 &
         &                                    +4.47d-3*redshift**2 &
         &                                    -2.54d-1*redshift**1 &
         &                                    -1.94d+0             &
         &                                   )                     &
         &                          *metallicitySolar
    return
  end function metallicityAstraeusVMetallicity
