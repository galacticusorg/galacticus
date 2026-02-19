!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  An intergalactic medium state decorator class which provides a fixed metallicity for the \gls{igm}.
  !!}

  !![
  <intergalacticMediumState name="intergalacticMediumStateMetallicityFixed">
   <description>
    An intergalactic medium state class which provides a fixed metallicity for the \gls{igm}, given by {\normalfont \ttfamily [metallicity]}.
   </description>
  </intergalacticMediumState>
  !!]
  type, extends(intergalacticMediumStateClass) :: intergalacticMediumStateMetallicityFixed
     !!{
     An intergalactic medium state class which provides a fixed metallicity for the \gls{igm}, given by {\normalfont \ttfamily [metallicity]}.
     !!}
     private
     class           (intergalacticMediumStateClass), pointer :: intergalacticMediumState_ => null()
     double precision                                         :: metallicity_
   contains
     final     ::                                metallicityFixedDestructor
     procedure :: electronFraction            => metallicityFixedElectronFraction
     procedure :: temperature                 => metallicityFixedTemperature
     procedure :: neutralHydrogenFraction     => metallicityFixedNeutralHydrogenFraction
     procedure :: neutralHeliumFraction       => metallicityFixedNeutralHeliumFraction
     procedure :: singlyIonizedHeliumFraction => metallicityFixedSinglyIonizedHeliumFraction
     procedure :: metallicity                 => metallicityFixedMetallicity
  end type intergalacticMediumStateMetallicityFixed

  interface intergalacticMediumStateMetallicityFixed
     !!{
     Constructors for the \refClass{intergalacticMediumStateMetallicityFixed} intergalactic medium state class.
     !!}
     module procedure metallicityFixedIGMConstructorParameters
     module procedure metallicityFixedIGMConstructorInternal
  end interface intergalacticMediumStateMetallicityFixed

contains

  function metallicityFixedIGMConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{intergalacticMediumStateMetallicityFixed} \gls{igm} state class which takes a parameter set as input.
    !!}
    use :: Input_Parameters                , only : inputParameter  , inputParameters
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    implicit none
    type            (intergalacticMediumStateMetallicityFixed)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (intergalacticMediumStateClass           ), pointer       :: intergalacticMediumState_
    double precision                                                          :: metallicity

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>metallicity</name>
      <source>parameters</source>
      <description>The metallicity (relative to Solar) of the \gls{igm}.</description>
    </inputParameter>
    <objectBuilder class="intergalacticMediumState" name="intergalacticMediumState_" source="parameters"/>
    !!]
    ! Convert metallicity to a mass fraction.
    metallicity=+metallicity      &
         &      *metallicitySolar
    ! Construct the object.
    self=intergalacticMediumStateMetallicityFixed(metallicity,intergalacticMediumState_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="intergalacticMediumState_"/>
    !!]
    return
  end function metallicityFixedIGMConstructorParameters

  function metallicityFixedIGMConstructorInternal(metallicity,intergalacticMediumState_) result(self)
    !!{
    Constructor for the \refClass{intergalacticMediumStateMetallicityFixed} \gls{igm} state class.
    !!}
    implicit none
    type            (intergalacticMediumStateMetallicityFixed)                        :: self
    double precision                                          , intent(in   )         :: metallicity
    class           (intergalacticMediumStateClass           ), intent(inout), target :: intergalacticMediumState_
    !![
    <constructorAssign variables="*intergalacticMediumState_"/>
    !!]

    self%metallicity_=metallicity
    return
  end function metallicityFixedIGMConstructorInternal

  subroutine metallicityFixedDestructor(self)
    !!{
    Destructor for the metallicityFixed \gls{igm} state class.
    !!}
    implicit none
    type(intergalacticMediumStateMetallicityFixed), intent(inout) :: self

    !![
    <objectDestructor name="self%intergalacticMediumState_"/>
    !!]
    return
  end subroutine metallicityFixedDestructor

  double precision function metallicityFixedElectronFraction(self,time)
    !!{
    Return the electron fraction of the \gls{igm}.
    !!}
    implicit none
    class           (intergalacticMediumStateMetallicityFixed), intent(inout) :: self
    double precision                                          , intent(in   ) :: time

    metallicityFixedElectronFraction=self%intergalacticMediumState_%electronFraction(time)
    return
  end function metallicityFixedElectronFraction

  double precision function metallicityFixedNeutralHydrogenFraction(self,time)
    !!{
    Return the neutral hydrogen fraction of the \gls{igm}.
    !!}
    implicit none
    class           (intergalacticMediumStateMetallicityFixed), intent(inout) :: self
    double precision                                          , intent(in   ) :: time

    metallicityFixedNeutralHydrogenFraction=self%intergalacticMediumState_%neutralHydrogenFraction(time)  
    return
  end function metallicityFixedNeutralHydrogenFraction

  double precision function metallicityFixedNeutralHeliumFraction(self,time)
    !!{
    Return the neutral helium fraction of the \gls{igm}.
    !!}
    implicit none
    class           (intergalacticMediumStateMetallicityFixed), intent(inout) :: self
    double precision                                          , intent(in   ) :: time

    metallicityFixedNeutralHeliumFraction=self%intergalacticMediumState_%neutralHeliumFraction(time)  
    return
  end function metallicityFixedNeutralHeliumFraction

  double precision function metallicityFixedSinglyIonizedHeliumFraction(self,time)
    !!{
    Return the singly-ionized helium fraction of the \gls{igm}.
    !!}
    implicit none
    class           (intergalacticMediumStateMetallicityFixed), intent(inout) :: self
    double precision                                          , intent(in   ) :: time

    metallicityFixedSinglyIonizedHeliumFraction=self%intergalacticMediumState_%singlyIonizedHeliumFraction(time)  
    return
  end function metallicityFixedSinglyIonizedHeliumFraction

  double precision function metallicityFixedTemperature(self,time)
    !!{
    Return the temperature of the \gls{igm}.
    !!}
    implicit none
    class           (intergalacticMediumStateMetallicityFixed), intent(inout) :: self
    double precision                                          , intent(in   ) :: time

    metallicityFixedTemperature=self%intergalacticMediumState_%temperature(time)  
    return
  end function metallicityFixedTemperature

  double precision function metallicityFixedMetallicity(self,time)
    !!{
    Return the metallicity of the \gls{igm}.
    !!}
    implicit none
    class           (intergalacticMediumStateMetallicityFixed), intent(inout) :: self
    double precision                                          , intent(in   ) :: time

    metallicityFixedMetallicity=self%metallicity_
    return
  end function metallicityFixedMetallicity
