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
  An intergalactic medium state decorator class which provides a fixed metallicity for the \gls{igm}.
  !!}

  !![
  <intergalacticMediumState name="intergalacticMediumStateMetallicityPolynomial">
   <description>
    An intergalactic medium state class which provides a fixed metallicity for the \gls{igm}, given by {\normalfont \ttfamily [metallicity]}.
   </description>
  </intergalacticMediumState>
  !!]
  type, extends(intergalacticMediumStateClass) :: intergalacticMediumStateMetallicityPolynomial
     !!{
     An intergalactic medium state class which provides a fixed metallicity for the \gls{igm}, given by {\normalfont \ttfamily [metallicity]}.
     !!}
     private
     class           (intergalacticMediumStateClass), pointer :: intergalacticMediumState_ => null()
     double precision, allocatable, dimension(:)              :: coefficients
   contains
     final     ::                                metallicityPolynomialDestructor
     procedure :: electronFraction            => metallicityPolynomialElectronFraction
     procedure :: temperature                 => metallicityPolynomialTemperature
     procedure :: neutralHydrogenFraction     => metallicityPolynomialNeutralHydrogenFraction
     procedure :: neutralHeliumFraction       => metallicityPolynomialNeutralHeliumFraction
     procedure :: singlyIonizedHeliumFraction => metallicityPolynomialSinglyIonizedHeliumFraction
     procedure :: metallicity                 => metallicityPolynomialMetallicity
  end type intergalacticMediumStateMetallicityPolynomial

  interface intergalacticMediumStateMetallicityPolynomial
     !!{
     Constructors for the \refClass{intergalacticMediumStateMetallicityPolynomial} intergalactic medium state class.
     !!}
     module procedure metallicityPolynomialIGMConstructorParameters
     module procedure metallicityPolynomialIGMConstructorInternal
  end interface intergalacticMediumStateMetallicityPolynomial

contains

  function metallicityPolynomialIGMConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{intergalacticMediumStateMetallicityPolynomial} \gls{igm} state class which takes a parameter set as input.
    !!}
    use :: Input_Parameters                , only : inputParameter  , inputParameters
    implicit none
    type            (intergalacticMediumStateMetallicityPolynomial)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                     ), pointer :: cosmologyFunctions_
    class           (intergalacticMediumStateClass           ), pointer       :: intergalacticMediumState_
    double precision, allocatable, dimension(:)                               :: coefficients

    ! Check and read parameters.
    allocate(coefficients(parameters%count('coefficients')))
    !![
    <inputParameter>
      <name>coefficients</name>
      <source>parameters</source>
      <description>The polynomial coefficients, $c_i$, in the function $Z(z)/Z_\odot = 10^{\sum_{i=0}^N c_i [\log_{10}(1+z)]^i}$.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="intergalacticMediumState" name="intergalacticMediumState_" source="parameters"/>
    !!]
    ! Construct the object.
    self=intergalacticMediumStateMetallicityPolynomial(coefficients,cosmologyFunctions_,intergalacticMediumState_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="intergalacticMediumState_"/>
    !!]
    return
  end function metallicityPolynomialIGMConstructorParameters

  function metallicityPolynomialIGMConstructorInternal(coefficients,cosmologyFunctions_,intergalacticMediumState_) result(self)
    !!{
    Constructor for the \refClass{intergalacticMediumStateMetallicityPolynomial} \gls{igm} state class.
    !!}
    implicit none
    type            (intergalacticMediumStateMetallicityPolynomial)                        :: self
    double precision                               , intent(in   ), dimension(:)          :: coefficients
    class(cosmologyFunctionsClass), intent(inout), target :: cosmologyFunctions_
    class           (intergalacticMediumStateClass           ), intent(inout), target :: intergalacticMediumState_
    !![
    <constructorAssign variables="coefficients,*cosmologyFunctions_,*intergalacticMediumState_"/>
    !!]

    return
  end function metallicityPolynomialIGMConstructorInternal

  subroutine metallicityPolynomialDestructor(self)
    !!{
    Destructor for the metallicityPolynomial \gls{igm} state class.
    !!}
    implicit none
    type(intergalacticMediumStateMetallicityPolynomial), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%intergalacticMediumState_"/>
    !!]
    return
  end subroutine metallicityPolynomialDestructor

  double precision function metallicityPolynomialElectronFraction(self,time)
    !!{
    Return the electron fraction of the \gls{igm}.
    !!}
    implicit none
    class           (intergalacticMediumStateMetallicityPolynomial), intent(inout) :: self
    double precision                                          , intent(in   ) :: time

    metallicityPolynomialElectronFraction=self%intergalacticMediumState_%electronFraction(time)
    return
  end function metallicityPolynomialElectronFraction

  double precision function metallicityPolynomialNeutralHydrogenFraction(self,time)
    !!{
    Return the neutral hydrogen fraction of the \gls{igm}.
    !!}
    implicit none
    class           (intergalacticMediumStateMetallicityPolynomial), intent(inout) :: self
    double precision                                          , intent(in   ) :: time

    metallicityPolynomialNeutralHydrogenFraction=self%intergalacticMediumState_%neutralHydrogenFraction(time)  
    return
  end function metallicityPolynomialNeutralHydrogenFraction

  double precision function metallicityPolynomialNeutralHeliumFraction(self,time)
    !!{
    Return the neutral helium fraction of the \gls{igm}.
    !!}
    implicit none
    class           (intergalacticMediumStateMetallicityPolynomial), intent(inout) :: self
    double precision                                          , intent(in   ) :: time

    metallicityPolynomialNeutralHeliumFraction=self%intergalacticMediumState_%neutralHeliumFraction(time)  
    return
  end function metallicityPolynomialNeutralHeliumFraction

  double precision function metallicityPolynomialSinglyIonizedHeliumFraction(self,time)
    !!{
    Return the singly-ionized helium fraction of the \gls{igm}.
    !!}
    implicit none
    class           (intergalacticMediumStateMetallicityPolynomial), intent(inout) :: self
    double precision                                          , intent(in   ) :: time

    metallicityPolynomialSinglyIonizedHeliumFraction=self%intergalacticMediumState_%singlyIonizedHeliumFraction(time)  
    return
  end function metallicityPolynomialSinglyIonizedHeliumFraction

  double precision function metallicityPolynomialTemperature(self,time)
    !!{
    Return the temperature of the \gls{igm}.
    !!}
    implicit none
    class           (intergalacticMediumStateMetallicityPolynomial), intent(inout) :: self
    double precision                                          , intent(in   ) :: time

    metallicityPolynomialTemperature=self%intergalacticMediumState_%temperature(time)  
    return
  end function metallicityPolynomialTemperature

  double precision function metallicityPolynomialMetallicity(self,time)
    !!{
    Return the metallicity of the \gls{igm}.
    !!}
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    implicit none
    class           (intergalacticMediumStateMetallicityPolynomial), intent(inout) :: self
    double precision                                          , intent(in   ) :: time
    integer :: i
    double precision :: redshift

    redshift =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))


    metallicityPolynomialMetallicity=0.0d0
    do i=1,size(self%coefficients)

    metallicityPolynomialMetallicity=+metallicityPolynomialMetallicity &
  &+self%coefficients(i)*log10(1.0d0+redshift)**(i-1)
    end do
    metallicityPolynomialMetallicity=metallicitySolar*10.0d0**metallicityPolynomialMetallicity

    return
  end function metallicityPolynomialMetallicity
