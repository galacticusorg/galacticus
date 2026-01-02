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

  !!{
  An implementation of Gaunt factors using the \cite{van_hoof_accurate_2014} fitting function.
  !!}

  use :: Atomic_Ionization_Potentials, only : atomicIonizationPotentialClass

  !![
  <gauntFactor name="gauntFactorVanHoof2014">
   <description>Gaunt factors are computed using the fitting function of \cite{van_hoof_accurate_2014}.</description>
  </gauntFactor>
  !!]
  type, extends(gauntFactorClass) :: gauntFactorVanHoof2014
     !!{
     A gaunt factor class implementing the fitting function of \cite{van_hoof_accurate_2014}.
     !!}
     private
     class(atomicIonizationPotentialClass), pointer :: atomicIonizationPotential_ => null()
   contains
     final     ::          vanHoof2014Destructor
     procedure :: total => vanHoof2014Total
  end type gauntFactorVanHoof2014

  interface gauntFactorVanHoof2014
     !!{
     Constructors for the \refClass{gauntFactorVanHoof2014} gaunt factor class.
     !!}
     module procedure vanHoof2014ConstructorParameters
     module procedure vanHoof2014ConstructorInternal
  end interface gauntFactorVanHoof2014

  ! Arrays to hold coefficients of the fitting function.
  double precision, dimension(0:4,2) :: fitA=reshape(                       &
       &                                             [                      &
       &                                              +1.43251926625281d+0, &
       &                                              +3.50626935257777d-1, &
       &                                              +4.36183448595035d-1, &
       &                                              +6.03536387105599d-2, &
       &                                              +3.66626405363100d-2, &
       &                                              +1.45481634667278d+0, &
       &                                              -9.55399384620923d-2, &
       &                                              +1.46327814151538d-1, &
       &                                              -1.41489406498468d-2, &
       &                                              +2.76891413242655d-3  &
       &                                             ]                    , &
       &                                             [5,2]                  &
       &                                            )
  double precision, dimension(0:4,2) :: fitB=reshape(                       &
       &                                             [                      &
       &                                              +1.00000000000000d+0, &
       &                                              +2.92525161994346d-1, &
       &                                              +4.05566949766954d-1, &
       &                                              +5.62573012783879d-2, &
       &                                              +3.33019373823972d-2, &
       &                                              +1.00000000000000d+0, &
       &                                              +3.31149751183539d-2, &
       &                                              +1.31127367293310d-1, &
       &                                              -1.32658217746618d-2, &
       &                                              +2.74809263365693d-3  &
       &                                             ]                    , &
       &                                             [5,2]                  &
       &                                            )

contains

  function vanHoof2014ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{gauntFactorVanHoof2014} gaunt factor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (gauntFactorVanHoof2014        )                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(atomicIonizationPotentialClass), pointer       :: atomicIonizationPotential_

    !![
    <objectBuilder class="atomicIonizationPotential" name="atomicIonizationPotential_" source="parameters"/>
    !!]
    self=gauntFactorVanHoof2014(atomicIonizationPotential_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="atomicIonizationPotential_"/>
    !!]
    return
  end function vanHoof2014ConstructorParameters

  function vanHoof2014ConstructorInternal(atomicIonizationPotential_) result(self)
    !!{
    Internal constructor for the \refClass{gauntFactorVanHoof2014} gaunt factor class.
    !!}
    implicit none
    type (gauntFactorVanHoof2014        )                        :: self
    class(atomicIonizationPotentialClass), intent(in   ), target :: atomicIonizationPotential_
    !![
    <constructorAssign variables="*atomicIonizationPotential_"/>
    !!]

    return
  end function vanHoof2014ConstructorInternal

  subroutine vanHoof2014Destructor(self)
    !!{
    Destructor for the \refClass{gauntFactorVanHoof2014} gaunt factor class.
    !!}
    implicit none
    type(gauntFactorVanHoof2014), intent(inout) :: self

    !![
    <objectDestructor name="self%atomicIonizationPotential_"/>
    !!]
    return
  end subroutine vanHoof2014Destructor

  double precision function vanHoof2014Total(self,atomicNumber,electronNumber,temperature)
    !!{
    Compute thermally averaged Gaunt factors for thermal electron distributions using the tabulations and fits of
    \cite{van_hoof_accurate_2014}.
    !!}
    use :: Error                       , only : Error_Report
    use :: Numerical_Constants_Physical, only : boltzmannsConstant
    use :: Numerical_Constants_Units   , only : rydberg
    implicit none
    class           (gauntFactorVanHoof2014), intent(inout) :: self
    integer                                 , intent(in   ) :: atomicNumber, electronNumber
    double precision                        , intent(in   ) :: temperature
    double precision                                        :: gammaSquared, g
    integer                                                 :: i

    ! Return zero for unphysical temperatures
    if (temperature <= 0.0d0) then
       vanHoof2014Total=0.0d0
       return
    end if
    ! Validate input.
    if (electronNumber > atomicNumber) call Error_Report('number of electrons exceeds atomic number'//{introspection:location})
    ! Return zero if ionization potential is not available for this ion.
    if (self%atomicIonizationPotential_%potential(atomicNumber,electronNumber) == 0.0d0) then
       vanHoof2014Total=0.0d0
    else
       ! Evaluate the γ parameter.
       gammaSquared=+dble(                &
            &             +atomicNumber   &
            &             -electronNumber &
            &             +1              &
            &            )**2             &
            &       *rydberg              &
            &       /boltzmannsConstant   &
            &       /temperature
       g           =log10(gammaSquared)
       ! Evaluate the Gaunt factor.
       if (g <= 0.8d0) then
          i=1
       else
          i=2
       end if
       vanHoof2014Total=+(                &
            &             +fitA(0,i)      &
            &             +fitA(1,i)*g    &
            &             +fitA(2,i)*g**2 &
            &             +fitA(3,i)*g**3 &
            &             +fitA(4,i)*g**4 &
            &            )                &
            &            /(               &
            &             +fitB(0,i)      &
            &             +fitB(1,i)*g    &
            &             +fitB(2,i)*g**2 &
            &             +fitB(3,i)*g**3 &
            &             +fitB(4,i)*g**4 &
            &            )
    end if
    return
  end function vanHoof2014Total
