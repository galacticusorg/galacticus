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

  !!{RST
  Implements a supernovae type Ia class which scales the rate of another class by a power of the metallicity.
  !!}

  !![
  <supernovaeTypeIa name="supernovaeTypeIaMetallicityDependentRate" docformat="rst">
   <description>
   A supernovae type Ia class which takes the rate of Type Ia supernovae from another ``supernovaeTypeIa`` class and multiplies it by a power of the metallicity,

   .. math::

      f(Z) = \left( \frac{\max(Z,Z_\mathrm{min})}{\mathrm{Z}_\odot} \right)^\beta,

   where :math:`\beta=` ``[exponent]`` and :math:`Z_\mathrm{min}=` ``[metallicityMinimum]``. Both the number of events and the yields are scaled, the latter being proportional to the former.

   :cite:t:`johnson_binaries_2023` argue that the anti-correlation of the close binary fraction with metallicity should enhance the Type Ia supernova rate at low metallicity, and find that a scaling of approximately :math:`Z^{-1/2}` on top of the differing star formation histories accounts for the high specific rates observed in dwarf galaxies. That is the default adopted here.

   The scaling is unity at :math:`Z=\mathrm{Z}_\odot`, so the normalization of the wrapped class retains its usual meaning as the rate at Solar metallicity. Because the scaling diverges as :math:`Z \rightarrow 0` for negative :math:`\beta`, the metallicity is floored at ``[metallicityMinimum]``; with the default values this caps the enhancement at a factor of ten. The floor should be chosen with the range over which the adopted scaling was calibrated in mind---it is not itself a physical statement.

   Separating the metallicity dependence of the rate from the shape of the delay time distribution in this way means it can be applied to any of the other ``supernovaeTypeIa`` classes.
   </description>
  </supernovaeTypeIa>
  !!]
  type, extends(supernovaeTypeIaClass) :: supernovaeTypeIaMetallicityDependentRate
     !!{RST
     A supernovae type Ia class which scales the rate of another class by a power of the metallicity.
     !!}
     private
     class           (supernovaeTypeIaClass), pointer :: supernovaeTypeIa_ => null()
     double precision                                 :: exponent                   , metallicityMinimum
   contains
     final     ::                     metallicityDependentRateDestructor
     procedure :: massInitialRange => metallicityDependentRateMassInitialRange
     procedure :: number           => metallicityDependentRateNumber
     procedure :: yield            => metallicityDependentRateYield
  end type supernovaeTypeIaMetallicityDependentRate

  interface supernovaeTypeIaMetallicityDependentRate
     !!{RST
     Constructors for the :galacticus-class:`supernovaeTypeIaMetallicityDependentRate` supernovae type Ia class.
     !!}
     module procedure metallicityDependentRateConstructorParameters
     module procedure metallicityDependentRateConstructorInternal
  end interface supernovaeTypeIaMetallicityDependentRate

contains

  function metallicityDependentRateConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`supernovaeTypeIaMetallicityDependentRate` supernovae type Ia class which takes a parameter list as input.
    !!}
    use :: Input_Parameters                , only : inputParameter, inputParameters
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    implicit none
    type            (supernovaeTypeIaMetallicityDependentRate)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (supernovaeTypeIaClass                   ), pointer       :: supernovaeTypeIa_
    double precision                                                          :: exponent         , metallicityMinimum

    !![
    <inputParameter docformat="rst">
      <name>exponent</name>
      <source>parameters</source>
      <defaultValue>-0.5d0</defaultValue>
      <defaultSource>
      :cite:p:`johnson_binaries_2023`
      </defaultSource>
      <description>
      The exponent :math:`\beta` in the scaling of the Type Ia supernova rate with metallicity, :math:`f(Z) \propto Z^\beta`; the default value of :math:`-1/2` follows :cite:t:`johnson_binaries_2023`.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>metallicityMinimum</name>
      <source>parameters</source>
      <defaultValue>1.0d-2*metallicitySolar</defaultValue>
      <description>
      The metallicity below which the scaling of the Type Ia supernova rate with metallicity is held constant. This is required because the scaling diverges as :math:`Z \rightarrow 0` for negative ``[exponent]``; with the default values it caps the enhancement of the rate at a factor of ten.
      </description>
    </inputParameter>
    <objectBuilder class="supernovaeTypeIa" name="supernovaeTypeIa_" source="parameters"/>
    !!]
    self=supernovaeTypeIaMetallicityDependentRate(exponent,metallicityMinimum,supernovaeTypeIa_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="supernovaeTypeIa_"/>
    !!]
    return
  end function metallicityDependentRateConstructorParameters

  function metallicityDependentRateConstructorInternal(exponent,metallicityMinimum,supernovaeTypeIa_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`supernovaeTypeIaMetallicityDependentRate` supernovae type Ia class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (supernovaeTypeIaMetallicityDependentRate)                        :: self
    class           (supernovaeTypeIaClass                   ), intent(in   ), target :: supernovaeTypeIa_
    double precision                                          , intent(in   )         :: exponent         , metallicityMinimum
    !![
    <constructorAssign variables="exponent, metallicityMinimum, *supernovaeTypeIa_"/>
    !!]

    if (metallicityMinimum <= 0.0d0) call Error_Report('[metallicityMinimum] must be positive, as the scaling diverges at zero metallicity for negative [exponent]'//{introspection:location})
    return
  end function metallicityDependentRateConstructorInternal

  subroutine metallicityDependentRateDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`supernovaeTypeIaMetallicityDependentRate` supernovae type Ia class.
    !!}
    implicit none
    type(supernovaeTypeIaMetallicityDependentRate), intent(inout) :: self

    !![
    <objectDestructor name="self%supernovaeTypeIa_"/>
    !!]
    return
  end subroutine metallicityDependentRateDestructor

  double precision function metallicityDependentRateFactor(self,metallicity) result(factor)
    !!{RST
    Return the factor by which the Type Ia supernova rate is scaled at the given ``metallicity``.
    !!}
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    implicit none
    class           (supernovaeTypeIaMetallicityDependentRate), intent(inout) :: self
    double precision                                          , intent(in   ) :: metallicity

    factor=(max(metallicity,self%metallicityMinimum)/metallicitySolar)**self%exponent
    return
  end function metallicityDependentRateFactor

  subroutine metallicityDependentRateMassInitialRange(self,initialMassFunction_,age,metallicity,massInitialMinimum,massInitialMaximum)
    !!{RST
    Return the range of initial stellar masses contributing to the Type Ia population. Scaling the rate does not change which stars contribute, so this simply defers to the wrapped class.
    !!}
    implicit none
    class           (supernovaeTypeIaMetallicityDependentRate), intent(inout) :: self
    class           (initialMassFunctionClass                ), intent(inout) :: initialMassFunction_
    double precision                                          , intent(in   ) :: age                 , metallicity
    double precision                                          , intent(  out) :: massInitialMinimum  , massInitialMaximum

    call self%supernovaeTypeIa_%massInitialRange(initialMassFunction_,age,metallicity,massInitialMinimum,massInitialMaximum)
    return
  end subroutine metallicityDependentRateMassInitialRange

  double precision function metallicityDependentRateNumber(self,initialMassFunction_,initialMass,age,metallicity) result(number)
    !!{RST
    Compute the cumulative number of Type Ia supernovae, scaled by a power of the metallicity.
    !!}
    implicit none
    class           (supernovaeTypeIaMetallicityDependentRate), intent(inout), target :: self
    class           (initialMassFunctionClass                ), intent(inout), target :: initialMassFunction_
    double precision                                          , intent(in   )         :: age                 , initialMass, &
         &                                                                               metallicity

    number=+metallicityDependentRateFactor(self,metallicity)                                          &
         & *self%supernovaeTypeIa_%number  (initialMassFunction_,initialMass,age,metallicity)
    return
  end function metallicityDependentRateNumber

  double precision function metallicityDependentRateYield(self,initialMassFunction_,initialMass,age,metallicity,atomIndex) result(yield)
    !!{RST
    Compute the cumulative yield from Type Ia supernovae, scaled by a power of the metallicity. The yield is proportional to the number of events, so it is scaled by the same factor.
    !!}
    implicit none
    class           (supernovaeTypeIaMetallicityDependentRate), intent(inout)           :: self
    class           (initialMassFunctionClass                ), intent(inout)           :: initialMassFunction_
    double precision                                          , intent(in   )           :: age                 , initialMass, &
         &                                                                                 metallicity
    integer                                                   , intent(in   ), optional :: atomIndex

    yield=+metallicityDependentRateFactor(self,metallicity)                                                     &
         & *self%supernovaeTypeIa_%yield  (initialMassFunction_,initialMass,age,metallicity,atomIndex)
    return
  end function metallicityDependentRateYield
