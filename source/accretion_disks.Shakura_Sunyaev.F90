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
  Implementation of a \cite{shakura_black_1973} accretion disk.
  !!}

  !![
  <accretionDisks name="accretionDisksShakuraSunyaev">
   <description>
    A circumnuclear accretion disk class, in which the accretion disks are always described by a radiatively efficient,
    geometrically thin accretion disk as described by \cite{shakura_black_1973}. The radiative efficiency of the flow is
    computed assuming that material falls into the black hole without further energy loss from the \gls{isco}, while the
    spin-up rate of the black hole is computed assuming that the material enters the black hole with the specific angular
    momentum of the \gls{isco} (i.e. there are no torques on the material once it begins to fall in from the \gls{isco};
    \citealt{bardeen_kerr_1970}). For these thin disks, jet power is computed, using the expressions from
    \citeauthor{meier_association_2001}~(\citeyear{meier_association_2001}; his equations 4 and 5).
   </description>
  </accretionDisks>
  !!]
  type, extends(accretionDisksClass) :: accretionDisksShakuraSunyaev
     !!{
     Implementation of a \cite{shakura_black_1973} accretion disk class.
     !!}
     private
   contains
     procedure :: efficiencyRadiative => shakuraSunyaevEfficiencyRadiative
     procedure :: powerJet            => shakuraSunyaevPowerJet
     procedure :: rateSpinUp          => shakuraSunyaevRateSpinUp
  end type accretionDisksShakuraSunyaev

  interface accretionDisksShakuraSunyaev
     !!{
     Constructors for the Eddington-limited accretion disk class.
     !!}
     module procedure shakuraSunyaevConstructorParameters
  end interface accretionDisksShakuraSunyaev

contains

  function shakuraSunyaevConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{shakura_black_1973} accretion disk class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(accretionDisksShakuraSunyaev)                :: self
    type(inputParameters             ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=accretionDisksShakuraSunyaev()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function shakuraSunyaevConstructorParameters

  double precision function shakuraSunyaevEfficiencyRadiative(self,blackHole,accretionRateMass)
    !!{
    Return the radiative efficiency of a \cite{shakura_black_1973} accretion disk.
    !!}
    use :: Black_Hole_Fundamentals, only : Black_Hole_ISCO_Specific_Energy, orbitPrograde, unitsGravitational
    implicit none
    class           (accretionDisksShakuraSunyaev), intent(inout) :: self
    class           (nodeComponentBlackHole      ), intent(inout) :: blackHole
    double precision                              , intent(in   ) :: accretionRateMass
    !$GLC attributes unused :: self, accretionRateMass

    shakuraSunyaevEfficiencyRadiative=+1.0d0                                                                                   &
         &                            -Black_Hole_ISCO_Specific_Energy(blackHole,units=unitsGravitational,orbit=orbitPrograde)
    return
  end function shakuraSunyaevEfficiencyRadiative

  double precision function shakuraSunyaevPowerJet(self,blackHole,accretionRateMass)
    !!{
    Computes the jet power for a \cite{shakura_black_1973} (thin) accretion disk, using the expressions from
    \citeauthor{meier_association_2001}~(\citeyear{meier_association_2001}; his equations 4 and 5).
    !!}
    use :: Black_Hole_Fundamentals         , only : Black_Hole_Eddington_Accretion_Rate
    use :: Numerical_Constants_Astronomical, only : gigaYear                           , massSolar
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Constants_Units       , only : ergs
    implicit none
    class           (accretionDisksShakuraSunyaev), intent(inout) :: self
    class           (nodeComponentBlackHole      ), intent(inout) :: blackHole
    double precision                              , intent(in   ) :: accretionRateMass
    double precision                              , parameter     :: alphaViscosity                =0.01d0
    double precision                              , parameter     :: alphaViscosityNormalized      =alphaViscosity/0.01d0
    double precision                              , parameter     :: powerNormalizationKerr        =(10.0d0**42.7)*ergs*gigaYear/massSolar/kilo**2
    double precision                              , parameter     :: powerNormalizationSchwarzchild=(10.0d0**41.7)*ergs*gigaYear/massSolar/kilo**2
    double precision                              , parameter     :: meierMassNormalization        =1.0d9
    double precision                                              :: accretionRateDimensionless                                                   , massBlackHoleDimensionless, &
         &                                                           spinBlackHole
    !$GLC attributes unused :: self

    ! Return immediately for non-positive accretion rates.
    if (accretionRateMass <= 0.0d0) then
       shakuraSunyaevPowerJet=0.0d0
    else
       ! Get the black hole spin and dimensionless accretion rate and mass as defined by Meier (2001).
       spinBlackHole             =+                                    blackHole%spin()
       accretionRateDimensionless=+accretionRateMass                                     &
            &                     /Black_Hole_Eddington_Accretion_Rate(blackHole       )
       massBlackHoleDimensionless=+                                    blackHole%mass()  &
            &                     /meierMassNormalization
       if (massBlackHoleDimensionless > 0.0d0 .and. accretionRateDimensionless > 0.0d0) then
          if (spinBlackHole > 0.8d0) then
             ! Use Meier's rapidly rotating (Kerr) solution for high spin black holes.
             shakuraSunyaevPowerJet=+powerNormalizationKerr            &
                  &                 *massBlackHoleDimensionless**0.9d0 &
                  &                 *accretionRateDimensionless**1.2d0 &
                  &                 /alphaViscosityNormalized  **0.1d0 &
                  &                 *(                                 &
                  &                   +1.00d0                          &
                  &                   +1.10d0*spinBlackHole            &
                  &                   +0.29d0*spinBlackHole    **2     &
                  &                  )
          else
             ! Use Meier's Schwarzchild solution for low spin black holes. We, somewhat arbitrarily, interpolate from the
             ! Schwarzchild solution assuming power grows exponentially with spin and matched to the Kerr solution at the
             ! transition spin.
             shakuraSunyaevPowerJet=+powerNormalizationSchwarzchild        &
                  &                 *massBlackHoleDimensionless    **0.9d0 &
                  &                 *accretionRateDimensionless    **1.2d0 &
                  &                 /alphaViscosityNormalized      **0.1d0 &
                  &                 *exp(                                  &
                  &                      +3.785d0                          &
                  &                      *spinBlackHole                    &
                  &                     )
          end if
       else
          shakuraSunyaevPowerJet=0.0d0
       end if
    end if
    return
  end function shakuraSunyaevPowerJet

  double precision function shakuraSunyaevRateSpinUp(self,blackHole,accretionRateMass) result(rateSpinUp)
    !!{
    Compute the rate of spin up of a black hole by a \cite{shakura_black_1973} accretion disk.
    !!}
    use :: Black_Hole_Fundamentals, only : Black_Hole_ISCO_Specific_Angular_Momentum, Black_Hole_ISCO_Specific_Energy, orbitPrograde, unitsGravitational
    implicit none
    class           (accretionDisksShakuraSunyaev), intent(inout) :: self
    class           (nodeComponentBlackHole      ), intent(inout) :: blackHole
    double precision                              , intent(in   ) :: accretionRateMass
    double precision                                              :: spinToMassRateOfChangeRatio
    !$GLC attributes unused :: self

    if (accretionRateMass /= 0.0d0) then
       spinToMassRateOfChangeRatio=+Black_Hole_ISCO_Specific_Angular_Momentum(blackHole       ,units=unitsGravitational,orbit=orbitPrograde) &
            &                      -2.0d0                                                                                                    &
            &                      *                                          blackHole%spin()                                               &
            &                      *Black_Hole_ISCO_Specific_Energy          (blackHole       ,units=unitsGravitational,orbit=orbitPrograde)
       rateSpinUp                 =+spinToMassRateOfChangeRatio &
            &                      *accretionRateMass           &
            &                      /blackHole%mass()
    else
       rateSpinUp                 =+0.0d0
    end if
    return
  end function shakuraSunyaevRateSpinUp
