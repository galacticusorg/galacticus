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

  !!{
  Implementation of an Eddington-limited accretion disk.
  !!}

  !![
  <accretionDisks name="accretionDisksEddingtonLimited">
   <description>
    A circumnuclear accretion disk class in which accretion is always Eddington-limited. This class does not assume any
    physical model for the accretion disk, but merely assumes that jets are powered at a fixed fraction {\normalfont \ttfamily
    [efficiencyJet]} of the Eddington luminosity. The radiative efficiency is similarly set at a fixed value of {\normalfont
    \ttfamily [efficiencyRadiation]}. Since no physical model for the disk is assumed, the black hole spin up rate is always
    set to zero.
   </description>
  </accretionDisks>
  !!]
  type, extends(accretionDisksClass) :: accretionDisksEddingtonLimited
     !!{
     Implementation of an accretion disk class in which accretion is always Eddington-limited.
     !!}
     private
     double precision :: efficiencyRadiation, efficiencyJet
   contains
     procedure :: efficiencyRadiative => eddingtonLimitedEfficiencyRadiative
     procedure :: powerJet            => eddingtonLimitedPowerJet
     procedure :: rateSpinUp          => eddingtonLimitedRateSpinUp
  end type accretionDisksEddingtonLimited

  interface accretionDisksEddingtonLimited
     !!{
     Constructors for the Eddington-limited accretion disk class.
     !!}
     module procedure eddingtonLimitedConstructorParameters
     module procedure eddingtonLimitedConstructorInternal
  end interface accretionDisksEddingtonLimited

contains

  function eddingtonLimitedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the Eddington-limited accretion disk class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (accretionDisksEddingtonLimited)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    double precision                                                :: efficiencyRadiation, efficiencyJet

    !![
    <inputParameter>
      <name>efficiencyRadiation</name>
      <source>parameters</source>
      <defaultValue>0.1d0</defaultValue>
      <description>The radiative efficiency of the Eddington-limited accretion disk.</description>
    </inputParameter>
    <inputParameter>
      <name>efficiencyJet</name>
      <source>parameters</source>
      <defaultValue>0.1d0</defaultValue>
      <description>The jet efficiency of the Eddington-limited accretion disk.</description>
    </inputParameter>
    !!]
    self=accretionDisksEddingtonLimited(efficiencyRadiation,efficiencyJet)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function eddingtonLimitedConstructorParameters

  function eddingtonLimitedConstructorInternal(efficiencyRadiation,efficiencyJet) result(self)
    !!{
    Internal constructor for the Eddington-limited accretion disk class.
    !!}
    implicit none
    type            (accretionDisksEddingtonLimited)                :: self
    double precision                                , intent(in   ) :: efficiencyRadiation, efficiencyJet

    !![
    <constructorAssign variables="efficiencyRadiation, efficiencyJet"/>
    !!]
    return
  end function eddingtonLimitedConstructorInternal

  double precision function eddingtonLimitedEfficiencyRadiative(self,blackHole,accretionRateMass)
    !!{
    Return the radiative efficiency of an Eddington-limited accretion disk.
    !!}
    implicit none
    class           (accretionDisksEddingtonLimited), intent(inout) :: self
    class           (nodeComponentBlackHole        ), intent(inout) :: blackHole
    double precision                                , intent(in   ) :: accretionRateMass
    !$GLC attributes unused :: blackHole, accretionRateMass

    eddingtonLimitedEfficiencyRadiative=self%efficiencyRadiation
    return
  end function eddingtonLimitedEfficiencyRadiative

  double precision function eddingtonLimitedPowerJet(self,blackHole,accretionRateMass) result(powerJet)
    !!{
    Return the jet power of an Eddington-limited accretion disk.
    !!}
    use :: Black_Hole_Fundamentals     , only : Black_Hole_Eddington_Accretion_Rate
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    class           (accretionDisksEddingtonLimited), intent(inout) :: self
    class           (nodeComponentBlackHole        ), intent(inout) :: blackHole
    double precision                                , intent(in   ) :: accretionRateMass
    double precision                                , parameter     :: rateFractionalCutOff                =0.01d0
    double precision                                                :: accretionRateEddingtonLimited              , &
         &                                                             accretionRateEddingtonLimitedReduced       , &
         &                                                             rateFractional

    if (accretionRateMass <= 0.0d0) then
       powerJet=0.0d0
       return
    end if
    accretionRateEddingtonLimited=+self%efficiencyJet                             &
         &                        *Black_Hole_Eddington_Accretion_Rate(blackHole)
    if (accretionRateMass > rateFractionalCutOff*accretionRateEddingtonLimited) then
       accretionRateEddingtonLimitedReduced=accretionRateEddingtonLimited
    else
       rateFractional                      =+accretionRateMass             &
            &                               /rateFractionalCutOff          &
            &                               /accretionRateEddingtonLimited
       accretionRateEddingtonLimitedReduced=+accretionRateEddingtonLimited &
            &                               *(                             &
            &                                 +3.0d0*rateFractional**2     &
            &                                 -2.0d0*rateFractional**3     &
            &                                )
    end if
    powerJet=+accretionRateEddingtonLimitedReduced &
         &   *(                                    &
         &     +speedLight                         &
         &     /kilo                               &
         &     )**2
    return
  end function eddingtonLimitedPowerJet

  double precision function eddingtonLimitedRateSpinUp(self,blackHole,accretionRateMass)
    !!{
    Computes the spin up rate of the given {\normalfont \ttfamily blackHole} due to accretion from an Eddington-limited
    accretion disk. This is always zero, as no physical model is specified for this accretion disk method.
    !!}
    implicit none
    class           (accretionDisksEddingtonLimited), intent(inout) :: self
    class           (nodeComponentBlackHole        ), intent(inout) :: blackHole
    double precision                                , intent(in   ) :: accretionRateMass
    !$GLC attributes unused :: self, blackHole, accretionRateMass

    eddingtonLimitedRateSpinUp=0.0d0
    return
  end function eddingtonLimitedRateSpinUp
