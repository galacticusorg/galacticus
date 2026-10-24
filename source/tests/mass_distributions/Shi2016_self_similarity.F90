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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program which tests that the :cite:t:`shi_outer_2016` accretion flow mass distribution is self-similar.
!!}

program Test_Mass_Distributions_Shi2016_Self_Similarity
  !!{RST
  Tests that the :galacticus-class:`massDistributionShi2016` accretion flow is self-similar. In a flat cosmology with a
  cosmological constant the scale-free solution of :cite:t:`shi_outer_2016` depends only on :math:`\Omega_\mathrm{m}(t)` at
  the epoch of observation and on the growth index :math:`s = \mathrm{d}\ln M/\mathrm{d}\ln a`, so that :math:`\rho(r) r_\mathrm{v}^3/M_\mathrm{v}`
  as a function of :math:`r/r_\mathrm{v}` must be the same for two halos with the same :math:`\Omega_\mathrm{m}(t)` and
  :math:`s`, whatever their redshift. A halo at :math:`z=0` in a cosmology with :math:`\Omega_\mathrm{m,0}=0.3` is compared
  with one at :math:`z=1` in a cosmology with :math:`\Omega_\mathrm{m,0}=0.3/5.9`, for which :math:`\Omega_\mathrm{m}(z=1)=0.3`.
  !!}
  use :: Coordinates             , only : coordinateSpherical           , assignment(=)
  use :: Cosmology_Functions     , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters    , only : cosmologyParametersSimple
  use :: Display                 , only : displayVerbositySet           , verbosityLevelStandard
  use :: Events_Hooks            , only : eventsHooksInitialize
  use :: Mass_Distributions      , only : massDistributionShi2016
  use :: Numerical_Constants_Math, only : Pi
  use :: Unit_Tests              , only : Assert                        , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (cosmologyParametersSimple     ), pointer               :: cosmologyParametersLate      , cosmologyParametersEarly
  type            (cosmologyFunctionsMatterLambda), pointer               :: cosmologyFunctionsLate       , cosmologyFunctionsEarly
  integer                                         , parameter             :: countRadii             =4
  double precision                                , dimension(countRadii) :: radiusScaleFree        =[1.5d0,2.0d0,3.0d0,5.0d0]
  double precision                                , dimension(countRadii) :: densityLate                  , densityEarly
  ! Halo mass [M☉], growth index, ratio of turnaround to virial radii, and density contrast (relative to the mean) defining the
  ! virial radius. The growth index avoids s=3/2, at which the initial estimate of the splashback radius changes form, so that
  ! rounding cannot place the two halos on different sides of that boundary.
  double precision                                , parameter             :: mass                   =1.0d12, growthIndex         =1.2d0, &
       &                                                                     ratioRadiusTurnaround  =2.0d0 , densityContrast     =2.0d2
  ! Present day matter density parameters of the two (flat, radiation-free) cosmologies, and the expansion factors at which the
  ! halos are observed.
  double precision                                , parameter             :: OmegaMatterLate        =0.3d0 , OmegaMatterEarly    =0.3d0/5.9d0, &
       &                                                                     expansionFactorLate    =1.0d0 , expansionFactorEarly=0.5d0

  call displayVerbositySet  (verbosityLevelStandard)
  call eventsHooksInitialize(                      )
  call Unit_Tests_Begin_Group("Shi (2016) accretion flow self-similarity")
  allocate(cosmologyParametersLate )
  allocate(cosmologyParametersEarly)
  allocate(cosmologyFunctionsLate  )
  allocate(cosmologyFunctionsEarly )
  !![
  <referenceConstruct object="cosmologyParametersLate"  constructor="cosmologyParametersSimple     (OmegaMatter=OmegaMatterLate ,OmegaBaryon=0.0d0,OmegaDarkEnergy=1.0d0-OmegaMatterLate ,temperatureCMB=0.0d0,HubbleConstant=70.0d0)"/>
  <referenceConstruct object="cosmologyParametersEarly" constructor="cosmologyParametersSimple     (OmegaMatter=OmegaMatterEarly,OmegaBaryon=0.0d0,OmegaDarkEnergy=1.0d0-OmegaMatterEarly,temperatureCMB=0.0d0,HubbleConstant=70.0d0)"/>
  <referenceConstruct object="cosmologyFunctionsLate"   constructor="cosmologyFunctionsMatterLambda(cosmologyParametersLate )"/>
  <referenceConstruct object="cosmologyFunctionsEarly"  constructor="cosmologyFunctionsMatterLambda(cosmologyParametersEarly)"/>
  !!]
  call Assert('Ωₘ(t) matches',cosmologyFunctionsEarly%OmegaMatterEpochal(expansionFactor=expansionFactorEarly),cosmologyFunctionsLate%OmegaMatterEpochal(expansionFactor=expansionFactorLate),relTol=1.0d-6)
  call evaluate(cosmologyFunctionsLate ,expansionFactorLate ,densityLate )
  call evaluate(cosmologyFunctionsEarly,expansionFactorEarly,densityEarly)
  call Assert('ρ(r) r_v³/M_v at z=1 equals that at z=0',densityEarly,densityLate,relTol=1.0d-4)
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  !![
  <objectDestructor name="cosmologyFunctionsLate"  />
  <objectDestructor name="cosmologyFunctionsEarly" />
  <objectDestructor name="cosmologyParametersLate" />
  <objectDestructor name="cosmologyParametersEarly"/>
  !!]

contains

  subroutine evaluate(cosmologyFunctions_,expansionFactor,density)
    !!{RST
    Construct the accretion flow for a halo observed at the given ``expansionFactor``, and evaluate its scaled density,
    :math:`\rho(r) r_\mathrm{v}^3/M_\mathrm{v}`, at the radii being tested.
    !!}
    implicit none
    type            (cosmologyFunctionsMatterLambda), intent(inout)                :: cosmologyFunctions_
    double precision                                , intent(in   )                :: expansionFactor
    double precision                                , intent(  out), dimension(:)  :: density
    type            (massDistributionShi2016       ), pointer                      :: massDistribution_
    type            (coordinateSpherical           )                               :: coordinates
    double precision                                                               :: time             , radiusVirial, &
         &                                                                            massAccretionRate
    integer                                                                        :: i

    time             =+cosmologyFunctions_%cosmicTime          (expansionFactor)
    massAccretionRate=+growthIndex                                                &
         &            *cosmologyFunctions_%expansionRate       (expansionFactor) &
         &            *mass
    radiusVirial     =+(                                                          &
         &              +3.0d0                                                    &
         &              *mass                                                     &
         &              /4.0d0                                                    &
         &              /Pi                                                       &
         &              /densityContrast                                          &
         &              /cosmologyFunctions_%matterDensityEpochal(time           ) &
         &             )**(1.0d0/3.0d0)
    allocate(massDistribution_)
    !![
    <referenceConstruct object="massDistribution_" constructor="massDistributionShi2016(mass,massAccretionRate,radiusVirial,ratioRadiusTurnaround,time,1.0d0,cosmologyFunctions_)"/>
    !!]
    do i=1,size(radiusScaleFree)
       coordinates=[radiusScaleFree(i)*radiusVirial,0.0d0,0.0d0]
       density(i) =+massDistribution_%density(coordinates) &
            &      *radiusVirial**3                       &
            &      /mass
    end do
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end subroutine evaluate

end program Test_Mass_Distributions_Shi2016_Self_Similarity
