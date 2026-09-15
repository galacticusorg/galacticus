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
Contains a program which tests black hole physics against independently computed reference values.
!!}

!+ Contributions to this file made by: Andrew Benson, Claude.

program Test_Black_Hole_Physics
  !!{RST
  Tests black hole physics (Kerr :term:`ISCO` properties and frame-dragging frequency, thin disk radiative efficiency, spin-up
  rate, and jet power, the
  Eddington accretion rate, Bondi-Hoyle-Lyttleton accretion, and binary merger remnant spins) against reference values computed
  independently from the cited papers and standard results by the ``blackHolePhysics.py`` script in the `galacticusDevTools
  <https://github.com/galacticusorg/galacticusDevTools>`_ repository. Tolerances are:

  * dimensionless closed-form quantities (:term:`ISCO` properties for :math:`j<0.99999`, merger remnant spins): :math:`10^{-8}`
    relative, limited only by double precision arithmetic;
  * :term:`ISCO` properties for :math:`j \ge 0.99999`, where Galacticus uses a first-order series in :math:`(1-j)^{1/3}`:
    :math:`10^{-3}` relative, the truncation error of that series being :math:`5\times 10^{-4}` at :math:`j=0.999995`;
  * dimensional quantities: the reference values use CODATA 2018 and IAU 2015 constants (via astropy), while Galacticus uses GSL
    constants, for which :math:`\mathrm{G}\mathrm{M}_\odot` differs by :math:`7\times 10^{-5}`, :math:`\mathrm{M}_\odot` by
    :math:`2.6\times 10^{-4}`, and the year definitions by :math:`2\times 10^{-5}`. Additionally, Galacticus uses the mass of the
    hydrogen atom in the Eddington accretion rate, where the reference uses the proton mass (a :math:`5.4\times 10^{-4}`
    difference). Tolerances of :math:`10^{-5}` (sound speed) to :math:`2\times 10^{-3}` (jet power, which scales as the Eddington
    rate to the power :math:`-1.2`) follow from these differences.
  !!}
  use :: Accretion_Disks                , only : accretionDisksShakuraSunyaev
  use :: Black_Hole_Binary_Mergers      , only : blackHoleBinaryMergerRezzolla2008
  use :: Black_Hole_Fundamentals        , only : Black_Hole_Eddington_Accretion_Rate , Black_Hole_ISCO_Radius                , Black_Hole_ISCO_Specific_Angular_Momentum, Black_Hole_ISCO_Specific_Energy, &
          &                                      Black_Hole_Frame_Dragging_Frequency , orbitPrograde                         , unitsGravitational
  use :: Bondi_Hoyle_Lyttleton_Accretion, only : Bondi_Hoyle_Lyttleton_Accretion_Rate, Bondi_Hoyle_Lyttleton_Accretion_Radius
  use :: Display                        , only : displayVerbositySet                 , verbosityLevelStandard
  use :: Galacticus_Nodes               , only : nodeComponentBlackHoleStandard
  use :: Ideal_Gases_Thermodynamics     , only : Ideal_Gas_Jeans_Length              , Ideal_Gas_Sound_Speed
  use :: Unit_Tests                     , only : Assert                              , Unit_Tests_Begin_Group                , Unit_Tests_End_Group                     , Unit_Tests_Finish
  implicit none
  ! Kerr ISCO properties (prograde orbits, gravitational units) from Bardeen, Press & Teukolsky (1972). The first six spins are
  ! evaluated using the closed-form expressions in Galacticus, the final two using its near-extremal series.
  double precision                                   , dimension(8) :: spinISCO                  =[                                                                        &
       &                                                                                           0.000000d+0     , 0.300000d+0     , 0.600000d+0     , 0.900000d+0     , &
       &                                                                                           0.990000d+0     , 0.999000d+0     , 0.999995d+0     , 1.000000d+0       &
       &                                                                                          ]
  double precision                                   , dimension(8) :: radiusISCOReference       =[                                                                        &
       &                                                                                           6.0000000000d+0 , 4.9786168306d+0 , 3.8290694188d+0 , 2.3208830418d+0 , &
       &                                                                                           1.4544979381d+0 , 1.1817646130d+0 , 1.0277936032d+0 , 1.0000000000d+0   &
       &                                                                                          ]
  double precision                                   , dimension(8) :: energyISCOReference       =[                                                                        &
       &                                                                                           9.4280904158d-1 , 9.3064171394d-1 , 9.0878671490d-1 , 8.4424700801d-1 , &
       &                                                                                           7.3596989988d-1 , 6.6020592652d-1 , 5.9275740308d-1 , 5.7735026919d-1   &
       &                                                                                          ]
  double precision                                   , dimension(8) :: angularMomentumISCOReference=[                                                                      &
       &                                                                                           3.4641016151d+0 , 3.1535982815d+0 , 2.7559862886d+0 , 2.0997847561d+0 , &
       &                                                                                           1.5683649770d+0 , 1.3418378381d+0 , 1.1861513528d+0 , 1.1547005384d+0   &
       &                                                                                          ]
  ! Thin disk spin-up function, s=L_ISCO-2jE_ISCO, for the closed-form spins.
  double precision                                   , dimension(6) :: spinUpReference           =[                                                                        &
       &                                                                                           3.4641016151d+0 , 2.5952132531d+0 , 1.6654422308d+0 , 5.8014014171d-1 , &
       &                                                                                           1.1114457519d-1 , 2.2746396856d-2                                       &
       &                                                                                          ]
  ! Rezzolla et al. (2008) remnant spins for aligned spins.
  double precision                                   , dimension(8) :: massA                     =[1.0d0,1.0d0,1.0d0,1.0d0,0.5d0,1.0d0,1.0d0,3.0d0]                      , &
       &                                                               massB                     =[1.0d0,1.0d0,1.0d0,0.5d0,1.0d0,0.1d0,1.0d-6,1.0d0]                     , &
       &                                                               spinA                     =[0.0d0,0.5d0,0.9d0,0.0d0,0.7d0,0.8d0,0.6d0,0.99d0]                     , &
       &                                                               spinB                     =[0.0d0,0.5d0,0.9d0,0.7d0,0.0d0,0.3d0,0.9d0,0.99d0]
  double precision                                   , dimension(8) :: spinMergerReference       =[                                                                        &
       &                                                                                           6.8691602878d-1 , 8.3110352878d-1 , 9.3484352878d-1 , 6.7827300528d-1 , &
       &                                                                                           6.7827300528d-1 , 8.7376377797d-1 , 6.0000180605d-1 , 9.9793994464d-1   &
       &                                                                                          ]
  ! Bondi-Hoyle-Lyttleton accretion (Edgar 2004) onto a 10⁸M☉ black hole.
  double precision                                   , dimension(4) :: temperatureBondi          =[1.0d2 ,1.0d2 ,1.0d7 ,1.0d7 ]                                          , &
       &                                                               densityBondi              =[1.0d18,1.0d18,1.0d15,1.0d15]                                          , &
       &                                                               velocityBondi             =[0.0d0 ,5.0d1 ,0.0d0 ,0.0d0 ]
  double precision                                   , dimension(4) :: rateBondiReference        =[                                                                        &
       &                                                                                           6.7078381708d+14, 1.8991969905d+10, 2.1212046796d+04, 6.1962782916d+09  &
       &                                                                                          ]
  ! Meier (2001) thin disk jet power for a 10⁹M☉ black hole accreting at 10⁸M☉/Gyr.
  double precision                                   , dimension(5) :: spinJet                   =[0.00d0,0.50d0,0.80d0,0.90d0,0.99d0]
  double precision                                   , dimension(5) :: powerJetReference         =[                                                                        &
       &                                                                                           1.9276168747d+13, 1.2791545693d+14, 3.9816621869d+14, 4.2887547845d+14, &
       &                                                                                           4.5746762679d+14                                                        &
       &                                                                                          ]
  double precision                                   , dimension(8) :: radiusISCO           , energyISCO, &
       &                                                               angularMomentumISCO  , spinMerger, &
       &                                                               massMerger
  double precision                                   , dimension(6) :: efficiencyRadiative  , spinUp
  double precision                                   , dimension(4) :: rateBondi
  double precision                                   , dimension(5) :: powerJet
  type            (nodeComponentBlackHoleStandard   )               :: blackHole
  type            (accretionDisksShakuraSunyaev     )               :: accretionDisk
  type            (blackHoleBinaryMergerRezzolla2008)               :: blackHoleBinaryMerger
  double precision                                                  :: spin
  integer                                                           :: i

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Black hole physics")

  ! Kerr ISCO properties.
  call Unit_Tests_Begin_Group("Kerr ISCO")
  do i=1,size(spinISCO)
     spin                  =spinISCO(i)
     call blackHole%massSet(1.0d0)
     call blackHole%spinSet(spin )
     radiusISCO         (i)=Black_Hole_ISCO_Radius                   (spin     ,orbit=orbitPrograde                         )
     energyISCO         (i)=Black_Hole_ISCO_Specific_Energy          (spin     ,orbit=orbitPrograde                         )
     angularMomentumISCO(i)=Black_Hole_ISCO_Specific_Angular_Momentum(blackHole,orbit=orbitPrograde,units=unitsGravitational)
  end do
  call Assert("radius {closed form}"                     ,radiusISCO         (1:6),radiusISCOReference         (1:6),relTol=1.0d-8)
  call Assert("specific energy {closed form}"            ,energyISCO         (1:6),energyISCOReference         (1:6),relTol=1.0d-8)
  call Assert("specific angular momentum {closed form}"  ,angularMomentumISCO(1:6),angularMomentumISCOReference(1:6),relTol=1.0d-8)
  call Assert("radius {near-extremal}"                   ,radiusISCO         (7:8),radiusISCOReference         (7:8),relTol=1.0d-8)
  call Assert("specific energy {near-extremal}"          ,energyISCO         (7:8),energyISCOReference         (7:8),relTol=1.0d-3)
  ! The standard black hole component limits spin to 0.9999, so the specific angular momentum (which takes a black hole
  ! component, not a spin) for spins above that limit must equal that at j=0.9999.
  call Assert("specific angular momentum {spin limited}" ,angularMomentumISCO(7:8),[1.2405842000d+0,1.2405842000d+0]  ,relTol=1.0d-8)
  call Unit_Tests_End_Group()

  ! Kerr frame-dragging frequency in the equatorial plane, ω=-g_tφ/g_φφ in Boyer-Lindquist coordinates. This tests the metric
  ! factor 𝒜 (which enters as ω=2j/𝒜r³).
  call Assert("frame-dragging frequency"                                 , &
       &      [                                                            &
       &       Black_Hole_Frame_Dragging_Frequency(0.50d0,2.0000000000d0), &
       &       Black_Hole_Frame_Dragging_Frequency(0.90d0,2.3208830418d0), &
       &       Black_Hole_Frame_Dragging_Frequency(0.99d0,1.5000000000d0)  &
       &      ]                                                          , &
       &      [1.1111111111d-1,1.1249052728d-1,2.9094756331d-1]          , &
       &      relTol=1.0d-8                                                &
       &     )

  ! Thin (Shakura-Sunyaev) accretion disk.
  call Unit_Tests_Begin_Group("Shakura-Sunyaev accretion disk")
  do i=1,6
     call blackHole%massSet(1.0d0      )
     call blackHole%spinSet(spinISCO(i))
     efficiencyRadiative(i)=accretionDisk%efficiencyRadiative(blackHole,accretionRateMass=1.0d0)
     spinUp             (i)=accretionDisk%rateSpinUp         (blackHole,accretionRateMass=1.0d0)
  end do
  call Assert("radiative efficiency",efficiencyRadiative,1.0d0-energyISCOReference(1:6),relTol=1.0d-8)
  call Assert("spin-up rate"        ,spinUp             ,spinUpReference               ,relTol=1.0d-8)
  do i=1,size(spinJet)
     call blackHole%massSet(1.0d9     )
     call blackHole%spinSet(spinJet(i))
     powerJet           (i)=accretionDisk%powerJet           (blackHole,accretionRateMass=1.0d8)
  end do
  call Assert("jet power {Meier 2001}",powerJet,powerJetReference,relTol=2.0d-3)
  call Unit_Tests_End_Group()

  ! Eddington accretion rate.
  call blackHole%massSet(1.0d8)
  call Assert("Eddington accretion rate",Black_Hole_Eddington_Accretion_Rate(blackHole),2.2198030698d+08,relTol=1.0d-3)

  ! Bondi-Hoyle-Lyttleton accretion.
  call Unit_Tests_Begin_Group("Bondi-Hoyle-Lyttleton accretion")
  do i=1,size(rateBondi)
     if (i < size(rateBondi)) then
        rateBondi(i)=Bondi_Hoyle_Lyttleton_Accretion_Rate(1.0d8,densityBondi(i),velocityBondi(i),temperatureBondi(i)              )
     else
        rateBondi(i)=Bondi_Hoyle_Lyttleton_Accretion_Rate(1.0d8,densityBondi(i),velocityBondi(i),temperatureBondi(i),radius=1.0d-3)
     end if
  end do
  call Assert("sound speed {T=10²K}"     ,Ideal_Gas_Sound_Speed                 (      1.0d2          ),1.5246411212d+00  ,relTol=1.0d-5)
  call Assert("sound speed {T=10⁷K}"     ,Ideal_Gas_Sound_Speed                 (      1.0d7          ),4.8213385574d+02  ,relTol=1.0d-5)
  call Assert("Jeans length"             ,Ideal_Gas_Jeans_Length                (      1.0d4 ,1.0d15  ),7.3516882867d-03  ,relTol=1.0d-4)
  call Assert("accretion radius {T=10²K}",Bondi_Hoyle_Lyttleton_Accretion_Radius(1.0d8,1.0d2          ),1.8502304790d-01  ,relTol=5.0d-4)
  call Assert("accretion radius {T=10⁷K}",Bondi_Hoyle_Lyttleton_Accretion_Radius(1.0d8,1.0d7          ),1.8502304790d-06  ,relTol=5.0d-4)
  call Assert("accretion rate"           ,rateBondi                                                    ,rateBondiReference,relTol=5.0d-4)
  call Unit_Tests_End_Group()

  ! Binary black hole mergers.
  call Unit_Tests_Begin_Group("Rezzolla et al. (2008) binary mergers")
  do i=1,size(massA)
     call blackHoleBinaryMerger%merge(massA(i),massB(i),spinA(i),spinB(i),massMerger(i),spinMerger(i))
  end do
  call Assert("remnant spin"                     ,spinMerger   ,spinMergerReference,relTol=1.0d-8)
  call Assert("remnant mass"                     ,massMerger   ,massA+massB        ,relTol=1.0d-12)
  call Assert("remnant spin symmetric in A and B",spinMerger(4),spinMerger(5)      ,relTol=1.0d-12)
  call Unit_Tests_End_Group()

  ! End the unit testing.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Black_Hole_Physics
