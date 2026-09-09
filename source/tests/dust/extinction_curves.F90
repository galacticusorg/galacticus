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
Contains a program to test dust extinction curves.
!!}

program Test_Dust_Extinction_Curves
  !!{RST
  Tests the ``dustExtinctionCurve`` class.

  Each curve is pinned to reference values recorded from the ``stellarSpectraDustAttenuation`` implementations it
  replaces, taken immediately before those were retired. That class returned ``vBandAttenuation`` multiplied by the
  shape of the extinction curve, so evaluated with a unit ``vBandAttenuation`` it returned exactly the quantity
  ``attenuationRelative`` returns, and the two were verified to agree to round-off at the time the values below were
  recorded. Pinning them here keeps that agreement enforced now that there is nothing left to compare against.
  !!}
  use :: Display                          , only : displayVerbositySet            , verbosityLevelStandard
  use :: Dust_Extinction_Curves           , only : dustExtinctionCurveCalzetti2000, dustExtinctionCurveCardelli1989 , dustExtinctionCurveGordon2003       , dustExtinctionCurveNoll2009    , &
       &                                           dustExtinctionCurvePowerLaw    , dustExtinctionCurvePrevotBouchet, dustExtinctionCurveWittGordon2000   , gordon2003SampleLMC            , &
       &                                           gordon2003SampleSMCBar         , wavelengthVBand                 , wittGordon2000ModelMilkyWayShellTau3, wittGordon2000ModelSMCShellTau3
  use :: Unit_Tests                       , only : Assert                         , Unit_Tests_Begin_Group          , Unit_Tests_End_Group                , Unit_Tests_Finish
  implicit none
  type            (dustExtinctionCurveCalzetti2000             )               :: curveCalzetti2000
  type            (dustExtinctionCurveCardelli1989             )               :: curveCardelli1989
  type            (dustExtinctionCurvePowerLaw                 )               :: curvePowerLaw
  type            (dustExtinctionCurveGordon2003               )               :: curveGordon2003SMC   , curveGordon2003LMC
  type            (dustExtinctionCurvePrevotBouchet            )               :: curvePrevotBouchet
  type            (dustExtinctionCurveWittGordon2000           )               :: curveWittGordon2000MW, curveWittGordon2000SMC
  type            (dustExtinctionCurveNoll2009                 )               :: curveNoll2009Plain   , curveNoll2009Bump
  ! Wavelengths spanning the ultraviolet to the near-infrared, in Angstroms, plus points outside the fitted ranges.
  double precision                                              , dimension(9) :: wavelengths      =[1.0d3,1.5d3,2.5d3,3.5d3,wavelengthVBand,7.0d3,1.0d4,2.0d4,3.0d4]
  ! Reference values of k(lambda)/k_V at each of those wavelengths, recorded from the `stellarSpectraDustAttenuation`
  ! implementations these curves replaced. Values outside a curve's fitted range are whatever that curve does there --
  ! zero for those which return no attenuation, and, where a tabulation is extrapolated, occasionally negative. Those
  ! are pinned along with the rest so that any change in out-of-range behavior is noticed.
  double precision                                              , dimension(9) :: referenceCalzetti2000     =[                          &
       & 0.00000000000000000d+00, 2.55158181984453547d+00, 1.92966518518518471d+00, 1.52238635712486015d+00, &
       & 9.98580621004780400d-01, 7.56234885361552034d-01, 4.63604197530864237d-01, 1.22201728395061743d-01, &
       & 0.00000000000000000d+00                                                                             &
       & ]
  double precision                                              , dimension(9) :: referenceCardelli1989     =[                          &
       & 0.00000000000000000d+00, 2.66387921471371580d+00, 2.31531710595831353d+00, 1.59162968321680975d+00, &
       & 9.97887601751360376d-01, 7.50173391182706428d-01, 4.03999999999999915d-01, 1.32349733789694668d-01, &
       & 6.88995118580363469d-02                                                                             &
       & ]
  double precision                                              , dimension(9) :: referencePowerLaw         =[                          &
       & 3.29996029053773388d+00, 2.48453336083206944d+00, 1.73760360058940511d+00, 1.37297011079780140d+00, &
       & 1.00000000000000000d+00, 8.45162240799131936d-01, 6.58428640860369674d-01, 4.05310371390765201d-01, &
       & 3.05157350559360385d-01                                                                             &
       & ]
  double precision                                              , dimension(9) :: referenceGordon2003SMC    =[                          &
       & 1.00595919268920291d+01, 4.81792727775525353d+00, 2.54692852318176799d+00, 1.74890399294187193d+00, &
       & 1.00000000000000000d+00, 7.17174639798813462d-01, 3.31776122143060959d-01, 6.16548150057808705d-02, &
       &-1.07380534166073002d-01                                                                             &
       & ]
  double precision                                              , dimension(9) :: referenceGordon2003LMC    =[                          &
       & 4.12754168958747680d+00, 2.80424591812278790d+00, 2.27222496310502553d+00, 1.58102177627214791d+00, &
       & 1.00000000000000000d+00, 7.16472268627521225d-01, 4.03367472592334519d-01, 7.65650387252765491d-02, &
       &-9.57891610026187917d-02                                                                             &
       & ]
  double precision                                              , dimension(9) :: referencePrevotBouchet    =[                          &
       & 9.27330585565893273d+00, 4.93187052997590403d+00, 2.67106808412647467d+00, 1.76987184474170056d+00, &
       & 1.00000000000000000d+00, 7.09132865210911345d-01, 3.87923033387325311d-01, 7.81458358622943161d-02, &
       & 1.30243059770491203d-02                                                                             &
       & ]
  double precision                                              , dimension(9) :: referenceWittGordon2000MW =[                          &
       & 5.20294159232237874d+00, 2.65872856288425385d+00, 2.31247089873966738d+00, 1.60490016388932810d+00, &
       & 1.00000000000000000d+00, 7.33475028664789752d-01, 4.03355436035962800d-01, 1.32391061531922011d-01, &
       & 6.85421403837320614d-02                                                                             &
       & ]
  double precision                                              , dimension(9) :: referenceWittGordon2000SMC=[                          &
       & 9.72433515012639837d+00, 4.85890663952229485d+00, 2.55658932074019596d+00, 1.74721034171619527d+00, &
       & 1.00000000000000000d+00, 7.07366960114033327d-01, 3.92576870040343917d-01, 1.02144728098181620d-01, &
       & 3.31737087535214983d-02                                                                             &
       & ]
  double precision                                                             :: new
  double precision                                                             :: excessBump           , excessWing
  integer                                                                      :: i
  character       (len=32                                      )               :: label

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Dust extinction curves")

  ! Every curve must be normalized to approximately unity in the V-band. The empirically fit curves are not exactly
  ! unity at the Buser V effective wavelength: each anchors its own normalization at a slightly different nominal
  ! wavelength (Cardelli et al. at x=1.82 inverse microns, i.e. 0.5495 microns) and the fitting functions do not
  ! evaluate to exactly one even there. The residual is a few parts in a thousand, and is deliberately not removed --
  ! see the class description -- so a loose tolerance here is checking normalization, not fit accuracy.
  call Unit_Tests_Begin_Group("normalization in the V-band")
  curveCalzetti2000=dustExtinctionCurveCalzetti2000(     )
  curveCardelli1989=dustExtinctionCurveCardelli1989(3.1d0)
  ! The reference wavelength is given explicitly: the reference values below come from
  ! `stellarSpectraDustAttenuationCharlotFall2000`, which normalized at the Buser V effective wavelength.
  curvePowerLaw    =dustExtinctionCurvePowerLaw    (0.7d0,wavelengthVBand)
  call Assert("calzetti2000",curveCalzetti2000%attenuationRelative(wavelengthVBand),1.0d0,relTol=3.0d-3)
  call Assert("cardelli1989",curveCardelli1989%attenuationRelative(wavelengthVBand),1.0d0,relTol=3.0d-3)
  ! The power law is exactly unity by construction.
  call Assert("powerLaw"    ,curvePowerLaw    %attenuationRelative(wavelengthVBand),1.0d0,relTol=1.0d-12)
  ! The tabulated curves normalize by their own interpolated V-band value, so they are exactly unity there.
  curveGordon2003SMC   =dustExtinctionCurveGordon2003    (gordon2003SampleSMCBar              )
  curvePrevotBouchet   =dustExtinctionCurvePrevotBouchet (2.7d0                               )
  curveWittGordon2000MW=dustExtinctionCurveWittGordon2000(wittGordon2000ModelMilkyWayShellTau3)
  call Assert("gordon2003"    ,curveGordon2003SMC   %attenuationRelative(wavelengthVBand),1.0d0,relTol=1.0d-12)
  call Assert("prevotBouchet" ,curvePrevotBouchet   %attenuationRelative(wavelengthVBand),1.0d0,relTol=1.0d-12)
  call Assert("wittGordon2000",curveWittGordon2000MW%attenuationRelative(wavelengthVBand),1.0d0,relTol=1.0d-12)
  call Unit_Tests_End_Group()

  ! Each curve must reproduce the reference values recorded from the class it replaced.
  call Unit_Tests_Begin_Group("reference values")
  curveGordon2003LMC    =dustExtinctionCurveGordon2003    (gordon2003SampleLMC            )
  curveWittGordon2000SMC=dustExtinctionCurveWittGordon2000(wittGordon2000ModelSMCShellTau3)
  do i=1,size(wavelengths)
     write (label,'(a,f9.1,a)') "lambda=",wavelengths(i)," A"
     new=curveCalzetti2000     %attenuationRelative(wavelengths(i))
     call Assert("calzetti2000 "             //trim(label),new,referenceCalzetti2000     (i),absTol=1.0d-12)
     new=curveCardelli1989     %attenuationRelative(wavelengths(i))
     call Assert("cardelli1989 "             //trim(label),new,referenceCardelli1989     (i),absTol=1.0d-12)
     new=curvePowerLaw         %attenuationRelative(wavelengths(i))
     call Assert("powerLaw "                 //trim(label),new,referencePowerLaw         (i),absTol=1.0d-12)
     new=curveGordon2003SMC    %attenuationRelative(wavelengths(i))
     call Assert("gordon2003 SMCBar "        //trim(label),new,referenceGordon2003SMC    (i),absTol=1.0d-12)
     new=curveGordon2003LMC    %attenuationRelative(wavelengths(i))
     call Assert("gordon2003 LMC "           //trim(label),new,referenceGordon2003LMC    (i),absTol=1.0d-12)
     new=curvePrevotBouchet    %attenuationRelative(wavelengths(i))
     call Assert("prevotBouchet "            //trim(label),new,referencePrevotBouchet    (i),absTol=1.0d-12)
     new=curveWittGordon2000MW %attenuationRelative(wavelengths(i))
     call Assert("wittGordon2000 MilkyWay "  //trim(label),new,referenceWittGordon2000MW (i),absTol=1.0d-12)
     new=curveWittGordon2000SMC%attenuationRelative(wavelengths(i))
     call Assert("wittGordon2000 SMC "       //trim(label),new,referenceWittGordon2000SMC(i),absTol=1.0d-12)
  end do
  call Unit_Tests_End_Group()

  ! With no bump and no slope modification, Noll et al. reduces exactly to Calzetti et al.; with a bump it must exceed
  ! it near 2175A and match it far from the bump.
  call Unit_Tests_Begin_Group("noll2009 reduces to calzetti2000")
  curveNoll2009Plain=dustExtinctionCurveNoll2009(0.0d0,0.0d0)
  curveNoll2009Bump =dustExtinctionCurveNoll2009(3.5d0,0.0d0)
  do i=1,size(wavelengths)
     write (label,'(a,f9.1,a)') "lambda=",wavelengths(i)," A"
     call Assert("no bump, no slope "//trim(label),curveNoll2009Plain%attenuationRelative(wavelengths(i)),curveCalzetti2000%attenuationRelative(wavelengths(i)),absTol=1.0d-12)
  end do
  ! The Drude profile is a Lorentzian, so it has power-law wings: the bump is strongly peaked at 2175A but does not
  ! vanish in the optical. Check that it dominates at the bump and is small -- but explicitly *not* zero -- far away.
  excessBump =+curveNoll2009Bump %attenuationRelative(2175.0d0) &
       &      /curveCalzetti2000 %attenuationRelative(2175.0d0) &
       &      -1.0d0
  excessWing =+curveNoll2009Bump %attenuationRelative(7000.0d0) &
       &      /curveCalzetti2000 %attenuationRelative(7000.0d0) &
       &      -1.0d0
  call Assert("bump dominates at 2175A"        ,excessBump > 2.0d-1                                     ,.true.)
  call Assert("wing present but small at 7000A",excessWing > 0.0d0             .and. excessWing < 1.0d-2,.true.)
  call Assert("wing far weaker than the bump " ,excessWing < 1.0d-1*excessBump                          ,.true.)
  call Unit_Tests_End_Group()

  ! Curves fit over a limited range must report that range, and must return zero outside it.
  call Unit_Tests_Begin_Group("wavelength ranges")
  call assertRange(curveCalzetti2000%attenuationRelative(1.0d3),"calzetti2000 below range")
  call assertRange(curveCalzetti2000%attenuationRelative(3.0d4),"calzetti2000 above range")
  call assertRange(curveCardelli1989%attenuationRelative(1.0d3),"cardelli1989 below range")
  call assertRange(curveCardelli1989%attenuationRelative(4.0d4),"cardelli1989 above range")
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

contains

  subroutine assertRange(value,label)
    !!{RST
    Assert that a curve returns zero outside of the range over which it is fit.
    !!}
    implicit none
    double precision       , intent(in   ) :: value
    character       (len=*), intent(in   ) :: label

    call Assert(label,value,0.0d0,absTol=1.0d-12)
    return
  end subroutine assertRange

end program Test_Dust_Extinction_Curves
