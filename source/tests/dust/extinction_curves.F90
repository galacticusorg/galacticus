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
  Tests the ``dustExtinctionCurve`` class against the ``stellarSpectraDustAttenuation`` implementations which it
  replaces.

  The old class returned ``vBandAttenuation`` multiplied by the shape of the extinction curve, so with a unit
  ``vBandAttenuation`` it returned exactly the quantity the new class returns from ``attenuationRelative``. Every
  curve must therefore agree to round-off, which is what pins the ports as faithful.
  !!}
  use :: Display                          , only : displayVerbositySet                      , verbosityLevelStandard
  use :: Dust_Extinction_Curves           , only : dustExtinctionCurveCalzetti2000          , dustExtinctionCurveCardelli1989          , dustExtinctionCurveGordon2003               , dustExtinctionCurveNoll2009              , &
       &                                           dustExtinctionCurvePowerLaw              , dustExtinctionCurvePrevotBouchet         , dustExtinctionCurveWittGordon2000           , gordon2003SampleLMC                      , &
       &                                           gordon2003SampleSMCBar                   , wavelengthVBand                          , wittGordon2000ModelMilkyWayShellTau3        , wittGordon2000ModelSMCShellTau3
  ! Both modules declare enumerations of the same name for the Gordon et al. samples and the Witt & Gordon models --
  ! they are distinct types in distinct modules. Rename the old module's on import so that each constructor is handed
  ! the enumeration belonging to its own module. The collision disappears when the old class is retired.
  use :: Stellar_Spectra_Dust_Attenuations, only : stellarSpectraDustAttenuationCalzetti2000   , stellarSpectraDustAttenuationCardelli1989          , stellarSpectraDustAttenuationCharlotFall2000, stellarSpectraDustAttenuationGordon2003, &
       &                                           stellarSpectraDustAttenuationPrevotBouchet  , stellarSpectraDustAttenuationWittGordon2000        , &
       &                                           oldGordon2003SampleSMCBar=>gordon2003SampleSMCBar                                                , oldGordon2003SampleLMC=>gordon2003SampleLMC , &
       &                                           oldWittGordon2000ModelMilkyWay=>wittGordon2000ModelMilkyWayShellTau3                             , oldWittGordon2000ModelSMC=>wittGordon2000ModelSMCShellTau3
  use :: Unit_Tests                       , only : Assert                                   , Unit_Tests_Begin_Group                    , Unit_Tests_End_Group       , Unit_Tests_Finish
  implicit none
  type            (dustExtinctionCurveCalzetti2000             )               :: curveCalzetti2000
  type            (dustExtinctionCurveCardelli1989             )               :: curveCardelli1989
  type            (dustExtinctionCurvePowerLaw                 )               :: curvePowerLaw
  type            (dustExtinctionCurveGordon2003               )               :: curveGordon2003SMC   , curveGordon2003LMC
  type            (dustExtinctionCurvePrevotBouchet            )               :: curvePrevotBouchet
  type            (dustExtinctionCurveWittGordon2000           )               :: curveWittGordon2000MW, curveWittGordon2000SMC
  type            (dustExtinctionCurveNoll2009                 )               :: curveNoll2009Plain   , curveNoll2009Bump
  type            (stellarSpectraDustAttenuationCalzetti2000   )               :: oldCalzetti2000
  type            (stellarSpectraDustAttenuationCardelli1989   )               :: oldCardelli1989
  type            (stellarSpectraDustAttenuationCharlotFall2000)               :: oldCharlotFall2000   , oldCharlotFall2000BirthClouds
  type            (stellarSpectraDustAttenuationGordon2003     )               :: oldGordon2003SMC     , oldGordon2003LMC
  type            (stellarSpectraDustAttenuationPrevotBouchet  )               :: oldPrevotBouchet
  type            (stellarSpectraDustAttenuationWittGordon2000 )               :: oldWittGordon2000MW  , oldWittGordon2000SMC
  ! Wavelengths spanning the ultraviolet to the near-infrared, in Angstroms, plus points outside the fitted ranges.
  double precision                                              , dimension(9) :: wavelengths      =[1.0d3,1.5d3,2.5d3,3.5d3,wavelengthVBand,7.0d3,1.0d4,2.0d4,3.0d4]
  double precision                                                             :: new                  , old
  double precision                                                             :: excessBump           , excessWing
  double precision                                                             :: attenuationVBand     , depthOpticalBirthClouds      , &
       &                                                                          factorBirthClouds    , age
  integer                                                                      :: i                    , j
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
  ! The reference wavelength is given explicitly: the curve is compared below against
  ! `stellarSpectraDustAttenuationCharlotFall2000`, which normalizes at the Buser V effective wavelength.
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

  ! Each port must reproduce the class it replaces, evaluated with a unit V-band attenuation.
  call Unit_Tests_Begin_Group("agreement with stellarSpectraDustAttenuation")
  oldCalzetti2000   =stellarSpectraDustAttenuationCalzetti2000   (                        )
  oldCardelli1989   =stellarSpectraDustAttenuationCardelli1989   (3.1d0                   )
  ! The Charlot & Fall class is a power law in wavelength; with the birth cloud optical depth equal to zero and an age
  ! beyond the birth cloud lifetime it reduces to the bare power-law curve.
  oldCharlotFall2000=stellarSpectraDustAttenuationCharlotFall2000(0.7d0,1.0d-2,1.0d0,1.0d0)
  curveGordon2003LMC    =dustExtinctionCurveGordon2003              (gordon2003SampleLMC            )
  curveWittGordon2000SMC=dustExtinctionCurveWittGordon2000          (wittGordon2000ModelSMCShellTau3)
  oldGordon2003SMC      =stellarSpectraDustAttenuationGordon2003    (oldGordon2003SampleSMCBar      )
  oldGordon2003LMC      =stellarSpectraDustAttenuationGordon2003    (oldGordon2003SampleLMC         )
  oldPrevotBouchet      =stellarSpectraDustAttenuationPrevotBouchet (2.7d0                          )
  oldWittGordon2000MW   =stellarSpectraDustAttenuationWittGordon2000(oldWittGordon2000ModelMilkyWay )
  oldWittGordon2000SMC  =stellarSpectraDustAttenuationWittGordon2000(oldWittGordon2000ModelSMC      )
  do i=1,size(wavelengths)
     write (label,'(a,f9.1,a)') "lambda=",wavelengths(i)," A"
     new=curveCalzetti2000%attenuationRelative(wavelengths(i))
     old=oldCalzetti2000  %attenuation        (wavelengths(i),1.0d0,1.0d0)
     call Assert("calzetti2000 "//trim(label),new,old,absTol=1.0d-12)
     new=curveCardelli1989%attenuationRelative(wavelengths(i))
     old=oldCardelli1989  %attenuation        (wavelengths(i),1.0d0,1.0d0)
     call Assert("cardelli1989 "//trim(label),new,old,absTol=1.0d-12)
     new=curvePowerLaw    %attenuationRelative(wavelengths(i))
     old=oldCharlotFall2000%attenuation       (wavelengths(i),1.0d0,1.0d0)
     call Assert("powerLaw "    //trim(label),new,old,absTol=1.0d-12)
     new=curveGordon2003SMC   %attenuationRelative(wavelengths(i))
     old=oldGordon2003SMC     %attenuation        (wavelengths(i),1.0d0,1.0d0)
     call Assert("gordon2003 SMCBar "        //trim(label),new,old,absTol=1.0d-12)
     new=curveGordon2003LMC   %attenuationRelative(wavelengths(i))
     old=oldGordon2003LMC     %attenuation        (wavelengths(i),1.0d0,1.0d0)
     call Assert("gordon2003 LMC "           //trim(label),new,old,absTol=1.0d-12)
     new=curvePrevotBouchet   %attenuationRelative(wavelengths(i))
     old=oldPrevotBouchet     %attenuation        (wavelengths(i),1.0d0,1.0d0)
     call Assert("prevotBouchet "            //trim(label),new,old,absTol=1.0d-12)
     new=curveWittGordon2000MW%attenuationRelative(wavelengths(i))
     old=oldWittGordon2000MW  %attenuation        (wavelengths(i),1.0d0,1.0d0)
     call Assert("wittGordon2000 MilkyWay "  //trim(label),new,old,absTol=1.0d-12)
     new=curveWittGordon2000SMC%attenuationRelative(wavelengths(i))
     old=oldWittGordon2000SMC  %attenuation       (wavelengths(i),1.0d0,1.0d0)
     call Assert("wittGordon2000 SMC "       //trim(label),new,old,absTol=1.0d-12)
  end do
  call Unit_Tests_End_Group()

  ! The SED fitting likelihood replaces the Charlot & Fall (2000) attenuation class with a power-law extinction curve
  ! scaled by the V-band attenuation and, for populations younger than the birth cloud lifetime, by a constant factor
  ! of one plus the birth cloud optical depth. Check that this reproduces the class it replaces for populations on
  ! both sides of that lifetime, and for an attenuation which is not unity.
  call Unit_Tests_Begin_Group("Charlot & Fall birth clouds as a factor on the curve")
  attenuationVBand             =1.3d0
  depthOpticalBirthClouds      =0.8d0
  oldCharlotFall2000BirthClouds=stellarSpectraDustAttenuationCharlotFall2000(0.7d0,1.0d-2,1.0d0,depthOpticalBirthClouds)
  do i=1,size(wavelengths)
     ! Ages either side of the 10 Myr birth cloud lifetime.
     do j=1,2
        if (j == 1) then
           age              =5.0d-3
           factorBirthClouds=1.0d0+depthOpticalBirthClouds
        else
           age              =5.0d-1
           factorBirthClouds=1.0d0
        end if
        write (label,'(a,f9.1,a,f6.3)') "lambda=",wavelengths(i)," age=",age
        new=+attenuationVBand                                   &
             & *curvePowerLaw%attenuationRelative(wavelengths(i)) &
             & *factorBirthClouds
        old=oldCharlotFall2000BirthClouds%attenuation(wavelengths(i),age,attenuationVBand)
        call Assert("charlotFall2000 "//trim(label),new,old,absTol=1.0d-12)
     end do
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
