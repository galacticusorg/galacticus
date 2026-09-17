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
Contains a program to test dust emission spectra.
!!}

program Test_Dust_Emission_Spectra
  !!{RST
  Tests the ``dustEmissionSpectrum`` class, and the far-infrared absorption opacity of the ``dustProperties`` class
  on which it rests.

  The central property is energy balance: dust must emit exactly the luminosity it absorbs. That is checked two ways
  for the modified blackbody---by integrating the returned spectrum numerically, and by integrating
  :math:`4\pi M_\mathrm{dust} \kappa_\mathrm{abs} B_\nu(T)` at the returned temperature using an independent Planck
  function---so that neither the closed-form temperature nor the analytic normalization can be wrong unnoticed.
  !!}
  use :: Cosmology_Functions             , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters            , only : cosmologyParametersSimple
  use :: Display                         , only : displayVerbositySet                         , verbosityLevelStandard
  use :: Dust_Emission_Spectra           , only : dustEmissionSpectrumBlackBodyModified       , dustEmissionSpectrumList , dustEmissionSpectrumSum
  use :: Dust_Properties                 , only : depthOpticalVPerSurfaceDensityMetalsMilkyWay, dustPropertiesSimple     , dustToMetalsRatioMilkyWay
  use :: Numerical_Constants_Astronomical, only : luminositySolar                             , massSolar
  use :: Numerical_Constants_Math        , only : Pi
  use :: Numerical_Constants_Physical    , only : speedLight
  use :: Numerical_Constants_Prefixes    , only : centi                                       , kilo
  use :: Numerical_Constants_Units       , only : metersToAngstroms
  use :: Thermodynamics_Radiation        , only : Blackbody_Emission                          , radianceTypeFrequency
  use :: Unit_Tests                      , only : Assert                                      , Unit_Tests_Begin_Group   , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  ! Wavelengths spanning 1 μm to 10 cm, finely enough, and far enough into the tails, that integrals of the spectra over
  ! frequency are accurate to much better than the tolerances asserted below.
  integer                                                , parameter                    :: countWavelengths      =20000
  double precision                                       , parameter                    :: wavelengthMinimum     =1.0d4 , wavelengthMaximum   =1.0d9
  ! Properties of the dust population used throughout.
  double precision                                       , parameter                    :: luminosityAbsorbed    =1.0d10, massDust            =1.0d7 , &
       &                                                                                   opacityReference      =0.77d0, wavelengthReference =8.5d6 , &
       &                                                                                   exponent              =2.0d0 , timePresent         =13.0d0
  double precision                                       , parameter                    :: toleranceIntegral     =1.0d-5
  type            (cosmologyParametersSimple            ), pointer                      :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda       ), pointer                      :: cosmologyFunctions_
  type            (dustPropertiesSimple                 ), pointer                      :: dustProperties_
  type            (dustEmissionSpectrumBlackBodyModified), pointer                      :: greybody_                     , greybodyCMB_              , &
       &                                                                                   greybodyFixed_                , greybodyWarm_             , &
       &                                                                                   greybodyCold_
  type            (dustEmissionSpectrumSum              ), pointer                      :: sum_
  type            (dustEmissionSpectrumList             ), pointer                      :: members
  double precision                                       , dimension(countWavelengths)  :: wavelengths                   , frequencies               , &
       &                                                                                   luminosity                    , luminosityCold            , &
       &                                                                                   luminosityWarm
  double precision                                                                      :: temperature                   , temperatureAbsorbed       , &
       &                                                                                   temperatureCMB                , timeHighRedshift          , &
       &                                                                                   emitted
  integer                                                                               :: i

  call displayVerbositySet(verbosityLevelStandard)
  ! Logarithmically spaced wavelengths, and the corresponding frequencies.
  do i=1,countWavelengths
     wavelengths(i)=exp(log(wavelengthMinimum)+log(wavelengthMaximum/wavelengthMinimum)*dble(i-1)/dble(countWavelengths-1))
  end do
  frequencies=speedLight*metersToAngstroms/wavelengths
  ! Construct the objects needed.
  allocate(cosmologyParameters_)
  allocate(cosmologyFunctions_ )
  allocate(dustProperties_     )
  allocate(greybody_           )
  allocate(greybodyCMB_        )
  allocate(greybodyFixed_      )
  allocate(greybodyWarm_       )
  allocate(greybodyCold_       )
  allocate(sum_                )
  !![
  <referenceConstruct object="cosmologyParameters_" >
   <constructor>
    cosmologyParametersSimple     (                                          &amp;
     &amp;                         OmegaMatter         = 0.30d0            , &amp;
     &amp;                         OmegaBaryon         = 0.05d0            , &amp;
     &amp;                         OmegaDarkEnergy     = 0.70d0            , &amp;
     &amp;                         temperatureCMB      = 2.72548d0         , &amp;
     &amp;                         HubbleConstant      =70.00d0              &amp;
     &amp;                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"  >
   <constructor>
    cosmologyFunctionsMatterLambda(                                          &amp;
     &amp;                         cosmologyParameters_=cosmologyParameters_ &amp;
     &amp;                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="dustProperties_" constructor="dustPropertiesSimple(dustToMetalsRatioMilkyWay,depthOpticalVPerSurfaceDensityMetalsMilkyWay/dustToMetalsRatioMilkyWay,opacityReference,wavelengthReference,exponent)"/>
  <referenceConstruct object="greybody_"       constructor="dustEmissionSpectrumBlackBodyModified(-1.0d0,.false.,dustProperties_,cosmologyFunctions_)"/>
  <referenceConstruct object="greybodyCMB_"    constructor="dustEmissionSpectrumBlackBodyModified(-1.0d0,.true. ,dustProperties_,cosmologyFunctions_)"/>
  <referenceConstruct object="greybodyFixed_"  constructor="dustEmissionSpectrumBlackBodyModified(25.0d0,.false.,dustProperties_,cosmologyFunctions_)"/>
  <referenceConstruct object="greybodyWarm_"   constructor="dustEmissionSpectrumBlackBodyModified(50.0d0,.false.,dustProperties_,cosmologyFunctions_)"/>
  <referenceConstruct object="greybodyCold_"   constructor="dustEmissionSpectrumBlackBodyModified(20.0d0,.false.,dustProperties_,cosmologyFunctions_)"/>
  !!]

  call Unit_Tests_Begin_Group("Far-infrared absorption opacity")
  call Assert("at the reference wavelength"     ,dustProperties_%opacityAbsorption(wavelengthReference      ),opacityReference                  ,relTol=1.0d-12)
  call Assert("at half the reference wavelength",dustProperties_%opacityAbsorption(wavelengthReference/2.0d0),opacityReference*2.0d0**exponent,relTol=1.0d-12)
  call Unit_Tests_End_Group()

  call Unit_Tests_Begin_Group("Modified blackbody: energy balance")
  temperature=greybody_%temperature(luminosityAbsorbed,massDust,timePresent)
  ! Independently integrate 4π M κ_abs(ν) B_ν(T) over frequency, in SI units, using the Planck function of
  ! `Thermodynamics_Radiation`, and convert to Solar luminosities.
  emitted=0.0d0
  do i=2,countWavelengths
     emitted=emitted+0.5d0*(emission(i)+emission(i-1))*(frequencies(i-1)-frequencies(i))
  end do
  call Assert("temperature is plausible for dust"                  ,temperature > 10.0d0 .and. temperature < 100.0d0                  ,.true.                                   )
  call Assert("emission at the temperature equals the absorption"  ,emitted/luminositySolar                                            ,luminosityAbsorbed,relTol=toleranceIntegral)
  luminosity=greybody_%luminosity(wavelengths,luminosityAbsorbed,massDust,timePresent)
  call Assert("spectrum integrates to the absorbed luminosity"     ,integral(luminosity)                                               ,luminosityAbsorbed,relTol=toleranceIntegral)
  call Assert("temperature scales as the (4+β)th root of luminosity",greybody_%temperature(2.0d0**(4.0d0+exponent)*luminosityAbsorbed,massDust,timePresent),2.0d0*temperature,relTol=1.0d-12)
  call Assert("no absorption, no emission"                         ,all(greybody_%luminosity(wavelengths,0.0d0,massDust,timePresent) == 0.0d0),.true.                  )
  call Unit_Tests_End_Group()

  call Unit_Tests_Begin_Group("Modified blackbody: fixed temperature")
  call Assert("temperature is that fixed"                          ,greybodyFixed_%temperature(luminosityAbsorbed,massDust,timePresent),25.0d0            ,relTol=1.0d-12          )
  luminosity=greybodyFixed_%luminosity(wavelengths,luminosityAbsorbed,massDust,timePresent)
  call Assert("spectrum integrates to the absorbed luminosity"     ,integral(luminosity)                                               ,luminosityAbsorbed,relTol=toleranceIntegral)
  call Unit_Tests_End_Group()

  call Unit_Tests_Begin_Group("Modified blackbody: heating by the CMB")
  ! At redshift six the CMB is at almost 20 K, comparable to the dust, so its heating is substantial.
  timeHighRedshift=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(6.0d0))
  temperatureCMB  =cosmologyFunctions_%temperatureCMBEpochal(time=timeHighRedshift)
  temperature     =greybodyCMB_%temperature(luminosityAbsorbed,massDust,timeHighRedshift,temperatureAbsorbed)
  call Assert("CMB temperature at z=6"                              ,temperatureCMB                                  ,2.72548d0*7.0d0                                                    ,relTol=1.0d-6           )
  call Assert("heated by starlight alone as without the CMB"        ,temperatureAbsorbed                             ,greybody_%temperature(luminosityAbsorbed,massDust,timeHighRedshift),relTol=1.0d-12          )
  call Assert("T^(4+β) = T★^(4+β) + T_CMB^(4+β)"                    ,temperature**(4.0d0+exponent)                   ,temperatureAbsorbed**(4.0d0+exponent)+temperatureCMB**(4.0d0+exponent),relTol=1.0d-10    )
  luminosity=greybodyCMB_%luminosity(wavelengths,luminosityAbsorbed,massDust,timeHighRedshift)
  call Assert("spectrum integrates to absorbed plus CMB luminosity" ,integral(luminosity)                            ,luminosityAbsorbed*(temperature/temperatureAbsorbed)**(4.0d0+exponent),relTol=toleranceIntegral)
  call Unit_Tests_End_Group()

  call Unit_Tests_Begin_Group("Sum of spectra")
  allocate(members     )
  allocate(members%next)
  members     %dustEmissionSpectrum_ => greybodyCold_
  members%next%dustEmissionSpectrum_ => greybodyWarm_
  !![
  <referenceConstruct object="sum_" constructor="dustEmissionSpectrumSum(members,[0.3d0,0.7d0])"/>
  !!]
  luminosity    =sum_         %luminosity(wavelengths,      luminosityAbsorbed,      massDust,timePresent)
  luminosityCold=greybodyCold_%luminosity(wavelengths,0.3d0*luminosityAbsorbed,0.3d0*massDust,timePresent)
  luminosityWarm=greybodyWarm_%luminosity(wavelengths,0.7d0*luminosityAbsorbed,0.7d0*massDust,timePresent)
  call Assert("sum is the weighted sum of its members"         ,all(abs(luminosity-(luminosityCold+luminosityWarm)) <= 1.0d-12*maxval(luminosity)),.true.                            )
  call Assert("spectrum integrates to the absorbed luminosity" ,integral(luminosity)                                                              ,luminosityAbsorbed,relTol=toleranceIntegral)
  call Unit_Tests_End_Group()

  !![
  <objectDestructor name="sum_"                />
  <objectDestructor name="greybodyCold_"       />
  <objectDestructor name="greybodyWarm_"       />
  <objectDestructor name="greybodyFixed_"      />
  <objectDestructor name="greybodyCMB_"        />
  <objectDestructor name="greybody_"           />
  <objectDestructor name="dustProperties_"     />
  <objectDestructor name="cosmologyFunctions_" />
  <objectDestructor name="cosmologyParameters_"/>
  !!]
  call Unit_Tests_Finish()

contains

  double precision function integral(luminosityPerFrequency)
    !!{RST
    Integrate a luminosity per unit frequency over frequency by the trapezoidal rule on the wavelength grid.
    !!}
    implicit none
    double precision, intent(in   ), dimension(countWavelengths) :: luminosityPerFrequency
    integer                                                      :: j

    integral=0.0d0
    do j=2,countWavelengths
       integral=integral+0.5d0*(luminosityPerFrequency(j)+luminosityPerFrequency(j-1))*(frequencies(j-1)-frequencies(j))
    end do
    return
  end function integral

  double precision function emission(j)
    !!{RST
    Return :math:`4\pi M_\mathrm{dust} \kappa_\mathrm{abs}(\nu) B_\nu(T)`, in W Hz⁻¹, at the ``j``-th frequency and
    the temperature found by energy balance. The opacity is converted from cm² g⁻¹ to m² kg⁻¹.
    !!}
    implicit none
    integer, intent(in   ) :: j

    emission=+4.0d0                                                                              &
         &   *Pi                                                                                 &
         &   *massDust                                                                           &
         &   *massSolar                                                                          &
         &   *opacityReference                                                                   &
         &   *centi**2                                                                           &
         &   *kilo                                                                               &
         &   *(frequencies(j)*wavelengthReference/speedLight/metersToAngstroms)**exponent        &
         &   *Blackbody_Emission(wavelengths(j),temperature,radianceTypeFrequency)
    return
  end function emission

end program Test_Dust_Emission_Spectra
