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

  The template spectra are checked against the tabulations themselves, read here independently of the classes: the
  shapes at tabulated wavelengths, the interpolation between templates, and, for :cite:t:`draine_infrared_2007`, the
  intensity of heating found by energy balance from the tabulated power.
  !!}
  use :: Cosmology_Functions             , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters            , only : cosmologyParametersSimple
  use :: Display                         , only : displayVerbositySet                         , verbosityLevelStandard
  use :: Dust_Emission_Spectra           , only : dustEmissionSpectrumBlackBodyModified       , dustEmissionSpectrumList        , dustEmissionSpectrumSum         , &
       &                                          dustEmissionSpectrumDaleHelou2002           , dustEmissionSpectrumDraineLi2007, draineLi2007GrainModelMilkyWay60
  use :: Dust_Properties                 , only : depthOpticalVPerSurfaceDensityMetalsMilkyWay, dustPropertiesSimple            , dustToMetalsRatioMilkyWay
  use :: IO_HDF5                         , only : hdf5File                                    , hdf5Group                       , ioHDF5AccessInitialize
  use :: Input_Paths                     , only : inputPath                                   , pathTypeDataStatic
  use :: ISO_Varying_String              , only : char                                        , operator(//)
  use :: Numerical_Constants_Astronomical, only : luminositySolar                             , massSolar
  use :: Numerical_Constants_Atomic      , only : massHydrogenAtom
  use :: Numerical_Constants_Math        , only : Pi
  use :: Numerical_Constants_Physical    , only : speedLight
  use :: Numerical_Constants_Prefixes    , only : centi                                       , kilo
  use :: Numerical_Constants_Units       , only : ergs                                        , metersToAngstroms
  use :: Thermodynamics_Radiation        , only : Blackbody_Emission                          , radianceTypeFrequency
  use :: Unit_Tests                      , only : Assert                                      , Unit_Tests_Begin_Group           , Unit_Tests_End_Group           , &
       &                                          Unit_Tests_Finish
  implicit none
  ! Wavelengths spanning 1 μm to 10 cm, finely enough, and far enough into the tails, that integrals of the spectra over
  ! frequency are accurate to much better than the tolerances asserted below.
  integer                                                , parameter                     :: countWavelengths      =20000
  double precision                                       , parameter                     :: wavelengthMinimum     =1.0d4 , wavelengthMaximum   =1.0d9
  ! Properties of the dust population used throughout .
  double precision                                       , parameter                     :: luminosityAbsorbed    =1.0d10, massDust            =1.0d7 , &
       &                                                                                    opacityReference      =0.77d0, wavelengthReference =8.5d6 , &
       &                                                                                    exponent              =2.0d0 , timePresent         =13.0d0
  double precision                                       , parameter                     :: toleranceIntegral     =1.0d-5
  ! Tolerance on the integrals of template spectra. These are truncated sharply at the ends of their tabulations, where
  ! the trapezoidal rule on the wavelength grid used here is accurate only to a fraction of a grid cell.
  double precision                                       , parameter                     :: toleranceTemplate     =1.0d-4
  ! The Draine & Li (2007) model tested: the fraction of dust mass heated by a power law of intensities, and the maximum
  ! intensity of that power law, which is the last one tabulated.
  double precision                                       , parameter                     :: fractionMassPowerLaw  =0.01d0, intensityMaximum    =1.0d6
  integer                                                , parameter                     :: indexIntensityMaximum =4
  type            (cosmologyParametersSimple            ), pointer                       :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda       ), pointer                       :: cosmologyFunctions_
  type            (dustPropertiesSimple                 ), pointer                       :: dustProperties_
  type            (dustEmissionSpectrumBlackBodyModified), pointer                       :: greybody_                     , greybodyCMB_              , &
       &                                                                                    greybodyFixed_                , greybodyWarm_             , &
       &                                                                                    greybodyCold_
  type            (dustEmissionSpectrumSum              ), pointer                       :: sum_
  type            (dustEmissionSpectrumDaleHelou2002    ), pointer                       :: daleHelou_                    , daleHelouMidpoint_
  type            (dustEmissionSpectrumDraineLi2007     ), pointer                       :: draineLi_                     , draineLiFixed_            , &
       &                                                                                    draineLiFixedNext_            , draineLiFixedBetween_     , &
       &                                                                                    draineLiFixedHighest_
  type            (hdf5File                             )                                :: fileDale                      , fileDraine
  type            (hdf5Group                            )                                :: groupDraine
  double precision                                       , allocatable, dimension(:    ) :: wavelengthsDale               , alphasDale                , &
       &                                                                                    wavelengthsDraine             , intensitiesDraine         , &
       &                                                                                    powerSingleDraine
  double precision                                       , allocatable, dimension(:,:  ) :: luminosityDale                , emissivitySingleDraine    , &
       &                                                                                    powerPowerLawDraine
  double precision                                       , allocatable, dimension(:,:,:) :: emissivityPowerLawDraine
  double precision                                       , dimension(2)                  :: wavelengthsNodes              , luminosityNodes           , &
       &                                                                                    luminosityNodesMidpoint
  double precision                                       , dimension(countWavelengths)   :: luminosityNext                , luminosityBetween
  double precision                                                                       :: massDustPerHydrogen           , luminositySpecific        , &
       &                                                                                    massDustBalance               , intensityMean             , &
       &                                                                                    emissivityNodeA               , emissivityNodeB
  integer                                                                                :: indexAlpha                    , indexIntensity            , &
       &                                                                                    indexNodeA                    , indexNodeB
  type            (dustEmissionSpectrumList             ), pointer                       :: members
  double precision                                       , dimension(countWavelengths)   :: wavelengths                   , frequencies               , &
       &                                                                                    luminosity                    , luminosityCold            , &
       &                                                                                    luminosityWarm
  double precision                                                                       :: temperature                   , temperatureAbsorbed       , &
       &                                                                                    temperatureCMB                , timeHighRedshift          , &
       &                                                                                    emitted
  integer                                                                                :: i

  call displayVerbositySet(verbosityLevelStandard)
  call ioHDF5AccessInitialize()
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

  call Unit_Tests_Begin_Group("Dale & Helou (2002) templates")
  allocate(daleHelou_        )
  allocate(daleHelouMidpoint_)
  !![
  <referenceConstruct object="daleHelou_"         constructor="dustEmissionSpectrumDaleHelou2002(2.0d0    )"/>
  <referenceConstruct object="daleHelouMidpoint_" constructor="dustEmissionSpectrumDaleHelou2002(2.03125d0)"/>
  !!]
  fileDale=hdf5File(char(inputPath(pathTypeDataStatic)//'dust/emission/daleHelou2002.hdf5'),readOnly=.true.)
  call fileDale%readDataset('wavelength',wavelengthsDale)
  call fileDale%readDataset('alpha'     ,alphasDale     )
  call fileDale%readDataset('luminosity',luminosityDale )
  ! α=2 is tabulated, and α=2.03125 lies midway between it and the next template. Compare the spectra at tabulated
  ! wavelengths of 10 and 100 μm, where L_ν ∝ λ νL_ν.
  indexAlpha      =findloc(alphasDale,2.0d0,dim=1)
  indexNodeA      =minloc (abs(wavelengthsDale-1.0d5),dim=1)
  indexNodeB      =minloc (abs(wavelengthsDale-1.0d6),dim=1)
  wavelengthsNodes=[wavelengthsDale(indexNodeA),wavelengthsDale(indexNodeB)]
  luminosity             =daleHelou_        %luminosity(wavelengths     ,luminosityAbsorbed,massDust,timePresent)
  luminosityNodes        =daleHelou_        %luminosity(wavelengthsNodes,luminosityAbsorbed,massDust,timePresent)
  luminosityNodesMidpoint=daleHelouMidpoint_%luminosity(wavelengthsNodes,luminosityAbsorbed,massDust,timePresent)
  call Assert("spectrum integrates to the absorbed luminosity"      ,integral(luminosity)                                                                                         ,luminosityAbsorbed                                                                                                                                         ,relTol=toleranceTemplate)
  call Assert("no emission shortward of 3 μm or longward of 1100 μm",all(daleHelou_%luminosity([2.0d4,2.0d7],luminosityAbsorbed,massDust,timePresent) == 0.0d0)                   ,.true.                                                                                                                                                                              )
  call Assert("shape is that of the tabulated template"             ,luminosityNodes(1)/luminosityNodes(2)                                                                        ,(wavelengthsNodes(1)*luminosityDale(indexNodeA,indexAlpha))/(wavelengthsNodes(2)*luminosityDale(indexNodeB,indexAlpha))                                    ,relTol=1.0d-9           )
  call Assert("templates are interpolated geometrically in α"       ,luminosityNodesMidpoint(1)/luminosityNodesMidpoint(2)                                                        ,(wavelengthsNodes(1)/wavelengthsNodes(2))*sqrt(luminosityDale(indexNodeA,indexAlpha)*luminosityDale(indexNodeA,indexAlpha+1)/luminosityDale(indexNodeB,indexAlpha)/luminosityDale(indexNodeB,indexAlpha+1)),relTol=1.0d-9)
  call Assert("mass of dust is not used"                            ,all(daleHelou_%luminosity(wavelengthsNodes,luminosityAbsorbed,2.0d0*massDust,timePresent) == luminosityNodes),.true.                                                                                                                                                                              )
  call Unit_Tests_End_Group()

  call Unit_Tests_Begin_Group("Draine & Li (2007) models")
  allocate(draineLi_            )
  allocate(draineLiFixed_       )
  allocate(draineLiFixedNext_   )
  allocate(draineLiFixedBetween_)
  allocate(draineLiFixedHighest_)
  !![
  <referenceConstruct object="draineLi_"             constructor="dustEmissionSpectrumDraineLi2007(draineLi2007GrainModelMilkyWay60,fractionMassPowerLaw,     -1.0d0 ,intensityMaximum)"/>
  <referenceConstruct object="draineLiFixed_"        constructor="dustEmissionSpectrumDraineLi2007(draineLi2007GrainModelMilkyWay60,fractionMassPowerLaw,      1.0d0 ,intensityMaximum)"/>
  <referenceConstruct object="draineLiFixedNext_"    constructor="dustEmissionSpectrumDraineLi2007(draineLi2007GrainModelMilkyWay60,fractionMassPowerLaw,      1.2d0 ,intensityMaximum)"/>
  <referenceConstruct object="draineLiFixedBetween_" constructor="dustEmissionSpectrumDraineLi2007(draineLi2007GrainModelMilkyWay60,fractionMassPowerLaw,sqrt( 1.2d0),intensityMaximum)"/>
  <referenceConstruct object="draineLiFixedHighest_" constructor="dustEmissionSpectrumDraineLi2007(draineLi2007GrainModelMilkyWay60,fractionMassPowerLaw,     25.0d0 ,intensityMaximum)"/>
  !!]
  fileDraine =hdf5File(char(inputPath(pathTypeDataStatic)//'dust/emission/draineLi2007.hdf5'),readOnly=.true.)
  call fileDraine %readDataset  ('wavelength'         ,wavelengthsDraine       )
  call fileDraine %readDataset  ('intensityMinimum'   ,intensitiesDraine       )
  groupDraine=fileDraine%openGroup('milkyWay60')
  call groupDraine%readDataset  ('emissivitySingle'   ,emissivitySingleDraine  )
  call groupDraine%readDataset  ('emissivityPowerLaw' ,emissivityPowerLawDraine)
  call groupDraine%readDataset  ('powerSingle'        ,powerSingleDraine       )
  call groupDraine%readDataset  ('powerPowerLaw'      ,powerPowerLawDraine     )
  call groupDraine%readAttribute('massDustPerHydrogen',massDustPerHydrogen     )
  ! Compare with the tabulation at U_min=1, at two tabulated wavelengths, where L_ν ∝ λ ν dP/dν of the mixture of
  ! single-intensity and power-law heating.
  indexIntensity  =findloc(intensitiesDraine,1.0d0,dim=1)
  indexNodeA      =200
  indexNodeB      =800
  wavelengthsNodes=[wavelengthsDraine(indexNodeA),wavelengthsDraine(indexNodeB)]
  emissivityNodeA =(1.0d0-fractionMassPowerLaw)*emissivitySingleDraine(indexNodeA,indexIntensity)+fractionMassPowerLaw*emissivityPowerLawDraine(indexNodeA,indexIntensity,indexIntensityMaximum)
  emissivityNodeB =(1.0d0-fractionMassPowerLaw)*emissivitySingleDraine(indexNodeB,indexIntensity)+fractionMassPowerLaw*emissivityPowerLawDraine(indexNodeB,indexIntensity,indexIntensityMaximum)
  luminosity      =draineLiFixed_%luminosity(wavelengths     ,luminosityAbsorbed,massDust,timePresent)
  luminosityNodes =draineLiFixed_%luminosity(wavelengthsNodes,luminosityAbsorbed,massDust,timePresent)
  call Assert("spectrum integrates to the absorbed luminosity"   ,integral(luminosity)                                                                          ,luminosityAbsorbed                                                        ,relTol=toleranceTemplate)
  call Assert("no emission shortward of 1 μm or longward of 1 cm",all(draineLiFixed_%luminosity([0.9d4,1.1d8],luminosityAbsorbed,massDust,timePresent) == 0.0d0),.true.                                                                                             )
  call Assert("shape is that of the tabulated model"             ,luminosityNodes(1)/luminosityNodes(2)                                                         ,(wavelengthsNodes(1)*emissivityNodeA)/(wavelengthsNodes(2)*emissivityNodeB),relTol=1.0d-9          )
  ! Energy balance. Find the luminosity per unit dust mass at U_min=1 from the tabulated power per H nucleon, and so the
  ! mass of dust which, absorbing the given luminosity, is heated with U_min=1.
  luminositySpecific=+(                                                                                        &
       &               +(1.0d0-fractionMassPowerLaw)*powerSingleDraine  (indexIntensity                      ) &
       &               +       fractionMassPowerLaw *powerPowerLawDraine(indexIntensity,indexIntensityMaximum) &
       &              )                                                                                        &
       &             *ergs                                                                                     &
       &             /massDustPerHydrogen                                                                      &
       &             /massHydrogenAtom                                                                         &
       &             *massSolar                                                                                &
       &             /luminositySolar
  intensityMean     =(1.0d0-fractionMassPowerLaw)*1.0d0+fractionMassPowerLaw*log(intensityMaximum)/(1.0d0-1.0d0/intensityMaximum)
  massDustBalance   =luminosityAbsorbed/luminositySpecific
  luminosityNext    =draineLiFixedNext_   %luminosity(wavelengths,luminosityAbsorbed,massDust,timePresent)
  luminosityBetween =draineLiFixedBetween_%luminosity(wavelengths,luminosityAbsorbed,massDust,timePresent)
  call Assert("luminosity per unit dust mass per unit mean intensity is about 136 L☉/M☉",luminositySpecific/intensityMean > 125.0d0 .and. luminositySpecific/intensityMean < 145.0d0                                                                                                                             ,.true.              )
  call Assert("energy balance recovers the minimum intensity"                            ,draineLi_%intensityMinimum(luminosityAbsorbed,massDustBalance)                                                                                                                                                          , 1.0d0,relTol=1.0d-9)
  call Assert("energy balance gives the spectrum at that intensity"                      ,maxval(abs(draineLi_%luminosity(wavelengths,luminosityAbsorbed,massDustBalance,timePresent)-luminosity)) <= 1.0d-9*maxval(luminosity)                                                                                   ,.true.              )
  call Assert("spectrum between tabulated intensities interpolates linearly in ln U_min" ,maxval(abs(luminosityBetween-0.5d0*(luminosity+luminosityNext))) <= 1.0d-9*maxval(luminosity)                                                                                                                           ,.true.              )
  call Assert("weak heating is limited to the lowest tabulated intensity"                ,draineLi_%intensityMinimum(luminosityAbsorbed,1.0d10*massDustBalance)                                                                                                                                                   , 0.1d0,relTol=1.0d-9)
  call Assert("intense heating is limited to the highest tabulated intensity"            ,draineLi_%intensityMinimum(luminosityAbsorbed,1.0d-10*massDustBalance)                                                                                                                                                  ,25.0d0,relTol=1.0d-9)
  call Assert("intense heating has the spectrum of the highest tabulated intensity"      ,maxval(abs(draineLi_%luminosity(wavelengths,luminosityAbsorbed,1.0d-10*massDustBalance,timePresent)-draineLiFixedHighest_%luminosity(wavelengths,luminosityAbsorbed,massDust,timePresent))) <= 1.0d-9*maxval(luminosity),.true.              )
  call Assert("no absorption, no emission"                                               ,all(draineLi_%luminosity(wavelengths,0.0d0,massDust,timePresent) == 0.0d0)                                                                                                                                              ,.true.              )
  call Unit_Tests_End_Group()

  !![
  <objectDestructor name="draineLiFixedHighest_"/>
  <objectDestructor name="draineLiFixedBetween_"/>
  <objectDestructor name="draineLiFixedNext_"   />
  <objectDestructor name="draineLiFixed_"       />
  <objectDestructor name="draineLi_"            />
  <objectDestructor name="daleHelouMidpoint_"   />
  <objectDestructor name="daleHelou_"           />
  <objectDestructor name="sum_"                 />
  <objectDestructor name="greybodyCold_"        />
  <objectDestructor name="greybodyWarm_"        />
  <objectDestructor name="greybodyFixed_"       />
  <objectDestructor name="greybodyCMB_"         />
  <objectDestructor name="greybody_"            />
  <objectDestructor name="dustProperties_"      />
  <objectDestructor name="cosmologyFunctions_"  />
  <objectDestructor name="cosmologyParameters_" />
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
