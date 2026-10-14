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
Contains a program which tests the AGN torus dust attenuation and emission classes.
!!}

program Test_Dust_Attenuation_AGN_Torus
  !!{RST
  Tests the :galacticus-class:`dustAttenuationAGNTorus` dust attenuation class and the
  :galacticus-class:`dustEmissionSpectrumAGNTorus` dust emission spectrum class.

  The checks are made against integrals computed here, independently of the classes:

  * the angular width found for a required covering fraction gives that covering fraction, integrating
    :math:`\{1-\exp[-N(\beta)]\}\cos\beta` over :math:`\beta` by the trapezoidal rule on a fine grid;
  * for clouds so opaque that a ray meeting one is absorbed completely, the absorption averaged over :math:`\cos i` is
    the covering fraction, again integrated on a fine grid, and the transmission at each inclination is
    :math:`\exp[-N(\beta)]`;
  * the absorption averaged over orientation by the 8-point rule used by default when finding the luminosity absorbed
    by dust is close to the exact covering fraction, which it would not be for a torus with a sharp edge;
  * only the continuum of the accretion disk is attenuated;
  * the emission spectrum integrates to the luminosity absorbed, and peaks between the wavelengths at which blackbodies
    at the minimum and sublimation temperatures peak.
  !!}
  use :: Black_Hole_Torus_Covering_Fractions, only : blackHoleTorusCoveringFractionFixed
  use :: Display                            , only : displayVerbositySet                , verbosityLevelStandard
  use :: Dust_Attenuation_Descriptors       , only : emissionDescriptor                 , emissionSourceAccretionDisk      , emissionSourceNebular, emissionSourceStellar
  use :: Dust_Attenuations                  , only : dustAttenuationAGNTorus            , gaussLegendreRule
  use :: Dust_Emission_Spectra              , only : dustEmissionSpectrumAGNTorus
  use :: Dust_Extinction_Curves             , only : dustExtinctionCurvePowerLaw
  use :: Error                              , only : Error_Handler_Register
  use :: Events_Hooks                       , only : eventsHooksInitialize
  use :: Functions_Global_Utilities         , only : Functions_Global_Set
  use :: Galactic_Inclinations              , only : galacticInclinationFixed
  use :: Galactic_Structure_Options         , only : componentTypeBlackHole             , componentTypeDisk
  use :: Galacticus_Nodes                   , only : nodeClassHierarchyInitialize       , treeNode
  use :: Input_Parameters                   , only : inputParameters
  use :: Node_Components                    , only : Node_Components_Initialize         , Node_Components_Thread_Initialize
  use :: Numerical_Constants_Math           , only : Pi
  use :: Numerical_Constants_Physical       , only : speedLight
  use :: Numerical_Constants_Units          , only : metersToAngstroms
  use :: Unit_Tests                         , only : Assert                             , Unit_Tests_Begin_Group           , Unit_Tests_End_Group , Unit_Tests_Finish
  implicit none
  ! Mean number of clouds on an equatorial ray, and the optical depths of the clouds: typical, and so opaque that a ray
  ! meeting one is absorbed completely.
  double precision                                     , parameter                         :: countClouds          =5.0d0  , depthOpticalTypical   =30.0d0, &
       &                                                                                      depthOpticalOpaque   =1.0d3
  ! The V-band wavelength, which is the reference wavelength of the extinction curve.
  double precision                                     , parameter                         :: wavelengthV          =5.5d3
  ! The covering fractions tested.
  double precision                                     , parameter, dimension(4)           :: coveringFractions    =[0.2d0,0.5d0,0.8d0,0.95d0]
  ! Number of points in the fine grids used for independent integrals.
  integer                                              , parameter                         :: countFine            =200001
  ! Properties of the emission spectrum tested.
  double precision                                     , parameter                         :: luminosityAbsorbed   =1.0d10  , temperatureSublimation=1.5d3, &
       &                                                                                      temperatureMinimum   =1.0d2
  ! The second radiation constant, hc/k in Å K, and the value of x=hν/kT at which νL_ν of a modified blackbody of
  ! emissivity index β=1 peaks: the root of x=(4+β)[1-exp(-x)].
  double precision                                     , parameter                         :: radiationConstant    =1.4388d8, argumentPeak          =4.965d0
  integer                                              , parameter                         :: countWavelengths     =20000
  double precision                                     , parameter                         :: wavelengthMinimum    =1.0d3   , wavelengthMaximum     =1.0d8
  type            (inputParameters                    )                                    :: parameters
  type            (treeNode                           ), pointer                           :: node
  type            (dustExtinctionCurvePowerLaw        ), pointer                           :: dustExtinctionCurve_
  type            (galacticInclinationFixed           ), pointer                           :: galacticInclination_
  type            (blackHoleTorusCoveringFractionFixed), pointer                           :: blackHoleTorusCoveringFraction_
  type            (dustAttenuationAGNTorus            ), pointer                           :: torus_                         , torusOpaque_
  type            (dustEmissionSpectrumAGNTorus       ), pointer                           :: emission_
  type            (emissionDescriptor                 ), dimension(4)                      :: descriptors
  double precision                                     , dimension(4)                      :: transmission
  double precision                                     , allocatable  , dimension(:  )     :: cosineInclination              , weight
  double precision                                     , allocatable  , dimension(:,:)     :: absorbed
  double precision                                     , dimension(countFine)              :: angle                          , integrand
  double precision                                     , dimension(countWavelengths)       :: wavelengths                    , frequencies         , &
       &                                                                                      luminosity
  double precision                                                                         :: widthAngular                   , coveringFraction    , &
       &                                                                                      absorbedAverage                , inclination         , &
       &                                                                                      wavelengthPeak
  integer                                                                                  :: i                              , j
  character       (len=64                             )                                    :: label

  call displayVerbositySet(verbosityLevelStandard)
  call Error_Handler_Register()
  parameters=inputParameters('testSuite/parameters/dustAttenuationAGNTorus.xml')
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  node => treeNode()
  ! Angles from the equatorial plane, from 0 to π/2, for independent integrals.
  do i=1,countFine
     angle(i)=0.5d0*Pi*dble(i-1)/dble(countFine-1)
  end do
  ! The 8-point Gauss-Legendre rule in cos(i), the default used when finding the luminosity absorbed by dust.
  call gaussLegendreRule(8,cosineInclination,weight)
  ! Parcels: accretion disk continuum at V, accretion disk continuum at 10 μm, AGN line emission, and stellar emission.
  descriptors%wavelength   =[wavelengthV                ,1.0d5                      ,wavelengthV           ,wavelengthV          ]
  descriptors%componentType=[componentTypeBlackHole     ,componentTypeBlackHole     ,componentTypeBlackHole,componentTypeDisk    ]
  descriptors%sourceType   =[emissionSourceAccretionDisk,emissionSourceAccretionDisk,emissionSourceNebular ,emissionSourceStellar]
  allocate(dustExtinctionCurve_)
  allocate(galacticInclination_)
  !![
  <referenceConstruct object="dustExtinctionCurve_" constructor="dustExtinctionCurvePowerLaw( 0.7d0,wavelengthV)"/>
  <referenceConstruct object="galacticInclination_" constructor="galacticInclinationFixed   (60.0d0            )"/>
  !!]

  do j=1,size(coveringFractions)
     write (label,'(a,f4.2)') "Dust attenuation: AGN torus, f_cov=",coveringFractions(j)
     call Unit_Tests_Begin_Group(trim(label))
     allocate(blackHoleTorusCoveringFraction_)
     allocate(torus_                         )
     allocate(torusOpaque_                   )
     !![
     <referenceConstruct object="blackHoleTorusCoveringFraction_" constructor="blackHoleTorusCoveringFractionFixed(coveringFractions(j))"                                                                         />
     <referenceConstruct object="torus_"                          constructor="dustAttenuationAGNTorus(countClouds,depthOpticalTypical,dustExtinctionCurve_,galacticInclination_,blackHoleTorusCoveringFraction_)"/>
     <referenceConstruct object="torusOpaque_"                    constructor="dustAttenuationAGNTorus(countClouds,depthOpticalOpaque ,dustExtinctionCurve_,galacticInclination_,blackHoleTorusCoveringFraction_)"/>
     !!]
     ! The angular width gives the required covering fraction.
     widthAngular    =torus_%widthAngular(node)
     integrand       =(1.0d0-exp(-countClouds*exp(-(angle/widthAngular)**2)))*cos(angle)
     coveringFraction=trapezoid(angle,integrand)
     call Assert("angular width gives the covering fraction"                        ,coveringFraction,coveringFractions(j),absTol=1.0d-4)
     ! Transmission of opaque clouds at several inclinations.
     do i=0,3
        inclination =dble(i)*Pi/6.0d0
        transmission=torusOpaque_%transmission(node,descriptors,inclination)
        write (label,'(a,i2,a)') "opaque clouds, i=",30*i,"°: accretion disk at V"
        call Assert(trim(label),transmission(1),exp(-countClouds*exp(-((0.5d0*Pi-inclination)/widthAngular)**2)),relTol=1.0d-6)
     end do
     ! Opaque clouds absorb, averaged over cos(i), the covering fraction. With cos(i)=sin(β), the average over cos(i)
     ! is the integral over β weighted by cos(β), so evaluate it on the fine grid of β.
     do i=1,countFine
        transmission=torusOpaque_%transmission(node,descriptors,0.5d0*Pi-angle(i))
        integrand(i)=(1.0d0-transmission(1))*cos(angle(i))
     end do
     absorbedAverage=trapezoid(angle,integrand)
     call Assert("opaque clouds absorb the covering fraction, averaged over cos(i)" ,absorbedAverage,coveringFractions(j),absTol=1.0d-4)
     ! The 8-point rule used for the absorbed luminosity is close to the exact average.
     absorbed=torusOpaque_%absorbedFractions(node,descriptors,cosineInclination,weight)
     call Assert("8-point orientation average of absorption is accurate"            ,absorbed(1,1),coveringFractions(j),absTol=5.0d-3)
     ! Only the accretion disk continuum is attenuated, and less so in the infrared.
     transmission=torus_%transmission(node,descriptors,0.5d0*Pi)
     call Assert("edge-on: accretion disk at V is attenuated"                       ,transmission(1) < 1.0d0          ,.true.               )
     call Assert("edge-on: accretion disk at 10μm is less attenuated than at V"     ,transmission(2) > transmission(1),.true.               )
     call Assert("edge-on: AGN line emission is not attenuated"                     ,transmission(3)                  ,1.0d0 ,relTol=1.0d-12)
     call Assert("edge-on: stellar emission is not attenuated"                      ,transmission(4)                  ,1.0d0 ,relTol=1.0d-12)
     ! Without an explicit inclination that of the galacticInclination object (60°) is used.
     call Assert("inclination is taken from the galacticInclination object"        ,torus_%transmission(node,descriptors),torus_%transmission(node,descriptors,Pi/3.0d0),relTol=1.0d-12)
     !![
     <objectDestructor name="torus_"                         />
     <objectDestructor name="torusOpaque_"                   />
     <objectDestructor name="blackHoleTorusCoveringFraction_"/>
     !!]
     call Unit_Tests_End_Group()
  end do

  call Unit_Tests_Begin_Group("Dust emission: AGN torus")
  do i=1,countWavelengths
     wavelengths(i)=exp(log(wavelengthMinimum)+log(wavelengthMaximum/wavelengthMinimum)*dble(i-1)/dble(countWavelengths-1))
  end do
  frequencies=speedLight*metersToAngstroms/wavelengths
  allocate(emission_)
  !![
  <referenceConstruct object="emission_" constructor="dustEmissionSpectrumAGNTorus(temperatureSublimation,temperatureMinimum,0.0d0,1.0d0)"/>
  !!]
  luminosity=emission_%luminosity(wavelengths,luminosityAbsorbed,0.0d0,13.8d0)
  call Assert("spectrum integrates to the absorbed luminosity"                      ,-trapezoid(frequencies,luminosity),luminosityAbsorbed,relTol=1.0d-5)
  wavelengthPeak=wavelengths(maxloc(luminosity*frequencies,dim=1))
  call Assert("νL_ν peaks between the hottest and coolest dust"                     ,wavelengthPeak > radiationConstant/argumentPeak/temperatureSublimation .and. wavelengthPeak < radiationConstant/argumentPeak/temperatureMinimum,.true.)
  call Assert("spectrum scales with the absorbed luminosity"                        ,emission_%luminosity(wavelengths,2.0d0*luminosityAbsorbed,0.0d0,13.8d0),2.0d0*luminosity,relTol=1.0d-12)
  !![
  <objectDestructor name="emission_"/>
  !!]
  call Unit_Tests_End_Group()
  !![
  <objectDestructor name="dustExtinctionCurve_"/>
  <objectDestructor name="galacticInclination_"/>
  !!]
  call Unit_Tests_Finish()

contains

  double precision function trapezoid(x,y)
    !!{RST
    Return the integral of ``y`` over ``x`` by the trapezoidal rule.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: x, y

    trapezoid=sum(0.5d0*(y(2:)+y(:size(y)-1))*(x(2:)-x(:size(x)-1)))
    return
  end function trapezoid

end program Test_Dust_Attenuation_AGN_Torus
