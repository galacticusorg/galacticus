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
Contains a program which tests the SKIRTOR AGN torus dust attenuation and emission classes.
!!}

program Test_Dust_Attenuation_AGN_Torus_SKIRTOR
  !!{RST
  Tests the :galacticus-class:`dustAttenuationAGNTorusSKIRTOR` dust attenuation class and the
  :galacticus-class:`dustEmissionSpectrumAGNTorusSKIRTOR` dust emission spectrum class, using a small synthetic
  tabulation written here in the format of the SKIRTOR tabulation, so that the expected results are known exactly. It
  holds three models, differing in the angle of the torus, with the parameters requested, and a decoy model with a
  different optical depth, whose transmission is zero, which would be noticed were it selected.

  The checks are made against quantities computed here, independently of the classes:

  * the transmission averaged over orientation in the :math:`V` band, integrated on a fine grid of inclination, is one
    minus the covering fraction requested, whether that falls between models or on one;
  * the transmission at a tabulated wavelength and inclination is the tabulated direct plus scattered transmission of
    the model whose covering fraction is requested, and between them is interpolated linearly in inclination and in
    the logarithm of wavelength;
  * where the scattered light would make the transmission averaged over orientation exceed one, that average is one;
  * a covering fraction beyond the largest of the models gives the transmission of that model;
  * only the continuum of the accretion disk is attenuated;
  * the emission spectrum integrates to the luminosity absorbed, has the tabulated shape, and is zero outside the range
    of wavelength over which the model emits.
  !!}
  use :: Black_Hole_Torus_Covering_Fractions, only : blackHoleTorusCoveringFractionFixed
  use :: Display                            , only : displayVerbositySet                , verbosityLevelStandard
  use :: Dust_Attenuation_Descriptors       , only : emissionDescriptor                 , emissionSourceAccretionDisk      , emissionSourceNebular, emissionSourceStellar
  use :: Dust_Attenuations                  , only : dustAttenuationAGNTorusSKIRTOR
  use :: Dust_Emission_Spectra              , only : dustEmissionSpectrumAGNTorusSKIRTOR
  use :: Error                              , only : Error_Handler_Register
  use :: Events_Hooks                       , only : eventsHooksInitialize
  use :: File_Utilities                     , only : Directory_Make
  use :: Functions_Global_Utilities         , only : Functions_Global_Set
  use :: Galactic_Inclinations              , only : galacticInclinationFixed
  use :: Galactic_Structure_Options         , only : componentTypeBlackHole             , componentTypeDisk
  use :: Galacticus_Nodes                   , only : nodeClassHierarchyInitialize       , treeNode
  use :: HDF5_Access                        , only : hdf5Access
  use :: Input_Parameters                   , only : inputParameters
  use :: Input_Paths                        , only : inputPath                          , pathTypeDataDynamic
  use :: IO_HDF5                            , only : hdf5File
  use :: ISO_Varying_String                 , only : assignment(=)                      , char                             , operator(//)         , var_str              , &
       &                                             varying_string
  use :: Node_Components                    , only : Node_Components_Initialize         , Node_Components_Thread_Initialize
  use :: Numerical_Constants_Math           , only : Pi
  use :: Numerical_Constants_Physical       , only : speedLight
  use :: Numerical_Constants_Units          , only : metersToAngstroms
  use :: Unit_Tests                         , only : Assert                             , Unit_Tests_Begin_Group           , Unit_Tests_End_Group , Unit_Tests_Finish
  implicit none
  ! The synthetic tabulation: wavelengths (in μm), inclinations and angles of the torus (in degrees).
  integer                                              , parameter                                             :: countWavelength  =6             , countInclination    =10   , &
       &                                                                                                          countModel       =4
  double precision                                     , parameter, dimension(countWavelength                ) :: wavelengthTable  =[0.01d0,0.1d0,0.55d0,2.0d0,10.0d0,100.0d0]
  ! Optical depth of the edge of the torus at each wavelength, per 10° beyond the edge.
  double precision                                     , parameter, dimension(countWavelength)                 :: depthOpticalTable=[6.0d0,4.0d0,3.0d0,1.0d0,0.5d0,0.1d0]
  ! Scattered transmission at each wavelength: zero in the V band, and so large at the shortest wavelength that it must
  ! be limited.
  double precision                                     , parameter, dimension(countWavelength)                 :: scatteredTable   =[5.0d0,0.02d0,0.0d0,0.02d0,0.02d0,0.02d0]
  ! Orientation-averaged dust emission, νL_ν per unit luminosity: zero at the shortest wavelength.
  double precision                                     , parameter, dimension(countWavelength)                 :: emissionTable    =[0.0d0,1.0d-3,1.0d-2,0.3d0,0.5d0,0.1d0]
  double precision                                     , parameter, dimension(countModel     )                 :: angleTorusTable  =[20.0d0,40.0d0,60.0d0,40.0d0], &
       &                                                                                                          depthOpticalModel=[ 7.0d0, 7.0d0, 7.0d0, 3.0d0]
  ! Index of the V band in the synthetic tabulation.
  integer                                              , parameter                                              :: indexV            =3
  ! Number of points in the fine grids used for independent integrals.
  integer                                              , parameter                                              :: countFine         =9001        , countWavelengthFine =200001
  double precision                                     , parameter                                              :: luminosityAbsorbed=1.0d10
  type            (inputParameters                    )                                                         :: parameters
  type            (treeNode                           ), pointer                                                :: node
  type            (galacticInclinationFixed           ), pointer                                                :: galacticInclination_
  type            (blackHoleTorusCoveringFractionFixed), pointer                                                :: blackHoleTorusCoveringFraction_
  type            (dustAttenuationAGNTorusSKIRTOR     ), pointer                                                :: torus_
  type            (dustEmissionSpectrumAGNTorusSKIRTOR), pointer                                                :: emission_
  type            (emissionDescriptor                 ), dimension(5)                                           :: descriptors
  type            (varying_string                     )                                                         :: fileName                       , pathDirectory
  double precision                                     , dimension(5)                                           :: transmission
  double precision                                     , dimension(countWavelength,countInclination,countModel) :: transmissionDirect             , transmissionScattered
  double precision                                     , dimension(countWavelength,                 countModel) :: emissionDustMean
  double precision                                     , dimension(                countInclination           ) :: inclinationTable
  double precision                                     , dimension(3)                                           :: coveringFractionModel
  double precision                                     , dimension(countFine)                                   :: inclinationFine                , integrand           , &
       &                                                                                                           transmissionFine
  double precision                                     , dimension(countWavelengthFine)                         :: wavelengthsFine                , luminosity
  double precision                                                                                              :: coveringFraction               , expected            , &
       &                                                                                                           weightModel                    , weightWavelength    , &
       &                                                                                                           transmissionLower              , transmissionUpper   , &
       &                                                                                                           luminosityRatio                , integral
  double precision                                     , dimension(1)                                           :: luminosityOne
  integer                                                                                                       :: i                              , j                   , &
       &                                                                                                           k                              , m
  character       (len=96                             )                                                         :: label

  call displayVerbositySet(verbosityLevelStandard)
  call Error_Handler_Register()
  parameters=inputParameters('testSuite/parameters/dustAttenuationAGNTorus.xml')
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  node => treeNode()
  ! Build the synthetic tabulation. The torus is transparent above its edge, at 90°-angleTorus from the axis, and
  ! beyond it the optical depth grows linearly with inclination.
  do i=1,countInclination
     inclinationTable(i)=10.0d0*dble(i-1)
  end do
  do m=1,countModel
     do j=1,countInclination
        do k=1,countWavelength
           if (depthOpticalModel(m) /= 7.0d0) then
              transmissionDirect(k,j,m)=0.0d0
           else
              transmissionDirect(k,j,m)=exp(-depthOpticalTable(k)*max(0.0d0,inclinationTable(j)-(90.0d0-angleTorusTable(m)))/10.0d0)
           end if
           transmissionScattered(k,j,m)=scatteredTable(k)
        end do
     end do
     emissionDustMean(:,m)=emissionTable
  end do
  fileName     ='SKIRTOR_unitTest.hdf5'
  pathDirectory=inputPath(pathTypeDataDynamic)//'dust/agnTorus'
  call Directory_Make(pathDirectory)
  !$ call hdf5Access%set()
  hdf5WriteScope : block
    type(hdf5File) :: file
    file=hdf5File(char(pathDirectory//'/'//fileName),overWrite=.true.)
    call file%writeDataset(wavelengthTable            ,'wavelength'           )
    call file%writeDataset(inclinationTable           ,'inclination'          )
    call file%writeDataset(depthOpticalModel          ,'opticalDepth9.7'      )
    call file%writeDataset(spread( 1.0d0,1,countModel),'exponentRadial'       )
    call file%writeDataset(spread( 1.0d0,1,countModel),'exponentPolar'        )
    call file%writeDataset(spread(20.0d0,1,countModel),'ratioRadii'           )
    call file%writeDataset(angleTorusTable            ,'angleTorus'           )
    call file%writeDataset(transmissionDirect         ,'transmissionDirect'   )
    call file%writeDataset(transmissionScattered      ,'transmissionScattered')
    call file%writeDataset(emissionDustMean           ,'emissionDustMean'     )
  end block hdf5WriteScope
  !$ call hdf5Access%unset()
  ! Fine grid of inclination, in radians, which includes every tabulated inclination.
  do i=1,countFine
     inclinationFine(i)=0.5d0*Pi*dble(i-1)/dble(countFine-1)
  end do
  ! The covering fraction of each model, found independently: one minus the V-band direct transmission, interpolated
  ! linearly in inclination, averaged over orientation. The integrand is smooth except at the tabulated inclinations,
  ! which fall on even-numbered points of the fine grid, so Simpson's rule integrates it to high accuracy.
  do m=1,3
     do i=1,countFine
        integrand(i)=interpolateInclination(transmissionDirect(indexV,:,m),inclinationFine(i))*sin(inclinationFine(i))
     end do
     coveringFractionModel(m)=1.0d0-simpson(inclinationFine,integrand)
  end do
  ! Parcels: accretion disk continuum at V, at 2 μm (tabulated), and at 1 μm (between tabulated wavelengths), AGN line
  ! emission, and stellar emission.
  descriptors%wavelength   =[5.5d3                      ,2.0d4                      ,1.0d4                      ,5.5d3                 ,5.5d3                ]
  descriptors%componentType=[componentTypeBlackHole     ,componentTypeBlackHole     ,componentTypeBlackHole     ,componentTypeBlackHole,componentTypeDisk    ]
  descriptors%sourceType   =[emissionSourceAccretionDisk,emissionSourceAccretionDisk,emissionSourceAccretionDisk,emissionSourceNebular ,emissionSourceStellar]
  allocate(galacticInclination_)
  !![
  <referenceConstruct object="galacticInclination_" constructor="galacticInclinationFixed(60.0d0)"/>
  !!]

  ! Covering fractions between models, and on one.
  do k=1,3
     select case (k)
     case (1)
        coveringFraction=0.5d0*(coveringFractionModel(1)+coveringFractionModel(2))
     case (2)
        coveringFraction=0.3d0*coveringFractionModel(2)+0.7d0*coveringFractionModel(3)
     case (3)
        coveringFraction=coveringFractionModel(2)
     end select
     write (label,'(a,f6.4)') "Dust attenuation: SKIRTOR AGN torus, f_cov=",coveringFraction
     call Unit_Tests_Begin_Group(trim(label))
     allocate(blackHoleTorusCoveringFraction_)
     allocate(torus_                         )
     !![
     <referenceConstruct object="blackHoleTorusCoveringFraction_" constructor="blackHoleTorusCoveringFractionFixed(coveringFraction)"                                                                                 />
     <referenceConstruct object="torus_"                          constructor="dustAttenuationAGNTorusSKIRTOR(fileName,var_str('none'),7.0d0,1.0d0,1.0d0,20.0d0,galacticInclination_,blackHoleTorusCoveringFraction_)"/>
     !!]
     ! The V-band transmission averaged over orientation is one minus the covering fraction.
     do i=1,countFine
        transmission       =torus_%transmission(node,descriptors,inclination=inclinationFine(i))
        transmissionFine(i)=transmission(1)
        integrand       (i)=transmission(1)*sin(inclinationFine(i))
     end do
     call Assert('V-band absorption averaged over orientation equals the covering fraction',1.0d0-simpson(inclinationFine,integrand),coveringFraction,absTol=1.0d-9)
     ! Interpolation between models, in inclination, and in wavelength, checked against values found here.
     m          =max(1,min(2,count(coveringFractionModel <= coveringFraction)))
     weightModel=(coveringFraction-coveringFractionModel(m))/(coveringFractionModel(m+1)-coveringFractionModel(m))
     transmission=torus_%transmission(node,descriptors,inclination=50.0d0*Pi/180.0d0)
     expected=+(1.0d0-weightModel)*(transmissionDirect(4,6,m  )+transmissionScattered(4,6,m  )) &
          &   +       weightModel *(transmissionDirect(4,6,m+1)+transmissionScattered(4,6,m+1))
     call Assert('transmission at a tabulated wavelength and inclination',transmission(2),expected,relTol=1.0d-9)
     transmission    =torus_%transmission(node,descriptors,inclination=45.0d0*Pi/180.0d0)
     weightWavelength=log(1.0d0/0.55d0)/log(2.0d0/0.55d0)
     transmissionLower=+(1.0d0-weightWavelength)*(0.5d0*(transmissionDirect(3,5,m  )+transmissionDirect(3,6,m  ))+scatteredTable(3)) &
          &            +       weightWavelength *(0.5d0*(transmissionDirect(4,5,m  )+transmissionDirect(4,6,m  ))+scatteredTable(4))
     transmissionUpper=+(1.0d0-weightWavelength)*(0.5d0*(transmissionDirect(3,5,m+1)+transmissionDirect(3,6,m+1))+scatteredTable(3)) &
          &            +       weightWavelength *(0.5d0*(transmissionDirect(4,5,m+1)+transmissionDirect(4,6,m+1))+scatteredTable(4))
     expected=(1.0d0-weightModel)*transmissionLower+weightModel*transmissionUpper
     call Assert('transmission between tabulated wavelengths and inclinations',transmission(3),expected,relTol=1.0d-9)
     ! Only the accretion disk continuum is attenuated.
     call Assert('AGN line emission is not attenuated',transmission(4),1.0d0,absTol=0.0d0)
     call Assert('stellar emission is not attenuated' ,transmission(5),1.0d0,absTol=0.0d0)
     ! Scattered light is limited so that the transmission averaged over orientation does not exceed one.
     descriptors(1)%wavelength=1.0d2
     do i=1,countFine
        transmission=torus_%transmission(node,descriptors,inclination=inclinationFine(i))
        integrand(i)=transmission(1)*sin(inclinationFine(i))
     end do
     call Assert('transmission averaged over orientation is limited to one',simpson(inclinationFine,integrand),1.0d0,absTol=1.0d-9)
     descriptors(1)%wavelength=5.5d3
     !![
     <objectDestructor name="torus_"                         />
     <objectDestructor name="blackHoleTorusCoveringFraction_"/>
     !!]
     call Unit_Tests_End_Group()
  end do

  ! A covering fraction beyond the largest of the models gives the transmission of that model.
  call Unit_Tests_Begin_Group("Dust attenuation: SKIRTOR AGN torus, covering fraction beyond the models")
  allocate(blackHoleTorusCoveringFraction_)
  allocate(torus_                         )
  !![
  <referenceConstruct object="blackHoleTorusCoveringFraction_" constructor="blackHoleTorusCoveringFractionFixed(0.999d0)"                                                                                          />
  <referenceConstruct object="torus_"                          constructor="dustAttenuationAGNTorusSKIRTOR(fileName,var_str('none'),7.0d0,1.0d0,1.0d0,20.0d0,galacticInclination_,blackHoleTorusCoveringFraction_)"/>
  !!]
  transmission=torus_%transmission(node,descriptors,inclination=50.0d0*Pi/180.0d0)
  call Assert('transmission is that of the model of largest covering fraction',transmission(2),transmissionDirect(4,6,3)+transmissionScattered(4,6,3),relTol=1.0d-9)
  !![
  <objectDestructor name="torus_"                         />
  <objectDestructor name="blackHoleTorusCoveringFraction_"/>
  !!]
  call Unit_Tests_End_Group()

  ! The emission spectrum.
  call Unit_Tests_Begin_Group("Dust emission: SKIRTOR AGN torus")
  allocate(emission_)
  !![
  <referenceConstruct object="emission_" constructor="dustEmissionSpectrumAGNTorusSKIRTOR(fileName,var_str('none'),7.0d0,1.0d0,1.0d0,20.0d0,40.0d0)"/>
  !!]
  ! Integrate L_ν over frequency on a fine grid spanning the range over which the model emits, as ∫ νL_ν dln λ.
  do i=1,countWavelengthFine
     wavelengthsFine(i)=1.0d3*exp(log(1.0d6/1.0d3)*dble(i-1)/dble(countWavelengthFine-1))
  end do
  luminosity=emission_%luminosity(wavelengthsFine,luminosityAbsorbed,0.0d0,13.8d0)
  integral  =trapezoid(log(wavelengthsFine),luminosity*speedLight*metersToAngstroms/wavelengthsFine)
  call Assert('emission integrates to the luminosity absorbed',integral,luminosityAbsorbed,relTol=1.0d-6)
  ! The shape: νL_ν at the tabulated wavelengths is in the ratio of the tabulated values.
  luminosityOne  =emission_%luminosity([1.0d5],luminosityAbsorbed,0.0d0,13.8d0)
  luminosityRatio=luminosityOne(1)*speedLight*metersToAngstroms/1.0d5
  luminosityOne  =emission_%luminosity([1.0d6],luminosityAbsorbed,0.0d0,13.8d0)
  call Assert('emission has the tabulated shape',luminosityRatio/(luminosityOne(1)*speedLight*metersToAngstroms/1.0d6),emissionTable(5)/emissionTable(6),relTol=1.0d-9)
  luminosityOne=emission_%luminosity([1.0d2],luminosityAbsorbed,0.0d0,13.8d0)
  call Assert('no emission where the model does not emit'   ,luminosityOne(1),0.0d0,absTol=0.0d0)
  luminosityOne=emission_%luminosity([1.0d7],luminosityAbsorbed,0.0d0,13.8d0)
  call Assert('no emission beyond the tabulated wavelengths',luminosityOne(1),0.0d0,absTol=0.0d0)
  !![
  <objectDestructor name="emission_"/>
  !!]
  call Unit_Tests_End_Group()
  !![
  <objectDestructor name="galacticInclination_"/>
  !!]
  call Unit_Tests_Finish()

contains

  double precision function interpolateInclination(values,inclination) result(value)
    !!{RST
    Interpolate ``values``, tabulated at the inclinations of the synthetic tabulation, linearly in inclination (given in
    radians).
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: values
    double precision, intent(in   )               :: inclination
    double precision                              :: inclinationDegrees, weight
    integer                                       :: j

    inclinationDegrees=inclination*180.0d0/Pi
    j                 =max(1,min(countInclination-1,int(inclinationDegrees/10.0d0)+1))
    weight            =(inclinationDegrees-inclinationTable(j))/10.0d0
    value             =(1.0d0-weight)*values(j)+weight*values(j+1)
    return
  end function interpolateInclination

  double precision function simpson(x,y) result(integral)
    !!{RST
    Integrate ``y`` over ``x``, which must be uniformly spaced and an odd number of points, by Simpson's rule.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: x, y
    integer                                       :: n

    n       =size(x)
    integral=+(x(2)-x(1))/3.0d0       &
         &   *(                       &
         &     +      y(1)            &
         &     +4.0d0*sum(y(2:n-1:2)) &
         &     +2.0d0*sum(y(3:n-2:2)) &
         &     +      y(n)            &
         &    )
    return
  end function simpson

  double precision function trapezoid(x,y) result(integral)
    !!{RST
    Integrate ``y`` over ``x`` by the trapezoidal rule.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: x, y

    integral=sum(0.5d0*(y(2:size(y))+y(1:size(y)-1))*(x(2:size(x))-x(1:size(x)-1)))
    return
  end function trapezoid

end program Test_Dust_Attenuation_AGN_Torus_SKIRTOR
