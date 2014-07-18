!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a program which computes the projected correlation function for galaxies in a given
!% mass range based on a halo occupation model.

program Projected_Correlation_Function
  !% Computes the projected correlation function for galaxies in a given mass range based on a
  !% halo occupation model.
  use, intrinsic :: ISO_C_Binding
  use               FGSL
  use               Command_Arguments
  use               ISO_Varying_String
  use               Memory_Management
  use               Input_Parameters
  use               Galacticus_Error
  use               Geometry_Surveys
  use               Galacticus_Nodes
  use               Cosmology_Functions
  use               Conditional_Mass_Functions
  use               Numerical_Integration
  use               Dark_Matter_Profiles_Concentration
  use               Node_Component_Dark_Matter_Profile_Scale
  use               Numerical_Ranges
  use               Numerical_Constants_Math
  use               Power_Spectra
  use               FFTLogs
  use               Tables
  use               IO_HDF5
  implicit none
  class           (cosmologyFunctionsClass            ), pointer                   :: cosmologyFunctions_
  class           (conditionalMassFunctionClass       ), pointer                   :: conditionalMassFunction_
  class           (surveyGeometryClass                ), pointer                   :: surveyGeometry_
  class           (darkMatterProfileConcentrationClass), pointer                   :: darkMatterProfileConcentration_
  type            (treeNode                           ), pointer                   :: thisNode
  class           (nodeComponentBasic                 ), pointer                   :: thisBasic
  class           (nodeComponentDarkMatterProfile     ), pointer                   :: thisDarkMatterProfile
  procedure       (projectionIntegrandWeight          ), pointer                   :: integrandWeightFunction
  double precision                                     , allocatable, dimension(:) :: powerSpectrumOneHalo                              , powerSpectrumTwoHalo                              , &
       &                                                                              wavenumber                                        , powerSpectrum                                     , &
       &                                                                              separation                                        , correlation                                       , &
       &                                                                              projectedCorrelation                              , projectedCorrelationBinned                        , &
       &                                                                              projectedSeparationBinned
  integer                                              , parameter                 :: wavenumberCountPerDecade                   =10
  double precision                                     , parameter                 :: wavenumberMinimum                          =1.0d-3, wavenumberMaximum                          =1.0d3
  type            (varying_string                     )                            :: parameterFile                                     , projectedCorrelationFunctionOutputFileName
  double precision                                                                 :: projectedCorrelationFunctionMassMinimum           , projectedCorrelationFunctionMassMaximum           , &
       &                                                                              projectedCorrelationFunctionHaloMassMinimum       , projectedCorrelationFunctionHaloMassMaximum       , &
       &                                                                              time                                              , timeMinimum                                       , &
       &                                                                              expansionFactor                                   , volume                                            , &
       &                                                                              timeMaximum                                       , galaxyDensity                                     , &
       &                                                                              projectedSeparation                               , projectedCorrelationFunctionLineOfSightDepth      , &
       &                                                                              projectedCorrelationFunctionSeparationMinimum     , projectedCorrelationFunctionSeparationMaximum     , &
       &                                                                              binWidthLogarithmic                               , binSeparationMinimum                              , &
       &                                                                              binSeparationMaximum
  integer                                                                          :: iField                                            , iWavenumber                                       , &
       &                                                                              wavenumberCount                                   , iSeparation                                       , &
       &                                                                              projectedCorrelationFunctionSeparationCount
  type            (fgsl_function                      )                            :: integrandFunction
  type            (fgsl_integration_workspace         )                            :: integrationWorkspace
  type            (c_ptr                              )                            :: parameterPointer
  logical                                                                          :: integrationReset
  type            (table1DLogarithmicLinear           )                            :: correlationTable
  type            (hdf5Object                         )                            :: outputFile

  ! Read in basic code memory usage.
  call Code_Memory_Usage('Projected_Correlation_Function.size')
  ! Check that correct number of arguments have been supplied.
  if (Command_Argument_Count() /= 1) call Galacticus_Error_Report(message="Usage: Projected_Correlation_Function.exe <parameterFile>")
  ! Get the name of the parameter file from the first command line argument.
  call Get_Argument              (1,parameterFile)
  ! Open the parameter file.
  call Input_Parameters_File_Open(  parameterFile)
  ! Read parameters controlling the calculation.
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionSeparationMinimum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The minimum separation at which to compute the projected correlation function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionSeparationMinimum',projectedCorrelationFunctionSeparationMinimum)
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionSeparationMaximum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The maximum separation at which to compute the projected correlation function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionSeparationMaximum',projectedCorrelationFunctionSeparationMaximum)
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionSeparationCount</name>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The number of separations at which to compute the projected correlation function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionSeparationCount',projectedCorrelationFunctionSeparationCount)
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionLineOfSightDepth</name>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The maximum line of sight depth to which to integrate when computing the projected correlation function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionLineOfSightDepth',projectedCorrelationFunctionLineOfSightDepth)
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionOutputFileName</name>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The name of the file to which to output projected correlation function data.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionOutputFileName',projectedCorrelationFunctionOutputFileName)
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionMassMinimum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$10^8M_\odot$</defaultValue>
  !@   <description>
  !@     The minimum mass of galaxies to include in the projected correlation function calculation.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionMassMinimum',projectedCorrelationFunctionMassMinimum,defaultValue=1.0d8)
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionMassMaximum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$10^{12}M_\odot$</defaultValue>
  !@   <description>
  !@     The maximum mass of galaxies to include in the projected correlation function calculation.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionMassMaximum',projectedCorrelationFunctionMassMaximum,defaultValue=1.0d12)
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionHaloMassMinimum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$10^6M_\odot$</defaultValue>
  !@   <description>
  !@     The minimum halo mass to use when integrating over the halo mass function.
  !@   </description>
  !@   <type>real</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionHaloMassMinimum',projectedCorrelationFunctionHaloMassMinimum,defaultValue=1.0d6)
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionHaloMassMaximum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$10^{16}M_\odot$</defaultValue>
  !@   <description>
  !@     The maximum halo mass to use when integrating over the halo mass function.
  !@   </description>
  !@   <type>real</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionHaloMassMaximum',projectedCorrelationFunctionHaloMassMaximum,defaultValue=1.0d16)  
  ! Ensure the nodes objects are initialized.
  call Galacticus_Nodes_Initialize                        ()
  call Node_Component_Dark_Matter_Profile_Scale_Initialize()
  ! Get the default cosmology functions, conditional mass function, and survey geometry objects.
  cosmologyFunctions_             => cosmologyFunctions            ()     
  conditionalMassFunction_        => conditionalMassFunction       ()
  surveyGeometry_                 => surveyGeometry                ()
  darkMatterProfileConcentration_ => darkMatterProfileConcentration()
  ! Create worker node.
  thisNode              => treeNode                  (                 )
  thisBasic             => thisNode%basic            (autoCreate=.true.)
  thisDarkMatterProfile => thisNode%darkMatterProfile(autoCreate=.true.)
  select type (thisDarkMatterProfile)
  type is (nodeComponentDarkMatterProfileScale)
     ! This is acceptable.
  class default
     ! This is not.
     call Galacticus_Error_Report('Projected_Correlation_Function','this code expects to use the "scale" dark matter profile component')
  end select
  ! Generate wavenumber range.
  wavenumberCount=int(log10(wavenumberMaximum/wavenumberMinimum)*dble(wavenumberCountPerDecade))+1
  wavenumber     =Make_Range(wavenumberMinimum,wavenumberMaximum,wavenumberCount,rangeTypeLogarithmic)
  call Alloc_Array(powerSpectrum       ,[wavenumberCount])
  call Alloc_Array(powerSpectrumOneHalo,[wavenumberCount])
  call Alloc_Array(powerSpectrumTwoHalo,[wavenumberCount])
  ! Iterate over wavenumbers.
  do iWavenumber=1,wavenumberCount
     ! Iterate over survey fields.
     volume                           =0.0d0
     galaxyDensity                    =0.0d0
     powerSpectrumOneHalo(iWavenumber)=0.0d0
     powerSpectrumTwoHalo(iWavenumber)=0.0d0
     do iField=1,surveyGeometry_%fieldCount()
        ! Find time range for volume-limited sample.
        timeMinimum=cosmologyFunctions_%timeAtDistanceComoving(surveyGeometry_%distanceMaximum(projectedCorrelationFunctionMassMinimum,iField))
        timeMaximum=cosmologyFunctions_%timeAtDistanceComoving(surveyGeometry_%distanceMinimum(projectedCorrelationFunctionMassMinimum,iField))
        ! Integrate the volume term.
        integrationReset    =.true.
        volume              =                                       &
             & +volume                                              &
             & +Integrate(                                          &
             &            timeMinimum                             , &
             &            timeMaximum                             , &
             &            volumeTimeIntegrand                     , &
             &            parameterPointer                        , &
             &            integrandFunction                       , &
             &            integrationWorkspace                    , &
             &            toleranceRelative                =1.0d-3, &
             &            reset=integrationReset                    &
             &           )                                          &
             & *surveyGeometry_%solidAngle(iField)
        call Integrate_Done(integrandFunction,integrationWorkspace)
        ! Integrate the normalization term.
        integrationReset    =.true.
        galaxyDensity       =                                       &
             & +galaxyDensity                                       &
             & +Integrate(                                          &
             &            timeMinimum                             , &
             &            timeMaximum                             , &
             &            normalizationTimeIntegrand              , &
             &            parameterPointer                        , &
             &            integrandFunction                       , &
             &            integrationWorkspace                    , &
             &            toleranceRelative                =1.0d-3, &
             &            reset=integrationReset                    &
             &           )                                          &
             & *surveyGeometry_%solidAngle(iField)
        call Integrate_Done(integrandFunction,integrationWorkspace)
        ! Integrate the one-halo term.
        integrationReset    =.true.
        powerSpectrumOneHalo        (iWavenumber)=                  &
             & +powerSpectrumOneHalo(iWavenumber)                   &
             & +Integrate(                                          &
             &            timeMinimum                             , &
             &            timeMaximum                             , &
             &            powerSpectrumOneHaloTimeIntegrand       , &
             &            parameterPointer                        , &
             &            integrandFunction                       , &
             &            integrationWorkspace                    , &
             &            toleranceRelative                =1.0d-3, &
             &            reset=integrationReset                    &
             &           )                                          &
             & *surveyGeometry_%solidAngle(iField)
        call Integrate_Done(integrandFunction,integrationWorkspace)
        ! Integrate the two-halo term.
        integrationReset    =.true.
        powerSpectrumTwoHalo        (iWavenumber)=                  &
             & +powerSpectrumTwoHalo(iWavenumber)                   &
             & +Integrate(                                          &
             &            timeMinimum                             , &
             &            timeMaximum                             , &
             &            powerSpectrumTwoHaloTimeIntegrand       , &
             &            parameterPointer                        , &
             &            integrandFunction                       , &
             &            integrationWorkspace                    , &
             &            toleranceRelative                =1.0d-3, &
             &            reset=integrationReset                    &
             &           )                                          &
             & *surveyGeometry_%solidAngle(iField)
        call Integrate_Done(integrandFunction,integrationWorkspace)
     end do
     ! Construct the net power spectrum.
     galaxyDensity                    =galaxyDensity                    /volume
     powerSpectrumOneHalo(iWavenumber)=powerSpectrumOneHalo(iWavenumber)/volume
     powerSpectrumTwoHalo(iWavenumber)=powerSpectrumTwoHalo(iWavenumber)/volume
     powerSpectrum(iWavenumber)=(                                                  &
          &                      +powerSpectrumOneHalo(           iWavenumber )    &
          &                      +Power_Spectrum      (wavenumber(iWavenumber))    &
          &                      *powerSpectrumTwoHalo(           iWavenumber )**2 &
          &                     )                                                  &
          &                     /galaxyDensity                                 **2
  end do
  ! Fourier transform to get the correlation function.
  call Alloc_Array(correlation,shape(wavenumber))
  call Alloc_Array(separation ,shape(wavenumber))
  call FFTLog(                                    &
       &      wavenumber, &
       &      separation                        , &
       &      +powerSpectrum                 &
       &      *wavenumber &
       &      * 4.0d0*Pi                          &
       &      /(2.0d0*Pi)**3                    , &
       &      correlation                  , &
       &      fftLogSine                        , &
       &      fftLogForward                       &
       &     )
  correlation=correlation/separation
  ! Project the correlation function.
  call Alloc_Array(projectedCorrelation,shape(wavenumber))
  call correlationTable%create(separation(1),separation(wavenumberCount),wavenumberCount,extrapolationTypeExtrapolate)
  integrandWeightFunction => projectionIntegrandWeight
  do iSeparation=1,wavenumberCount
     projectedSeparation=separation(iSeparation)
     projectedCorrelation          (iSeparation)=sum(                                                                                            &
          &                                           correlationTable%integrationWeights(                                                       &
          &                                                                               projectedSeparation                                  , &
          &                                                                               sqrt(                                                  &
          &                                                                                    +projectedSeparation                         **2  &
          &                                                                                    +projectedCorrelationFunctionLineOfSightDepth**2  &
          &                                                                                   )                                                , &
          &                                                                               integrandWeightFunction                                &
          &                                                                              )                                                       &
          &                                          *correlation                                                                                &
          &                                         )
  end do
  ! Average the projected correlation function into bins.
  call Alloc_Array(projectedSeparationBinned ,[projectedCorrelationFunctionSeparationCount])
  call Alloc_Array(projectedCorrelationBinned,[projectedCorrelationFunctionSeparationCount])
  projectedSeparationBinned=Make_Range(projectedCorrelationFunctionSeparationMinimum,projectedCorrelationFunctionSeparationMaximum,projectedCorrelationFunctionSeparationCount,rangeTypeLogarithmic)
  binWidthLogarithmic=log(projectedSeparationBinned(2)/projectedSeparationBinned(1))
  integrandWeightFunction => binningIntegrandWeight
  do iSeparation=1,projectedCorrelationFunctionSeparationCount
     binSeparationMinimum         =projectedSeparationBinned(iSeparation)*exp(-0.5d0*binWidthLogarithmic)
     binSeparationMaximum         =projectedSeparationBinned(iSeparation)*exp(+0.5d0*binWidthLogarithmic)
     projectedCorrelationBinned(iSeparation)=sum(                                                              &
          &                                       correlationTable%integrationWeights(                         &
          &                                                                           binSeparationMinimum   , &
          &                                                                           binSeparationMaximum   , &
          &                                                                           integrandWeightFunction  &
          &                                                                          )                         &
          &                                       /Pi                                                          &
          &                                       /(                                                           &
          &                                         +binSeparationMaximum**2                                   &
          &                                         -binSeparationMinimum**2                                   &
          &                                       )                                                            &
          &                                      *projectedCorrelation                                         &
          &                                     )
  end do
  ! Write the data to file.
  call outputFile%openFile(char(projectedCorrelationFunctionOutputFileName))
  call outputFile%writeDataset(projectedSeparationBinned ,"separation"          ,commentText="Projected separation in units of Mpc." )
  call outputFile%writeDataset(projectedCorrelationBinned,"projectedCorrelation",commentText="Projected correlation in units of Mpc.")
  call outputFile%close()
  ! Close the parameter file.
  call Input_Parameters_File_Close

contains  

  double precision function projectionIntegrandWeight(separation)
    !% The weight function applied to the correlation function when integrating to get the projected correlation function.
    implicit none
    double precision, intent(in   ) :: separation
    
    if (separation > projectedSeparation) then
       projectionIntegrandWeight=2.0d0*separation/sqrt(separation**2-projectedSeparation**2)
    else
       projectionIntegrandWeight=0.0d0
    end if
    return
  end function projectionIntegrandWeight
    
    double precision function binningIntegrandWeight(separation)
      !% The weight function applied to the projected correlation function when integrating into bins.
      implicit none
      double precision, intent(in   ) :: separation

      binningIntegrandWeight=2.0d0*Pi*separation
      return
    end function binningIntegrandWeight

  function powerSpectrumOneHaloTimeIntegrand(timePrime,parameterPointer) bind(c)
    !% Time integrand for the one-halo term in the power spectrum.
    implicit none
    real   (c_double                  )        :: powerSpectrumOneHaloTimeIntegrand
    real   (c_double                  ), value :: timePrime
    type   (c_ptr                     ), value :: parameterPointer
    type   (fgsl_function             )        :: integrandFunctionTime
    type   (fgsl_integration_workspace)        :: integrationWorkspaceTime
    logical                                    :: integrationResetTime
    
    time           =timePrime
    expansionFactor=cosmologyFunctions_%expansionFactor(time)
    call thisBasic%timeSet                             (time)
    call thisBasic%timeLastIsolatedSet                 (time)
    integrationResetTime             =.true.
    powerSpectrumOneHaloTimeIntegrand=                                    &
         & +Integrate(                                                    &
         &            projectedCorrelationFunctionHaloMassMinimum       , &
         &            projectedCorrelationFunctionHaloMassMaximum       , &
         &            powerSpectrumOneHaloIntegrand                     , &
         &            parameterPointer                                  , &
         &            integrandFunctionTime                             , &
         &            integrationWorkspaceTime                          , &
         &            toleranceRelative                          =1.0d-3, &
         &            reset=integrationResetTime                          &
         &           )                                                    &
         & *cosmologyFunctions_%comovingVolumeElementTime(time)
    call Integrate_Done(integrandFunctionTime,integrationWorkspaceTime)
    return
  end function powerSpectrumOneHaloTimeIntegrand
  
  function powerSpectrumOneHaloIntegrand(massHalo,parameterPointer) bind(c)
    !% Integrand for the one-halo term in the power spectrum.
    use Dark_Matter_Halo_Biases
    use Dark_Matter_Profiles
    use Halo_Mass_Function
    use Dark_Matter_Halo_Scales
    use Galacticus_Calculations_Resets
    implicit none
    real(c_double)        :: powerSpectrumOneHaloIntegrand
    real(c_double), value :: massHalo
    type(c_ptr   ), value :: parameterPointer

    call Galacticus_Calculations_Reset(thisNode)
    call thisBasic            % massSet(                                                         &
         &                              massHalo                                                 &
         &                             )
    call thisDarkMatterProfile%scaleSet(                                                         &
         &                               Dark_Matter_Halo_Virial_Radius               (thisNode) &
         &                              /darkMatterProfileConcentration_%concentration(thisNode) &
         &                             )
    powerSpectrumOneHaloIntegrand=                                                                                                       &
         & +Halo_Mass_Function_Differential(time    ,massHalo                               )                                            &
         & *(                                                                                                                            &
         &   +Dark_Matter_Profile_kSpace   (thisNode,waveNumber(iWavenumber)/expansionFactor)                                            &
         &   *max(                                                                                                                       &
         &        +0.0d0                                                                                                               , &
         &        +conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMinimum,haloModelGalaxyTypeCentral  )  &
         &        -conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMaximum,haloModelGalaxyTypeCentral  )  &
         &       )                                                                                                                       &
         &   *max(                                                                                                                       &
         &        +0.0d0                                                                                                               , &
         &        +conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMinimum,haloModelGalaxyTypeSatellite)  &
         &        -conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMaximum,haloModelGalaxyTypeSatellite)  &
         &       )                                                                                                                       &
         &   +Dark_Matter_Profile_kSpace   (thisNode,waveNumber(iWavenumber)/expansionFactor)**2                                         &
         &   *max(                                                                                                                       &
         &        +0.0d0                                                                                                               , &
         &        +conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMinimum,haloModelGalaxyTypeSatellite)  &
         &        -conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMaximum,haloModelGalaxyTypeSatellite)  &
         &       )**2                                                                                                                    &
         &  )
    return
  end function powerSpectrumOneHaloIntegrand

  function powerSpectrumTwoHaloTimeIntegrand(timePrime,parameterPointer) bind(c)
    !% Time integrand for the two-halo term in the power spectrum.
    implicit none
    real   (c_double                  )        :: powerSpectrumTwoHaloTimeIntegrand
    real   (c_double                  ), value :: timePrime
    type   (c_ptr                     ), value :: parameterPointer
    type   (fgsl_function             )        :: integrandFunctionTime
    type   (fgsl_integration_workspace)        :: integrationWorkspaceTime
    logical                                    :: integrationResetTime
    
    time           =timePrime
    expansionFactor=cosmologyFunctions_%expansionFactor(time)
    call thisBasic%timeSet                             (time)
    call thisBasic%timeLastIsolatedSet                 (time)
    integrationResetTime             =.true.
    powerSpectrumTwoHaloTimeIntegrand=                                    &
         & +Integrate(                                                    &
         &            projectedCorrelationFunctionHaloMassMinimum       , &
         &            projectedCorrelationFunctionHaloMassMaximum       , &
         &            powerSpectrumTwoHaloIntegrand                     , &
         &            parameterPointer                                  , &
         &            integrandFunctionTime                             , &
         &            integrationWorkspaceTime                          , &
         &            toleranceRelative                          =1.0d-3, &
         &            reset=integrationResetTime                          &
         &           )                                                    &
         & *cosmologyFunctions_%comovingVolumeElementTime(time)
    call Integrate_Done(integrandFunctionTime,integrationWorkspaceTime)
    return
  end function powerSpectrumTwoHaloTimeIntegrand
  
  function powerSpectrumTwoHaloIntegrand(massHalo,parameterPointer) bind(c)
    !% Integrand for the two-halo term in the power spectrum.
    use Dark_Matter_Halo_Biases
    use Dark_Matter_Profiles
    use Halo_Mass_Function
    use Dark_Matter_Halo_Scales
    use Galacticus_Calculations_Resets
    use Linear_Growth
    implicit none
    real(c_double)        :: powerSpectrumTwoHaloIntegrand
    real(c_double), value :: massHalo
    type(c_ptr   ), value :: parameterPointer

    call Galacticus_Calculations_Reset(thisNode)
    call thisBasic            % massSet(                                                         &
         &                              massHalo                                                 &
         &                             )
    call thisDarkMatterProfile%scaleSet(                                                         &
         &                               Dark_Matter_Halo_Virial_Radius               (thisNode) &
         &                              /darkMatterProfileConcentration_%concentration(thisNode) &
         &                             )
    powerSpectrumTwoHaloIntegrand=                                                                        &
         & +Halo_Mass_Function_Differential(time    ,massHalo                               )             &
         & *Dark_Matter_Halo_Bias          (thisNode                                        )             &
         & *Dark_Matter_Profile_kSpace     (thisNode,waveNumber(iWavenumber)/expansionFactor)             &
         & *Linear_Growth_Factor           (time                                            )             &
         & *max(                                                                                          &
         &      +0.0d0                                                                                  , &
         &      +conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMinimum)  &
         &      -conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMaximum)  &
         &     )
    return
  end function powerSpectrumTwoHaloIntegrand

  function normalizationTimeIntegrand(timePrime,parameterPointer) bind(c)
    !% Time integrand for the normalization term in the power spectrum.
    implicit none
    real   (c_double                  )        :: normalizationTimeIntegrand
    real   (c_double                  ), value :: timePrime
    type   (c_ptr                     ), value :: parameterPointer
    type   (fgsl_function             )        :: integrandFunctionTime
    type   (fgsl_integration_workspace)        :: integrationWorkspaceTime
    logical                                    :: integrationResetTime

    time           =timePrime
    expansionFactor=cosmologyFunctions_%expansionFactor(time)
    call thisBasic%timeSet                             (time)
    call thisBasic%timeLastIsolatedSet                 (time)
    integrationResetTime=.true.
    normalizationTimeIntegrand=                                           &
         & +Integrate(                                                    &
         &            projectedCorrelationFunctionHaloMassMinimum       , &
         &            projectedCorrelationFunctionHaloMassMaximum       , &
         &            normalizationIntegrand                            , &
         &            parameterPointer                                  , &
         &            integrandFunctionTime                             , &
         &            integrationWorkspaceTime                          , &
         &            toleranceRelative                          =1.0d-3, &
         &            reset=integrationResetTime                          &
         &           )                                                    &
         & *cosmologyFunctions_%comovingVolumeElementTime(time)
    call Integrate_Done(integrandFunctionTime,integrationWorkspaceTime)
    return
  end function normalizationTimeIntegrand
  
  function normalizationIntegrand(massHalo,parameterPointer) bind(c)
    !% Integrand for the normalization term in the power spectrum.
    use Halo_Mass_Function
    implicit none
    real(c_double)        :: normalizationIntegrand
    real(c_double), value :: massHalo
    type(c_ptr   ), value :: parameterPointer

    normalizationIntegrand =                                                                              &
         & +Halo_Mass_Function_Differential(time,massHalo)                                                &
         & *max(                                                                                          &
         &      +0.0d0                                                                                  , &
         &      +conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMinimum)  &
         &      -conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMaximum)  &
         &     )
    return
  end function normalizationIntegrand

  function volumeTimeIntegrand(timePrime,parameterPointer) bind(c)
    !% Volume integrand for the normalization term in the power spectrum.
    implicit none
    real   (c_double                  )        :: volumeTimeIntegrand
    real   (c_double                  ), value :: timePrime
    type   (c_ptr                     ), value :: parameterPointer
 
    time               =timePrime
    volumeTimeIntegrand=cosmologyFunctions_%comovingVolumeElementTime(time)
    return
  end function volumeTimeIntegrand
  
end program Projected_Correlation_Function
