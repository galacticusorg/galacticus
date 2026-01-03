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

!!{
Contains a module which implements calculations of projected correlation functions using the halo model.
!!}

module Halo_Model_Projected_Correlations
  !!{
  Implements calculations of projected correlation functions using the halo model.
  !!}
  private
  public :: Halo_Model_Projected_Correlation

  abstract interface
     double precision function integrandWeight(x)
       double precision, intent(in   ) :: x
     end function integrandWeight
  end interface

contains

  subroutine Halo_Model_Projected_Correlation(                                              &
       &                                      conditionalMassFunction_                    , &
       &                                      powerSpectrum_                              , &
       &                                      cosmologyFunctions_                         , &
       &                                      surveyGeometry_                             , &
       &                                      darkMatterHaloScale_                        , &
       &                                      haloMassFunction_                           , &
       &                                      darkMatterProfileDMO_                       , &
       &                                      darkMatterHaloBias_                         , &
       &                                      darkMatterProfileScaleRadius_               , &
       &                                      projectedSeparationBinned                   , &
       &                                      projectedCorrelationFunctionMassMinimum     , &
       &                                      projectedCorrelationFunctionMassMaximum     , &
       &                                      projectedCorrelationFunctionHaloMassMinimum , &
       &                                      projectedCorrelationFunctionHaloMassMaximum , &
       &                                      projectedCorrelationFunctionLineOfSightDepth, &
       &                                      projectedCorrelationFunctionHalfIntegral    , &
       &                                      projectedCorrelationBinned                    &
       &                                     )
    !!{
    Compute the projected correlation function of galaxies above a specified mass using the halo model.
    !!}
    use :: Conditional_Mass_Functions, only : conditionalMassFunctionClass
    use :: Cosmology_Functions       , only : cosmologyFunctionsClass
    use :: Dark_Matter_Halo_Biases   , only : darkMatterHaloBiasClass
    use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleClass
    use :: Dark_Matter_Profile_Scales, only : darkMatterProfileScaleRadiusClass
    use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOClass
    use :: FFTLogs                   , only : FFTLogSineTransform              , fftLogForward
    use :: Error                     , only : Error_Report
    use :: Galacticus_Nodes          , only : nodeComponentBasic               , nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScale, treeNode
    use :: Geometry_Surveys          , only : surveyGeometryClass
    use :: Halo_Mass_Functions       , only : haloMassFunctionClass
    use :: Linear_Growth             , only : linearGrowthClass
    use :: Numerical_Constants_Math  , only : Pi
    use :: Numerical_Integration     , only : GSL_Integ_Gauss61                , integrator
    use :: Numerical_Ranges          , only : Make_Range                       , rangeTypeLogarithmic
    use :: Power_Spectra             , only : powerSpectrumClass
    use :: Table_Labels              , only : extrapolationTypeExtrapolate
    use :: Tables                    , only : table1DLogarithmicLinear
    implicit none
    class           (conditionalMassFunctionClass       ), intent(inout)                                             :: conditionalMassFunction_
    class           (powerSpectrumClass                 ), intent(inout)                                             :: powerSpectrum_
    class           (cosmologyFunctionsClass            ), intent(inout)                                             :: cosmologyFunctions_
    class           (surveyGeometryClass                ), intent(inout)                                             :: surveyGeometry_
    class           (darkMatterHaloScaleClass           ), intent(inout)                                             :: darkMatterHaloScale_
    class           (haloMassFunctionClass              ), intent(inout)                                             :: haloMassFunction_
    class           (darkMatterProfileDMOClass          ), intent(inout)                                             :: darkMatterProfileDMO_
    class           (darkMatterHaloBiasClass            ), intent(inout)                                             :: darkMatterHaloBias_
    class           (darkMatterProfileScaleRadiusClass  ), intent(inout)                                             :: darkMatterProfileScaleRadius_
    double precision                                     , intent(in   ), dimension(                             : ) :: projectedSeparationBinned
    double precision                                     , intent(in   )                                             :: projectedCorrelationFunctionMassMinimum     , projectedCorrelationFunctionMassMaximum    , &
         &                                                                                                              projectedCorrelationFunctionHaloMassMinimum , projectedCorrelationFunctionHaloMassMaximum, &
         &                                                                                                              projectedCorrelationFunctionLineOfSightDepth
    logical                                              , intent(in   )                                             :: projectedCorrelationFunctionHalfIntegral
    double precision                                     , intent(  out), dimension(size(projectedSeparationBinned)) :: projectedCorrelationBinned
    type            (treeNode                           ), pointer                                                   :: node
    class           (nodeComponentBasic                 ), pointer                                                   :: basic
    class           (nodeComponentDarkMatterProfile     ), pointer                                                   :: darkMatterProfileHalo
    procedure       (integrandWeight                    ), pointer                                                   :: integrandWeightFunction
    double precision                                     , allocatable  , dimension(                             : ) :: powerSpectrumOneHalo                              , powerSpectrumTwoHalo                              , &
         &                                                                                                              wavenumber                                        , powerSpectrumTotal                                , &
         &                                                                                                              separation                                        , correlation                                       , &
         &                                                                                                              projectedCorrelation
    integer                                              , parameter                                                 :: wavenumberCountPerDecade                   =10
    double precision                                     , parameter                                                 :: wavenumberMinimum                          =1.0d-3, wavenumberMaximum                          =1.0d4, &
         &                                                                                                              virialWavenumberMultiplierMaximum          =1.0d+3
    logical                                                                                                          :: integrationWarningIssued=.false.
    double precision                                                                                                 :: time                                              , timeMinimum                                       , &
         &                                                                                                              expansionFactor                                   , volume                                            , &
         &                                                                                                              timeMaximum                                       , galaxyDensity                                     , &
         &                                                                                                              projectedSeparation                               , binSeparationMaximum                              , &
         &                                                                                                              binWidthLogarithmic                               , binSeparationMinimum
    integer                                                                                                          :: iField                                            , iWavenumber                                       , &
         &                                                                                                              wavenumberCount                                   , iSeparation                                       , &
         &                                                                                                              projectedCorrelationFunctionSeparationCount
    type            (integrator                         )                                                            :: integratorVolumeTime                              , integratorNormalizationTime                       , &
         &                                                                                                              integratorPowerSpectrumOneHaloTime                , integratorPowerSpectrumTwoHaloTime
    type            (table1DLogarithmicLinear           )                                                            :: correlationTable

    ! Create worker node.
    node                  => treeNode              (                 )
    basic                 => node%basic            (autoCreate=.true.)
    darkMatterProfileHalo => node%darkMatterProfile(autoCreate=.true.)
    select type (darkMatterProfileHalo)
    type is (nodeComponentDarkMatterProfileScale)
       ! This is acceptable.
    class default
       ! This is not.
       call Error_Report('this code expects to use the "scale" dark matter profile component'//{introspection:location})
    end select
    ! Generate wavenumber range.
    wavenumberCount=int(log10(wavenumberMaximum/wavenumberMinimum)*dble(wavenumberCountPerDecade))+1
    wavenumber     =Make_Range(wavenumberMinimum,wavenumberMaximum,wavenumberCount,rangeTypeLogarithmic)
    allocate(powerSpectrumTotal  (wavenumberCount))
    allocate(powerSpectrumOneHalo(wavenumberCount))
    allocate(powerSpectrumTwoHalo(wavenumberCount))
    ! Initialize.
    volume              =0.0d0
    galaxyDensity       =0.0d0
    powerSpectrumOneHalo=0.0d0
    powerSpectrumTwoHalo=0.0d0
    ! Iterate over survey fields.
    do iField=1,surveyGeometry_%fieldCount()
       ! Find time range for volume-limited sample.
       timeMinimum=cosmologyFunctions_%timeAtDistanceComoving(surveyGeometry_%distanceMaximum(projectedCorrelationFunctionMassMinimum,field=iField))
       timeMaximum=cosmologyFunctions_%timeAtDistanceComoving(surveyGeometry_%distanceMinimum(projectedCorrelationFunctionMassMinimum,field=iField))
       ! Integrate the volume term.
       integratorVolumeTime= integrator                      (volumeTimeIntegrand,toleranceRelative=1.0d-2     )
       volume              =+volume                                                                              &
            &               +integratorVolumeTime%integrate  (timeMinimum        ,                  timeMaximum) &
            &               *surveyGeometry_      %solidAngle(iField                                           )
       ! Integrate the normalization term.
       integratorNormalizationTime=integrator                             (normalizationTimeIntegrand,toleranceRelative=1.0d-2     )
       galaxyDensity              =+galaxyDensity                                                                                    &
            &                      +integratorNormalizationTime%integrate (timeMinimum               ,                  timeMaximum) &
            &                      *surveyGeometry_            %solidAngle(iField                                                  )
       ! Iterate over wavenumbers.
       integratorPowerSpectrumOneHaloTime=integrator(powerSpectrumOneHaloTimeIntegrand,toleranceRelative=1.0d-2,integrationRule=GSL_Integ_Gauss61)
       integratorPowerSpectrumTwoHaloTime=integrator(powerSpectrumTwoHaloTimeIntegrand,toleranceRelative=1.0d-2,integrationRule=GSL_Integ_Gauss61)
       do iWavenumber=1,wavenumberCount
          ! Integrate the one-halo term.
          powerSpectrumOneHalo(iWavenumber)=+powerSpectrumOneHalo                         (iWavenumber            ) &
               &                            +integratorPowerSpectrumOneHaloTime%integrate (timeMinimum,timeMaximum) &
               &                            *surveyGeometry_                   %solidAngle(iField                 )
          ! Integrate the two-halo term.
          powerSpectrumTwoHalo(iWavenumber)=+powerSpectrumTwoHalo(iWavenumber)                                      &
               &                            +integratorPowerSpectrumTwoHaloTime%integrate (timeMinimum,timeMaximum) &
               &                            *surveyGeometry_                   %solidAngle(iField                 )
       end do
    end do
    ! Check for non-zero galaxy density.
    if (galaxyDensity > 0.0d0) then
       ! Construct the net power spectrum.
       galaxyDensity       =galaxyDensity       /volume
       powerSpectrumOneHalo=powerSpectrumOneHalo/volume
       powerSpectrumTwoHalo=powerSpectrumTwoHalo/volume
       do iWavenumber=1,wavenumberCount
          powerSpectrumTotal(iWavenumber)=(                                       &
               &                           +powerSpectrumOneHalo(iWavenumber )    &
               &                           +powerSpectrumTwoHalo(iWavenumber )**2 &
               &                          )                                       &
               &                          /galaxyDensity                      **2
       end do
       ! Fourier transform to get the correlation function.
       allocate(correlation,mold=wavenumber)
       allocate(separation,mold=wavenumber)
       call FFTLogSineTransform(                     &
            &                   wavenumber         , &
            &                   separation         , &
            &                   +powerSpectrumTotal  &
            &                   *wavenumber          &
            &                   * 4.0d0*Pi           &
            &                   /(2.0d0*Pi)**3     , &
            &                   correlation        , &
            &                   fftLogForward        &
            &                  )
       correlation=correlation/separation
       ! Project the correlation function.
       allocate(projectedCorrelation,mold=wavenumber)
       call correlationTable%create(separation(1),separation(wavenumberCount),wavenumberCount,extrapolationType=[extrapolationTypeExtrapolate,extrapolationTypeExtrapolate])
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
       ! If the integral was taken over the half range, 0<pi<pi_max, rather than the full range, -pi_max<pi<pi_max, then divide
       ! the projected correlation function by two.
       if (projectedCorrelationFunctionHalfIntegral) projectedCorrelation=projectedCorrelation/2.0d0
       ! Average the projected correlation function into bins.
       projectedCorrelationFunctionSeparationCount=size(projectedSeparationBinned)
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
               &                                        )                                                           &
               &                                      *projectedCorrelation                                         &
               &                                     )
       end do
    else
       ! Galaxy density is zero. Return a zero correlation function.
       projectedCorrelationBinned=0.0d0
    end if

  contains

    double precision function projectionIntegrandWeight(separation)
      !!{
      The weight function applied to the correlation function when integrating to get the projected correlation function.
      !!}
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
      !!{
      The weight function applied to the projected correlation function when integrating into bins.
      !!}
      implicit none
      double precision, intent(in   ) :: separation

      binningIntegrandWeight=2.0d0*Pi*separation
      return
    end function binningIntegrandWeight

    double precision function powerSpectrumOneHaloTimeIntegrand(timePrime)
      !!{
      Time integrand for the one-halo term in the power spectrum.
      !!}
      use :: Display         , only : displayMessage    , verbosityLevelWarn, displayMagenta, displayReset

      use :: Error, only : errorStatusSuccess
      implicit none
      double precision            , intent(in   ) :: timePrime
      type            (integrator)                :: integratorTime
      integer                                     :: errorStatus

      time           =timePrime
      expansionFactor=cosmologyFunctions_%expansionFactor(time)
      call basic%timeSet                             (time)
      call basic%timeLastIsolatedSet                 (time)
      integratorTime                   = integrator(powerSpectrumOneHaloIntegrand,toleranceRelative=1.0d-2,integrationRule=GSL_Integ_Gauss61)
      powerSpectrumOneHaloTimeIntegrand=+integratorTime     %integrate                (                                                    &
           &                                                                                  projectedCorrelationFunctionHaloMassMinimum, &
           &                                                                                  projectedCorrelationFunctionHaloMassMaximum, &
           &                                                                           status=errorStatus                                  &
           &                                                                          )                                                    &
           &                            *cosmologyFunctions_%comovingVolumeElementTime(time)
      if (errorStatus /= errorStatusSuccess .and. .not.integrationWarningIssued) then
         call displayMessage(displayMagenta()//'WARNING:'//displayReset()//' [powerSpectrumOneHaloTimeIntegrand] integration failed - likely due to oscillatory nature of integrand - proceeding anyway',verbosity=verbosityLevelWarn)
         integrationWarningIssued=.true.
      end if
      return
    end function powerSpectrumOneHaloTimeIntegrand

    double precision function powerSpectrumOneHaloIntegrand(massHalo)
      !!{
      Integrand for the one-halo term in the power spectrum.
      !!}
      use :: Conditional_Mass_Functions, only : haloModelGalaxyTypeCentral, haloModelGalaxyTypeSatellite
      use :: Calculations_Resets       , only : Calculations_Reset
      use :: Mass_Distributions        , only : massDistributionClass
     implicit none
      double precision                       , intent(in   ) :: massHalo
      class           (massDistributionClass), pointer       :: massDistribution_
      double precision                                       :: darkMatterProfileKSpace, numberCentrals   , &
           &                                                    numberSatellites       , wavenumberMaximum, &
           &                                                    radiusVirial

      call Calculations_Reset(node)
      call basic                % massSet(massHalo                                  )
      call Calculations_Reset(node)
      call darkMatterProfileHalo%scaleSet(darkMatterProfileScaleRadius_%radius(node))
      ! Return zero if we're more than some maximum factor above the virial wavenumber for this halo. This avoids attempting to
      ! integrate rapidly oscillating Fourier profiles.
      wavenumberMaximum=virialWavenumberMultiplierMaximum/(darkMatterHaloScale_%radiusVirial(node)/expansionFactor)
      if (waveNumber(iWavenumber) > wavenumberMaximum) then
         powerSpectrumOneHaloIntegrand=0.0d0
      else
         massDistribution_       => darkMatterProfileDMO_%get             (node                                                )
         radiusVirial            =  darkMatterHaloScale_ %radiusVirial    (node                                                )
         darkMatterProfileKSpace =  massDistribution_    %fourierTransform(radiusVirial,waveNumber(iWavenumber)/expansionFactor)
         numberCentrals          =  max(                                                                                                                       &
              &                         +0.0d0                                                                                                               , &
              &                         +conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMinimum,haloModelGalaxyTypeCentral  )  &
              &                         -conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMaximum,haloModelGalaxyTypeCentral  )  &
              &                        )
         numberSatellites        =  max(                                                                                                                       &
              &                         +0.0d0                                                                                                               , &
              &                         +conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMinimum,haloModelGalaxyTypeSatellite)  &
              &                         -conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMaximum,haloModelGalaxyTypeSatellite)  &
              &                        )
         !![
	 <objectDestructor name="massDistribution_"/>
         !!]
         ! Note that we include 2 times the central-satellite term since we want to count each pair twice (i.e. central-satellite and
         ! then satellite-central). This is consistent with the N(N-1) counting for the satellite-satellite term, and with the
         ! counting in the two-halo term.
         powerSpectrumOneHaloIntegrand=                        &
              & +haloMassFunction_%differential(time,massHalo) &
              & *(                                             &
              &   +darkMatterProfileKSpace                     &
              &   *2.0d0                                       &
              &   *numberCentrals                              &
              &   *numberSatellites                            &
              &   +darkMatterProfileKSpace**2                  &
              &   *numberSatellites       **2                  &
              &  )
      end if
      return
    end function powerSpectrumOneHaloIntegrand

    double precision function powerSpectrumTwoHaloTimeIntegrand(timePrime)
      !!{
      Time integrand for the two-halo term in the power spectrum.
      !!}
      use :: Display, only : displayMessage    , verbosityLevelWarn, displayMagenta, displayReset
      use :: Error  , only : errorStatusSuccess
      implicit none
      double precision            , intent(in   ) :: timePrime
      type            (integrator)                :: integratorTime
      integer                                     :: errorStatus
      
      time           =timePrime
      expansionFactor=cosmologyFunctions_%expansionFactor(time)
      call basic%timeSet                             (time)
      call basic%timeLastIsolatedSet                 (time)
      integratorTime                  =       integrator                                   (                                                               &
           &                                                                                                  powerSpectrumTwoHaloIntegrand              , &
           &                                                                                toleranceRelative=1.0d-2                                     , &
           &                                                                                integrationRule  =GSL_Integ_Gauss61                            &
           &                                                                               )
      powerSpectrumTwoHaloTimeIntegrand=+     integratorTime     %integrate                (                                                               &
           &                                                                                                  projectedCorrelationFunctionHaloMassMinimum, &
           &                                                                                                  projectedCorrelationFunctionHaloMassMaximum, &
           &                                                                                status           =errorStatus                                  &
           &                                                                               )                                                               &
           &                            *sqrt(powerSpectrum_     %power                    (wavenumber(iWavenumber),time))                                 &
           &                            *     cosmologyFunctions_%comovingVolumeElementTime(                        time)
      if (errorStatus /= errorStatusSuccess .and. .not.integrationWarningIssued) then
         call displayMessage(displayMagenta()//'WARNING:'//displayReset()//' [powerSpectrumTwoHaloTimeIntegrand] integration failed - likely due to oscillatory nature of integrand - proceeding anyway',verbosity=verbosityLevelWarn)
         integrationWarningIssued=.true.
      end if
      return
    end function powerSpectrumTwoHaloTimeIntegrand

    double precision function powerSpectrumTwoHaloIntegrand(massHalo)
      !!{
      Integrand for the two-halo term in the power spectrum.
      !!}
      use :: Calculations_Resets, only : Calculations_Reset
      use :: Mass_Distributions , only : massDistributionClass
     implicit none
      double precision                       , intent(in   ) :: massHalo
      class           (massDistributionClass), pointer       :: massDistribution_
      double precision                                       :: wavenumberMaximum, radiusVirial

      call Calculations_Reset(node)
      call basic                % massSet(massHalo                                  )
      call Calculations_Reset(node)
      call darkMatterProfileHalo%scaleSet(darkMatterProfileScaleRadius_%radius(node))
      ! Return zero if we're more than some maximum factor above the virial wavenumber for this halo. This avoids attempting to
      ! integrate rapidly oscillating Fourier profiles.
      wavenumberMaximum=virialWavenumberMultiplierMaximum/(darkMatterHaloScale_%radiusVirial(node)/expansionFactor)
      if (waveNumber(iWavenumber) > wavenumberMaximum) then
         powerSpectrumTwoHaloIntegrand=0.0d0
      else
         massDistribution_             =>  darkMatterProfileDMO_%get             (node                                   )
         radiusVirial                  =   darkMatterHaloScale_ %radiusVirial    (node                                   )
         powerSpectrumTwoHaloIntegrand =  +haloMassFunction_    %differential    (time,massHalo                          )               &
              &                           *darkMatterHaloBias_  %bias            (node                                   )               &
              &                           *massDistribution_    %fourierTransform(radiusVirial,waveNumber(iWavenumber)/expansionFactor)               &
              &                           *max(                                                                                          &
              &                                +0.0d0                                                                                  , &
              &                                +conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMinimum)  &
              &                                -conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMaximum)  &
              &                               )
         !![
	 <objectDestructor name="massDistribution_"/>
         !!]
      end if
      return
    end function powerSpectrumTwoHaloIntegrand

    double precision function normalizationTimeIntegrand(timePrime)
      !!{
      Time integrand for the normalization term in the power spectrum.
      !!}
      implicit none
      double precision            , intent(in   ) :: timePrime
      type            (integrator)                :: integratorTime

      time           =timePrime
      integratorTime=integrator(normalizationIntegrand,toleranceRelative=1.0d-2)
      normalizationTimeIntegrand=+integratorTime     %integrate                (                                             &
           &                                                                    projectedCorrelationFunctionHaloMassMinimum, &
           &                                                                    projectedCorrelationFunctionHaloMassMaximum  &
           &                                                                   )                                             &
           &                     *cosmologyFunctions_%comovingVolumeElementTime(time)
      return
    end function normalizationTimeIntegrand

    double precision function normalizationIntegrand(massHalo)
      !!{
      Integrand for the normalization term in the power spectrum.
      !!}
      implicit none
      double precision, intent(in   ) :: massHalo

      normalizationIntegrand =                                                                              &
           & +haloMassFunction_%differential(time,massHalo)                                                 &
           & *max(                                                                                          &
           &      +0.0d0                                                                                  , &
           &      +conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMinimum)  &
           &      -conditionalMassFunction_%massFunction(massHalo,projectedCorrelationFunctionMassMaximum)  &
           &     )
      return
    end function normalizationIntegrand

    double precision function volumeTimeIntegrand(timePrime)
      !!{
      Volume integrand for the normalization term in the power spectrum.
      !!}
      implicit none
      double precision, intent(in   ) :: timePrime

      time               =timePrime
      volumeTimeIntegrand=cosmologyFunctions_%comovingVolumeElementTime(time)
      return
    end function volumeTimeIntegrand

  end subroutine Halo_Model_Projected_Correlation

end module Halo_Model_Projected_Correlations
