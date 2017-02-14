!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which computes mass function covariances.

module Statistics_Mass_Function_Covariance
  !% Implements calculations of mass function covariances.
  use Cosmology_Functions
  use Geometry_Surveys
  private
  public :: Mass_Function_Covariance_Matrix

  ! Record of whether this module is intialized.
  logical          :: moduleInitialized=.false.

  ! The cosmic time at which the calculation is to be performed.
  double precision :: time
  !$omp threadprivate(time)

  ! The wavenumber for which LSS integrations are currently being performed.
  double precision :: waveNumberGlobal
  !$omp threadprivate(waveNumberGlobal)

  ! Integration limits.
  double precision :: logMassLower,logMassUpper

  ! Minimum and maximum masses for the bins being considered.
  double precision, dimension(:), allocatable :: log10MassBinWidth,logMassBinWidth
  integer          :: binI,binJ,lssBin
  !$omp threadprivate(lssBin,binI,binJ)
  double precision :: massBinCenterI,massBinMinimumI,massBinMaximumI
  double precision :: massBinCenterJ,massBinMinimumJ,massBinMaximumJ

  ! Table of biases.
  integer          :: timeBinCount
  double precision, dimension(:  ), allocatable :: timeTable
  double precision, dimension(:,:), allocatable :: biasTable
  
  ! Cosmological functions.
  class           (cosmologyFunctionsClass), pointer                    :: cosmologyFunctions_

  ! Survey geometry.
  class           (surveyGeometryClass    ), pointer                    :: surveyGeometry_
  double precision                                                      :: surveyRedshiftMinimum, surveyRedshiftMaximum
  double precision                         , allocatable, dimension(:)  :: volumeNormalizationI , volumeNormalizationJ , &
       &                                                                   timeMinimumI         , timeMinimumJ         , &
       &                                                                   timeMaximumI         , timeMaximumJ         , &
       &                                                                   logMassBinCenter
  !$omp threadprivate(timeMinimumI,timeMinimumJ,timeMaximumI,timeMaximumJ,volumeNormalizationI,volumeNormalizationJ)

contains

  subroutine Mass_Function_Covariance_Matrix(redshiftMinimum,redshiftMaximum,massBinCount,massMinimum,massMaximum,massObserved,massWidthObserved,massFunctionObserved,completenessObserved,numberObserved,completenessErrorObserved,includePoisson,includeHalo,includeLSS&
       &,mass,massFunction,covariance ,covariancePoisson,covarianceHalo,covarianceLSS,correlation)
    !% Compute the mass function covariance matrix.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use FFTW3
    use Memory_Management
    use Input_Parameters
    use Numerical_Ranges
    use Numerical_Integration
    use Numerical_Constants_Math
    use Power_Spectra_Nonlinear
    use Galacticus_Error
    use Galacticus_Display
    implicit none
    integer                            , intent(in   )                                        :: massBinCount
    double precision                   , intent(in   )                                        :: redshiftMinimum,redshiftMaximum&
         &,massMinimum,massMaximum,completenessErrorObserved
    logical                            , intent(in   )                                        :: includePoisson,includeHalo&
         &,includeLSS
    double precision                   , intent(inout), allocatable, dimension(:    )         :: mass
    double precision                   , intent(inout), allocatable, dimension(:    ), target :: massFunction,massFunctionObserved, completenessObserved, numberObserved, massObserved, massWidthObserved
    double precision                   , intent(inout), allocatable, dimension(:,:  )         :: covariance,covariancePoisson &
         &,covarianceHalo,covarianceLSS,correlation
    double precision                   ,                allocatable, dimension(:,:  )         :: varianceLSS    , volume
    double precision                   ,                pointer    , dimension(:    )         :: massFunctionUse
    double precision                   , parameter                                            :: timePointsPerDecade=100
    double precision                   , parameter                                            :: massFunctionMinimum=1.0d-50
    logical                                                                                   :: integrationReset, useCompleteness, useNumber
    integer                                                                                   :: i,j &
         &,iTime,iField,fieldCount
    double precision                                                                          :: logMassMinimum ,logMassMaximum&
         &,normalization ,massFunctionCovarianceHaloMassMinimum&
         &,massFunctionCovarianceHaloMassMaximum,timeMinimum,timeMaximum,volumeNormalization
    double precision                                                                          :: binCompleteness
    type   (fgsl_function             )                                                       :: integrandFunction
    type   (fgsl_integration_workspace)                                                       :: integrationWorkspace
    
    ! Initialize the module if necessary.
    if (.not.moduleInitialized) then
       ! Read controlling parameters.
       !@ <inputParameter>
       !@   <name>massFunctionCovarianceHaloMassMinimum</name>
       !@   <defaultValue>$10^{10}M_\odot$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The minimum halo mass to use when computing mass function covariance matrices.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('massFunctionCovarianceHaloMassMinimum',massFunctionCovarianceHaloMassMinimum,defaultValue=1.0d10)
       !@ <inputParameter>
       !@   <name>massFunctionCovarianceHaloMassMaximum</name>
       !@   <defaultValue>$10^{10}M_\odot$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The minimum halo mass to use when computing mass function covariance matrices.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('massFunctionCovarianceHaloMassMaximum',massFunctionCovarianceHaloMassMaximum,defaultValue=1.0d15)

       ! Record that the module is initialized.
       moduleInitialized=.true.
    end if
    
    ! Get the default cosmology functions object.
    cosmologyFunctions_   => cosmologyFunctions        ()

    ! Get the default survey geometry.
    surveyGeometry_       => surveyGeometry            ()
    fieldCount            =  surveyGeometry_%fieldCount()
    surveyRedshiftMinimum = redshiftMinimum
    surveyRedshiftMaximum = redshiftMaximum

    ! Determine number of times over which to tabulate bias.
    timeMaximum =cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum))
    timeMinimum =cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum))
    timeBinCount=int(log10(timeMaximum/timeMinimum)*dble(timePointsPerDecade))+1

    ! Allocate arrays.
    call allocateArray(mass             ,[massBinCount             ])
    call allocateArray(logMassBinCenter ,[massBinCount             ])
    call allocateArray(log10MassBinWidth,[massBinCount             ])
    call allocateArray(logMassBinWidth  ,[massBinCount             ])
    call allocateArray(massFunction     ,[massBinCount             ])
    call allocateArray(volume           ,[massBinCount,fieldCount  ])
    call allocateArray(covariance       ,[massBinCount,massBinCount])
    call allocateArray(covariancePoisson,[massBinCount,massBinCount])
    call allocateArray(covarianceHalo   ,[massBinCount,massBinCount])
    call allocateArray(covarianceLSS    ,[massBinCount,massBinCount])
    call allocateArray(correlation      ,[massBinCount,massBinCount])
    call allocateArray(varianceLSS      ,[massBinCount,massBinCount])
    call allocateArray(timeTable        ,[timeBinCount             ])
    call allocateArray(biasTable        ,[timeBinCount,massBinCount])

    ! Create time bins.
    timeTable=Make_Range(timeMinimum,timeMaximum,timeBinCount,rangeType=rangeTypeLogarithmic)

    ! Create mass bins.
    if (allocated(massWidthObserved)) then
       mass             =massObserved
       log10MassBinWidth=log10(massWidthObserved)
       logMassBinWidth  =log(10.0d0)*log10MassBinWidth
       logMassBinCenter =log10(mass)
    else
       logMassMinimum   =log10(massMinimum)
       logMassMaximum   =log10(massMaximum)
       logMassBinCenter =Make_Range(logMassMinimum,logMassMaximum,massBinCount,rangeType=rangeTypeLinear)
       log10MassBinWidth=logMassBinCenter(2)-logMassBinCenter(1)
       logMassBinWidth  =log(10.0d0)*log10MassBinWidth
       mass             =10.0d0**logMassBinCenter
    end if

    ! Halo mass limits for integrations.
    logMassLower    =log10(massFunctionCovarianceHaloMassMinimum)
    logMassUpper    =log10(massFunctionCovarianceHaloMassMaximum)

    ! Determine which mass function to use.
    if (allocated(massFunctionObserved)) then
       if (size(massFunctionObserved) /= massBinCount) call Galacticus_Error_Report('Mass_Function_Covariance_Matrix','observed mass function has incorrect number of bins')
       massFunctionUse => massFunctionObserved
    else
       massFunctionUse => massFunction
    end if

    ! Determine if completeness and/or number is available.
    useCompleteness=allocated(completenessObserved)
    if (useCompleteness .and. size(completenessObserved) /= massBinCount) &
         & call Galacticus_Error_Report('Mass_Function_Covariance_Matrix','observed completeness has incorrect number of bins')
    useNumber      =allocated(      numberObserved)
    if (useNumber       .and. size(      numberObserved) /= massBinCount) &
         & call Galacticus_Error_Report('Mass_Function_Covariance_Matrix','observed number has incorrect number of bins'      )

    ! Compute the mass function and bias averaged over each bin.
    massFunction=0.0d0
    do i=1,massBinCount

       ! Find limits on mass for this bin.
       massBinCenterI =10.0** logMassBinCenter(i)
       massBinMinimumI=10.0**(logMassBinCenter(i)-0.5d0*log10MassBinWidth(i))
       massBinMaximumI=10.0**(logMassBinCenter(i)+0.5d0*log10MassBinWidth(i))

       ! Iterate over fields.
       volumeNormalization=0.0d0
       do iField=1,fieldCount

          ! Find integration limits for this bin.
          timeMaximum=    cosmologyFunctions_%cosmicTime            (cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum       ))
          timeMinimum=max(                                                                                                                                 &
               &          cosmologyFunctions_%cosmicTime            (cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum       )), &
               &          cosmologyFunctions_%timeAtDistanceComoving(surveyGeometry_%distanceMaximum                      (massBinCenterI ,iField))  &
               &         )
          ! Get the normalizing volume integral.
          integrationReset=.true.
          volumeNormalization=                                     &
               &              +volumeNormalization                 &
               &              +Integrate(                          &
               &                         timeMinimum             , &
               &                         timeMaximum             , &
               &                         Volume_Integrand        , &
               &                         integrandFunction       , &
               &                         integrationWorkspace    , &
               &                         toleranceRelative=1.0d-3, &
               &                         reset=integrationReset    &
               &                        ) 

          ! Integrate mass function over the bin.
          integrationReset=.true.
          massFunction(i)=                                            &
               &          +massFunction(i)                            &
               &          +Integrate(                                 &
               &                     timeMinimum                    , &
               &                     timeMaximum                    , &
               &                     Mass_Function_Time_Integrand_I , &
               &                     integrandFunction              , &
               &                     integrationWorkspace           , &
               &                     toleranceRelative=1.0d-3       , &
               &                     reset=integrationReset           &
               &                    )
          call Integrate_Done(integrandFunction,integrationWorkspace)
          
          ! Find the effective volume of the survey at this mass.
          volume(i,iField)=surveyGeometry_%volumeMaximum(massBinCenterI,iField)
       end do

       ! Normalize the mass function.
       massFunction(i)=massFunction(i)/logMassBinWidth(i)/volumeNormalization
       ! Tabulate the bias as a function of time in this bin.
       integrationReset=.true.
       do iTime=1,timeBinCount
          time=timeTable(iTime)
          biasTable(iTime,i)=  Integrate(                  &
               &                 logMassLower            , &
               &                 logMassUpper            , &
               &                 Bias_Integrand_I        , &
               &                 integrandFunction       , &
               &                 integrationWorkspace    , &
               &                 toleranceRelative=1.0d-2, &
               &                 reset=integrationReset    &
               &                )                          &
               & /logMassBinWidth(i)
       end do
       call Integrate_Done(integrandFunction,integrationWorkspace)
    end do

    ! Compute LSS variance if necessary.
    if (includeLSS) then
       ! Compute large-scale structure variance for each cell pair.
       ! If angular power spectrum of survey window function is available, use it to compute LSS contribution to variance.
       if (surveyGeometry_%angularPowerAvailable()) then
          call Variance_LSS_Angular_Spectrum(massBinCount,redshiftMinimum,redshiftMaximum,varianceLSS)
       ! If survey function function is available, use it to compute LSS contribution to variance.
       else if (surveyGeometry_%windowFunctionAvailable()) then
          call Variance_LSS_Window_Function(massBinCount,redshiftMinimum,redshiftMaximum,varianceLSS)
       ! No method exists to compute the LSS contribution to variance. Abort.
       else
          call Galacticus_Error_Report('Mass_Function_Covariance_Matrix','no method exists to compute LSS contribution to covariance matrix')
       end if       
    end if

    ! Construct the covariance matrix.
    covariancePoisson=0.0d0
    covarianceHalo   =0.0d0
    covarianceLSS    =0.0d0
    do i   =1,massBinCount
       massBinCenterI    =10.0d0** logMassBinCenter(i)
       massBinMinimumI   =10.0d0**(logMassBinCenter(i)-0.5d0*log10MassBinWidth(i))
       massBinMaximumI   =10.0d0**(logMassBinCenter(i)+0.5d0*log10MassBinWidth(i))
       do j=i,massBinCount
          massBinCenterJ =10.0d0** logMassBinCenter(j)
          massBinMinimumJ=10.0d0**(logMassBinCenter(j)-0.5d0*log10MassBinWidth(j))
          massBinMaximumJ=10.0d0**(logMassBinCenter(j)+0.5d0*log10MassBinWidth(j))
          ! Poisson term.
          if (includePoisson .and. i == j) then
             if      (useCompleteness) then
                binCompleteness=completenessObserved(i)
             else if (useNumber      ) then
                if (numberObserved(i) > 0.0d0) then
                   binCompleteness=     numberObserved     (i  )  &
                        &          /    massFunctionUse    (i  )  &
                        &          /sum(volume             (i,:)) &
                        &          /    logMassBinWidth(i)
                else
                   binCompleteness=1.0d0
                end if
             else
                binCompleteness=1.0d0
             end if
             if (massFunctionUse(i) > 0.0d0) then
                covariancePoisson(i,j)=massFunctionUse(i)/binCompleteness/(sum(volume(i,:))*logMassBinWidth(i))
             else
                covariancePoisson(i,j)=1.0d0                             /(sum(volume(i,:))*logMassBinWidth(i))**2
             end if
          end if
          
          ! Halo occupancy covariance.
          if (includeHalo) then
             ! Iterate over fields.
             do iField=1,fieldCount
                integrationReset=.true.
                ! Find integration limits for this bin.
                timeMaximum=    cosmologyFunctions_%cosmicTime            (cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum))
                timeMinimum=max(                                                                                                              &
                     &          cosmologyFunctions_%cosmicTime            (cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum)), &
                     &          cosmologyFunctions_%timeAtDistanceComoving(                                                                   &
                     &                                                     min(                                                               &
                     &                                                         surveyGeometry_%distanceMaximum(massBinCenterI,iField),        &
                     &                                                         surveyGeometry_%distanceMaximum(massBinCenterJ,iField)         &
                     &                                                        )                                                               &
                     &                                                    )                                                                   &
                     &         )
                if (timeMaximum > timeMinimum) then
                   ! Integrate over the volume. Note that the following expression is multiplied through by the volumes of both
                   ! fields such that we accumulate a volume-weighted covariance, which will be normalized below.
                   covarianceHalo(i,j)=                                          &
                        &             + covarianceHalo(i,j)                      &
                        &             + Integrate(                               &
                        &                         timeMinimum                  , &
                        &                         timeMaximum                  , &
                        &                         Halo_Occupancy_Time_Integrand, &
                        &                         integrandFunction            , &
                        &                         integrationWorkspace         , &
                        &                         toleranceRelative=1.0d-3     , &
                        &                         reset=integrationReset         &
                        &                        )                               &
                        &              *surveyGeometry_%solidAngle(iField)       &
                        &              /logMassBinWidth(i)                       &
                        &              /logMassBinWidth(j)
                   call Integrate_Done(integrandFunction,integrationWorkspace)
                end if
             end do
             ! Normalize the covariance for the total field volume.
             if (sum(volume(i,:)) > 0.0d0 .and. sum(volume(j,:)) > 0.0d0) covarianceHalo(i,j)=covarianceHalo(i,j)/sum(volume(i,:))/sum(volume(j,:))
             ! Renormalize to actual mass function. Accounts for any difference between model and data. Including incompleteness.            
             if     (                                                                                         &
                  &   massFunctionUse(i) > massFunctionMinimum .and. massFunctionUse(j) > massFunctionMinimum &
                  &  .and.                                                                                    &
                  &   massFunction   (i) > massFunctionMinimum .and. massFunction   (j) > massFunctionMinimum &
                  & ) then
                covarianceHalo(i,j)= covarianceHalo( i,j) &
                     &              *massFunctionUse(i  ) &
                     &              /massFunction   (i  ) &
                     &              *massFunctionUse(  j) & 
                     &              /massFunction   (  j)
             end if
          end if
         
          ! Large-scale structure term.
          if (includeLSS) then
             covarianceLSS(i,j)=varianceLSS(i,j)
             if     (                                                                                         &
                  &   massFunctionUse(i) > massFunctionMinimum .and. massFunctionUse(j) > massFunctionMinimum &
                  &  .and.                                                                                    &
                  &   massFunction   (i) > massFunctionMinimum .and. massFunction   (j) > massFunctionMinimum &
                  & ) then
                ! Renormalize to actual mass function. Accounts for any difference between model and data. Including incompleteness.
                covarianceLSS(i,j)= covarianceLSS  (i,j) &
                     &             *massFunctionUse(i  ) &
                     &             /massFunction   (i  ) &
                     &             *massFunctionUse(  j) &
                     &             /massFunction   (  j)
             end if
          end if
       end do
    end do

    ! Symmetrize the covariance matrices.
    do i   =1,massBinCount
       do j=1,massBinCount
          if (j < i) then
             covariancePoisson(i,j)=covariancePoisson(j,i)
             covarianceHalo   (i,j)=covarianceHalo   (j,i)
             covarianceLSS    (i,j)=covarianceLSS    (j,i)
          end if
       end do
    end do

    ! Sum covariances.
    covariance=covariancePoisson+covarianceHalo+covarianceLSS      
    
    ! Add in any covariance arising from uncertainty in the incompleteness.
    do i   =1,massBinCount
       do j=1,massBinCount
          covariance(i,j)=covariance(i,j)+completenessErrorObserved**2*massFunctionUse(i)*massFunctionUse(j)
       end do
    end do

    ! Compute the corresponding correlation matrix.
    do i   =1,massBinCount
       do j=1,massBinCount
          normalization=sqrt(covariance(i,i)*covariance(j,j))
          if (normalization > 0.0d0) then
             correlation(i,j)=covariance(i,j)/normalization
          else
             correlation(i,j)=0.0d0
          end if
       end do
    end do

    ! Deallocate arrays.
    call deallocateArray(logMassBinCenter    )
    call deallocateArray(volume              )
    call deallocateArray(varianceLSS         )

    return
  end subroutine Mass_Function_Covariance_Matrix
  
  double precision function Galaxy_Root_Power_Spectrum(iBin,timeMinimum,timeMaximum)
    !% Computes the quantity $\int_{t_{\mathrm min}}^{t_{\mathrm max}} {\mathrm d} t b(t) \sqrt{P(k,t)} {\mathrm d} V / {\mathrm d}t$, where $b(t)$ is
    !% galaxy bias, and $P(k,t)$ is the non-linear galaxy power spectrum.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Numerical_Integration
    implicit none
    integer                                     , intent(in   ) :: iBin
    double precision                            , intent(in   ) :: timeMinimum         , timeMaximum
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    logical                                                     :: integrationReset

    lssBin                    =iBin
    integrationReset          =.true.
    Galaxy_Root_Power_Spectrum=Integrate(                          &
         &                               timeMinimum             , &
         &                               timeMaximum             , &
         &                               LSS_Integrand           , &
         &                               integrandFunction       , &
         &                               integrationWorkspace    , &
         &                               toleranceRelative=1.0d-2, &
         &                               reset=integrationReset    &
         &                              )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return
  end function Galaxy_Root_Power_Spectrum

  double precision function Angular_Power_Integrand(wavenumber)
    !% Integrand for large scale structure variance computed using survey mask angular power spectrum.
    implicit none
    double precision, intent(in   ) :: wavenumber
    integer                         :: iField                 , jField               , &
         &                             l
    double precision                :: x0i                    , x1i                  , &
         &                             x0j                    , x1j                  , &
         &                             powerSpectrumI         , powerSpectrumJ       , &
         &                             surveyDistanceMinimum  , surveyDistanceMaximum, &
         &                             angularFactor

    Angular_Power_Integrand=0.0d0
    if (wavenumber <= 0.0d0) return
    wavenumberGlobal=wavenumber
    surveyDistanceMinimum                                                                 &
         &     =cosmologyFunctions_  %distanceComoving           (                        &
         &       cosmologyFunctions_ %cosmicTime                  (                       &
         &        cosmologyFunctions_%expansionFactorFromRedshift  (                      &
         &                                                          surveyRedshiftMinimum &
         &                                                         )                      &
         &                                                        )                       &
         &                                                       )
    surveyDistanceMaximum                                                                 &
         &     =cosmologyFunctions_  %distanceComoving           (                        &
         &       cosmologyFunctions_ %cosmicTime                  (                       &
         &        cosmologyFunctions_%expansionFactorFromRedshift  (                      &
         &                                                          surveyRedshiftMaximum &
         &                                                         )                      &
         &                                                        )                       &
         &                                                       )
    do iField=1,surveyGeometry_%fieldCount()
       if (timeMinimumI(iField) >= timeMaximumI(iField)) cycle
       powerSpectrumI=+surveyGeometry_%solidAngle(             iField)   &
            &         *Galaxy_Root_Power_Spectrum(                       &
            &                                                  binI    , &
            &                                     timeMinimumI(iField) , &
            &                                     timeMaximumI(iField)   &
            &                                    )
       if (surveyRedshiftMinimum <= 0.0d0) then
          x0i=   +0.0d0
       else
          x0i=                                                                                &
               & +wavenumberGlobal                                                            &
               & *    surveyDistanceMinimum  
       end if
       x1i=                                                                                   &
            &    +wavenumberGlobal                                                            &
            &    *min(                                                                        &
            &         surveyDistanceMaximum                                                 , &
            &         surveyGeometry_%distanceMaximum(10.0d0**logMassBinCenter(binI),iField)  &
            &        )
       do jField=1,surveyGeometry_%fieldCount()
          if (timeMinimumJ(jField) >= timeMaximumJ(jField)) cycle
          powerSpectrumJ=+surveyGeometry_%solidAngle(             jField)   &
               &         *Galaxy_Root_Power_Spectrum(                       &
               &                                                  binJ    , &
               &                                     timeMinimumJ(jField) , &
               &                                     timeMaximumJ(jField)   &
               &                                    )
          if (surveyRedshiftMinimum <= 0.0d0) then
             x0j=   +0.0d0
          else
             x0j=                                                                                &
                  & +wavenumberGlobal                                                            &
                  & *    surveyDistanceMinimum
          end if
          x1j=                                                                                   &
               &    +wavenumberGlobal                                                            &
               &    *min(                                                                        &
               &         surveyDistanceMaximum                                                 , &
               &         surveyGeometry_%distanceMaximum(10.0d0**logMassBinCenter(binJ),jField)  &
               &       )
          angularFactor=0.0d0
          !$omp parallel do reduction(+:angularFactor)
          do l=0,surveyGeometry_%angularPowerMaximumDegree()
             angularFactor=                                               &
                  &        +angularFactor                                 &
                  &        +dble(2*l+1)                                   &
                  &        *surveyGeometry_%angularPower(iField,jField,l) &
                  &        *Angular_Power_Radial_Term   (x0i   ,x1i   ,l) &
                  &        *Angular_Power_Radial_Term   (x0j   ,x1j   ,l)
          end do
          !$omp end parallel do
          Angular_Power_Integrand=Angular_Power_Integrand+powerSpectrumI*powerSpectrumJ*angularFactor
       end do
    end do
    Angular_Power_Integrand=Angular_Power_Integrand/wavenumberGlobal**4
    return
  end function Angular_Power_Integrand

  double precision function Angular_Power_Radial_Term(x0,x1,l)
    !% Computes the radial term in the expression for large scale structure variance.
    use Numerical_Constants_Math
    use Gamma_Functions
    use Hypergeometric_Functions
    implicit none
    double precision, intent(in   ) :: x0                , x1
    integer         , intent(in   ) :: l
    double precision, parameter     :: xMaximum  = 512.0d0
    double precision, parameter     :: aMinimum  =-750.0d0
    integer         , save          :: lPrevious =-1
    double precision, save          :: x0Previous=-1.0d0 , x1Previous=-1.0d0
    !$omp threadprivate(lPrevious,x0Previous,x1Previous)
    double precision, save          :: h0                , h1, &
         &                             logGammas
    !$omp threadprivate(h0,h1,logGammas)
    double precision                :: a0                , a1

    ! Evaluate combination of logarithms of Gamma functions.
    if (l /= lPrevious) then
       logGammas=+Gamma_Function_Logarithmic(0.5d0*(3.0d0+dble(l))) &
            &    -Gamma_Function_Logarithmic(       1.5d0+dble(l) ) &
            &    -Gamma_Function_Logarithmic(0.5d0*(5.0d0+dble(l)))
    end if
    ! Evaluate hypergeometric terms and power-law terms, catching the x=0 special case.
    if (l /= lPrevious .or. x0 /= x0Previous) then
       if (x0 <= 0.0d0 .or. x0 > xMaximum) then 
          h0     =0.0d0
       else
          a0=                                                               &
               &      +logGammas                                            &
               &      +dble(3+l)                                            &
               &      *log (x0 )                                            &
               &      -dble(2+l)                                            &
               &      *ln2
          if (a0 > aMinimum) then
             h0     =                                                          &
                  & +Hypergeometric_pFq(                                       &
                  &                     [              0.5d0*(3.0d0+dble(l))], &
                  &                     [1.5d0+dble(l),0.5d0*(5.0d0+dble(l))], &
                  &                     -x0**2/4.0d0                           &
                  &                    )                                       &
                  & *exp(a0)
          else
             h0=0.0d0
          end if
       end if
       x0Previous=x0
    end if
    if (l /= lPrevious .or. x1 /= x1Previous) then
       if (x1 <= 0.0d0 .or. x1 > xMaximum) then 
          h1     =0.0d0
       else
          a1=                                                               &
               &      +logGammas                                            &
               &      +dble(3+l)                                            &
               &      *log (x1 )                                            &
               &      -dble(2+l)                                            &
               &      *ln2
          if (a1 > aMinimum) then
             h1     =                                                          &
                  & +Hypergeometric_pFq(                                       &
                  &                     [              0.5d0*(3.0d0+dble(l))], &
                  &                     [1.5d0+dble(l),0.5d0*(5.0d0+dble(l))], &
                  &                     -x1**2/4.0d0                           &
                  &                    )                                       &
                  & *exp(a1)
          else
             h1=0.0d0
          end if
       end if
       x1Previous=x1
    end if
    Angular_Power_Radial_Term=          &
         &                    +sqrt(Pi) &
         &                    *(        &
         &                      +h1     &
         &                      -h0     &
         &                     )
    lPrevious=l
    return
  end function Angular_Power_Radial_Term
  
  double precision function Volume_Integrand(time)
    !% Integral for comoving volume.
    use Cosmology_Functions
    implicit none
    double precision, intent(in   ) :: time

    Volume_Integrand=cosmologyFunctions_%comovingVolumeElementTime(time)
    return
  end function Volume_Integrand

  double precision function Mass_Function_Time_Integrand_I(timePrime)
    !% Integral for comoving volume.
    use Cosmology_Functions
    use FGSL
    use Numerical_Integration
    implicit none
    double precision                            , intent(in   ) :: timePrime
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    logical                                                     :: integrationReset
    double precision                                            :: massFunction

    time=timePrime
    integrationReset=.true.
    massFunction=Integrate(                            &
         &                 logMassLower              , &
         &                 logMassUpper              , &
         &                 Mass_Function_Integrand_I , &
         &                 integrandFunction         , &
         &                 integrationWorkspace      , &
         &                 toleranceRelative=1.0d-3  , &
         &                 reset=integrationReset      &
         &                )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    Mass_Function_Time_Integrand_I=massFunction*cosmologyFunctions_%comovingVolumeElementTime(time)
    return
  end function Mass_Function_Time_Integrand_I

  double precision function Mass_Function_Integrand_I(logMass)
    !% Integral for mass function.
    use Halo_Mass_Functions
    use Conditional_Mass_Functions
    implicit none
    double precision                               , intent(in   ) :: logMass
    class            (conditionalMassFunctionClass), pointer       :: conditionalMassFunction_
    class            (haloMassFunctionClass       ), pointer       :: haloMassFunction_
    double precision                                               :: mass

    conditionalMassFunction_ => conditionalMassFunction()
    haloMassFunction_        => haloMassFunction       ()
    mass=10.0d0**logMass
    Mass_Function_Integrand_I=+haloMassFunction_%differential(time,mass)                         &
         &                    *                                    mass                          &
         &                    *log(10.0d0)                                                       &
         &                    *max(                                                              &
         &                         +conditionalMassFunction_%massFunction(mass,massBinMinimumI)  &
         &                         -conditionalMassFunction_%massFunction(mass,massBinMaximumI), &
         &                          0.0d0                                                        &
         &                        )
    return
  end function Mass_Function_Integrand_I

  double precision function LSS_Integrand(timePrime)
    !% Integral for LSS contribution to the covariance matrix.
    use Cosmology_Functions
    use FGSL
    use Power_Spectra_Nonlinear
    use Numerical_Interpolation
    implicit none
    double precision                             , intent(in   ) :: timePrime
    class           (powerSpectrumNonlinearClass), pointer       :: powerSpectrumNonlinear_
    type            (fgsl_interp                )                :: interpolationObject
    type            (fgsl_interp_accel          )                :: interpolationAccelerator
    logical                                                      :: interpolationReset
    double precision                                             :: bias                    , powerSpectrumValue

    ! Get required objects.
    powerSpectrumNonLinear_ => powerSpectrumNonLinear()
    ! Copy the time to module scope.
    time=timePrime
    ! Get the bias-mass function product for the I bin.
    interpolationReset=.true.
    bias=Interpolate(timeTable,biasTable(:,lssBin),interpolationObject,interpolationAccelerator,time,reset=interpolationReset)
    call Interpolate_Done(interpolationObject,interpolationAccelerator,interpolationReset)
    ! Get the nonlinear power spectrum for the current wavenumber and time.
    powerSpectrumValue=powerSpectrumNonlinear_%value(waveNumberGlobal,time)
    ! Return the cross-correlation biased power spectrum multiplied by the volume element.
    LSS_Integrand=bias*sqrt(powerSpectrumValue)*cosmologyFunctions_%comovingVolumeElementTime(time)
    return
  end function LSS_Integrand

  double precision function Bias_Integrand_I(logMass)
    !% Integral for bias.
    use Dark_Matter_Halo_Biases
    implicit none
    double precision, intent(in   ) :: logMass
    double precision                :: mass

    mass=10.0d0**logMass
    Bias_Integrand_I=Mass_Function_Integrand_I(logMass)*Dark_Matter_Halo_Bias(mass,time)
  return
  end function Bias_Integrand_I

  double precision function Halo_Occupancy_Time_Integrand(timePrime)
    !% Integral for comoving volume.
    use Cosmology_Functions
    use FGSL
    use Numerical_Integration
    implicit none
    double precision                            , intent(in   ) :: timePrime
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    logical                                                     :: integrationReset
    double precision                                            :: massFunction

    time=timePrime
    integrationReset=.true.
    massFunction=Integrate(                          &
         &                 logMassLower            , &
         &                 logMassUpper            , &
         &                 Halo_Occupancy_Integrand, &
         &                 integrandFunction       , &
         &                 integrationWorkspace    , &
         &                 toleranceRelative=1.0d-3, &
         &                 reset=integrationReset    &
         &                )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    Halo_Occupancy_Time_Integrand=massFunction*cosmologyFunctions_%comovingVolumeElementTime(time)
    return
  end function Halo_Occupancy_Time_Integrand

  double precision function Halo_Occupancy_Integrand(logMass)
    !% Integral for mass function.
    use Halo_Mass_Functions
    use Conditional_Mass_Functions
    implicit none
    double precision                               , intent(in   ) :: logMass
    class            (conditionalMassFunctionClass), pointer       :: conditionalMassFunction_
    class            (haloMassFunctionClass       ), pointer       :: haloMassFunction_
    double precision                                               :: mass
    
    conditionalMassFunction_ => conditionalMassFunction()
    haloMassFunction_        => haloMassFunction       ()
    mass=10.0d0**logMass
    Halo_Occupancy_Integrand=+haloMassFunction_%differential(time,mass)                         &
         &                   *                                    mass                          &
         &                   *log(10.0d0)                                                       &
         &                   *max(                                                              &
         &                        +conditionalMassFunction_%massFunction(mass,massBinMinimumI)  &
         &                        -conditionalMassFunction_%massFunction(mass,massBinMaximumI), &
         &                         0.0d0                                                        &
         &                       )                                                              &
         &                   *max(                                                              &
         &                        +conditionalMassFunction_%massFunction(mass,massBinMinimumJ)  &
         &                        -conditionalMassFunction_%massFunction(mass,massBinMaximumJ), &
         &                         0.0d0                                                        &
         &                       )
    return
  end function Halo_Occupancy_Integrand

  subroutine Compute_Volume_Normalizations(logMass,surveyGeometry_,redshiftMinimum,redshiftMaximum,timeMinimum,timeMaximum,volumeNormalization)
    !% Compute volume normalization factors for LSS covariance calculations.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Numerical_Integration
    use Geometry_Surveys
    implicit none
    double precision                            , intent(in   )               :: logMass             , redshiftMinimum, &
         &                                                                       redshiftMaximum
    class           (surveyGeometryClass       ), intent(inout)               :: surveyGeometry_
    double precision                            , intent(  out), dimension(:) :: timeMinimum         , timeMaximum    , &
         &                                                                       volumeNormalization
    integer                                                                   :: iField
    type            (fgsl_function             )                              :: integrandFunction
    type            (fgsl_integration_workspace)                              :: integrationWorkspace
    logical                                                                   :: integrationReset

    do iField=1,surveyGeometry_%fieldCount() 
       ! Find integration limits for this bin. We want the maximum of the volumes associated with the two bins.
       timeMaximum(iField)=    cosmologyFunctions_%cosmicTime            (cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum       ))
       timeMinimum(iField)=min(                                                                                                                         &
            &                  max(                                                                                                                     &
            &                      cosmologyFunctions_%cosmicTime            (cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum       )), &
            &                      cosmologyFunctions_%timeAtDistanceComoving(surveyGeometry_    %distanceMaximum            (10.0d0**logMass,iField))  &
            &                     )                                                                                                                   , &
            &                  timeMaximum(iField)                                                                                                      &
            &                 )
       ! Get the normalizing volume integral for bin i.
       integrationReset=.true.
       volumeNormalization(iField)= Integrate(                          &
            &                                 timeMinimum(iField)     , &
            &                                 timeMaximum(iField)     , &
            &                                 Volume_Integrand        , &
            &                                 integrandFunction       , &
            &                                 integrationWorkspace    , &
            &                                 toleranceRelative=1.0d-3, &
            &                                 reset=integrationReset    &
            &                                )                          &
            &                      *surveyGeometry_%solidAngle(iField)
      call Integrate_Done(integrandFunction,integrationWorkspace)
    end do
    return
  end subroutine Compute_Volume_Normalizations

  subroutine Variance_LSS_Window_Function(massBinCount,redshiftMinimum,redshiftMaximum,varianceLSS)
    !% Compute variance due to large scale structure by directly summing over the Fourier transform
    !% of the survey selection function.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Memory_Management
    use FFTW3
    use Numerical_Constants_Math
    use Galacticus_Display
    use Input_Parameters
    implicit none
    integer                                     , intent(in   )                   :: massBinCount
    double precision                            , intent(in   )                   :: redshiftMinimum,redshiftMaximum
    double precision                            , intent(  out), dimension(:,:  ) :: varianceLSS
    complex         (c_double_complex          ), allocatable  , dimension(:,:,:) :: windowFunctionI,windowFunctionJ
    logical                                     , save                            :: functionInitialized=.false.
    integer                                                                       :: i,j,u,w,v,taskCount,taskTotal,fieldCount,iField,massFunctionCovarianceFFTGridSize
    double precision                                                              :: waveNumberU,waveNumberV,waveNumberW, variance,powerSpectrumI,powerSpectrumJ,powerSpectrum,normalizationI,normalizationJ,multiplier,boxLength
    
    if (.not.functionInitialized) then
       ! Read controlling parameters.
       !@ <inputParameter>
       !@   <name>massFunctionCovarianceFFTGridSize</name>
       !@   <defaultValue>$64$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The size of the FFT grid to use in computing window functions for mass function covariance matrices.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('massFunctionCovarianceFFTGridSize',massFunctionCovarianceFFTGridSize,defaultValue=64)
       functionInitialized=.true.
    end if
    ! Allocate arrays for times and volume normalizations.
    fieldCount=surveyGeometry_%fieldCount()
    ! Allocate arrays for survey window functions if these will be used.
    allocate(windowFunctionI(                                   &
         &                   massFunctionCovarianceFFTGridSize, &
         &                   massFunctionCovarianceFFTGridSize, &
         &                   massFunctionCovarianceFFTGridSize  &
         &                  )                                   &
         &  )
    allocate(windowFunctionJ(                                   &
         &                   massFunctionCovarianceFFTGridSize, &
         &                   massFunctionCovarianceFFTGridSize, &
         &                   massFunctionCovarianceFFTGridSize  &
         &                  )                                   &
         &  )
    taskTotal  =massBinCount*(massBinCount+1)/2
    taskCount  =0
    varianceLSS=0.0d0
    do i   =1,massBinCount
       binI=i
       massBinCenterI =10.0d0** logMassBinCenter(i)
       massBinMinimumI=10.0d0**(logMassBinCenter(i)-0.5d0*log10MassBinWidth(i))
       massBinMaximumI=10.0d0**(logMassBinCenter(i)+0.5d0*log10MassBinWidth(i))
       do j=i,massBinCount
          binJ=j
          massBinCenterJ =10.0d0** logMassBinCenter(j)
          massBinMinimumJ=10.0d0**(logMassBinCenter(j)-0.5d0*log10MassBinWidth(j))
          massBinMaximumJ=10.0d0**(logMassBinCenter(j)+0.5d0*log10MassBinWidth(j))
          ! Update progress.
          call Galacticus_Display_Counter(                                              &
               &                          int(100.0d0*dble(taskCount)/dble(taskTotal)), &
               &                          isNew=(taskCount==0)                          &
               &                         )
          taskCount=taskCount+1
          ! Compute window functions for this pair of cells.
          call surveyGeometry_%windowFunctions(                                   &
               &                               massBinCenterI                   , &
               &                               massBinCenterJ                   , &
               &                               massFunctionCovarianceFFTGridSize, &
               &                               boxLength                        , &
               &                               windowFunctionI                  , &
               &                               windowFunctionJ                    &
               &                              )
          ! Integrate the large-scale structure variance over the window functions. Note that
          ! FFTW3 works in terms of inverse wavelengths, not wavenumbers (as we want to use
          ! here). According to FFTW3 documentation
          ! (http://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029)
          ! for a 1-D FFT, the k^th output corresponds to "frequency" k/T where T is the total
          ! period of the box. This then is the side length of each cell of the FFT in terms of
          ! inverse wavelengths, if we associated T with the total box length, L. In terms of
          ! wavenumber, that means that each cell of the FFT has side of length 2π/L.
          variance=0.0d0
          !$omp parallel private (u,v,w,waveNumberU,waveNumberV,waveNumberW,multiplier,normalizationI,normalizationJ,powerSpectrumI,powerSpectrumJ,powerSpectrum,iField)
          call allocateArray(volumeNormalizationI,[fieldCount])
          call allocateArray(volumeNormalizationJ,[fieldCount])
          call allocateArray(timeMinimumI        ,[fieldCount])
          call allocateArray(timeMinimumJ        ,[fieldCount])
          call allocateArray(timeMaximumI        ,[fieldCount])
          call allocateArray(timeMaximumJ        ,[fieldCount])
          call Compute_Volume_Normalizations(                         &
               &                             logMassBinCenter    (i), &
               &                             surveyGeometry_        , &
               &                             redshiftMinimum        , &
               &                             redshiftMaximum        , &
               &                             timeMinimumI           , &
               &                             timeMaximumI           , &
               &                             volumeNormalizationI     &
               &                            )   
          call Compute_Volume_Normalizations(                         &
               &                             logMassBinCenter    (j), &
               &                             surveyGeometry_        , &
               &                             redshiftMinimum        , &
               &                             redshiftMaximum        , &
               &                             timeMinimumJ           , &
               &                             timeMaximumJ           , &
               &                             volumeNormalizationJ     &
               &                            )
          !$omp do reduction(+:variance)
          do u      =1,massFunctionCovarianceFFTGridSize/2+1
             waveNumberU      =FFTW_Wavenumber(u,massFunctionCovarianceFFTGridSize)*2.0d0*Pi/boxLength
             do v   =1,massFunctionCovarianceFFTGridSize/2+1
                waveNumberV   =FFTW_Wavenumber(v,massFunctionCovarianceFFTGridSize)*2.0d0*Pi/boxLength
                do w=1,massFunctionCovarianceFFTGridSize/2+1
                   waveNumberW=FFTW_Wavenumber(w,massFunctionCovarianceFFTGridSize)*2.0d0*Pi/boxLength
                   ! Compute the wavenumber for this cell.
                   waveNumberGlobal=sqrt(waveNumberU**2+waveNumberV**2+waveNumberW**2)
                   ! Find the power spectrum for this wavenumber.
                   if (waveNumberGlobal > 0.0d0) then      
                      ! Integrate the power spectrum, weighted by the galaxy bias, over the
                      ! volume of interest. Then normalize by that volume.
                      normalizationI=0.0d0
                      normalizationJ=0.0d0
                      powerSpectrumI=0.0d0
                      powerSpectrumJ=0.0d0
                      do iField=1,fieldCount
                         powerSpectrumI  =                                                   &
                              &           +powerSpectrumI                                    &
                              &           +surveyGeometry_%solidAngle(             iField)   &
                              &           *Galaxy_Root_Power_Spectrum(                       &
                              &                                                    binI    , &
                              &                                       timeMinimumI(iField) , &
                              &                                       timeMaximumI(iField)   &
                              &                                      )
                         normalizationI  =                                                   &
                              &           +normalizationI                                    &
                              &           +volumeNormalizationI      (             iField)
                         powerSpectrumJ  =                                                   &
                              &           +powerSpectrumJ                                    &
                              &           +surveyGeometry_%solidAngle(             iField)   &
                              &           *Galaxy_Root_Power_Spectrum(                       &
                              &                                                    binJ    , &
                              &                                       timeMinimumJ(iField) , &
                              &                                       timeMaximumJ(iField)   &
                              &                                      )
                         normalizationJ  =                                                   &
                              &           +normalizationJ                                    &
                              &           +volumeNormalizationJ      (             iField)
                      end do
                      powerSpectrum=                &
                           &        +powerSpectrumI &
                           &        /normalizationI &
                           &        *powerSpectrumJ &
                           &        /normalizationJ
                   else
                      powerSpectrum=0.0d0
                   end if
                   ! Add the contribution from this cell to the total variance.
                   multiplier=2.0d0
                   if     (                                                &
                        &  u == massFunctionCovarianceFFTGridSize/2+1 .or. &
                        &  v == massFunctionCovarianceFFTGridSize/2+1 .or. &
                        &  w == massFunctionCovarianceFFTGridSize/2+1      &
                        & ) multiplier=1.0d0
                   variance=+variance                            &
                        &   +multiplier                          &
                        &   *powerSpectrum                       &
                        &   *real(                               &
                        &                windowFunctionI(u,v,w)  &
                        &         *conjg(windowFunctionJ(u,v,w)) &
                        &        )
                end do
             end do
          end do
          !$omp end do
          call deallocateArray(volumeNormalizationI)
          call deallocateArray(volumeNormalizationJ)
          call deallocateArray(timeMinimumI        )
          call deallocateArray(timeMinimumJ        )
          call deallocateArray(timeMaximumI        )
          call deallocateArray(timeMaximumJ        )
          !$omp end parallel
          ! Normalize the variance. We multiply by (2π/L)³ to account for the volume of each FFT
          ! cell, and divide by (2π)³ as defined in eqn. (66) of Smith (2012; MNRAS; 426; 531).
          varianceLSS(i,j)=dble(variance)/boxLength**3
       end do
    end do
    if (allocated(windowFunctionI)) deallocate(windowFunctionI)
    if (allocated(windowFunctionJ)) deallocate(windowFunctionJ)
    call Galacticus_Display_Counter_Clear()
    return
  end subroutine Variance_LSS_Window_Function

  subroutine Variance_LSS_Angular_Spectrum(massBinCount,redshiftMinimum,redshiftMaximum,varianceLSS)
    !% Compute variance due to large scale structure by integration over the angular power spectrum.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Galacticus_Display
    use Numerical_Constants_Math
    use Numerical_Integration
    use Memory_Management
    implicit none
    integer                                     , intent(in   )                   :: massBinCount
    double precision                            , intent(in   )                   :: redshiftMinimum,redshiftMaximum
    double precision                            , intent(  out), dimension(:,:  ) :: varianceLSS
    ! Dimensionless factor controlling the highest wavenumber to be used when integrating over
    ! angular power spectra.
    double precision                            , parameter                       :: wavenumberMaximumFactor=1.0d0
    integer                                                                       :: i,j,fieldCount,iField,taskCount,taskTotal
    double precision                                                              :: wavenumberMinimum,wavenumberMaximum,distanceMaximum
    logical                                                                       :: integrationReset
    type            (fgsl_function             )                                  :: integrandFunction
    type            (fgsl_integration_workspace)                                  :: integrationWorkspace

    fieldCount=surveyGeometry_%fieldCount()
    call omp_set_nested(.true.)
    taskTotal  =massBinCount*(massBinCount+1)/2
    taskCount  =0
    !$omp parallel do private (i,j,wavenumberMinimum,wavenumberMaximum,integrationReset,integrandFunction,integrationWorkspace) schedule (dynamic)
    do i=1,massBinCount
       ! Allocate arrays for times and volume normalizations.
       if (.not.allocated(volumeNormalizationI)) then
          call allocateArray(volumeNormalizationI,[fieldCount])
          call allocateArray(volumeNormalizationJ,[fieldCount])
          call allocateArray(timeMinimumI        ,[fieldCount])
          call allocateArray(timeMinimumJ        ,[fieldCount])
          call allocateArray(timeMaximumI        ,[fieldCount])
          call allocateArray(timeMaximumJ        ,[fieldCount])
       end if
       binI=i
       call Compute_Volume_Normalizations(                        &
            &                             logMassBinCenter(binI), &
            &                             surveyGeometry_       , &
            &                             redshiftMinimum       , &
            &                             redshiftMaximum       , &
            &                             timeMinimumI          , &
            &                             timeMaximumI          , &
            &                             volumeNormalizationI    &
            &                            )
       do j=binI,massBinCount
          ! Update progress.
          call Galacticus_Display_Counter(                                              &
               &                          int(100.0d0*dble(taskCount)/dble(taskTotal)), &
               &                          isNew=(taskCount==0)                          &
               &                         )
          binJ=j
          call Compute_Volume_Normalizations(                        &
               &                             logMassBinCenter(binJ), &
               &                             surveyGeometry_       , &
               &                             redshiftMinimum       , &
               &                             redshiftMaximum       , &
               &                             timeMinimumJ          , &
               &                             timeMaximumJ          , &
               &                             volumeNormalizationJ    &
               &                            )
          distanceMaximum=huge(1.0d0)
          do iField=1,surveyGeometry_%fieldCount()
             distanceMaximum=min(distanceMaximum,surveyGeometry_%distanceMaximum(10.0d0**logMassBinCenter(binI),iField),surveyGeometry_%distanceMaximum(10.0d0**logMassBinCenter(binJ),iField))
          end do
          wavenumberMinimum=0.0d0
          wavenumberMaximum=wavenumberMaximumFactor                                               &
               &            *max(                                                                 &
               &                  1.0d0                                                        ,  &
               &                  surveyGeometry_%angularPowerMaximumDegree()                     &
               &                 /2.0d0                                                           &
               &                 /Pi                                                              &
               &                )                                                                 &
               &            /distanceMaximum
          integrationReset=.true.
          varianceLSS(binI,binJ)=                                               &
               &           +2.0d0                                               &
               &           /Pi                                                  &
               &           /sum(volumeNormalizationI)**2                        &
               &           /sum(volumeNormalizationJ)**2                        &
               &           *Integrate(                                          &
               &                      wavenumberMinimum                       , &
               &                      wavenumberMaximum                       , &
               &                      Angular_Power_Integrand                 , &
               &                      integrandFunction                       , &
               &                      integrationWorkspace                    , &
               &                      toleranceRelative      =1.0d-2          , &
               &                      reset                  =integrationReset  &
               &                     )
          call Integrate_Done(integrandFunction,integrationWorkspace)
          !$omp atomic
          taskCount=taskCount+1
       end do
       ! Allocate arrays for times and volume normalizations.
       if (allocated(volumeNormalizationI)) then
          call deallocateArray(volumeNormalizationI)
          call deallocateArray(volumeNormalizationJ)
          call deallocateArray(timeMinimumI        )
          call deallocateArray(timeMinimumJ        )
          call deallocateArray(timeMaximumI        )
          call deallocateArray(timeMaximumJ        )
       end if
    end do
    !$omp end parallel do
    return
  end subroutine Variance_LSS_Angular_Spectrum

end module Statistics_Mass_Function_Covariance
