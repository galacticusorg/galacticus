!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which computes mass function covariances.

module Statistics_Mass_Function_Covariance
  !% Implements calculations of mass function covariances.
  private
  public :: Mass_Function_Covariance_Matrix

  ! Record of whether this module is intialized.
  logical          :: moduleInitialized=.false.

  ! The cosmic time at which the calculation is to be performed.
  double precision :: time
  !$omp threadprivate(time)

  ! The wavenumber for which LSS integrations are currently being performed.
  double precision :: waveNumber
  !$omp threadprivate(waveNumber)

  ! Integration limits.
  double precision :: logMassLower,logMassUpper

  ! Minimum and maximum masses for the bins being considered.
  double precision :: log10MassBinWidth,logMassBinWidth
  integer          :: binI,binJ
  double precision :: massBinCenterI,massBinMinimumI,massBinMaximumI
  double precision :: massBinCenterJ,massBinMinimumJ,massBinMaximumJ

  ! Table of biases.
  integer          :: timeBinCount
double precision, dimension(:), allocatable :: timeTable
double precision, dimension(:,:), allocatable :: biasTable

contains

  subroutine Mass_Function_Covariance_Matrix(redshiftMinimum,redshiftMaximum,massBinCount,massMinimum,massMaximum,massFunctionObserved,includePoisson,includeHalo,includeLSS&
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
    use Cosmology_Functions
    use Power_Spectra_Nonlinear
    use Galacticus_Error
    use Galacticus_Display
    use Geometry_Surveys
    implicit none
    integer                            , intent(in   )                                        :: massBinCount
    double precision                   , intent(in   )                                        :: redshiftMinimum,redshiftMaximum&
         &,massMinimum,massMaximum
    logical                            , intent(in   )                                        :: includePoisson,includeHalo&
         &,includeLSS
    double precision                   , intent(inout), allocatable, dimension(:    )         :: mass
    double precision                   , intent(inout), allocatable, dimension(:    ), target :: massFunction,massFunctionObserved
    double precision                   , intent(inout), allocatable, dimension(:,:  )         :: covariance,covariancePoisson &
         &,covarianceHalo,covarianceLSS,correlation
    double precision                   ,                allocatable, dimension(:    )         :: logMassBinCenter,volume
    double precision                   ,                allocatable, dimension(:,:  )         :: varianceLSS
    complex(c_double_complex          ),                allocatable, dimension(:,:,:)         :: windowFunctionI,windowFunctionJ
    double precision                   ,                pointer    , dimension(:    )         :: massFunctionUse
    double precision                   , parameter                                            :: timePointsPerDecade=100
    logical                                                                                   :: integrationReset
    integer                                                                                   :: i,j,u,w,v,taskCount,taskTotal &
         &,massFunctionCovarianceFFTGridSize,iTime
    double precision                                                                          :: logMassMinimum ,logMassMaximum&
         &,normalization,boxLength,waveNumberU,waveNumberV,waveNumberW ,powerSpectrum ,massFunctionCovarianceHaloMassMinimum&
         &,massFunctionCovarianceHaloMassMaximum,timeMinimum,timeMaximum,volumeNormalization
    double precision                                                                          :: variance,multiplier
    type   (fgsl_function             )                                                       :: integrandFunction
    type   (fgsl_integration_workspace)                                                       :: integrationWorkspace
    type   (c_ptr                     )                                                       :: parameterPointer
    
    ! Initialize the module if necessary.
    if (.not.moduleInitialized) then
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

    ! Determine number of times over which to tabulate bias.
    timeMaximum =Cosmology_Age(Expansion_Factor_from_Redshift(redshiftMinimum))
    timeMinimum =Cosmology_Age(Expansion_Factor_from_Redshift(redshiftMaximum))
    timeBinCount=int(log10(timeMaximum/timeMinimum)*dble(timePointsPerDecade))+1

    ! Allocate arrays.
    call Alloc_Array(mass             ,[massBinCount             ])
    call Alloc_Array(logMassBinCenter ,[massBinCount             ])
    call Alloc_Array(massFunction     ,[massBinCount             ])
    call Alloc_Array(volume           ,[massBinCount             ])
    call Alloc_Array(covariance       ,[massBinCount,massBinCount])
    call Alloc_Array(covariancePoisson,[massBinCount,massBinCount])
    call Alloc_Array(covarianceHalo   ,[massBinCount,massBinCount])
    call Alloc_Array(covarianceLSS    ,[massBinCount,massBinCount])
    call Alloc_Array(correlation      ,[massBinCount,massBinCount])
    call Alloc_Array(varianceLSS      ,[massBinCount,massBinCount])
    call Alloc_Array(timeTable        ,[timeBinCount             ])
    call Alloc_Array(biasTable        ,[timeBinCount,massBinCount])

    ! Create time bins.
    timeTable=Make_Range(timeMinimum,timeMaximum,timeBinCount,rangeType=rangeTypeLogarithmic)

    ! Create mass bins.
    logMassMinimum   =log10(massMinimum)
    logMassMaximum   =log10(massMaximum)
    logMassBinCenter =Make_Range(logMassMinimum,logMassMaximum,massBinCount,rangeType=rangeTypeLinear)
    log10MassBinWidth=logMassBinCenter(2)-logMassBinCenter(1)
    logMassBinWidth  =log(10.0d0)*log10MassBinWidth
    mass             =10.0d0**logMassBinCenter

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

    ! Compute the mass function and bias averaged over each bin.
    do i=1,massBinCount

       ! Find limits on mass for this bin.
       massBinCenterI =10.0** logMassBinCenter(i)
       massBinMinimumI=10.0**(logMassBinCenter(i)-0.5d0*log10MassBinWidth)
       massBinMaximumI=10.0**(logMassBinCenter(i)+0.5d0*log10MassBinWidth)

       ! Find integration limits for this bin.
       timeMaximum=    Cosmology_Age(              Expansion_Factor_from_Redshift  (redshiftMinimum))
       timeMinimum=max(                                                                                &
            &          Cosmology_Age              (Expansion_Factor_from_Redshift  (redshiftMaximum)), &
            &          Time_From_Comoving_Distance(Geometry_Survey_Distance_Maximum(massBinCenterI ))  &
            &         )
       ! Get the normalizing volume integral.
       integrationReset=.true.
       volumeNormalization=Integrate(                          &
            &                        timeMinimum             , &
            &                        timeMaximum             , &
            &                        Volume_Integrand        , &
            &                        parameterPointer        , &
            &                        integrandFunction       , &
            &                        integrationWorkspace    , &
            &                        toleranceRelative=1.0d-3, &
            &                        reset=integrationReset    &
            &                       ) 

       ! Integrate mass function over the bin.
       integrationReset=.true.
       massFunction(i)=Integrate(                                 &
            &                    timeMinimum                    , &
            &                    timeMaximum                    , &
            &                    Mass_Function_Time_Integrand_I , &
            &                    parameterPointer               , &
            &                    integrandFunction              , &
            &                    integrationWorkspace           , &
            &                    toleranceRelative=1.0d-3       , &
            &                    reset=integrationReset           &
            &                   )                                 &
            &                   /logMassBinWidth                  &
            &                   /volumeNormalization
       call Integrate_Done(integrandFunction,integrationWorkspace)

       ! Find the effective volume of the survey at this mass.
       volume(i)=Geometry_Survey_Volume_Maximum(massBinCenterI)

       ! Tabulate the bias as a function of time in this bin.
       integrationReset=.true.
       do iTime=1,timeBinCount
          time=timeTable(iTime)
          biasTable(iTime,i)=  Integrate(                  &
               &                 logMassLower            , &
               &                 logMassUpper            , &
               &                 Bias_Integrand_I        , &
               &                 parameterPointer        , &
               &                 integrandFunction       , &
               &                 integrationWorkspace    , &
               &                 toleranceRelative=1.0d-2, &
               &                 reset=integrationReset    &
               &                )                          &
               & /logMassBinWidth
       end do
       call Integrate_Done(integrandFunction,integrationWorkspace)

    end do

    ! Compute LSS variance if necessary.
    if (includeLSS) then
       ! Allocate arrays for survey window functions.
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
       ! Compute large-scale structure variance for each cell pair.
       taskTotal=massBinCount*(massBinCount+1)/2
       taskCount=0
       do i   =1,massBinCount
          binI=i
          massBinCenterI =10.0** logMassBinCenter(i)
          massBinMinimumI=10.0**(logMassBinCenter(i)-0.5d0*log10MassBinWidth)
          massBinMaximumI=10.0**(logMassBinCenter(i)+0.5d0*log10MassBinWidth)
          do j=i,massBinCount
             binJ=j
             massBinCenterJ =10.0** logMassBinCenter(j)
             massBinMinimumJ=10.0**(logMassBinCenter(j)-0.5d0*log10MassBinWidth)
             massBinMaximumJ=10.0**(logMassBinCenter(j)+0.5d0*log10MassBinWidth)

             ! Update progress.
             call Galacticus_Display_Counter(int(100.0d0*dble(taskCount)/dble(taskTotal)),isNew=(taskCount==0))
             taskCount=taskCount+1
             
             ! Compute window functions for this pair of cells.
             call Geometry_Survey_Window_Functions(                                   &
                  &                                massBinCenterI                   , &
                  &                                massBinCenterJ                   , &
                  &                                boxLength                        , &
                  &                                massFunctionCovarianceFFTGridSize, &
                  &                                windowFunctionI                  , &
                  &                                windowFunctionJ                    &
                  &                               )

             ! Find integration limits for this bin. We want the maximum of the volumes associated with the two bins.
             timeMaximum=    Cosmology_Age(                  Expansion_Factor_from_Redshift  (redshiftMinimum))
             timeMinimum=max(                                                                                    &
                  &          Cosmology_Age              (    Expansion_Factor_from_Redshift  (redshiftMaximum)), &
                  &          Time_From_Comoving_Distance(                                                        &
                  &                                      max(                                                    &
                  &                                          Geometry_Survey_Distance_Maximum(massBinCenterI ),  &
                  &                                          Geometry_Survey_Distance_Maximum(massBinCenterJ )   &
                  &                                         )                                                    &
                  &                                     )                                                        &
                  &         )
             
             ! Get the normalizing volume integral.
             integrationReset=.true.
             volumeNormalization=Integrate(                          &
                  &                        timeMinimum             , &
                  &                        timeMaximum             , &
                  &                        Volume_Integrand        , &
                  &                        parameterPointer        , &
                  &                        integrandFunction       , &
                  &                        integrationWorkspace    , &
                  &                        toleranceRelative=1.0d-3, &
                  &                        reset=integrationReset    &
                  &                       ) 
             call Integrate_Done(integrandFunction,integrationWorkspace)

             ! Integrate the large-scale structure variance over the window functions.
             variance=0.0d0
             !$omp parallel do private (u,v,w,waveNumberU,waveNumberV,waveNumberW,integrationReset,integrandFunction,integrationWorkspace,powerSpectrum,multiplier), reduction (+:variance)
             do u      =1,massFunctionCovarianceFFTGridSize/2+1
                waveNumberU      =FFTW_Wavenumber(u,massFunctionCovarianceFFTGridSize)/boxLength
                do v   =1,massFunctionCovarianceFFTGridSize/2+1
                   waveNumberV   =FFTW_Wavenumber(v,massFunctionCovarianceFFTGridSize)/boxLength
                   do w=1,massFunctionCovarianceFFTGridSize/2+1
                      waveNumberW=FFTW_Wavenumber(w,massFunctionCovarianceFFTGridSize)/boxLength

                      ! Compute the wavenumber for this cell.
                      waveNumber=sqrt(waveNumberU**2+waveNumberV**2+waveNumberW**2)
                      
                      ! Find the power spectrum for this wavenumber.
                      if (waveNumber > 0.0d0) then      
                         ! Integrate the power spectrum, weighted by the galaxy bias, over the volume of interest. Then normalize
                         ! by that volume.
                         integrationReset=.true.
                         powerSpectrum= Integrate(                          &
                              &                   timeMinimum             , &
                              &                   timeMaximum             , &
                              &                   LSS_Integrand           , &
                              &                   parameterPointer        , &
                              &                   integrandFunction       , &
                              &                   integrationWorkspace    , &
                              &                   toleranceRelative=1.0d-2, &
                              &                   reset=integrationReset    &
                              &                  )                          &
                              &        /volumeNormalization                         
                         call Integrate_Done(integrandFunction,integrationWorkspace)
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
                      variance=variance+multiplier*powerSpectrum*real(windowFunctionI(u,v,w)*conjg(windowFunctionJ(u,v,w)))
                   end do
                end do
             end do
             !$omp end parallel do

             ! Normalize the variance.
             varianceLSS(i,j)=real(variance)/(2.0d0*Pi*boxLength)**3
             
          end do
       end do
       deallocate(windowFunctionI)
       deallocate(windowFunctionJ)
       call Galacticus_Display_Counter_Clear()
    end if

    ! Construct the covariance matrix.
    covariancePoisson=0.0d0
    covarianceHalo   =0.0d0
    covarianceLSS    =0.0d0
    do i   =1,massBinCount
       massBinCenterI    =10.0** logMassBinCenter(i)
       massBinMinimumI   =10.0**(logMassBinCenter(i)-0.5d0*log10MassBinWidth)
       massBinMaximumI   =10.0**(logMassBinCenter(i)+0.5d0*log10MassBinWidth)
       do j=i,massBinCount
          massBinCenterJ =10.0** logMassBinCenter(j)
          massBinMinimumJ=10.0**(logMassBinCenter(j)-0.5d0*log10MassBinWidth)
          massBinMaximumJ=10.0**(logMassBinCenter(j)+0.5d0*log10MassBinWidth)
  
          ! Poisson term.
          if (includePoisson .and. i == j) covariancePoisson(i,j)=massFunctionUse(i)/volume(i)/logMassBinWidth

          ! Halo occupancy covariance.
          if (includeHalo) then
             integrationReset=.true.
             ! Find integration limits for this bin.
             timeMaximum=    Cosmology_Age(                  Expansion_Factor_from_Redshift  (redshiftMinimum))
             timeMinimum=max(                                                                                    &
                  &          Cosmology_Age              (    Expansion_Factor_from_Redshift  (redshiftMaximum)), &
                  &          Time_From_Comoving_Distance(                                                        &
                  &                                      min(                                                    &
                  &                                          Geometry_Survey_Distance_Maximum(massBinCenterI ),  &
                  &                                          Geometry_Survey_Distance_Maximum(massBinCenterJ )   &
                  &                                         )                                                    &
                  &                                     )                                                        &
                  &         )
             ! Integrate over the volume.
             covarianceHalo(i,j)= Integrate(                               &
                  &                         timeMinimum                  , &
                  &                         timeMaximum                  , &
                  &                         Halo_Occupancy_Time_Integrand, &
                  &                         parameterPointer             , &
                  &                         integrandFunction            , &
                  &                         integrationWorkspace         , &
                  &                         toleranceRelative=1.0d-3     , &
                  &                         reset=integrationReset         &
                  &                        )                               &
                  &              *Geometry_Survey_Solid_Angle()            &
                  &              /logMassBinWidth**2                       &
                  &              /volume(i)                                &
                  &              /volume(j)
             call Integrate_Done(integrandFunction,integrationWorkspace)
          end if

          ! Large-scale structure term.
          if (includeLSS) covarianceLSS(i,j)=varianceLSS(i,j)

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
    call Dealloc_Array(logMassBinCenter)
    call Dealloc_Array(volume          )
    call Dealloc_Array(varianceLSS     )

    return
  end subroutine Mass_Function_Covariance_Matrix
  
  function Volume_Integrand(time,parameterPointer) bind(c)
    !% Integral for comoving volume.
    use, intrinsic :: ISO_C_Binding
    use Cosmology_Functions
    implicit none
    real(c_double)        :: Volume_Integrand
    real(c_double), value :: time
    type(c_ptr),    value :: parameterPointer

    Volume_Integrand=Comoving_Volume_Element_Time(time)
    return
  end function Volume_Integrand

  function Mass_Function_Time_Integrand_I(timePrime,parameterPointer) bind(c)
    !% Integral for comoving volume.
    use, intrinsic :: ISO_C_Binding
    use Cosmology_Functions
    use FGSL
    use Numerical_Integration
    implicit none
    real(c_double                  )        :: Mass_Function_Time_Integrand_I
    real(c_double                  ), value :: timePrime
    type(c_ptr                     ), value :: parameterPointer
    type(fgsl_function             )        :: integrandFunction
    type(fgsl_integration_workspace)        :: integrationWorkspace
    logical                                 :: integrationReset
    double precision                        :: massFunction

    time=timePrime
    integrationReset=.true.
    massFunction=Integrate(                            &
         &                 logMassLower              , &
         &                 logMassUpper              , &
         &                 Mass_Function_Integrand_I , &
         &                 parameterPointer          , &
         &                 integrandFunction         , &
         &                 integrationWorkspace      , &
         &                 toleranceRelative=1.0d-3  , &
         &                 reset=integrationReset      &
         &                )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    Mass_Function_Time_Integrand_I=massFunction*Comoving_Volume_Element_Time(time)
    return
  end function Mass_Function_Time_Integrand_I

  function Mass_Function_Integrand_I(logMass,parameterPointer) bind(c)
    !% Integral for mass function.
    use, intrinsic :: ISO_C_Binding
    use Halo_Mass_Function
    use Conditional_Stellar_Mass_Functions
    implicit none
    real(c_double)        :: Mass_Function_Integrand_I
    real(c_double), value :: logMass
    type(c_ptr),    value :: parameterPointer
    double precision      :: mass

    mass=10.0d0**logMass
    Mass_Function_Integrand_I= Halo_Mass_Function_Differential(time,mass)                       &
         &                *                                         mass                        &
         &                *log(10.0d0)                                                          &
         &                *(                                                                    &
         &                   Cumulative_Conditional_Stellar_Mass_Function(mass,massBinMinimumI) &
         &                  -Cumulative_Conditional_Stellar_Mass_Function(mass,massBinMaximumI) &
         &                 )
    return
  end function Mass_Function_Integrand_I

  function LSS_Integrand(timePrime,parameterPointer) bind(c)
    !% Integral for LSS contribution to the covariance matrix.
    use, intrinsic :: ISO_C_Binding
    use Cosmology_Functions
    use FGSL
    use Power_Spectra_Nonlinear
    use Numerical_Interpolation
    implicit none
    real(c_double                  )        :: LSS_Integrand
    real(c_double                  ), value :: timePrime
    type(c_ptr                     ), value :: parameterPointer
    type(fgsl_interp)                       :: interpolationObject
    type(fgsl_interp_accel)                 :: interpolationAccelerator
    logical                                 :: interpolationReset
    logical                                 :: integrationReset
    double precision                        :: biasI,biasJ,powerSpectrum

    ! Copy the time to module scope.
    time=timePrime
    ! Get the bias-mass function product for the I bin.
    interpolationReset=.true.
    biasI=Interpolate(timeBinCount,timeTable,biasTable(:,binI),interpolationObject,interpolationAccelerator,time,reset=interpolationReset)
    call Interpolate_Done(interpolationObject,interpolationAccelerator,interpolationReset)
    ! Get the bias-mass function product for the J bin.
    if (binJ == binI) then
       biasJ=biasI
    else
       interpolationReset=.true.
       biasJ=Interpolate(timeBinCount,timeTable,biasTable(:,binJ),interpolationObject,interpolationAccelerator,time,reset=interpolationReset)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,interpolationReset)
    end if
    ! Get the nonlinear power spectrum for the current wavenumber and time.
    powerSpectrum=Power_Spectrum_Nonlinear(waveNumber,time)
    ! Return the cross-correlation biased power spectrum multiplied by the volume element.
    LSS_Integrand=biasI*biasJ*powerSpectrum*Comoving_Volume_Element_Time(time)
    return
  end function LSS_Integrand

  function Bias_Integrand_I(logMass,parameterPointer) bind(c)
    !% Integral for bias.
    use, intrinsic :: ISO_C_Binding
    use Dark_Matter_Halo_Biases
    implicit none
    real(c_double)        :: Bias_Integrand_I
    real(c_double), value :: logMass
    type(c_ptr),    value :: parameterPointer
    double precision      :: mass

    mass=10.0d0**logMass
    Bias_Integrand_I=Mass_Function_Integrand_I(logMass,parameterPointer)*Dark_Matter_Halo_Bias(mass,time)
  return
  end function Bias_Integrand_I

  function Halo_Occupancy_Time_Integrand(timePrime,parameterPointer) bind(c)
    !% Integral for comoving volume.
    use, intrinsic :: ISO_C_Binding
    use Cosmology_Functions
    use FGSL
    use Numerical_Integration
    implicit none
    real(c_double                  )        :: Halo_Occupancy_Time_Integrand
    real(c_double                  ), value :: timePrime
    type(c_ptr                     ), value :: parameterPointer
    type(fgsl_function             )        :: integrandFunction
    type(fgsl_integration_workspace)        :: integrationWorkspace
    logical                                 :: integrationReset
    double precision                        :: massFunction

    time=timePrime
    integrationReset=.true.
    massFunction=Integrate(                          &
         &                 logMassLower            , &
         &                 logMassUpper            , &
         &                 Halo_Occupancy_Integrand, &
         &                 parameterPointer        , &
         &                 integrandFunction       , &
         &                 integrationWorkspace    , &
         &                 toleranceRelative=1.0d-3, &
         &                 reset=integrationReset    &
         &                )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    Halo_Occupancy_Time_Integrand=massFunction*Comoving_Volume_Element_Time(time)
    return
  end function Halo_Occupancy_Time_Integrand

  function Halo_Occupancy_Integrand(logMass,parameterPointer) bind(c)
    !% Integral for mass function.
    use, intrinsic :: ISO_C_Binding
    use Halo_Mass_Function
    use Conditional_Stellar_Mass_Functions
    implicit none
    real(c_double)        :: Halo_Occupancy_Integrand
    real(c_double), value :: logMass
    type(c_ptr),    value :: parameterPointer
    double precision      :: mass

    mass=10.0d0**logMass
    Halo_Occupancy_Integrand= Halo_Mass_Function_Differential(time,mass)                           &
         &                   *                                     mass                            &
         &                   *log(10.0d0)                                                          &
         &                   *(                                                                    &
         &                      Cumulative_Conditional_Stellar_Mass_Function(mass,massBinMinimumI) &
         &                     -Cumulative_Conditional_Stellar_Mass_Function(mass,massBinMaximumI) &
         &                    )                                                                    &
         &                   *(                                                                    &
         &                      Cumulative_Conditional_Stellar_Mass_Function(mass,massBinMinimumJ) &
         &                     -Cumulative_Conditional_Stellar_Mass_Function(mass,massBinMaximumJ) &
         &                    )
    return
  end function Halo_Occupancy_Integrand

end module Statistics_Mass_Function_Covariance
