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
Implements the gravitational lensing distributions of \cite{takahashi_probability_2011}.
!!}

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Locks                  , only : ompReadWriteLock
  use :: Power_Spectra_Nonlinear, only : powerSpectrumNonlinearClass
  use :: Tables                 , only : table1DGeneric             , table1DLogarithmicLinear

  !![
  <gravitationalLensing name="gravitationalLensingTakahashi2011">
   <description>
    A gravitational lensing distribution class utilizing the fitting functions of \cite{takahashi_probability_2011} to compute
    the effects of gravitational lensing. Specifically, eqn.~11 of \cite{takahashi_probability_2011} is used. The parameters
    $\kappa_\mathrm{empty}$ and $\langle \kappa^2 \rangle$ are computed from the assumed cosmology and non-linear power
    spectrum as described by \cite[][eqns.~5 and 2 respectively]{takahashi_probability_2011}. The parameters, $N_\kappa$,
    $A_\kappa$, and $\omega_\kappa$ of the lensing convergence distribution are determined using the conditions given by
    \cite[][eqn.~9]{takahashi_probability_2011}.
   </description>
  </gravitationalLensing>
  !!]
  type, extends(gravitationalLensingClass) :: gravitationalLensingTakahashi2011
     private
     class           (cosmologyParametersClass   ), pointer :: cosmologyParameters_    => null()
     class           (cosmologyFunctionsClass    ), pointer :: cosmologyFunctions_     => null()
     class           (powerSpectrumNonlinearClass), pointer :: powerSpectrumNonlinear_ => null()
     logical                                                :: tableInitialized        =  .false., cdfInitialized                      =.false.
     !$ type         (ompReadWriteLock           )          :: lock
     type            (table1DGeneric             )          :: convergencePDF
     type            (table1DLogarithmicLinear   )          :: magnificationCDFTable
     double precision                                       :: redshiftPrevious                  , convergenceEmptyBeam                        , &
          &                                                    convergenceVariance               , convergenceScale                            , &
          &                                                    convergenceVarianceScaled         , convergenceDistributionNormalization        , &
          &                                                    aConvergence                      , omegaConvergence                            , &
          &                                                    scaleSourcePrevious
   contains
     !![
     <methods>
       <method description="Construct the gravitational lensing distribution functions for the specified redshift and source scale." method="lensingDistributionConstruct" />
       <method description="Returns the gravitational lensing convergence probability density function at the given convergence." method="convergenceDistribution" />
     </methods>
     !!]
     final     ::                                 takahashi2011Destructor
     procedure :: magnificationPDF             => takahashi2011MagnificationPDF
     procedure :: magnificationCDF             => takahashi2011MagnificationCDF
     procedure :: convergenceDistribution      => takahashi2011ConvergenceDistribution
     procedure :: lensingDistributionConstruct => takahashi2011LensingDistributionConstruct
   end type gravitationalLensingTakahashi2011

  interface gravitationalLensingTakahashi2011
     !!{
     Constructors for the \cite{takahashi_probability_2011} gravitational lensing class.
     !!}
     module procedure takahashi2011ConstructorParameters
     module procedure takahashi2011ConstructorInternal
  end interface gravitationalLensingTakahashi2011

  ! Smallest variance for which calculations are stable.
  double precision, parameter :: convergenceVarianceSmall=1.0d-5

  ! Smallest redshift for which to compute lensing.
  double precision, parameter :: redshiftTiny            =1.5d-2

  ! CDF tabulation range.
  double precision, parameter :: magnificationMinimum    =   1.0d-2
  double precision, parameter :: magnificationMaximum    =1000.0d+0

contains

  function takahashi2011ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{takahashi_probability_2011} gravitational lensing class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (gravitationalLensingTakahashi2011)                :: self
    type (inputParameters                  ), intent(inout) :: parameters
    class(cosmologyParametersClass         ), pointer       :: cosmologyParameters_
    class(cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class(powerSpectrumNonlinearClass      ), pointer       :: powerSpectrumNonlinear_

    !![
    <objectBuilder class="cosmologyParameters"    name="cosmologyParameters_"    source="parameters"/>
    <objectBuilder class="cosmologyFunctions"     name="cosmologyFunctions_"     source="parameters"/>
    <objectBuilder class="powerSpectrumNonlinear" name="powerSpectrumNonlinear_" source="parameters"/>
    !!]
    self=gravitationalLensingTakahashi2011(cosmologyParameters_,cosmologyFunctions_,powerSpectrumNonlinear_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"   />
    <objectDestructor name="cosmologyFunctions_"    />
    <objectDestructor name="powerSpectrumNonlinear_"/>
    !!]
    return
  end function takahashi2011ConstructorParameters

  function takahashi2011ConstructorInternal(cosmologyParameters_,cosmologyFunctions_,powerSpectrumNonlinear_) result(self)
    !!{
    Internal for the \cite{takahashi_probability_2011} gravitational lensing class.
    !!}
    implicit none
    type (gravitationalLensingTakahashi2011)                        :: self
    class(cosmologyParametersClass         ), intent(in   ), target :: cosmologyParameters_
    class(cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
    class(powerSpectrumNonlinearClass      ), intent(in   ), target :: powerSpectrumNonlinear_
    !![
    <constructorAssign variables="*cosmologyParameters_, *cosmologyFunctions_, *powerSpectrumNonlinear_"/>
    !!]

    self   %tableInitialized=.false.
    self   %cdfInitialized  =.false.
    self   %redshiftPrevious=-2.0d0
    !$ self%lock            =ompReadWriteLock()
   return
  end function takahashi2011ConstructorInternal

  subroutine takahashi2011Destructor(self)
    !!{
    Destructor for the \refClass{gravitationalLensingTakahashi2011} gravitational lensing class.
    !!}
    implicit none
    type(gravitationalLensingTakahashi2011), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"   />
    <objectDestructor name="self%cosmologyFunctions_"    />
    <objectDestructor name="self%powerSpectrumNonlinear_"/>
    !!]
    return
  end subroutine takahashi2011Destructor

  double precision function takahashi2011MagnificationPDF(self,magnification,redshift,scaleSource)
    !!{
    Compute the magnification probability density function at the given {\normalfont \ttfamily magnification} and {\normalfont \ttfamily redshift} using the
    \cite{takahashi_probability_2011} formalism.
    !!}
    implicit none
    class           (gravitationalLensingTakahashi2011), intent(inout) :: self
    double precision                                   , intent(in   ) :: magnification, redshift, &
         &                                                                scaleSource

    ! Handle redshift zero case.
    if (redshift <= redshiftTiny) then
       if (magnification == 1.0d0) then
          takahashi2011MagnificationPDF=1.0d0
       else
          takahashi2011MagnificationPDF=0.0d0
       end if
       return
    else
       !$ call self%lock%setRead()
       ! Construct the distribution.
       call self%lensingDistributionConstruct(redshift,scaleSource)
       ! Approximate a δ-function for small redshifts.
       if (self%convergenceVariance < convergenceVarianceSmall) then
          if (magnification == 1.0d0) then
             takahashi2011MagnificationPDF=0.0d0
          else
             takahashi2011MagnificationPDF=1.0d0
          end if
          !$ call self%lock%unsetRead()
          return
       end if
       ! Evaluate the magnification PDF (eqn. 11 of Takahashi et al.).
       takahashi2011MagnificationPDF=takahashi2011MagnificationDistribution(self,magnification)
       !$ call self%lock%unsetRead()
    end if
    return
  end function takahashi2011MagnificationPDF

  double precision function takahashi2011MagnificationCDF(self,magnification,redshift,scaleSource)
    !!{
    Compute the magnification probability density function at the given {\normalfont \ttfamily magnification} and {\normalfont \ttfamily redshift} using the
    \cite{takahashi_probability_2011} formalism.
    !!}
    implicit none
    class           (gravitationalLensingTakahashi2011), intent(inout) :: self
    double precision                                   , intent(in   ) :: magnification, redshift, &
         &                                                                scaleSource

    ! Handle redshift zero case.
    if (redshift <= redshiftTiny) then
       if (magnification < 1.0d0) then
          takahashi2011MagnificationCDF=0.0d0
       else
          takahashi2011MagnificationCDF=1.0d0
       end if
       return
    else
       !$ call self%lock%setRead()
       ! Construct the distribution.
       call self%lensingDistributionConstruct(redshift,scaleSource)
       ! Approximate a δ-function for small redshifts.
       if (self%convergenceVariance < convergenceVarianceSmall) then
          if (magnification < 1.0d0) then
             takahashi2011MagnificationCDF=0.0d0
          else
             takahashi2011MagnificationCDF=1.0d0
          end if
          !$ call self%lock%unsetRead()
          return
       end if
       ! Interpolate in the tabulated cumulative distribution function.
       if      (magnification < magnificationMinimum) then
          takahashi2011MagnificationCDF=0.0d0
       else if (magnification > magnificationMaximum) then
          takahashi2011MagnificationCDF=1.0d0
       else
          takahashi2011MagnificationCDF=self%magnificationCDFTable%interpolate(magnification)/self%magnificationCDFTable%y(-1)
       end if
       !$ call self%lock%unsetRead()
    end if
    return

  contains

    double precision function magnificationPDFIntegrand(magnification)
      !!{
      Integral for the magnification probability distribution function.
      !!}
      implicit none
      double precision, intent(in   ) :: magnification

      magnificationPDFIntegrand=takahashi2011MagnificationDistribution(self,magnification)
      return
    end function magnificationPDFIntegrand

  end function takahashi2011MagnificationCDF

  double precision function takahashi2011MagnificationDistribution(self,magnification)
    !!{
    The gravitational lensing magnification distribution from \cite[][eq.~11]{takahashi_probability_2011}.
    !!}
    implicit none
    class           (gravitationalLensingTakahashi2011), intent(inout) :: self
    double precision                                   , intent(in   ) :: magnification
    double precision                                   , parameter     :: magnificationZeroPoint=3.0d0
    double precision                                   , parameter     :: convergenceZeroPoint  =1.0d0-1.0d0/sqrt(magnificationZeroPoint)
    double precision                                                   :: convergence

    if (magnification <= 0.0d0) then
       takahashi2011MagnificationDistribution=0.0d0
    else
       convergence                           =          &
            & +1.0d0                                    &
            & -1.0d0                                    &
            & /sqrt(magnification)
       takahashi2011MagnificationDistribution=          &
            & +0.5d0                                    &
            & *(                                        &
            &   +1.0d0                                  &
            &   -convergence                            &
            & )**3                                      &
            & *self%convergenceDistribution(convergence)
       if (magnification > 1.0d0)                                  &
            & takahashi2011MagnificationDistribution=              &
            &  +takahashi2011MagnificationDistribution             &
            &  +exp(                                               &
            &       -0.25d0                                        &
            &       /(                                             &
            &         +magnification                               &
            &         -1.0d0                                       &
            &        )**4                                          &
            &      )                                               &
            &  *0.5d0                                              &
            &  *(                                                  &
            &    +1.0d0                                            &
            &    -convergenceZeroPoint                             &
            &   )**3                                               &
            &  *self%convergenceDistribution(convergenceZeroPoint) &
            &  /(                                                  &
            &     magnification                                    &
            &    /magnificationZeroPoint                           &
            &   )**3
    end if
    return
  end function takahashi2011MagnificationDistribution

  double precision function takahashi2011ConvergenceDistribution(self,convergence)
    !!{
    The distribution function for gravitational lensing convergence \citep[][eqn.~8]{takahashi_probability_2011}.
    !!}
    implicit none
    class           (gravitationalLensingTakahashi2011), intent(inout) :: self
    double precision                                   , intent(in   ) :: convergence
    double precision                                                   :: scaledConvergence

    scaledConvergence=convergence/self%convergenceScale
    if (scaledConvergence <= -1.0d0) then
       takahashi2011ConvergenceDistribution=+0.0d0
    else
       takahashi2011ConvergenceDistribution=                                           &
            &                               +self%convergenceDistributionNormalization &
            &                               *exp(                                      &
            &                                    -0.5d0                                &
            &                                    /self%omegaConvergence**2             &
            &                                    *(                                    &
            &                                      +log(                               &
            &                                           +1.0d0                         &
            &                                           +scaledConvergence             &
            &                                          )                               &
            &                                      +self%omegaConvergence**2           &
            &                                      /2.0d0                              &
            &                                     )**2                                 &
            &                                    *(                                    &
            &                                      +1.0d0                              &
            &                                      +self%aConvergence                  &
            &                                      /(                                  &
            &                                        +1.0d0                            &
            &                                        +scaledConvergence                &
            &                                       )                                  &
            &                                     )                                    &
            &                                   )                                      &
            &                               /(                                         &
            &                                 +1.0d0                                   &
            &                                 +scaledConvergence                       &
            &                                )                                         &
            &                               /self%convergenceScale
    end if
    return
  end function takahashi2011ConvergenceDistribution

  subroutine takahashi2011LensingDistributionConstruct(self,redshift,scaleSource)
    !!{
    Construct the lensing distribution function for the \cite{takahashi_probability_2011} formalism.
    !!}
    use :: File_Utilities       , only : Directory_Make              , File_Exists                 , File_Lock            , File_Unlock, &
         &                               lockDescriptor
    use :: Error                , only : Error_Report
    use :: Input_Paths          , only : inputPath                   , pathTypeDataDynamic
    use :: HDF5_Access          , only : hdf5Access
    use :: IO_HDF5              , only : hdf5Object
    use :: Numerical_Comparison , only : Values_Differ
    use :: Numerical_Integration, only : integrator
    use :: Numerical_Ranges     , only : Make_Range                  , rangeTypeLogarithmic
    use :: Root_Finder          , only : rangeExpandMultiplicative   , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    use :: Table_Labels         , only : extrapolationTypeExtrapolate, extrapolationTypeFix
    implicit none
    class           (gravitationalLensingTakahashi2011), intent(inout)               :: self
    double precision                                   , intent(in   )               :: redshift                                        , scaleSource
    double precision                                   , parameter                   :: redshiftZero                          =   0.0d+0
    double precision                                   , parameter                   :: convergenceParametersToleranceAbsolute=   0.0d+0
    double precision                                   , parameter                   :: convergenceParametersToleranceRelative=   1.0d-6
    double precision                                   , parameter                   :: convergenceMaximumFactor              =   1.0d+1
    double precision                                   , parameter                   :: omegaKappaMinimum                     =   1.0d-2
    double precision                                   , parameter                   :: omegaKappaMaximum                     =   1.0d+0
    double precision                                   , parameter                   :: magnificationMinimum                  =   1.0d-3
    double precision                                   , parameter                   :: magnificationMaximum                  =   1.0d+3
    integer                                            , parameter                   :: omegaKappaCount                       =1000
    integer                                            , parameter                   :: cdfMagnificationCount                 =1000
    double precision                                   , allocatable  , dimension(:) :: tableOmegaKappa                           , tableAKappa                   , &
         &                                                                              tableNKappa                               , tableConvergenceVariance
    integer                                                                          :: i
    type            (rootFinder                       )                              :: finder
    type            (integrator                       )                              :: integratorMoment0                         , integratorMoment1             , &
         &                                                                              integratorMoment2                         , integratorEmptyBeamConvergence, &
         &                                                                              integratorConvergenceVariance             , integratorMagnificationPdf
    double precision                                                                 :: distanceComovingSource                    , timeLens                      , &
         &                                                                              convergencePdfMoment0                     , convergencePdfMoment1         , &
         &                                                                              convergencePdfMoment2                     , convergenceMinimum            , &
         &                                                                              convergenceMaximum                        , magnificationPdfMoment0       , &
         &                                                                              magnificationLower                        , magnificationUpper            , &
         &                                                                              cdfPrevious                               , cdf
    type            (hdf5Object                       )                              :: parametersFile
    type            (lockDescriptor                   )                              :: fileLock
    type            (varying_string                   )                              :: fileName

    ! Construct tabulation of convergence distribution function parameters if necessary.
    if (.not.self%tableInitialized) then
       !$ call self%lock%setWrite(haveReadLock=.true.)
       if (.not.self%tableInitialized) then
          ! Determine the parameters, A_κ and ω_κ, of the convergence distribution (eq. 8
          ! of Takahashi et al.). To do this, we use a look-up table of precomputed values.
          ! Check if a precomputed file exists.
          fileName=inputPath(pathTypeDataDynamic)//"largeScaleStructure/"//self%objectType()//'_'//self%hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)//".hdf5"
          call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
          if (File_Exists(fileName)) then
             ! Read the results from file.
             !$ call hdf5Access%set()
             parametersFile=hdf5Object(char(fileName),readOnly=.true.)
             call parametersFile%readDataset("convergenceVariance",tableConvergenceVariance)
             call parametersFile%readDataset(             "NKappa",tableNKappa             )
             call parametersFile%readDataset(             "AKappa",tableAKappa             )
             call parametersFile%readDataset(         "omegaKappa",tableOmegaKappa         )
             !$ call hdf5Access%unset()
          else
             call File_Unlock(fileLock)
             call File_Lock(char(fileName),fileLock,lockIsShared=.false.)
             ! Create a root-finder to solve the parameters.
             finder=rootFinder(&
                  &            rootFunction                 =convergencePdfParameterSolver         , &
                  &            toleranceAbsolute            =convergenceParametersToleranceAbsolute, &
                  &            toleranceRelative            =convergenceParametersToleranceRelative, &
                  &            rangeExpandUpward            =2.0d0                                 , &
                  &            rangeExpandDownward          =0.5d0                                 , &
                  &            rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive         , &
                  &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative         , &
                  &            rangeExpandType              =rangeExpandMultiplicative               &
                  &           )
             ! Construct a range of ω_κ values.
             tableOmegaKappa=Make_Range(                      &
                  &                     omegaKappaMinimum   , &
                  &                     omegaKappaMaximum   , &
                  &                     omegaKappaCount     , &
                  &                     rangeTypeLogarithmic  &
                  &                    )
             allocate(tableAKappa             (omegaKappaCount))
             allocate(tableNKappa             (omegaKappaCount))
             allocate(tableConvergenceVariance(omegaKappaCount))
             ! Iterate over values of ω_κ.
             integratorMoment0=integrator(convergenceDistributionMoment0Integrand,toleranceRelative=1.0d-4)
             integratorMoment1=integrator(convergenceDistributionMoment1Integrand,toleranceRelative=1.0d-4)
             integratorMoment2=integrator(convergenceDistributionMoment2Integrand,toleranceRelative=1.0d-4)
             do i=1,omegaKappaCount
                ! Set ω_κ parameter.
                self%omegaConvergence                    =tableOmegaKappa(i)
                ! Set convergence scale to unity for table building. This implies a minimum convergence of -1.
                self%convergenceScale                    =+1.0d0
                convergenceMinimum                       =-self%convergenceScale
                convergenceMaximum                       =+convergenceMaximumFactor*sqrt(self%omegaConvergence)
                ! Set normalization to unity.
                self%convergenceDistributionNormalization=1.0d0
                ! Solve for a value of A_k that gives consistent first and second moments of the
                ! convergence distribution (eq. 9 of Takahashi et al.).
                self%aConvergence                        =finder%find(rootGuess=1.0d0)
                ! Evaluate zeroth and second moments of the convergence distribution.
                convergencePdfMoment0=integratorMoment0%integrate(convergenceMinimum,convergenceMaximum)
                convergencePdfMoment1=integratorMoment1%integrate(convergenceMinimum,convergenceMaximum)
                convergencePdfMoment2=integratorMoment2%integrate(convergenceMinimum,convergenceMaximum)
                ! Store the A parameter.
                tableAKappa             (i)=self%aConvergence
                ! Compute the normalization of the convergence distribution.
                tableNKappa             (i)=1.0d0                /convergencePdfMoment0
                ! Determine what convergence variance this parameter combination corresponds to.
                tableConvergenceVariance(i)=convergencePdfMoment2/convergencePdfMoment0
                ! Check that the ratio of first to second moments equals -2.
                if (Values_Differ(convergencePdfMoment1/convergencePdfMoment2,-2.0d0,absTol=2.0d-3)) &
                     & call Error_Report('convergence PDF does not satisfy consistency criterion'//{introspection:location})
             end do
             ! Store the results to file.
             call Directory_Make(inputPath(pathTypeDataDynamic)//'largeScaleStructure')
             !$ call hdf5Access%set()
             parametersFile=hdf5Object(char(fileName))
             call parametersFile%writeDataset(tableConvergenceVariance,"convergenceVariance","Dimensionless variance of lensing convergence"     )
             call parametersFile%writeDataset(tableNKappa             ,"NKappa"             ,"Parameter N_kappa from Takahashi et al. (2011)"    )
             call parametersFile%writeDataset(tableAKappa             ,"AKappa"             ,"Parameter A_kappa from Takahashi et al. (2011)"    )
             call parametersFile%writeDataset(tableOmegaKappa         ,"omegaKappa"         ,"Parameter omega_kappa from Takahashi et al. (2011)")
             !$ call hdf5Access%unset()
          end if
          call File_Unlock(fileLock)
          ! Create a table. We fix the extrapolation for large scaled variances. In these cases, the convergence distribution is
          ! extremely narrow (effectively a delta-function) so it does not matter much what we do. If we do allow extrapolation, the
          ! values obtained for the parameter result in convergence distributions which have secondary peaks at κ≅2, which results in
          ! unrealistic magnification distributions.
          call self%convergencePDF%create  (tableConvergenceVariance,tableCount=3,extrapolationType=[extrapolationTypeExtrapolate,extrapolationTypeFix])
          call self%convergencePDF%populate(tableNKappa             ,           1                                                                      )
          call self%convergencePDF%populate(tableAKappa             ,           2                                                                      )
          call self%convergencePDF%populate(tableOmegaKappa         ,           3                                                                      )
          deallocate(tableConvergenceVariance)
          deallocate(tableNKappa             )
          deallocate(tableAKappa             )
          deallocate(tableOmegaKappa         )
          ! Record that the table is now initialized.
          self%tableInitialized=.true.
       end if
       !$ call self%lock%unsetWrite(haveReadLock=.true.)
    end if
    ! Check if redshift has changed since previous call.
    if (redshift /= self%redshiftPrevious .or. scaleSource /= self%scaleSourcePrevious) then
       !$ call self%lock%setWrite(haveReadLock=.true.)
       if (redshift /= self%redshiftPrevious .or. scaleSource /= self%scaleSourcePrevious) then
          ! Update convergences.
          self%redshiftPrevious   =redshift
          self%scaleSourcePrevious=scaleSource
          ! Find the comoving distance to the source.
          distanceComovingSource=self%cosmologyFunctions_%distanceComoving           (           &
               &                 self%cosmologyFunctions_%cosmicTime                  (          &
               &                 self%cosmologyFunctions_%expansionFactorFromRedshift  (         &
               &                                                                        redshift &
               &                                                                       )         &
               &                                                                      )          &
               &                                                                     )
          ! Find the convergence of an empty beam.
          integratorEmptyBeamConvergence=integrator                              (emptyBeamConvergenceIntegrand,toleranceRelative=1.0d-3)
          self%convergenceEmptyBeam     =integratorEmptyBeamConvergence%integrate(redshiftZero                 ,redshift                )
          ! Find the variance of the convergence.
          integratorConvergenceVariance=integrator                               (convergenceVarianceIntegrand ,toleranceRelative=1.0d-3)
          self%convergenceVariance     =integratorConvergenceVariance%integrate  (redshiftZero                 ,redshift                )
          ! Determine the parameters of the convergence distribution.
          self%convergenceScale                    =abs(self%convergenceEmptyBeam)
          self%convergenceVarianceScaled           =self%convergenceVariance/self%convergenceScale**2
          self%convergenceDistributionNormalization=self%convergencePDF%interpolate(self%convergenceVarianceScaled,1)
          self%aConvergence                        =self%convergencePDF%interpolate(self%convergenceVarianceScaled,2)
          self%omegaConvergence                    =self%convergencePDF%interpolate(self%convergenceVarianceScaled,3)
          ! Integrate the modified magnification distribution in order to find the normalization.
          integratorMagnificationPdf=integrator                          (magnificationPDFIntegrand,toleranceRelative   =1.0d-6)
          magnificationPdfMoment0   =integratorMagnificationPdf%integrate(magnificationMinimum     ,magnificationMaximum       )
          self%convergenceDistributionNormalization         &
               & =self%convergenceDistributionNormalization &
               & /magnificationPdfMoment0
          ! Tabulate the cumulative distribution function if table does not yet exist.
          if (self%cdfInitialized) call self%magnificationCDFTable%destroy()
          call self%magnificationCDFTable%create(magnificationMinimum,magnificationMaximum,cdfMagnificationCount,tableCount=1)
          do i=1,cdfMagnificationCount
             if (i == 1 ) then
                magnificationLower=0.0d0
                cdfPrevious       =0.0d0
             else
                magnificationLower=self%magnificationCDFTable%x(i-1)
                cdfPrevious       =self%magnificationCDFTable%y(i-1)
             end if
             magnificationUpper   =self%magnificationCDFTable%x(i  )
             cdf=integratorMagnificationPdf%integrate(magnificationLower,magnificationUpper)
             call self%magnificationCDFTable%populate(cdf+cdfPrevious,i)
          end do
          self%cdfInitialized=.true.
       end if
       !$ call self%lock%unsetWrite(haveReadLock=.true.)
    end if
    return

  contains

    double precision function magnificationPDFIntegrand(magnification)
      !!{
      Integral for the magnification probability distribution function.
      !!}
      implicit none
      double precision, intent(in   ) :: magnification

      magnificationPDFIntegrand=takahashi2011MagnificationDistribution(self,magnification)
      return
    end function magnificationPDFIntegrand

    double precision function convergencePdfParameterSolver(a)
      !!{
      Root function used in finding equivalent circular orbits.
      !!}
      implicit none
      double precision            , intent(in   ) :: a
      type            (integrator)                :: integratorConvergenceDistributionMoment1, integratorConvergenceDistributionMoment2

      self%aConvergence    =a
      integratorConvergenceDistributionMoment1=integrator                                        (convergenceDistributionMoment1Integrand,toleranceRelative =1.0d-3)
      integratorConvergenceDistributionMoment2=integrator                                        (convergenceDistributionMoment2Integrand,toleranceRelative =1.0d-3)
      convergencePdfMoment1                   =integratorConvergenceDistributionMoment1%integrate(convergenceMinimum                     ,convergenceMaximum       )
      convergencePdfMoment2                   =integratorConvergenceDistributionMoment2%integrate(convergenceMinimum                     ,convergenceMaximum       )
      if (convergencePdfMoment2 <= 0.0d0) then
         convergencePdfParameterSolver=0.0d0
      else
         convergencePdfParameterSolver=2.0d0+convergencePdfMoment1/convergencePdfMoment2*self%convergenceScale
      end if
      return
    end function convergencePdfParameterSolver

    double precision function emptyBeamConvergenceIntegrand(redshiftLens)
      !!{
      Integral for gravitational lensing convergence in an empty beam.
      !!}
      use :: Numerical_Constants_Physical, only : speedLight
      use :: Numerical_Constants_Prefixes, only : kilo
      implicit none
      double precision, intent(in   ) :: redshiftLens
      double precision                :: distanceComovingLens

      ! Find cosmic time at this redshift.
      timeLens            =self%cosmologyFunctions_%cosmicTime                 (              &
           &               self%cosmologyFunctions_%expansionFactorFromRedshift (             &
           &                                                                     redshiftLens &
           &                                                                    )             &
           &                                                                   )
      ! Find comoving distance to the lens.
      distanceComovingLens=self%cosmologyFunctions_%distanceComoving           (              &
           &                                                                     timeLens     &
           &                                                                   )
      ! Evaluate the integrand.
      emptyBeamConvergenceIntegrand=                                                               &
           &                        -1.5d0                                                         &
           &                        /speedLight                                                    &
           &                        *kilo                                                          &
           &                        *self%cosmologyParameters_%HubbleConstant        (        )**2 &
           &                        /self%cosmologyFunctions_ %hubbleParameterEpochal(timeLens)    &
           &                        *self%cosmologyParameters_%OmegaMatter           (        )    &
           &                        *(                                                             &
           &                          +1.0d0                                                       &
           &                          +redshiftLens                                                &
           &                         )                                                             &
           &                        *  distanceComovingLens                                        &
           &                        *(                                                             &
           &                          +distanceComovingSource                                      &
           &                          -distanceComovingLens                                        &
           &                        )                                                              &
           &                        /  distanceComovingSource
      return
    end function emptyBeamConvergenceIntegrand

    double precision function convergenceVarianceIntegrand(redshiftLens)
      !!{
      Integral for variance in the gravitational lensing convergence.
      !!}
      use :: Numerical_Constants_Math    , only : Pi
      use :: Numerical_Constants_Physical, only : speedLight
      use :: Numerical_Constants_Prefixes, only : kilo
      implicit none
      double precision            , intent(in   ) :: redshiftLens
      double precision            , parameter     :: wavenumberDynamicRange =1.0d-6
      double precision                            :: distanceComovingLens          , logWavenumberMaximum, &
           &                                         lensingPower                  , logWavenumberMinimum
      type            (integrator)                :: integratorPowerSpectrum

      ! Find cosmic time at this redshift.
      timeLens            =self%cosmologyFunctions_%cosmicTime                 (              &
           &               self%cosmologyFunctions_%expansionFactorFromRedshift (             &
           &                                                                     redshiftLens &
           &                                                                    )             &
           &                                                                   )
      ! Find comoving distance to the lens.
      distanceComovingLens=self%cosmologyFunctions_%distanceComoving           (              &
           &                                                                     timeLens     &
           &                                                                   )
      ! Integrate over the power spectrum.
      logWavenumberMaximum   =                    -log(scaleSource           )
      logWavenumberMinimum   =logWavenumberMaximum+log(wavenumberDynamicRange)
      integratorPowerSpectrum=integrator                       (convergenceVariancePowerSpectrumIntegrand,toleranceRelative   =1.0d-3)
      lensingPower           =integratorPowerSpectrum%integrate(logWavenumberMinimum                     ,logWavenumberMaximum       )
      ! Evaluate the integrand.
      convergenceVarianceIntegrand =                                                               &
           &                        +9.0d0                                                         &
           &                        /8.0d0                                                         &
           &                        /Pi                                                            &
           &                        /speedLight                                                **3 &
           &                        *kilo                                                      **3 &
           &                        *self%cosmologyParameters_%HubbleConstant        (        )**4 &
           &                        /self%cosmologyFunctions_ %hubbleParameterEpochal(timeLens)    &
           &                        *self%cosmologyParameters_%OmegaMatter           (        )**2 &
           &                        *(                                                             &
           &                          +1.0d0                                                       &
           &                          +redshiftLens                                                &
           &                         )                                                         **2 &
           &                        *(                                                             &
           &                          +  distanceComovingLens                                      &
           &                          *(                                                           &
           &                            +distanceComovingSource                                    &
           &                            -distanceComovingLens                                      &
           &                           )                                                           &
           &                          / distanceComovingSource                                     &
           &                         )                                                         **2 &
           &                        *lensingPower
      return
    end function convergenceVarianceIntegrand

    double precision function convergenceVariancePowerSpectrumIntegrand(logWavenumber)
      !!{
      Integral over power spectrum used in computing the variance in the gravitational lensing convergence.
      !!}
      implicit none
      double precision, intent(in   ) :: logWavenumber
      double precision                :: wavenumber

      ! Compute integrand.
      wavenumber=exp(logWavenumber)
      convergenceVariancePowerSpectrumIntegrand=wavenumber**2*self%powerSpectrumNonLinear_%value(wavenumber,timeLens)
      return
    end function convergenceVariancePowerSpectrumIntegrand

    double precision function convergenceDistributionMoment0Integrand(scaledConvergence)
      !!{
      Integral over scaled convergence distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: scaledConvergence

      convergenceDistributionMoment0Integrand=self%convergenceDistribution(scaledConvergence)
      return
    end function convergenceDistributionMoment0Integrand

    double precision function convergenceDistributionMoment1Integrand(scaledConvergence)
      !!{
      Integral of first moment over scaled convergence distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: scaledConvergence

      convergenceDistributionMoment1Integrand=scaledConvergence*self%convergenceDistribution(scaledConvergence)
      return
    end function convergenceDistributionMoment1Integrand

    double precision function convergenceDistributionMoment2Integrand(scaledConvergence)
      !!{
      Integral of second moment over scaled convergence distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: scaledConvergence

      convergenceDistributionMoment2Integrand=scaledConvergence**2*self%convergenceDistribution(scaledConvergence)
      return
    end function convergenceDistributionMoment2Integrand

  end subroutine takahashi2011LensingDistributionConstruct
