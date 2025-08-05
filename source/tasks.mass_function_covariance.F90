!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  use :: Conditional_Mass_Functions, only : conditionalMassFunction, conditionalMassFunctionClass
  use :: Cosmology_Functions       , only : cosmologyFunctions     , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Biases   , only : darkMatterHaloBias     , darkMatterHaloBiasClass
  use :: Geometry_Surveys          , only : surveyGeometry         , surveyGeometryClass
  use :: Halo_Mass_Functions       , only : haloMassFunction       , haloMassFunctionClass
  use :: Power_Spectra_Nonlinear   , only : powerSpectrumNonlinear , powerSpectrumNonlinearClass

  !![
  <task name="taskMassFunctionCovariance">
   <description>
    A task class which computes and stores covariance matrices for mass functions. In general, for constraints corresponding to
    mass functions (whether stellar mass or HI mass), the covariance matrix of the observational data is determined using the
    analytic model of \cite{smith_how_2012}. This requires knowledge of both the survey geometry (angular mask and radial extent
    as a function of mass) and of the \gls{hod} of the observed galaxies.
  
    Details of the survey geometry and depth are given for each individual constraints. Computing the large-scale structure
    contribution to the covariance function requires integration of the non-linear matter power spectrum over the Fourier
    transform of the survey window function. We use the method of \cite{peacock_non-linear_1996} to determine the non-linear
    matter power spectrum, because of its simplicity and speed. We have checked that using a more accurate non-linear matter
    power spectrum (e.g. \citealt{lawrence_coyote_2010}) makes negligible difference to our results. If the angular power
    spectrum of the survey mask is available\footnote{Typically if the survey geometry is defined by \protect\gls{mangle}
    polygons, allowing the angular power spectrum to be found using the \protect\gls{mangle} {\normalfont \ttfamily harmonize}
    utility.}, this is used to compute the relation
    \begin{equation}
     \sigma^2(M_\mu,M_\nu) = {2 \over \pi V_\mu V_\nu}\int_0^\infty \mathrm{d} k\, k^{-4} P(k) \sum_i \sum_j \sum_{\ell=0}^\infty (2\ell+1) C^{ij}_\ell R^i_{\ell}(kr_{\mu 0},kr_{\mu 1}) R^j_{\ell}(kr_{\nu 0},kr_{\nu 1}),
    \end{equation}
    where $(2\ell+1) C^{ij}_\ell = \sum_{m=-\ell}^{+\ell} \Psi^i_{\ell m} \Psi^{j*}_{\ell m}$, $\Psi^i_{\ell m}$ are the
    spherical harmonic coefficients of the $i^\mathrm{th}$ field of the survey, $V$ is the maximum distance to which a galaxy of
    mass $M$ can be seen, $P(k)$ is the nonlinear power spectrum and
    \begin{equation}
     R_{\ell}(x_0,x_1) \equiv \int_{x_0}^{x_1} x^2 j_\ell(x) \mathrm{d}x = \sqrt{\pi} 2^{-2-\ell} \Gamma\left({1\over 2}[3+\ell]\right) \left[ x^{3+\ell} \tensor*[_1]{\stackrel{\sim}{F}}{_2} \left({1\over 2}[3+\ell]; \ell+{3\over 2},{1\over 2}(5+\ell);-{x^2\over 4}\right)\right]_{x_0}^{x_1},
   \end{equation}
   where $\tensor*[_1]{\stackrel{\sim}{F}}{_2}$ is the regularized generalized hypergeometric function. In other cases, where
   the angular power spectrum is not available, the survey geometry is realized on a grid which is when Fourier transformed to
   obtain the appropriate window function.
  
   To find a suitable \gls{hod} to describe the observed galaxies we adopt the model of \cite{behroozi_comprehensive_2010}. This
   is an 11 parameter model which describes separately the numbers of satellite and central galaxies occupying a halo of given
   mass---the reader is referred to \cite{behroozi_comprehensive_2010} for a complete description of the functional form of this
   parametric \gls{hod}. An \gls{mcmc} approach is used to to constrain the \gls{hod} parameters to fit the observational
   data. We use a likelihood
   \begin{equation}
    \ln \mathcal{L} = -{1\over 2} \Delta\cdot \mathcal{C}^{-1}\cdot \Delta^\mathrm{T} - {N \over 2} \ln(2\pi) - {\ln |\mathcal{C}| \over 2},
   \end{equation}
   where $N$ is the number of bins in the mass function, $\mathcal{C}$ is the covariance matrix of the observed mass function,
   and $\Delta_i = \phi_i^\mathrm{(HOD)} - \phi_i^\mathrm{(observed)}$. Of course, it is precisely this covariance matrix,
   $\mathcal{C}$, that we are trying to compute. We therefore adopt an iterative approach as follows:
   \begin{enumerate}
    \item make an initial estimate of the covariance matrix, assuming that only Poisson errors contribute (the covariance
    matrix is therefore diagonal, and the terms are easily computed from the measured mass function and the survey volume as a
    function of stellar mass);
    \item find the maximum likelihood parameters of the \gls{hod} given the observed mass function and the current estimate of
    the covariance matrix;
    \item using this \gls{hod} and the framework of \cite{smith_how_2012}, compute a new estimate of the covariance matrix,
    including all three contributions;
    \item repeat steps 2 and 3 until convergence in the covariance matrix is achieved.       
   \end{enumerate}
  
   In practice we find that this procedure often leads to an \gls{hod} and covariance matrix which oscillate between two states
   in successive iterations. The differences in the covariance matrix are relatively small however, so we choose to
   conservatively adopt the covariance matrix with the larger values.
   </description>
  </task>
  !!]
  type, extends(taskClass) :: taskMassFunctionCovariance
     !!{
     Implementation of a task which computes and stores covariance matrices for mass functions.
     !!}
     private
     class           (cosmologyFunctionsClass     ), pointer                     :: cosmologyFunctions_      => null()
     class           (surveyGeometryClass         ), pointer                     :: surveyGeometry_          => null()
     class           (powerSpectrumNonlinearClass ), pointer                     :: powerSpectrumNonlinear_  => null()
     class           (darkMatterHaloBiasClass     ), pointer                     :: darkMatterHaloBias_      => null()
     class           (conditionalMassFunctionClass), pointer                     :: conditionalMassFunction_ => null()
     class           (haloMassFunctionClass       ), pointer                     :: haloMassFunction_        => null()
     integer                                                                     :: countMassBins                     , sizeGridFFT
     double precision                                                            :: massMinimum                       , massMaximum          , &
          &                                                                         massHaloMinimum                   , massHaloMaximum      , &
          &                                                                         completenessErrorObserved
     logical                                                                     :: includePoisson                    , includeHalo          , &
          &                                                                         includeLSS
     type            (varying_string              )                              :: massFunctionFileName
     ! State used in integrations etc.
     double precision                                                            :: time
     double precision                                                            :: waveNumberGlobal
     double precision                                                            :: logMassLower                      , logMassUpper
     double precision                              , dimension(:  ), allocatable :: log10MassBinWidth                 , logMassBinWidth
     integer                                                                     :: binI                              , binJ                 , &
          &                                                                         lssBin
     double precision                                                            :: massBinCenterI                    , massBinMinimumI      , &
          &                                                                         massBinMaximumI
     double precision                                                            :: massBinCenterJ                    , massBinMinimumJ      , &
          &                                                                         massBinMaximumJ
     integer                                                                     :: countTimeBins
     double precision                              , dimension(:  ), allocatable :: timeTable
     double precision                              , dimension(:,:), allocatable :: biasTable
     double precision                                                            :: surveyRedshiftMinimum             , surveyRedshiftMaximum
     double precision                              , dimension(:  ), allocatable :: volumeNormalizationI              , volumeNormalizationJ , &
          &                                                                         timeMinimumI                      , timeMinimumJ         , &
          &                                                                         timeMaximumI                      , timeMaximumJ         , &
          &                                                                         logMassBinCenter
   contains
     final     ::                       massFunctionCovarianceDestructor
     procedure :: perform            => massFunctionCovariancePerform
     procedure :: requiresOutputFile => massFunctionCovarianceRequiresOutputFile
  end type taskMassFunctionCovariance

  interface taskMassFunctionCovariance
     !!{
     Constructors for the \refClass{taskMassFunctionCovariance} task.
     !!}
     module procedure massFunctionCovarianceConstructorParameters
     module procedure massFunctionCovarianceConstructorInternal
  end interface taskMassFunctionCovariance

  ! Module-scope pointer to self used in integrands.
  class(taskMassFunctionCovariance), pointer :: self_, selfCopy
  !$omp threadprivate(self_,selfCopy)

contains

  function massFunctionCovarianceConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskMassFunctionCovariance} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (taskMassFunctionCovariance  )                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass     ), pointer       :: cosmologyFunctions_
    class           (surveyGeometryClass         ), pointer       :: surveyGeometry_
    class           (powerSpectrumNonlinearClass ), pointer       :: powerSpectrumNonlinear_
    class           (darkMatterHaloBiasClass     ), pointer       :: darkMatterHaloBias_
    class           (conditionalMassFunctionClass), pointer       :: conditionalMassFunction_
    class           (haloMassFunctionClass       ), pointer       :: haloMassFunction_
    double precision                                              :: surveyRedshiftMinimum   , surveyRedshiftMaximum, &
         &                                                           massMinimum             , massMaximum          , &
         &                                                           massHaloMinimum         , massHaloMaximum
    integer                                                       :: countMassBins           , sizeGridFFT
    logical                                                       :: includePoisson          , includeHalo          , &
         &                                                           includeLSS
    type            (varying_string              )                :: massFunctionFileName

    !![
    <inputParameter>
      <name>massFunctionFileName</name>
      <description>The name of the file to which the covariance matrix should be written.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>surveyRedshiftMinimum</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum redshift at which calculations of the mass function covariance should be carried out.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>surveyRedshiftMaximum</name>
      <defaultValue>0.1d0</defaultValue>
      <description>The maximum redshift at which calculations of the mass function covariance should be carried out.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>countMassBins</name>
      <defaultValue>10</defaultValue>
      <description>The number of bins in the mass function for covariance calculations.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMinimum</name>
      <defaultValue>1.0d08</defaultValue>
      <description>The minimum mass in the mass function for covariance calculations.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <defaultValue>1.0d13</defaultValue>
      <description>The maximum mass in the mass function for covariance calculations.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includePoisson</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether or not to include the Poisson contribution to mass function covariance matrices.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeHalo</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether or not to include the halo contribution to mass function covariance matrices.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeLSS</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether or not to include the large-scale structure contribution to mass function covariance matrices.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massHaloMinimum</name>
      <defaultValue>1.0d10</defaultValue>
      <description>The minimum halo mass to use when computing mass function covariance matrices.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massHaloMaximum</name>
      <defaultValue>1.0d15</defaultValue>
      <description>The minimum halo mass to use when computing mass function covariance matrices.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>sizeGridFFT</name>
      <defaultValue>64</defaultValue>
      <description>The size of the FFT grid to use in computing window functions for mass function covariance matrices.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"      name="cosmologyFunctions_"      source="parameters"/>
    <objectBuilder class="surveyGeometry"          name="surveyGeometry_"          source="parameters"/>
    <objectBuilder class="powerSpectrumNonlinear"  name="powerSpectrumNonlinear_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloBias"      name="darkMatterHaloBias_"      source="parameters"/>
    <objectBuilder class="conditionalMassFunction" name="conditionalMassFunction_" source="parameters"/>
    <objectBuilder class="haloMassFunction"        name="haloMassFunction_"        source="parameters"/>
    !!]
    self=taskMassFunctionCovariance(massFunctionFileName,surveyredshiftMinimum, surveyRedshiftMaximum, massMinimum, massMaximum, massHaloMinimum, massHaloMaximum,countMassBins,sizeGridFFT,includePoisson, includeHalo, includeLSS, cosmologyFunctions_,surveyGeometry_,powerSpectrumNonlinear_,darkMatterHaloBias_,conditionalMassFunction_,haloMassFunction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"     />
    <objectDestructor name="surveyGeometry_"         />
    <objectDestructor name="powerSpectrumNonlinear_" />
    <objectDestructor name="darkMatterHaloBias_"     />
    <objectDestructor name="conditionalMassFunction_"/>
    <objectDestructor name="haloMassFunction_"       />
    !!]
    return
  end function massFunctionCovarianceConstructorParameters

  function massFunctionCovarianceConstructorInternal(massFunctionFileName,surveyRedshiftMinimum,surveyRedshiftMaximum,massMinimum,massMaximum,massHaloMinimum,massHaloMaximum,countMassBins,sizeGridFFT,includePoisson,includeHalo,includeLSS,cosmologyFunctions_,surveyGeometry_,powerSpectrumNonlinear_,darkMatterHaloBias_,conditionalMassFunction_,haloMassFunction_) result(self)
    !!{
    Internal constructor for the \refClass{taskMassFunctionCovariance} task class.
    !!}
    implicit none
    type            (taskMassFunctionCovariance  )                        :: self
    class           (cosmologyFunctionsClass     ), intent(in   ), target :: cosmologyFunctions_
    class           (surveyGeometryClass         ), intent(in   ), target :: surveyGeometry_
    class           (powerSpectrumNonlinearClass ), intent(in   ), target :: powerSpectrumNonlinear_
    class           (darkMatterHaloBiasClass     ), intent(in   ), target :: darkMatterHaloBias_
    class           (conditionalMassFunctionClass), intent(in   ), target :: conditionalMassFunction_
    class           (haloMassFunctionClass       ), intent(in   ), target :: haloMassFunction_
    type            (varying_string              ), intent(in   )         :: massFunctionFileName
    double precision                              , intent(in   )         :: surveyRedshiftMinimum   , surveyRedshiftMaximum, &
         &                                                                   massMinimum             , massMaximum          , &
         &                                                                   massHaloMinimum         , massHaloMaximum
    integer                                       , intent(in   )         :: countMassBins           , sizeGridFFT
    logical                                       , intent(in   )         :: includePoisson          , includeHalo          , &
         &                                                                   includeLSS
    !![
    <constructorAssign variables="massFunctionFileName, surveyRedshiftMinimum, surveyRedshiftMaximum, massMinimum, massMaximum, massHaloMinimum, massHaloMaximum,countMassBins,sizeGridFFT,includePoisson, includeHalo, includeLSS, *cosmologyFunctions_,*surveyGeometry_,*powerSpectrumNonlinear_,*darkMatterHaloBias_,*conditionalMassFunction_,*haloMassFunction_"/>
    !!]

    return
  end function massFunctionCovarianceConstructorInternal

  subroutine massFunctionCovarianceDestructor(self)
    !!{
    Destructor for the \refClass{taskMassFunctionCovariance} task class.
    !!}
    implicit none
    type(taskMassFunctionCovariance), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"     />
    <objectDestructor name="self%surveyGeometry_"         />
    <objectDestructor name="self%powerSpectrumNonlinear_" />
    <objectDestructor name="self%darkMatterHaloBias_"     />
    <objectDestructor name="self%conditionalMassFunction_"/>
    <objectDestructor name="self%haloMassFunction_"       />
    !!]
    return
  end subroutine massFunctionCovarianceDestructor

  subroutine massFunctionCovariancePerform(self,status)
    !!{
    Compute and output the halo mass function.
    !!}
    use :: Display                         , only : displayIndent          , displayUnindent
    use :: Error                , only : Error_Report, errorStatusSuccess
    use :: IO_HDF5                         , only : hdf5Object
    use :: Numerical_Constants_Astronomical, only : massSolar              , megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Integration           , only : integrator
    use :: Numerical_Ranges                , only : Make_Range             , rangeTypeLinear   , rangeTypeLogarithmic
    implicit none
    class           (taskMassFunctionCovariance), intent(inout), target                   :: self
    integer                                     , intent(  out), optional                 :: status
    double precision                            , allocatable  , dimension(:    )         :: mass                        , massWidthObserved          , &
         &                                                                                   completenessObserved        , numberObserved             , &
         &                                                                                   massObserved
    double precision                            , allocatable  , dimension(:    ), target :: massFunction                , massFunctionObserved
    double precision                            , allocatable  , dimension(:,:  )         :: covariance                  , covariancePoisson          , &
         &                                                                                   covarianceHalo              , covarianceLSS              , &
         &                                                                                   correlation
    double precision                            , allocatable  , dimension(:,:  )         :: varianceLSS                 , volume
    double precision                            , pointer      , dimension(:    )         :: massFunctionUse
    double precision                            , parameter                               :: timePointsPerDecade =100
    double precision                            , parameter                               :: massFunctionMinimum =1.0d-50
    logical                                                                               :: useCompleteness             ,  useNumber
    integer                                                                               :: i                           , j                          , &
         &                                                                                   iTime                       , iField                     , &
         &                                                                                   countFields
    double precision                                                                      :: logMassMinimum              , logMassMaximum             , &
         &                                                                                   normalization               , volumeNormalization        , &
         &                                                                                   timeMinimum                 , timeMaximum
    double precision                                                                      :: binCompleteness
    type            (integrator                )                                          :: integratorVolume            , integratorTimeI            , &
         &                                                                                   integratorBiasI             , integratorHaloOccupancyTime
    type            (hdf5Object                )                                          :: massFunctionFile            , dataset

    call displayIndent('Begin task: mass function covariance' )
    ! Set a module-scope pointer to our self.
    self_    => self
    selfCopy => self
    ! Open the mass function file.
    massFunctionFile=hdf5Object(char(self%massFunctionFileName),overWrite=.true.)
    ! Read the observed mass function if available.
    self%completenessErrorObserved=0.0d0
    if (massFunctionFile%hasDataset  ("massFunctionObserved")) call massFunctionFile%readDataset  ("massFunctionObserved",     massFunctionObserved     )
    if (massFunctionFile%hasDataset  ("completenessObserved")) call massFunctionFile%readDataset  ("completenessObserved",     completenessObserved     )
    if (massFunctionFile%hasDataset  ("numberObserved"      )) call massFunctionFile%readDataset  ("numberObserved"      ,     numberObserved           )
    if (massFunctionFile%hasAttribute("completenessError"   )) call massFunctionFile%readAttribute("completenessError"   ,self%completenessErrorObserved)
    if (massFunctionFile%hasDataset  ("massObserved"        )) call massFunctionFile%readDataset  ("massObserved"        ,     massObserved             )
    if (massFunctionFile%hasDataset  ("massWidthObserved"   )) call massFunctionFile%readDataset  ("massWidthObserved"   ,     massWidthObserved        )
    ! Determine number of times over which to tabulate bias.
    timeMaximum =self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%surveyRedshiftMinimum))
    timeMinimum =self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%surveyRedshiftMaximum))
    self%countTimeBins=int(log10(timeMaximum/timeMinimum)*dble(timePointsPerDecade))+1
    ! Allocate arrays.
    countFields=self%surveyGeometry_%fieldCount()
    allocate(     mass             (self%countMassBins                   ))
    allocate(self%logMassBinCenter (self%countMassBins                   ))
    allocate(self%log10MassBinWidth(self%countMassBins                   ))
    allocate(self%logMassBinWidth  (self%countMassBins                   ))
    allocate(     massFunction     (self%countMassBins                   ))
    allocate(     volume           (self%countMassBins,     countFields  ))
    allocate(     covariance       (self%countMassBins,self%countMassBins))
    allocate(     covariancePoisson(self%countMassBins,self%countMassBins))
    allocate(     covarianceHalo   (self%countMassBins,self%countMassBins))
    allocate(     covarianceLSS    (self%countMassBins,self%countMassBins))
    allocate(     correlation      (self%countMassBins,self%countMassBins))
    allocate(     varianceLSS      (self%countMassBins,self%countMassBins))
    allocate(self%timeTable        (self%countTimeBins                   ))
    allocate(self%biasTable        (self%countTimeBins,self%countMassBins))
    ! Build integrators.
    integratorVolume           =integrator(massFunctionCovarianceVolumeIntegrand           ,toleranceRelative=1.0d-3)
    integratorTimeI            =integrator(massFunctionCovarianceMassFunctionTimeIntegrandI,toleranceRelative=1.0d-3)
    integratorBiasI            =integrator(massFunctionCovarianceBiasIntegrandI            ,toleranceRelative=1.0d-2)
    integratorHaloOccupancyTime=integrator(massFunctionCovarianceHaloOccupancyTimeInterand ,toleranceRelative=1.0d-3)
    ! Create time bins.
    self%timeTable=Make_Range(timeMinimum,timeMaximum,self%countTimeBins,rangeType=rangeTypeLogarithmic)
    ! Create mass bins.
    if (allocated(massWidthObserved)) then
       mass                  =      massObserved
       self%log10MassBinWidth=log10(massWidthObserved)
       self%logMassBinWidth  =log  (10.0d0           )*self%log10MassBinWidth
       self%logMassBinCenter =log10(mass)
    else
       logMassMinimum        =log10(self%massMinimum)
       logMassMaximum        =log10(self%massMaximum)
       self%logMassBinCenter =Make_Range(logMassMinimum,logMassMaximum,self%countMassBins,rangeType=rangeTypeLinear)
       self%log10MassBinWidth=+self%logMassBinCenter(2) &
            &                 -self%logMassBinCenter(1)
       self%logMassBinWidth  =log(10.0d0)*self%log10MassBinWidth
       mass                  =10.0d0**self%logMassBinCenter
    end if
    ! Halo mass limits for integrations.
    self%logMassLower        =log10(self%massHaloMinimum)
    self%logMassUpper        =log10(self%massHaloMaximum)
    ! Determine which mass function to use.
    if (allocated(massFunctionObserved)) then
       if (size(massFunctionObserved) /= self%countMassBins) call Error_Report('observed mass function has incorrect number of bins'//{introspection:location})
       massFunctionUse => massFunctionObserved
    else
       massFunctionUse => massFunction
    end if
    ! Determine if completeness and/or number is available.
    useCompleteness=allocated(completenessObserved)
    if (useCompleteness .and. size(completenessObserved) /= self%countMassBins) &
         & call Error_Report('observed completeness has incorrect number of bins'//{introspection:location})
    useNumber      =allocated(      numberObserved)
    if (useNumber       .and. size(      numberObserved) /= self%countMassBins) &
         & call Error_Report('observed number has incorrect number of bins'      //{introspection:location})
    ! Compute the mass function and bias averaged over each bin.
    massFunction=0.0d0
    do i=1,self%countMassBins
       ! Find limits on mass for this bin.
       self%massBinCenterI =10.0d0** self%logMassBinCenter(i)
       self%massBinMinimumI=10.0d0**(self%logMassBinCenter(i)-0.5d0*self%log10MassBinWidth(i))
       self%massBinMaximumI=10.0d0**(self%logMassBinCenter(i)+0.5d0*self%log10MassBinWidth(i))
       ! Iterate over fields.
       volumeNormalization=0.0d0
       do iField=1,countFields
          ! Find integration limits for this bin.
          timeMaximum=    self%cosmologyFunctions_%cosmicTime            (self%cosmologyFunctions_%expansionFactorFromRedshift(self%surveyRedshiftMinimum             ))
          timeMinimum=max(                                                                                                                                                &
               &          self%cosmologyFunctions_%cosmicTime            (self%cosmologyFunctions_%expansionFactorFromRedshift(self%surveyRedshiftMaximum             )), &
               &          self%cosmologyFunctions_%timeAtDistanceComoving(self%surveyGeometry_    %distanceMaximum            (self%massBinCenterI       ,field=iField))  &
               &         )
          ! Get the normalizing volume integral.
          volumeNormalization=+volumeNormalization                                 &
               &              +integratorVolume%integrate(timeMinimum,timeMaximum)

          ! Integrate mass function over the bin.
          massFunction(i)=+massFunction             (i                      ) &
               &          +integratorTimeI%integrate(timeMinimum,timeMaximum)
          ! Find the effective volume of the survey at this mass.
          volume(i,iField)=self%surveyGeometry_%volumeMaximum(self%massBinCenterI,iField)
       end do
       ! Normalize the mass function.
       massFunction(i)=+     massFunction      (i) &
            &          /self%logMassBinWidth   (i) &
            &          /     volumeNormalization
       ! Tabulate the bias as a function of time in this bin.
       do iTime=1,self%countTimeBins
          self%time              =+self           %timeTable      (iTime                              )
          self%biasTable(iTime,i)=+integratorBiasI%integrate      (self%logMassLower,self%logMassUpper) &
               &                  /self           %logMassBinWidth(i                                  )
       end do
    end do
    ! Compute LSS variance if necessary.
    if (self%includeLSS) then
       ! Compute large-scale structure variance for each cell pair.
       ! If angular power spectrum of survey window function is available, use it to compute LSS contribution to variance.
       if (self%surveyGeometry_%angularPowerAvailable()) then
          call massFunctionCovarianceLSSAngularSpectrum(self,self%countMassBins,self%surveyRedshiftMinimum,self%surveyRedshiftMaximum,varianceLSS)
       ! If survey function function is available, use it to compute LSS contribution to variance.
       else if (self%surveyGeometry_%windowFunctionAvailable()) then
          call massFunctionCovarianceLSSWindowFunction (self,self%countMassBins,self%surveyRedshiftMinimum,self%surveyRedshiftMaximum,varianceLSS)
       ! No method exists to compute the LSS contribution to variance. Abort.
       else
          call Error_Report('no method exists to compute LSS contribution to covariance matrix'//{introspection:location})
       end if
    end if
    ! Construct the covariance matrix.
    covariancePoisson=0.0d0
    covarianceHalo   =0.0d0
    covarianceLSS    =0.0d0
    do i   =1,self%countMassBins
       self%massBinCenterI    =10.0d0** self%logMassBinCenter(i)
       self%massBinMinimumI   =10.0d0**(self%logMassBinCenter(i)-0.5d0*self%log10MassBinWidth(i))
       self%massBinMaximumI   =10.0d0**(self%logMassBinCenter(i)+0.5d0*self%log10MassBinWidth(i))
       do j=i,self%countMassBins
          self%massBinCenterJ =10.0d0** self%logMassBinCenter(j)
          self%massBinMinimumJ=10.0d0**(self%logMassBinCenter(j)-0.5d0*self%log10MassBinWidth(j))
          self%massBinMaximumJ=10.0d0**(self%logMassBinCenter(j)+0.5d0*self%log10MassBinWidth(j))
          ! Poisson term.
          if (self%includePoisson .and. i == j) then
             if      (useCompleteness) then
                binCompleteness=completenessObserved(i)
             else if (useNumber      ) then
                if (numberObserved(i) > 0.0d0) then
                   binCompleteness=     numberObserved     (i  )  &
                        &          /    massFunctionUse    (i  )  &
                        &          /sum(volume             (i,:)) &
                        &          /    self%logMassBinWidth(i)
                else
                   binCompleteness=1.0d0
                end if
             else
                binCompleteness=1.0d0
             end if
             if (massFunctionUse(i) > 0.0d0) then
                covariancePoisson(i,j)=massFunctionUse(i)/binCompleteness/(sum(volume(i,:))*self%logMassBinWidth(i))
             else
                covariancePoisson(i,j)=1.0d0                             /(sum(volume(i,:))*self%logMassBinWidth(i))**2
             end if
          end if
          ! Halo occupancy covariance.
          if (self%includeHalo) then
             ! Iterate over fields.
             do iField=1,countFields
                ! Find integration limits for this bin.
                timeMaximum=    self%cosmologyFunctions_%cosmicTime            (self%cosmologyFunctions_%expansionFactorFromRedshift(self%surveyRedshiftMinimum))
                timeMinimum=max(                                                                                                                                   &
                     &          self%cosmologyFunctions_%cosmicTime            (self%cosmologyFunctions_%expansionFactorFromRedshift(self%surveyRedshiftMaximum)), &
                     &          self%cosmologyFunctions_%timeAtDistanceComoving(                                                                                   &
                     &                                                          min(                                                                               &
                     &                                                              self%surveyGeometry_%distanceMaximum(self%massBinCenterI,field=iField)       , &
                     &                                                              self%surveyGeometry_%distanceMaximum(self%massBinCenterJ,field=iField)         &
                     &                                                             )                                                                               &
                     &                                                         )                                                                                   &
                     &         )
                if (timeMaximum > timeMinimum) then
                   ! Integrate over the volume. Note that the following expression is multiplied through by the volumes of both
                   ! fields such that we accumulate a volume-weighted covariance, which will be normalized below.
                   covarianceHalo(i,j)=+covarianceHalo                                             (i          ,j          ) &
                        &              +integratorHaloOccupancyTime                %integrate      (timeMinimum,timeMaximum) &
                        &              *self                       %surveyGeometry_%solidAngle     (iField                 ) &
                        &              /self                                       %logMassBinWidth(i                      ) &
                        &              /self                                       %logMassBinWidth(            j          )
                end if
             end do
             ! Normalize the covariance for the total field volume.
             if     (                                                 &
                  &   sum(volume(i,:)) > 0.0d0                        &
                  &  .and.                                            &
                  &   sum(volume(j,:)) > 0.0d0                        &
                  & ) covarianceHalo(i,j)=+    covarianceHalo(i,j  )  &
                  &                       /sum(volume        (i  ,:)) &
                  &                       /sum(volume        (  j,:))
             ! Renormalize to actual mass function. Accounts for any difference between model and data. Including incompleteness.
             if     (                                                                                         &
                  &   massFunctionUse(i) > massFunctionMinimum .and. massFunctionUse(j) > massFunctionMinimum &
                  &  .and.                                                                                    &
                  &   massFunction   (i) > massFunctionMinimum .and. massFunction   (j) > massFunctionMinimum &
                  & ) then
                covarianceHalo(i,j)=+covarianceHalo( i,j) &
                     &              *massFunctionUse(i  ) &
                     &              /massFunction   (i  ) &
                     &              *massFunctionUse(  j) &
                     &              /massFunction   (  j)
             end if
          end if
          ! Large-scale structure term.
          if (self%includeLSS) then
             covarianceLSS(i,j)=varianceLSS(i,j)
             if     (                                                                                         &
                  &   massFunctionUse(i) > massFunctionMinimum .and. massFunctionUse(j) > massFunctionMinimum &
                  &  .and.                                                                                    &
                  &   massFunction   (i) > massFunctionMinimum .and. massFunction   (j) > massFunctionMinimum &
                  & ) then
                ! Renormalize to actual mass function. Accounts for any difference between model and data. Including incompleteness.
                covarianceLSS(i,j)=+covarianceLSS  (i,j) &
                     &             *massFunctionUse(i  ) &
                     &             /massFunction   (i  ) &
                     &             *massFunctionUse(  j) &
                     &             /massFunction   (  j)
             end if
          end if
       end do
    end do
    ! Symmetrize the covariance matrices.
    do i   =1,self%countMassBins
       do j=1,self%countMassBins
          if (j < i) then
             covariancePoisson(i,j)=covariancePoisson(j,i)
             covarianceHalo   (i,j)=covarianceHalo   (j,i)
             covarianceLSS    (i,j)=covarianceLSS    (j,i)
          end if
       end do
    end do
    ! Sum covariances.
    covariance=+covariancePoisson &
         &     +covarianceHalo    &
         &     +covarianceLSS
    ! Add in any covariance arising from uncertainty in the incompleteness.
    do i   =1,self%countMassBins
       do j=1,self%countMassBins
          covariance(i,j)=+     covariance               (i,j)    &
               &          +self%completenessErrorObserved     **2 &
               &          *     massFunctionUse          (i  )    &
               &          *     massFunctionUse          (  j)
       end do
    end do
    ! Compute the corresponding correlation matrix.
    do    i=1,self%countMassBins
       do j=1,self%countMassBins
          normalization=sqrt(                 &
               &             +covariance(i,i) &
               &             *covariance(j,j) &
               &            )
          if (normalization > 0.0d0) then
             correlation(i,j)=covariance(i,j)/normalization
          else
             correlation(i,j)=0.0d0
          end if
       end do
    end do
    ! Deallocate arrays.
    deallocate(self%logMassBinCenter)
    deallocate(     volume          )
    deallocate(     varianceLSS     )
    ! Write out the covariance matrix.
    call massFunctionFile %writeDataset  (mass               ,"mass"             ,"Mass; M [M☉]"                                    ,datasetReturned=dataset)
    call dataset%writeAttribute(massSolar          ,"unitsInSI"                                    )
    call massFunctionFile %writeDataset  (massFunction       ,"massFunction"     ,"Mass function; dn/dln(M) [Mpc⁻³]"                ,datasetReturned=dataset)
    call dataset%writeAttribute(1.0d0/megaParsec**3,"unitsInSI"                                    )
    call massFunctionFile %writeDataset  (covariance         ,"covariance"       ,"Covariance of mass function; [Mpc⁻⁶]"            ,datasetReturned=dataset)
    call dataset%writeAttribute(1.0d0/megaParsec**3,"unitsInSI"                                    )
    call massFunctionFile %writeDataset  (covariancePoisson  ,"covariancePoisson","Covariance due to Poisson noise; [Mpc⁻⁶]"        ,datasetReturned=dataset)
    call dataset%writeAttribute(1.0d0/megaParsec**3,"unitsInSI"                                    )
    call massFunctionFile %writeDataset  (covarianceHalo     ,"covarianceHalo"   ,"Covariance due to halo effect; [Mpc⁻⁶]"          ,datasetReturned=dataset)
    call dataset%writeAttribute(1.0d0/megaParsec**3,"unitsInSI"                                    )
    call massFunctionFile %writeDataset  (covarianceLSS      ,"covarianceLSS"    ,"Covariance due to large scale structure; [Mpc⁻⁶]",datasetReturned=dataset)
    call dataset%writeAttribute(1.0d0/megaParsec**3,"unitsInSI"                                    )
    call massFunctionFile %writeDataset  (correlation        ,"correlation"      ,"Correlation matrix for stellar mass function; []"                        )
    ! Done.
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: mass function covariance' )
    return
  end subroutine massFunctionCovariancePerform

  double precision function massFunctionCovarianceGalaxyRootPowerSpectrum(iBin,timeMinimum,timeMaximum)
    !!{
    Computes the quantity $\int_{t_\mathrm{min}}^{t_\mathrm{max}} \mathrm{d} t b(t) \sqrt{P(k,t)} \mathrm{d} V / \mathrm{d}t$, where $b(t)$ is
    galaxy bias, and $P(k,t)$ is the non-linear galaxy power spectrum.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    integer                     , intent(in   ) :: iBin
    double precision            , intent(in   ) :: timeMinimum, timeMaximum
    type            (integrator)                :: integrator_

    selfCopy%lssBin        =iBin
    integrator_                                  =integrator           (massFunctionCovarianceLargeScaleStructureIntegrand,toleranceRelative=1.0d-2     )
    massFunctionCovarianceGalaxyRootPowerSpectrum=integrator_%integrate(timeMinimum                                       ,                  timeMaximum)
    return
  end function massFunctionCovarianceGalaxyRootPowerSpectrum

  double precision function massFunctionCovarianceAngularPowerIntegrand(wavenumber)
    !!{
    Integrand for large scale structure variance computed using survey mask angular power spectrum.
    !!}
    implicit none
    double precision, intent(in   ) :: wavenumber
    integer                         :: iField               , jField               , &
         &                             l
    double precision                :: x0i                  , x1i                  , &
         &                             x0j                  , x1j                  , &
         &                             powerSpectrumI       , powerSpectrumJ       , &
         &                             surveyDistanceMinimum, surveyDistanceMaximum, &
         &                             angularFactor

    massFunctionCovarianceAngularPowerIntegrand=0.0d0
    if (wavenumber <= 0.0d0) return
    selfCopy%wavenumberGlobal=wavenumber
    surveyDistanceMinimum=selfCopy%cosmologyFunctions_ %distanceComoving          (                         &
         &                selfCopy%cosmologyFunctions_%cosmicTime                  (                        &
         &                selfCopy%cosmologyFunctions_%expansionFactorFromRedshift  (                       &
         &                selfCopy%                                                  surveyRedshiftMinimum  &
         &                                                                                                ) &
         &                                                                                               )  &
         &                                                                                              )
    surveyDistanceMaximum=selfCopy%cosmologyFunctions_%distanceComoving           (                         &
         &                selfCopy%cosmologyFunctions_%cosmicTime                  (                        &
         &                selfCopy%cosmologyFunctions_%expansionFactorFromRedshift  (                       &
         &                selfCopy%                                                  surveyRedshiftMaximum  &
         &                                                                                                ) &
         &                                                                                               )  &
         &                                                                                              )
    do iField=1,selfCopy%surveyGeometry_%fieldCount()
       if (selfCopy%timeMinimumI(iField) >= selfCopy%timeMaximumI(iField)) cycle
       powerSpectrumI=+selfCopy%surveyGeometry_                     %solidAngle                      (iField)  &
            &         *massFunctionCovarianceGalaxyRootPowerSpectrum           (                               &
            &                                                                   selfCopy%binI                , &
            &                                                                   selfCopy%timeMinimumI(iField), &
            &                                                                   selfCopy%timeMaximumI(iField)  &
            &                                                                  )
       if (selfCopy%surveyRedshiftMinimum <= 0.0d0) then
          x0i=   +0.0d0
       else
          x0i=                                                                                                 &
               & +selfCopy%wavenumberGlobal                                                                    &
               & *    surveyDistanceMinimum
       end if
       x1i=                                                                                                    &
            &    +selfCopy%wavenumberGlobal                                                                    &
            &    *min(                                                                                                         &
            &         surveyDistanceMaximum                                                                                  , &
            &         selfCopy%surveyGeometry_%distanceMaximum(10.0d0**selfCopy%logMassBinCenter(selfCopy%binI),field=iField)  &
            &        )
       do jField=1,selfCopy%surveyGeometry_%fieldCount()
          if (selfCopy%timeMinimumJ(jField) >= selfCopy%timeMaximumJ(jField)) cycle
          powerSpectrumJ=+selfCopy                                     %surveyGeometry_%solidAngle                      (jField)   &
               &         *massFunctionCovarianceGalaxyRootPowerSpectrum                           (                                &
               &                                                                                   selfCopy%binJ                 , &
               &                                                                                   selfCopy%timeMinimumJ(jField) , &
               &                                                                                   selfCopy%timeMaximumJ(jField)   &
               &                                                                                  )
          if (selfCopy%surveyRedshiftMinimum <= 0.0d0) then
             x0j=   +0.0d0
          else
             x0j=                                                                                                                  &
                  & +selfCopy%wavenumberGlobal                                                                                     &
                  & *    surveyDistanceMinimum
          end if
          x1j=                                                                                                                     &
               &    +selfCopy%wavenumberGlobal                                                                                     &
               &    *min(                                                                                                          &
               &         surveyDistanceMaximum                                                                                   , &
               &         selfCopy%surveyGeometry_%distanceMaximum(10.0d0**selfCopy%logMassBinCenter(selfCopy%binJ),field=jField)   &
               &       )
          angularFactor=0.0d0
          !$omp parallel do reduction(+:angularFactor) copyin(selfCopy)
          do l=0,selfCopy%surveyGeometry_%angularPowerMaximumDegree()
             angularFactor=+angularFactor                                                                              &
                  &        +dble(2*l+1)                                                                                &
                  &        *selfCopy                                    %surveyGeometry_%angularPower(iField,jField,l) &
                  &        *massFunctionCovarianceAngularPowerRadialTerm                             (x0i   ,x1i   ,l) &
                  &        *massFunctionCovarianceAngularPowerRadialTerm                             (x0j   ,x1j   ,l)
          end do
          !$omp end parallel do
          massFunctionCovarianceAngularPowerIntegrand=massFunctionCovarianceAngularPowerIntegrand+powerSpectrumI*powerSpectrumJ*angularFactor
       end do
    end do
    massFunctionCovarianceAngularPowerIntegrand=+massFunctionCovarianceAngularPowerIntegrand                     &
         &                                      /selfCopy                                   %wavenumberGlobal**4
    return
  end function massFunctionCovarianceAngularPowerIntegrand

  double precision function massFunctionCovarianceAngularPowerRadialTerm(x0,x1,l)
    !!{
    Computes the radial term in the expression for large scale structure variance.
    !!}
    use :: Gamma_Functions         , only : Gamma_Function_Logarithmic
    use :: Hypergeometric_Functions, only : Hypergeometric_pFq
    use :: Numerical_Constants_Math, only : Pi                        , ln2
    implicit none
    double precision, intent(in   ) :: x0                , x1
    integer         , intent(in   ) :: l
    double precision, parameter     :: xMaximum  = 512.0d0
    double precision, parameter     :: aMinimum  =-750.0d0
    integer         , save          :: lPrevious =-1
    double precision, save          :: x0Previous=-1.0d0 , x1Previous=-1.0d0
    !$omp threadprivate(lPrevious,x0Previous,x1Previous)
    double precision, save          :: h0                , h1               , &
         &                             logGammas
    !$omp threadprivate(h0,h1,logGammas)
    double precision                :: a0                , a1

    ! Evaluate combination of logarithms of Gamma functions.
    if (l /= lPrevious)                                                 &
         & logGammas=+Gamma_Function_Logarithmic(0.5d0*(3.0d0+dble(l))) &
         &           -Gamma_Function_Logarithmic(       1.5d0+dble(l) ) &
         &           -Gamma_Function_Logarithmic(0.5d0*(5.0d0+dble(l)))
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
    massFunctionCovarianceAngularPowerRadialTerm=+sqrt(Pi) &
         &                                       *(        &
         &                                         +h1     &
         &                                         -h0     &
         &                                        )
    lPrevious=l
    return
  end function massFunctionCovarianceAngularPowerRadialTerm

  double precision function massFunctionCovarianceVolumeIntegrand(time)
    !!{
    Integral for comoving volume.
    !!}
    implicit none
    double precision, intent(in   ) :: time

    massFunctionCovarianceVolumeIntegrand=selfCopy%cosmologyFunctions_%comovingVolumeElementTime(time)
    return
  end function massFunctionCovarianceVolumeIntegrand

  double precision function massFunctionCovarianceMassFunctionTimeIntegrandI(timePrime)
    !!{
    Integral for comoving volume.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    double precision            , intent(in   ) :: timePrime
    type            (integrator)                :: integrator_
    double precision                            :: massFunction

    self_                                           %time=+timePrime
    integrator_                                          = integrator                                         (massFunctionCovarianceMassFunctionIntegrandI,toleranceRelative=1.0d-3            )
    massFunction                                         =+integrator_              %integrate                (self_%logMassLower                          ,                  self_%logMassUpper)
    massFunctionCovarianceMassFunctionTimeIntegrandI     =+massFunction                                                                                                                           &
         &                                                *self_%cosmologyFunctions_%comovingVolumeElementTime(self_%time                                                                       )
    return
  end function massFunctionCovarianceMassFunctionTimeIntegrandI

  double precision function massFunctionCovarianceMassFunctionIntegrandI(logMass)
    !!{
    Integral for mass function.
    !!}
    implicit none
    double precision, intent(in   ) :: logMass
    double precision                :: mass

    mass=10.0d0**logMass
    massFunctionCovarianceMassFunctionIntegrandI=+     self_%haloMassFunction_       %differential(self_%time,mass)             &
         &                                       *                                                                      mass    &
         &                                       *log(10.0d0)                                                                   &
         &                                       *max(                                                                          &
         &                                            +self_%conditionalMassFunction_%massFunction(mass,self_%massBinMinimumI)  &
         &                                            -self_%conditionalMassFunction_%massFunction(mass,self_%massBinMaximumI), &
         &                                            +0.0d0                                                                    &
         &                                           )
    return
  end function massFunctionCovarianceMassFunctionIntegrandI

  double precision function massFunctionCovarianceLargeScaleStructureIntegrand(timePrime)
    !!{
    Integral for LSS contribution to the covariance matrix.
    !!}
    use :: Numerical_Interpolation, only : interpolator
    implicit none
    double precision              , intent(in   ) :: timePrime
    type            (interpolator)                :: interpolator_
    double precision                              :: bias         , powerSpectrumValue

    ! Copy the time to module scope.
    selfCopy%time=timePrime
    ! Get the bias-mass function product for the I bin.
    interpolator_=interpolator             (selfCopy%timeTable,selfCopy%biasTable(:,selfCopy%lssBin))
    bias         =interpolator_%interpolate(selfCopy%time                                           )
    ! Get the nonlinear power spectrum for the current wavenumber and time.
    powerSpectrumValue=selfCopy%powerSpectrumNonlinear_%value(selfCopy%waveNumberGlobal,selfCopy%time)
    ! Return the cross-correlation biased power spectrum multiplied by the volume element.
    massFunctionCovarianceLargeScaleStructureIntegrand=+bias                                                                  &
         &                                             *sqrt(powerSpectrumValue)                                              &
         &                                             *selfCopy%cosmologyFunctions_%comovingVolumeElementTime(selfCopy%time)
    return
  end function massFunctionCovarianceLargeScaleStructureIntegrand

  double precision function massFunctionCovarianceBiasIntegrandI(logMass)
    !!{
    Integral for bias.
    !!}
    implicit none
    double precision, intent(in   ) :: logMass
    double precision                :: mass

    mass                                =+10.0d0**logMass
    massFunctionCovarianceBiasIntegrandI=+massFunctionCovarianceMassFunctionIntegrandI                         (logMass           ) &
         &                               *self_                                       %darkMatterHaloBias_%bias(   mass,self_%time)
    return
  end function massFunctionCovarianceBiasIntegrandI

  double precision function massFunctionCovarianceHaloOccupancyTimeInterand(timePrime)
    !!{
    Integral for comoving volume.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    double precision            , intent(in   ) :: timePrime
    type            (integrator)                :: integrator_
    double precision                            :: massFunction

    self_%time=timePrime
    integrator_                                    =integrator           (                                                                &
         &                                                                                  massFunctionCovarianceHaloOccupancyIntegrand, &
         &                                                                toleranceRelative=1.0d-3                                        &
         &                                                               )
    massFunction                                   =integrator_%integrate(                                                                &
         &                                                                                  self_%logMassLower                          , &
         &                                                                                  self_%logMassUpper                            &
         &                                                               )
    massFunctionCovarianceHaloOccupancyTimeInterand=+massFunction                                                                         &
         &                                          *self_%cosmologyFunctions_%comovingVolumeElementTime(self_%time)
      return
  end function massFunctionCovarianceHaloOccupancyTimeInterand

  double precision function massFunctionCovarianceHaloOccupancyIntegrand(logMass)
    !!{
    Integral for mass function.
    !!}
    implicit none
    double precision, intent(in   ) :: logMass
    double precision                :: mass

    mass                                        =+10.0d0**logMass
    massFunctionCovarianceHaloOccupancyIntegrand=+     self_%haloMassFunction_       %differential(self_%time,mass)             &
         &                                       *                                                                      mass    &
         &                                       *log(10.0d0)                                                                   &
         &                                       *max(                                                                          &
         &                                            +self_%conditionalMassFunction_%massFunction(mass,self_%massBinMinimumI)  &
         &                                            -self_%conditionalMassFunction_%massFunction(mass,self_%massBinMaximumI), &
         &                                            +0.0d0                                                                    &
         &                                           )                                                                          &
         &                                       *max(                                                                          &
         &                                            +self_%conditionalMassFunction_%massFunction(mass,self_%massBinMinimumJ)  &
         &                                            -self_%conditionalMassFunction_%massFunction(mass,self_%massBinMaximumJ), &
         &                                            +0.0d0                                                                    &
         &                                           )
    return
  end function massFunctionCovarianceHaloOccupancyIntegrand

  subroutine massFunctionCovarianceComputeVolumeNormalizations(logMass,redshiftMinimum,redshiftMaximum,timeMinimum,timeMaximum,volumeNormalization)
    !!{
    Compute volume normalization factors for LSS covariance calculations.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    double precision            , intent(in   )               :: logMass            , redshiftMinimum, &
         &                                                       redshiftMaximum
    double precision            , intent(  out), dimension(:) :: timeMinimum        , timeMaximum    , &
         &                                                       volumeNormalization
    integer                                                   :: iField
    type            (integrator)                              :: integrator_

    integrator_=integrator(massFunctionCovarianceVolumeIntegrand,toleranceRelative=1.0d-3)
    do iField=1,selfCopy%surveyGeometry_%fieldCount()
       ! Find integration limits for this bin. We want the maximum of the volumes associated with the two bins.
       timeMaximum(iField)=        selfCopy%cosmologyFunctions_%cosmicTime            (selfCopy%cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum            ))
       timeMinimum(iField)=min(                                                                                                                                                 &
            &                  max(                                                                                                                                             &
            &                      selfCopy%cosmologyFunctions_%cosmicTime            (selfCopy%cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum             )), &
            &                      selfCopy%cosmologyFunctions_%timeAtDistanceComoving(selfCopy%surveyGeometry_    %distanceMaximum            (10.0d0**logMass,field=iField))  &
            &                     )                                                                                                                                           , &
            &                  timeMaximum(iField)                                                                                                                              &
            &                 )
       ! Get the normalizing volume integral for bin i.
       volumeNormalization(iField)=+integrator_             %integrate (                     &
            &                                                           timeMinimum(iField), &
            &                                                           timeMaximum(iField)  &
            &                                )                                               &
            &                      *selfCopy%surveyGeometry_%solidAngle(                     &
            &                                                                       iField   &
            &                                                          )
    end do
    return
  end subroutine massFunctionCovarianceComputeVolumeNormalizations

  subroutine massFunctionCovarianceLSSWindowFunction(self,massBinCount,redshiftMinimum,redshiftMaximum,varianceLSS)
    !!{
    Compute variance due to large scale structure by directly summing over the Fourier transform
    of the survey selection function.
    !!}
    use            :: Display                 , only : displayCounter  , displayCounterClear
    use            :: FFTW3                   , only : FFTW_Wavenumber
    use, intrinsic :: ISO_C_Binding           , only : c_double_complex
    use            :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (taskMassFunctionCovariance), intent(inout), target           :: self
    integer                                     , intent(in   )                   :: massBinCount
    double precision                            , intent(in   )                   :: redshiftMinimum, redshiftMaximum
    double precision                            , intent(  out), dimension(:,:  ) :: varianceLSS
    complex         (c_double_complex          ), allocatable  , dimension(:,:,:) :: windowFunctionI, windowFunctionJ
    integer                                                                       :: i              , j              , &
         &                                                                           u              , w              , &
         &                                                                           v              , taskCount      , &
         &                                                                           taskTotal      , countFields    , &
         &                                                                           iField
    double precision                                                              :: waveNumberU    , waveNumberV    , &
         &                                                                           waveNumberW    , variance       , &
         &                                                                           powerSpectrumI , powerSpectrumJ , &
         &                                                                           powerSpectrum  , normalizationI , &
         &                                                                           normalizationJ , multiplier     , &
         &                                                                           boxLength

    ! Allocate arrays for times and volume normalizations.
    countFields=self_%surveyGeometry_%fieldCount()
    ! Allocate arrays for survey window functions if these will be used.
    allocate(windowFunctionI(                      &
         &                   selfCopy%sizeGridFFT, &
         &                   selfCopy%sizeGridFFT, &
         &                   selfCopy%sizeGridFFT  &
         &                  )                      &
         &  )
    allocate(windowFunctionJ(                      &
         &                   selfCopy%sizeGridFFT, &
         &                   selfCopy%sizeGridFFT, &
         &                   selfCopy%sizeGridFFT  &
         &                  )                      &
         &  )
    taskTotal  =massBinCount*(massBinCount+1)/2
    taskCount  =0
    varianceLSS=0.0d0
    do i   =1,massBinCount
       selfCopy%binI=i
       selfCopy%massBinCenterI =10.0d0** selfCopy%logMassBinCenter(i)
       selfCopy%massBinMinimumI=10.0d0**(selfCopy%logMassBinCenter(i)-0.5d0*selfCopy%log10MassBinWidth(i))
       selfCopy%massBinMaximumI=10.0d0**(selfCopy%logMassBinCenter(i)+0.5d0*selfCopy%log10MassBinWidth(i))
       do j=i,massBinCount
          selfCopy%binJ=j
          selfCopy%massBinCenterJ =10.0d0** selfCopy%logMassBinCenter(j)
          selfCopy%massBinMinimumJ=10.0d0**(selfCopy%logMassBinCenter(j)-0.5d0*selfCopy%log10MassBinWidth(j))
          selfCopy%massBinMaximumJ=10.0d0**(selfCopy%logMassBinCenter(j)+0.5d0*selfCopy%log10MassBinWidth(j))
          ! Update progress.
          call displayCounter(                                              &
               &              int(100.0d0*dble(taskCount)/dble(taskTotal)), &
               &              isNew=(taskCount==0)                          &
               &             )
          taskCount=taskCount+1
          ! Compute window functions for this pair of cells.
          call selfCopy%surveyGeometry_%windowFunctions(                                           &
               &                                                          selfCopy%massBinCenterI, &
               &                                                          selfCopy%massBinCenterJ, &
               &                                                          selfCopy%sizeGridFFT   , &
               &                                                          boxLength              , &
               &                                                          windowFunctionI        , &
               &                                                          windowFunctionJ          &
               &                                                         )
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
          allocate(selfCopy,mold=self)
          !$omp critical(massFunctionCovarianceDeepCopy)
          !![
          <deepCopyReset variables="self"/>
          <deepCopy source="self" destination="selfCopy"/>
          <deepCopyFinalize variables="selfCopy"/>
          !!]
          !$omp end critical(massFunctionCovarianceDeepCopy)
          allocate(selfCopy%volumeNormalizationI(countFields))
          allocate(selfCopy%volumeNormalizationJ(countFields))
          allocate(selfCopy%timeMinimumI        (countFields))
          allocate(selfCopy%timeMinimumJ        (countFields))
          allocate(selfCopy%timeMaximumI        (countFields))
          allocate(selfCopy%timeMaximumJ        (countFields))
          call massFunctionCovarianceComputeVolumeNormalizations(                                  &
               &                                                 selfCopy%logMassBinCenter    (i), &
               &                                                 redshiftMinimum                 , &
               &                                                 redshiftMaximum                 , &
               &                                                 selfCopy%timeMinimumI           , &
               &                                                 selfCopy%timeMaximumI           , &
               &                                                 selfCopy%volumeNormalizationI     &
               &                                                )
          call massFunctionCovarianceComputeVolumeNormalizations(                                  &
               &                                                 selfCopy%logMassBinCenter    (j), &
               &                                                 redshiftMinimum                 , &
               &                                                 redshiftMaximum                 , &
               &                                                 selfCopy%timeMinimumJ           , &
               &                                                 selfCopy%timeMaximumJ           , &
               &                                                 selfCopy%volumeNormalizationJ     &
               &                                                )
          !$omp do reduction(+:variance)
          do u      =1,selfCopy%sizeGridFFT/2+1
             waveNumberU      =FFTW_Wavenumber(u,selfCopy%sizeGridFFT)*2.0d0*Pi/boxLength
             do v   =1,selfCopy%sizeGridFFT/2+1
                waveNumberV   =FFTW_Wavenumber(v,selfCopy%sizeGridFFT)*2.0d0*Pi/boxLength
                do w=1,selfCopy%sizeGridFFT/2+1
                   waveNumberW=FFTW_Wavenumber(w,selfCopy%sizeGridFFT)*2.0d0*Pi/boxLength
                   ! Compute the wavenumber for this cell.
                   selfCopy%waveNumberGlobal=sqrt(waveNumberU**2+waveNumberV**2+waveNumberW**2)
                   ! Find the power spectrum for this wavenumber.
                   if (selfCopy%waveNumberGlobal > 0.0d0) then
                      ! Integrate the power spectrum, weighted by the galaxy bias, over the
                      ! volume of interest. Then normalize by that volume.
                      normalizationI=0.0d0
                      normalizationJ=0.0d0
                      powerSpectrumI=0.0d0
                      powerSpectrumJ=0.0d0
                      do iField=1,countFields
                         powerSpectrumI  =+powerSpectrumI                                                               &
                              &           +selfCopy%surveyGeometry_%solidAngle                                (iField)  &
                              &           *massFunctionCovarianceGalaxyRootPowerSpectrum(                               &
                              &                                                          selfCopy%binI                , &
                              &                                                          selfCopy%timeMinimumI(iField), &
                              &                                                          selfCopy%timeMaximumI(iField)  &
                              &                                                         )
                         normalizationI  =+normalizationI                                                               &
                              &           +selfCopy%volumeNormalizationI                                      (iField)
                         powerSpectrumJ  =+powerSpectrumJ                                                               &
                              &           +selfCopy%surveyGeometry_%solidAngle                                (iField)  &
                              &           *massFunctionCovarianceGalaxyRootPowerSpectrum(                               &
                              &                                                          selfCopy%binJ                , &
                              &                                                          selfCopy%timeMinimumJ(iField), &
                              &                                                          selfCopy%timeMaximumJ(iField)  &
                              &                                                         )
                         normalizationJ  =+normalizationJ                                                               &
                              &           +selfCopy%volumeNormalizationJ                                      (iField)
                      end do
                      powerSpectrum=+powerSpectrumI &
                           &        /normalizationI &
                           &        *powerSpectrumJ &
                           &        /normalizationJ
                   else
                      powerSpectrum=0.0d0
                   end if
                   ! Add the contribution from this cell to the total variance.
                   multiplier=2.0d0
                   if     (                                   &
                        &  u == selfCopy%sizeGridFFT/2+1 .or. &
                        &  v == selfCopy%sizeGridFFT/2+1 .or. &
                        &  w == selfCopy%sizeGridFFT/2+1      &
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
          !![
          <objectDestructor name="selfCopy"/>
          !!]
          !$omp end parallel
          selfCopy => self
          ! Normalize the variance. We multiply by (2π/L)³ to account for the volume of each FFT
          ! cell, and divide by (2π)³ as defined in eqn. (66) of Smith (2012; MNRAS; 426; 531).
          varianceLSS(i,j)=dble(variance)/boxLength**3
       end do
    end do
    if (allocated(windowFunctionI)) deallocate(windowFunctionI)
    if (allocated(windowFunctionJ)) deallocate(windowFunctionJ)
    call displayCounterClear()
    return
  end subroutine massFunctionCovarianceLSSWindowFunction

  subroutine massFunctionCovarianceLSSAngularSpectrum(self,massBinCount,redshiftMinimum,redshiftMaximum,varianceLSS)
    !!{
    Compute variance due to large scale structure by integration over the angular power spectrum.
    !!}
    use :: Display                 , only : displayCounter, displayCounterClear
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    implicit none
    class           (taskMassFunctionCovariance), intent(inout), target           :: self
    integer                                     , intent(in   )                   :: massBinCount
    double precision                            , intent(in   )                   :: redshiftMinimum               , redshiftMaximum
    double precision                            , intent(  out), dimension(:,:  ) :: varianceLSS
    ! Dimensionless factor controlling the highest wavenumber to be used when integrating over
    ! angular power spectra.
    double precision                            , parameter                       :: wavenumberMaximumFactor=1.0d0
    integer                                                                       :: i                            , j                , &
         &                                                                           countFields                  , iField           , &
         &                                                                           taskCount                    , taskTotal
    double precision                                                              :: wavenumberMinimum            , wavenumberMaximum, &
         &                                                                           distanceMaximum
    type            (integrator                )                                  :: integrator_

    countFields=self_%surveyGeometry_%fieldCount()
    !$ call OMP_Set_Nested(.true.)
    taskTotal  =massBinCount*(massBinCount+1)/2
    taskCount  =0
    !$omp parallel private (i,j,wavenumberMinimum,wavenumberMaximum,integrator_)
    allocate(selfCopy,mold=self)
    !$omp critical(massFunctionCovarianceDeepCopy)
    !![
    <deepCopyReset variables="self"/>
    <deepCopy source="self" destination="selfCopy"/>
    <deepCopyFinalize variables="selfCopy"/>
    !!]
    !$omp end critical(massFunctionCovarianceDeepCopy)
    allocate(selfCopy%volumeNormalizationI(countFields))
    allocate(selfCopy%volumeNormalizationJ(countFields))
    allocate(selfCopy%timeMinimumI        (countFields))
    allocate(selfCopy%timeMinimumJ        (countFields))
    allocate(selfCopy%timeMaximumI        (countFields))
    allocate(selfCopy%timeMaximumJ        (countFields))
    integrator_=integrator(massFunctionCovarianceAngularPowerIntegrand,toleranceRelative=1.0d-2)
    !$omp do schedule (dynamic)
    do i=1,massBinCount
       selfCopy%binI=i
       call massFunctionCovarianceComputeVolumeNormalizations(                                              &
            &                                                 selfCopy%logMassBinCenter    (selfCopy%binI), &
            &                                                 redshiftMinimum                             , &
            &                                                 redshiftMaximum                             , &
            &                                                 selfCopy%timeMinimumI                       , &
            &                                                 selfCopy%timeMaximumI                       , &
            &                                                 selfCopy%volumeNormalizationI                 &
            &                                                )
       do j=selfCopy%binI,massBinCount
          ! Update progress.
          call displayCounter(                                              &
               &              int(100.0d0*dble(taskCount)/dble(taskTotal)), &
               &              isNew=(taskCount==0)                          &
               &             )
          selfCopy%binJ=j
          call massFunctionCovarianceComputeVolumeNormalizations(                                              &
               &                                                 selfCopy%logMassBinCenter    (selfCopy%binJ), &
               &                                                 redshiftMinimum                             , &
               &                                                 redshiftMaximum                             , &
               &                                                 selfCopy%timeMinimumJ                       , &
               &                                                 selfCopy%timeMaximumJ                       , &
               &                                                 selfCopy%volumeNormalizationJ                 &
               &                                                )
          distanceMaximum=huge(1.0d0)
          do iField=1,selfCopy%surveyGeometry_%fieldCount()
             distanceMaximum=min(                                                                                                         &
                  &                                       distanceMaximum                                                               , &
                  &              selfCopy%surveyGeometry_%distanceMaximum(10.0d0**selfCopy%logMassBinCenter(selfCopy%binI),field=iField), &
                  &              selfCopy%surveyGeometry_%distanceMaximum(10.0d0**selfCopy%logMassBinCenter(selfCopy%binJ),field=iField)  &
                  &             )
          end do
          wavenumberMinimum=0.0d0
          wavenumberMaximum=wavenumberMaximumFactor                                     &
               &            *max(                                                       &
               &                  1.0d0                                               , &
               &                  selfCopy%surveyGeometry_%angularPowerMaximumDegree()  &
               &                 /2.0d0                                                 &
               &                 /Pi                                                    &
               &                )                                                       &
               &            /distanceMaximum
          varianceLSS(selfCopy%binI,selfCopy%binJ)=+2.0d0                                        &
               &                                   /Pi                                           &
               &                                   /sum(selfCopy   %volumeNormalizationI)**2     &
               &                                   /sum(selfCopy   %volumeNormalizationJ)**2     &
               &                                   *    integrator_%integrate(                   &
               &                                                              wavenumberMinimum, &
               &                                                              wavenumberMaximum  &
               &                                                             )
          !$omp atomic
          taskCount=taskCount+1
       end do
    end do
    !$omp end do
    !![
    <objectDestructor name="selfCopy"/>
    !!]
    !$omp end parallel
    selfCopy => self
    return
  end subroutine massFunctionCovarianceLSSAngularSpectrum

  logical function massFunctionCovarianceRequiresOutputFile(self)
    !!{
    Specifies that this task does not requires the main output file.
    !!}
    implicit none
    class(taskMassFunctionCovariance), intent(inout) :: self
    !$GLC attributes unused :: self

    massFunctionCovarianceRequiresOutputFile=.false.
    return
  end function massFunctionCovarianceRequiresOutputFile
