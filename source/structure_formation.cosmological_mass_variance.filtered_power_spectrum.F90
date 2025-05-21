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

  !+    Contributions to this file made by: Andrew Benson, Christoph Behrens, Xiaolong Du.

  !!{
  An implementation of cosmological density field mass variance computed using a filtered power spectrum.
  !!}

  !![
  <cosmologicalMassVariance name="cosmologicalMassVarianceFilteredPower">
   <description>
    Mass variance of cosmological density fields computed from a filtered power spectrum:
    \begin{equation}
     \sigma^2(M) = {1 \over 2 \pi^2} \int_0^\infty P(k) T^2(k) W^2(k) k^2 \mathrm{d}k
    \end{equation}
    where $P(k)$ is the primordial power spectrum (see \refPhysics{powerSpectrumPrimordial}), $T(k)$ is the transfer function
    (see \refPhysics{transferFunction}), and $W(k)$ is the power spectrum variance window function (see
    \refPhysics{powerSpectrumWindowFunction}).
  
    The normalization of the mass variance is specified via the {\normalfont \ttfamily [sigma\_8]} parameter, which defines the
    linear theory root-variance of the density field in spheres of radii $8h^{-1}$Mpc. Note that when computing the
    normalization of the power spectrum to match the specified value of $\sigma_8$ a top-hat real-space window function is
    used (as per the definition of $\sigma_8$), unless a different window function is explicitly defined via the {\normalfont
    \ttfamily [powerSpectrumWindowFunctionTopHat]} parameter.
  
    The mass variance, $\sigma(M)$, is found by integration over the linear theory power spectrum, with the specified power
    spectrum window function. The fractional tolerance for this integration can be set via the {\normalfont \ttfamily
    [tolerance]} parameter. (The normalization of $\sigma(M)$ to give the desired $\sigma_8$ always uses a top-hat window
    function. For this integration the tolerance can be set via the {\normalfont \ttfamily [toleranceTopHat]} parameter.) This
    is tabulated across the required range.
  
    Cubic spline interpolation is then used to interpolate in this table to give $\sigma(M)$ at any required value of $M$. The
    tabulation is always forced to be monotonically decreasing with $M$. However, the interpolation is not necessarily
    monotonic---for example in cases where $\sigma(M)$ becomes constant or close to constant as a function of $M$ the
    interpolation can become non-monotonic over some ranges of $M$. If strict monotonicity is required set {\normalfont
    \ttfamily [monotonicInterpolation]}={\normalfont \ttfamily true}. This causes a monotonic spline interpolator to be used
    instead which guarantees monotonicity.
   </description>
  </cosmologicalMassVariance>
  !!]
  use :: Cosmology_Functions                 , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters                , only : cosmologyParametersClass
  use :: Linear_Growth                       , only : linearGrowthClass
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredClass
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionClass
  use :: Transfer_Functions                  , only : transferFunctionClass
  use :: Tables                              , only : table1DLinearCSpline
  use :: File_Utilities                      , only : lockDescriptor

  !![
  <stateStorable class="uniqueTable"/>
  !!]

  type :: uniqueTable
     !!{
     Type used to store unique values of the mass variance.
     !!}
     private
     double precision, allocatable, dimension(:) :: rootVariance
     integer         , allocatable, dimension(:) :: index
  end type uniqueTable

  type, extends(cosmologicalMassVarianceClass) :: cosmologicalMassVarianceFilteredPower
     !!{
     A cosmological mass variance class computing variance from a filtered power spectrum.
     !!}
     private
     class           (cosmologyParametersClass               ), pointer                   :: cosmologyParameters_                => null()
     class           (cosmologyFunctionsClass                ), pointer                   :: cosmologyFunctions_                 => null()
     class           (powerSpectrumPrimordialTransferredClass), pointer                   :: powerSpectrumPrimordialTransferred_ => null() , powerSpectrumPrimordialTransferredReference => null()
     class           (cosmologicalMassVarianceClass          ), pointer                   :: cosmologicalMassVarianceReference   => null()
     class           (linearGrowthClass                      ), pointer                   :: linearGrowth_                       => null()
     class           (transferFunctionClass                  ), pointer                   :: transferFunction_                   => null()
     class           (powerSpectrumWindowFunctionClass       ), pointer                   :: powerSpectrumWindowFunction_        => null() , powerSpectrumWindowFunctionTopHat_          => null()
     logical                                                                              :: initialized                         =  .false., nonMonotonicIsFatal                                  , &
          &                                                                                  integrationFailureIsFatal
     double precision                                                                     :: tolerance                                     , toleranceTopHat                                      , &
          &                                                                                  sigma8Value                                   , sigmaNormalization                                   , &
          &                                                                                  massMinimum                                   , massMaximum                                          , &
          &                                                                                  timeMinimum                                   , timeMaximum                                          , &
          &                                                                                  timeMinimumLogarithmic                        , timeLogarithmicDeltaInverse                          , &
          &                                                                                  wavenumberReference                           , wavenumberHalfMode                                   , &
          &                                                                                  rootVarianceLogarithmicGradientTolerance      , amplitudeScalar
     double precision                                         , allocatable, dimension(:) :: times
     class           (table1DLinearCSpline                   ), allocatable, dimension(:) :: rootVarianceTable
     type            (varying_string                         )                            :: fileName
     type            (lockDescriptor                         )                            :: fileLock
     ! Unique values in the variance table and their corresponding indices.
     type            (uniqueTable                            ), allocatable, dimension(:) :: rootVarianceUniqueTable
     logical                                                                              :: monotonicInterpolation                             , growthIsMassDependent_                               , &
          &                                                                                  normalizationSigma8                        =.false., truncateAtParticleHorizon                            , &
          &                                                                                  storeTabulations                                   , warnedNonIncreasing
   contains
     !![
     <methods>
       <method description="Tabulate cosmological mass variance."        method="retabulate"      />
       <method description="Compute the interpolating factors in time."  method="interpolantsTime"/>
       <method description="Write the tabulated mass variance to file."  method="fileWrite"       />
       <method description="Read the tabulated mass variance from file." method="fileRead"        />
       <method description="Return true if the table must be remade."    method="remakeTable"     />
     </methods>
     !!]
     final     ::                                        filteredPowerDestructor
     procedure :: sigma8                              => filteredPowerSigma8
     procedure :: powerNormalization                  => filteredPowerPowerNormalization
     procedure :: rootVariance                        => filteredPowerRootVariance
     procedure :: rootVarianceLogarithmicGradient     => filteredPowerRootVarianceLogarithmicGradient
     procedure :: rootVarianceLogarithmicGradientTime => filteredPowerRootVarianceLogarithmicGradientTime
     procedure :: rootVarianceAndLogarithmicGradient  => filteredPowerRootVarianceAndLogarithmicGradient
     procedure :: mass                                => filteredPowerMass
     procedure :: retabulate                          => filteredPowerRetabulate
     procedure :: interpolantsTime                    => filteredPowerInterpolantsTime
     procedure :: growthIsMassDependent               => filteredPowerGrowthIsMassDependent
     procedure :: fileWrite                           => filteredPowerFileWrite
     procedure :: fileRead                            => filteredPowerFileRead
     procedure :: remakeTable                         => filteredPowerRemakeTable
     procedure :: descriptor                          => filteredPowerDescriptor
     procedure :: descriptorNormalizationOnly         => filteredPowerDescriptorNormalizationOnly
  end type cosmologicalMassVarianceFilteredPower

  interface cosmologicalMassVarianceFilteredPower
     !!{
     Constructors for the \refClass{cosmologicalMassVarianceFilteredPower} cosmological mass variance class.
     !!}
     module procedure filteredPowerConstructorParameters
     module procedure filteredPowerConstructorInternal
  end interface cosmologicalMassVarianceFilteredPower

  ! Number of points per decade to use in tabulation of σ(M).
  integer                         , parameter :: pointsPerDecade=10, timePointsPerDecade=100

  ! Module-scope time used in integrals.
  double precision                            :: time__
  !$omp threadprivate(time__)

  ! Cached copies of tabulated solutions. These are used to avoid re-reading from file if the same variance is requested multiple times.
  type :: cachedVariance
     type            (varying_string)                              :: fileName
     double precision                                              :: sigma8Value           , sigmaNormalization         , &
          &                                                           massMinimum           , massMaximum                , &
          &                                                           timeMinimum           , timeMaximum                , &
          &                                                           timeMinimumLogarithmic, timeLogarithmicDeltaInverse
     double precision                , dimension(:  ), allocatable :: massTmp               , timesTmp
     double precision                , dimension(:,:), allocatable :: rootVarianceTmp       , rootVarianceUniqueTmp
     integer                         , dimension(:  ), allocatable :: uniqueSizeTmp
     integer                         , dimension(:,:), allocatable :: indexTmp
  end type cachedVariance
  
  integer                , parameter            :: sizeCache      =25
  integer                                       :: countCache     = 0, lastCache=0
  type   (cachedVariance), dimension(sizeCache) :: cachedVariances
  
contains

  function filteredPowerConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{cosmologicalMassVarianceFilteredPower} cosmological mass variance class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    use :: Error           , only : Error_Report
    implicit none
    type            (cosmologicalMassVarianceFilteredPower  )                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    type            (inputParameters                        )                :: referenceParameters
    class           (cosmologicalMassVarianceClass          ), pointer       :: cosmologicalMassVarianceReference
    class           (cosmologyParametersClass               ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    class           (powerSpectrumPrimordialTransferredClass), pointer       :: powerSpectrumPrimordialTransferred_     , powerSpectrumPrimordialTransferredReference
    class           (powerSpectrumWindowFunctionClass       ), pointer       :: powerSpectrumWindowFunction_            , powerSpectrumWindowFunctionTopHat_
    class           (linearGrowthClass                      ), pointer       :: linearGrowth_
    class           (transferFunctionClass                  ), pointer       :: transferFunction_
    double precision                                                         :: sigma8Value                             , tolerance                                  , &
         &                                                                      toleranceTopHat                         , wavenumberReference                        , &
         &                                                                      rootVarianceLogarithmicGradientTolerance, amplitudeScalar
    logical                                                                  :: monotonicInterpolation                  , nonMonotonicIsFatal                        , &
         &                                                                      truncateAtParticleHorizon               , storeTabulations                           , &
         &                                                                      integrationFailureIsFatal
    
    !![
    <objectBuilder    class="cosmologyParameters"                name="cosmologyParameters_"                        source="parameters"                                                  />
    <objectBuilder    class="cosmologyFunctions"                 name="cosmologyFunctions_"                         source="parameters"                                                  />
    <objectBuilder    class="powerSpectrumPrimordialTransferred" name="powerSpectrumPrimordialTransferred_"         source="parameters"                                                  />
    <objectBuilder    class="powerSpectrumWindowFunction"        name="powerSpectrumWindowFunction_"                source="parameters"                                                  />
    <objectBuilder    class="linearGrowth"                       name="linearGrowth_"                               source="parameters"                                                  />
    <objectBuilder    class="transferFunction"                   name="transferFunction_"                           source="parameters"                                                  />
    !!]
    if (parameters%isPresent('powerSpectrumWindowFunctionTopHat')) then
       !![
       <objectBuilder class="powerSpectrumWindowFunction"        name="powerSpectrumWindowFunctionTopHat_"          source="parameters" parameterName="powerSpectrumWindowFunctionTopHat"/>
       !!]
    else
       nullify(powerSpectrumWindowFunctionTopHat_)
    end if
    if      (parameters%isPresent('wavenumberReference')) then
       if (              parameters%isPresent('sigma_8'                                                )) call Error_Report('sigma_8 must not be specified if a reference wavenumber is specified'                //{introspection:location})
       if (              parameters%isPresent('amplitudeScalar'                                        )) call Error_Report('amplitudeScalar must not be specified if a reference wavenumber is specified'        //{introspection:location})
       if (.not.         parameters%isPresent('reference'                         ,requireValue=.false.)) call Error_Report('parameters must contain a "reference" section if a reference wavenumber is specified'//{introspection:location})
       referenceParameters=parameters%subParameters('reference',requireValue=.false.)
       if (.not.referenceParameters%isPresent('cosmologicalMassVariance'                               )) call Error_Report('"reference" section must explicitly defined a "cosmologicalMassVariance"'            //{introspection:location})
       if (.not.referenceParameters%isPresent('powerSpectrumPrimordialTransferred'                     )) call Error_Report('"reference" section must explicitly defined a "powerSpectrumPrimordialTransferred"'  //{introspection:location})
       !![
       <objectBuilder class="cosmologicalMassVariance"           name="cosmologicalMassVarianceReference"           source="referenceParameters"                                         />
       <objectBuilder class="powerSpectrumPrimordialTransferred" name="powerSpectrumPrimordialTransferredReference" source="referenceParameters"                                         />
       <inputParameter>
         <name>wavenumberReference</name>
         <source>parameters</source>
         <description>The reference wavenumber at which the amplitude of the power spectrum is matched to that in the reference model.</description>
       </inputParameter>
       !!]
    else if (parameters%isPresent('amplitudeScalar')) then
       if (parameters%isPresent('sigma_8')) call Error_Report('sigma_8 must not be specified if a power spectrum amplitude is specified'//{introspection:location})
       !![
       <inputParameter>
         <name>amplitudeScalar</name>
         <source>parameters</source>
         <description>The amplitude of the primordial scalar power spectrum, $A_\mathrm{s}$, such that $P_\chi(k) = A_\mathrm{s} (k/k_\mathrm{s0})^{n_\mathrm{s}-1}$ with $k_\mathrm{s0}=0.05$~Mpc$^{-1}$.</description>
       </inputParameter>
       !!]
    else       
       !![
       <inputParameter>
         <name>sigma_8</name>
         <source>parameters</source>
         <variable>sigma8Value</variable>
         <defaultValue>0.8111d0</defaultValue>
         <defaultSource>(\citealt{planck_collaboration_planck_2018}; TT,TE,EE$+$lowE$+$lensing)</defaultSource>
         <description>The fractional mass fluctuation in the linear density field at the present day in spheres of radius 8~Mpc/h.</description>
       </inputParameter>
       !!]
    end if
    !![
    <inputParameter>
      <name>toleranceTopHat</name>
      <source>parameters</source>
      <defaultValue>1.0d-6</defaultValue>
      <description>The relative tolerance to use in integrating over the linear power spectrum using a top-hat (real space) window function to compute the cosmological mass variance.</description>
    </inputParameter>
    <inputParameter>
      <name>tolerance</name>
      <source>parameters</source>
      <defaultValue>4.0d-6</defaultValue>
      <description>The relative tolerance to use in integrating over the linear power spectrum to compute the cosmological mass variance.</description>
    </inputParameter>
    <inputParameter>
      <name>nonMonotonicIsFatal</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true any non-monotonicity in the tabulated $\sigma(M)$ is treated as a fatal error. Otherwise a only a warning is issued.</description>
    </inputParameter>
    <inputParameter>
      <name>integrationFailureIsFatal</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true any failed integrals when evaluating $\sigma(M)$ are treated as fatal errors. Otherwise a only a warning is issued.</description>
    </inputParameter>
    <inputParameter>
      <name>monotonicInterpolation</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true use a monotonic cubic spline interpolator to interpolate in the $\sigma(M)$ table. Otherwise use a standard cubic spline interpolator. Use of the monotonic interpolator can be helpful is $\sigma(M)$ must be strictly monotonic but becomes a very weak function of $M$ at low masses.</description>
    </inputParameter>
    <inputParameter>
      <name>rootVarianceLogarithmicGradientTolerance</name>
      <source>parameters</source>
      <defaultValue>1.0d-12</defaultValue>
      <description>The tolerance in $\mathrm{d}\log\sigma/\mathrm{d}\log M$ to allow before reporting errors when monotonic interpolation is used..</description>
    </inputParameter>
    <inputParameter>
      <name>truncateAtParticleHorizon</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true then integration over the power spectrum is truncated at a wavenumber $k=1/H_\mathrm{p}(t_0)$, where $H_\mathrm{p}(t_0)$ is the comoving distance to the particle horizon at the present epoch. Otherwise, integration continues to $k=0$.</description>
    </inputParameter>
    <inputParameter>
      <name>storeTabulations</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true then tabulated $\sigma(M)$ results are stored to file for future re-use.</description>
    </inputParameter>
    <conditionalCall>
     <call>
      self=filteredPowerConstructorInternal(                                                                                   &amp;
       &amp;                                tolerance                               =tolerance                               , &amp;
       &amp;                                toleranceTopHat                         =toleranceTopHat                         , &amp;
       &amp;                                nonMonotonicIsFatal                     =nonMonotonicIsFatal                     , &amp;
       &amp;                                integrationFailureIsFatal               =integrationFailureIsFatal               , &amp;
       &amp;                                monotonicInterpolation                  =monotonicInterpolation                  , &amp;
       &amp;                                rootVarianceLogarithmicGradientTolerance=rootVarianceLogarithmicGradientTolerance, &amp;
       &amp;                                truncateAtParticleHorizon               =truncateAtParticleHorizon               , &amp;
       &amp;                                storeTabulations                        =storeTabulations                        , &amp;
       &amp;                                cosmologyParameters_                    =cosmologyParameters_                    , &amp;
       &amp;                                cosmologyFunctions_                     =cosmologyFunctions_                     , &amp;
       &amp;                                linearGrowth_                           =linearGrowth_                           , &amp;
       &amp;                                transferFunction_                       =transferFunction_                       , &amp;
       &amp;                                powerSpectrumPrimordialTransferred_     =powerSpectrumPrimordialTransferred_     , &amp;
       &amp;                                powerSpectrumWindowFunction_            =powerSpectrumWindowFunction_              &amp;
       &amp;                                {conditions}                                                                       &amp;
       &amp;                               )
     </call>
     <argument name="sigma8"                                      value="sigma8Value"                                 condition=".not.parameters%isPresent('wavenumberReference').and..not.parameters%isPresent('amplitudeScalar')"/>
     <argument name="amplitudeScalar"                             value="amplitudeScalar"                             condition="                                                          parameters%isPresent('amplitudeScalar')"/>
     <argument name="cosmologicalMassVarianceReference"           value="cosmologicalMassVarianceReference"           condition="     parameters%isPresent('wavenumberReference')"                                                 />
     <argument name="wavenumberReference"                         value="wavenumberReference"                         condition="     parameters%isPresent('wavenumberReference')"                                                 />
     <argument name="powerSpectrumPrimordialTransferredReference" value="powerSpectrumPrimordialTransferredReference" condition="     parameters%isPresent('wavenumberReference')"                                                 />
     <argument name="powerSpectrumWindowFunctionTopHat_"          value="powerSpectrumWindowFunctionTopHat_"          parameterPresent="parameters" parameterName="powerSpectrumWindowFunctionTopHat"                              />
    </conditionalCall>
    <inputParametersValidate source="parameters" extraAllowedNames="reference"/>
    !!]
    if (parameters%isPresent('wavenumberReference')) then
       !![
       <objectDestructor name="powerSpectrumPrimordialTransferredReference"/>
       <objectDestructor name="cosmologicalMassVarianceReference"          />
       !!]
    end if
    !![
    <objectDestructor    name="cosmologyParameters_"                       />
    <objectDestructor    name="cosmologyFunctions_"                        />
    <objectDestructor    name="linearGrowth_"                              />
    <objectDestructor    name="powerSpectrumPrimordialTransferred_"        />
    <objectDestructor    name="powerSpectrumWindowFunction_"               />
    <objectDestructor    name="transferFunction_"                          />
    !!]
    if (parameters%isPresent('powerSpectrumWindowFunctionTopHat')) then
       !![
       <objectDestructor name="powerSpectrumWindowFunctionTopHat_"         />
       !!]
    end if
    return
  end function filteredPowerConstructorParameters

  function filteredPowerConstructorInternal(sigma8,amplitudeScalar,cosmologicalMassVarianceReference,powerSpectrumPrimordialTransferredReference,wavenumberReference,tolerance,toleranceTopHat,nonMonotonicIsFatal,integrationFailureIsFatal,monotonicInterpolation,rootVarianceLogarithmicGradientTolerance,truncateAtParticleHorizon,storeTabulations,cosmologyParameters_,cosmologyFunctions_,linearGrowth_,transferFunction_,powerSpectrumPrimordialTransferred_,powerSpectrumWindowFunction_,powerSpectrumWindowFunctionTopHat_) result(self)
    !!{
    Internal constructor for the \refClass{cosmologicalMassVarianceFilteredPower} linear growth class.
    !!}
    use :: File_Utilities                 , only : Directory_Make                   , File_Path
    use :: Error                          , only : Error_Report
    use :: Input_Paths                    , only : inputPath                        , pathTypeDataDynamic
    use :: Power_Spectrum_Window_Functions, only : powerSpectrumWindowFunctionTopHat
    use :: Interfaces_CLASS               , only : Interface_CLASS_Normalization
    use :: Error                          , only : errorStatusSuccess
    use :: Numerical_Constants_Math       , only : Pi
    implicit none
    type            (cosmologicalMassVarianceFilteredPower  )                                  :: self
    double precision                                         , intent(in   )                   :: tolerance                                  , toleranceTopHat       , &
         &                                                                                        rootVarianceLogarithmicGradientTolerance
    double precision                                         , intent(in   )        , optional :: wavenumberReference                        , sigma8                , &
         &                                                                                        amplitudeScalar
    logical                                                  , intent(in   )                   :: nonMonotonicIsFatal                        , monotonicInterpolation, &
         &                                                                                        truncateAtParticleHorizon                  , storeTabulations      , &
         &                                                                                        integrationFailureIsFatal
    class           (cosmologyParametersClass               ), intent(in   ), target           :: cosmologyParameters_
    class           (cosmologyFunctionsClass                ), intent(in   ), target           :: cosmologyFunctions_
    class           (powerSpectrumPrimordialTransferredClass), intent(in   ), target           :: powerSpectrumPrimordialTransferred_
    class           (powerSpectrumWindowFunctionClass       ), intent(in   ), target           :: powerSpectrumWindowFunction_
    class           (powerSpectrumPrimordialTransferredClass), intent(in   ), target, optional :: powerSpectrumPrimordialTransferredReference
    class           (cosmologicalMassVarianceClass          ), intent(in   ), target, optional :: cosmologicalMassVarianceReference
    class           (powerSpectrumWindowFunctionClass       ), intent(in   ), target, optional :: powerSpectrumWindowFunctionTopHat_
    class           (linearGrowthClass                      ), intent(in   ), target           :: linearGrowth_
    class           (transferFunctionClass                  ), intent(in   ), target, optional :: transferFunction_ 
    double precision                                                                           :: halfModeMass
    integer                                                                                    :: status
    !![
    <constructorAssign variables="tolerance, toleranceTopHat, nonMonotonicIsFatal, integrationFailureIsFatal, monotonicInterpolation, rootVarianceLogarithmicGradientTolerance, truncateAtParticleHorizon, storeTabulations, *cosmologyParameters_, *cosmologyFunctions_, *linearGrowth_, *transferFunction_, *powerSpectrumPrimordialTransferred_, *powerSpectrumWindowFunction_, *powerSpectrumWindowFunctionTopHat_"/>
    !!]

    if (.not.present(powerSpectrumWindowFunctionTopHat_)) then
       allocate(powerSpectrumWindowFunctionTopHat :: self%powerSpectrumWindowFunctionTopHat_)
       select type (powerSpectrumWindowFunctionTopHat__ => self%powerSpectrumWindowFunctionTopHat_)
       type is (powerSpectrumWindowFunctionTopHat)
          !![
          <referenceConstruct isResult="yes" object="powerSpectrumWindowFunctionTopHat__" constructor="powerSpectrumWindowFunctionTopHat(cosmologyParameters_)"/>  
          !!]
       end select
    end if
    if      (present(sigma8         )) then
       if (     present(wavenumberReference).or.     present(cosmologicalMassVarianceReference).or.     present(powerSpectrumPrimordialTransferredReference)) call Error_Report('sigma8 is specified, can not also specify matched power spectrum'//{introspection:location})
       if (     present(amplitudeScalar    )                                                                                                                ) call Error_Report('sigma8 is specified, can not also specify scalar amplitude'      //{introspection:location})
       self%sigma8Value                                 =  sigma8
       self%normalizationSigma8                         =  .true.
       self%cosmologicalMassVarianceReference           => null()
       self%powerSpectrumPrimordialTransferredReference => null()
    else if (present(amplitudeScalar)) then
       if (     present(wavenumberReference).or.     present(cosmologicalMassVarianceReference).or.     present(powerSpectrumPrimordialTransferredReference)) call Error_Report('sigma8 is specified, can not also specify matched power spectrum'//{introspection:location})
       self%amplitudeScalar                             =  amplitudeScalar
       self%sigma8Value                                 =  sqrt(amplitudeScalar*Interface_CLASS_Normalization(self%cosmologyParameters_))
       self%normalizationSigma8                         =  .true.
       self%cosmologicalMassVarianceReference           => null()
       self%powerSpectrumPrimordialTransferredReference => null()
    else
       if (.not.present(wavenumberReference).or..not.present(cosmologicalMassVarianceReference).or..not.present(powerSpectrumPrimordialTransferredReference)) call Error_Report('sigma8 is not specified, must specify matched power spectrum'    //{introspection:location})
       self%sigma8Value                                 =  -1.0d0 
       self%normalizationSigma8                         =  .false.
       self%wavenumberReference                         =  wavenumberReference
       !![
       <referenceAcquire isResult="yes" owner="self" target="cosmologicalMassVarianceReference"           source="cosmologicalMassVarianceReference"          />
       <referenceAcquire isResult="yes" owner="self" target="powerSpectrumPrimordialTransferredReference" source="powerSpectrumPrimordialTransferredReference"/>
       !!]
    end if
    self%initialized           =.false.
    self%warnedNonIncreasing   =.false.
    self%growthIsMassDependent_=self%powerSpectrumPrimordialTransferred_%growthIsWavenumberDependent()
    self%fileName              =inputPath(pathTypeDataDynamic)                                                       // &
         &                      'largeScaleStructure/'                                                               // &
         &                      self%objectType      (                                                              )// &
         &                      '_'                                                                                  // &
         &                      self%hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)// &
         &                      '.hdf5'
    call Directory_Make(File_Path(self%fileName))
    if (present(transferFunction_)) then
       halfModeMass              =self%transferFunction_%halfModeMass(status)
       if (status == errorStatusSuccess) then
          self%wavenumberHalfMode=+(                                             &
               &                    +4.0d0                                       &
               &                    *Pi                                          &
               &                    /3.0d0                                       &
               &                    *self%cosmologyParameters_%OmegaMatter    () &
               &                    *self%cosmologyParameters_%densityCritical() &
               &                    /halfModeMass                                &
               &                   )**(1.0d0/3.0d0)                              &
               &                  *Pi
       else
          self%wavenumberHalfMode=-1.0d0
       end if
    else
       self%wavenumberHalfMode=-1.0d0
    end if
    return
  end function filteredPowerConstructorInternal

  subroutine filteredPowerDestructor(self)
    !!{
    Destructor for the \refClass{cosmologicalMassVarianceFilteredPower} linear growth class.
    !!}
    implicit none
    type   (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    integer                                                       :: i

    !![
    <objectDestructor name="self%cosmologyParameters_"               />
    <objectDestructor name="self%cosmologyFunctions_"                />
    <objectDestructor name="self%linearGrowth_"                      />
    <objectDestructor name="self%powerSpectrumPrimordialTransferred_"/>
    <objectDestructor name="self%powerSpectrumWindowFunction_"       />
    <objectDestructor name="self%powerSpectrumWindowFunctionTopHat_" />
    <objectDestructor name="self%transferFunction_"                  />
    !!]
    if (.not.self%normalizationSigma8) then
       !![
       <objectDestructor name="self%powerSpectrumPrimordialTransferredReference"/>
       <objectDestructor name="self%cosmologicalMassVarianceReference"          />
       !!]
    end if
    if (self%initialized) then
       do i=1,size(self%rootVarianceTable)
          call self%rootVarianceTable(i)%destroy()
       end do
    end if
    return
  end subroutine filteredPowerDestructor

  double precision function filteredPowerPowerNormalization(self)
    !!{
    Return the normalization of the power spectrum.
    !!}
    implicit none
    class(cosmologicalMassVarianceFilteredPower), intent(inout) :: self

    call self%retabulate()
    filteredPowerPowerNormalization=self%sigmaNormalization**2
    return
  end function filteredPowerPowerNormalization

  double precision function filteredPowerSigma8(self)
    !!{
    Return the value of $\sigma_8$.
    !!}
    implicit none
    class(cosmologicalMassVarianceFilteredPower), intent(inout) :: self

    filteredPowerSigma8=self%sigma8Value
    return
  end function filteredPowerSigma8

  double precision function filteredPowerRootVariance(self,mass,time)
    !!{
    Return the root-variance of the cosmological density field in a spherical region containing the given {\normalfont
    \ttfamily mass} on average.
    !!}
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: mass, time
    double precision                                                       :: h
    integer                                                                :: i

    call self%retabulate      (mass,time)
    if (self%growthIsMassDependent_) then
       call self%interpolantsTime(time,i,h)
       filteredPowerRootVariance=+self%rootVarianceTable(i  )%interpolate(mass)*(1.0d0-h) &
            &                    +self%rootVarianceTable(i+1)%interpolate(mass)*       h
    else
       filteredPowerRootVariance=+self%rootVarianceTable(1  )%interpolate(mass)           &
            &                    *self%linearGrowth_         %value      (time)
    end if
    return
  end function filteredPowerRootVariance

  double precision function filteredPowerRootVarianceLogarithmicGradient(self,mass,time)
    !!{
    Return the logarithmic gradient with respect to mass of the root-variance of the cosmological density field in a spherical
    region containing the given {\normalfont \ttfamily mass} on average.
    !!}
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: mass        , time
    double precision                                                       :: rootVariance

    call self%rootVarianceAndLogarithmicGradient(mass,time,rootVariance,filteredPowerRootVarianceLogarithmicGradient)
    return
  end function filteredPowerRootVarianceLogarithmicGradient

  double precision function filteredPowerRootVarianceLogarithmicGradientTime(self,mass,time)
    !!{
    Return the logarithmic gradient with respect to time of the root-variance of the cosmological density field in a spherical
    region containing the given {\normalfont \ttfamily mass} on average.
    !!}
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: mass, time
    double precision                                                       :: h
    integer                                                                :: i

    call self%retabulate(mass,time)
    if (self%growthIsMassDependent_) then
       call self%interpolantsTime(time,i,h)
       filteredPowerRootVarianceLogarithmicGradientTime=+(                                                         &
            &                                             -self%rootVarianceTable(i  )%interpolate(mass)           &
            &                                             +self%rootVarianceTable(i+1)%interpolate(mass)           &
            &                                            )                                                         &
            &                                           /(                                                         &
            &                                             +self%rootVarianceTable(i  )%interpolate(mass)*(1.0d0-h) &
            &                                             +self%rootVarianceTable(i+1)%interpolate(mass)*       h  &
            &                                            )                                                         &
            &                                           *self%timeLogarithmicDeltaInverse
    else
       filteredPowerRootVarianceLogarithmicGradientTime=+self%linearGrowth_       %logarithmicDerivativeExpansionFactor (time) &
            &                                           *self%cosmologyFunctions_ %expansionRate                       (       &
            &                                            self%cosmologyFunctions_ %expansionFactor                      (time) &
            &                                                                                                          )       &
            &                                           *                                                                time
    end if
    return
  end function filteredPowerRootVarianceLogarithmicGradientTime

  subroutine filteredPowerRootVarianceAndLogarithmicGradient(self,mass,time,rootVariance,rootVarianceLogarithmicGradient)
    !!{
    Return the value and logarithmic gradient with respect to mass of the root-variance of the cosmological density field in a
    spherical region containing the given {\normalfont \ttfamily mass} on average.
    !!}
    use :: Display                 , only : displayGreen, displayBlue, displayYellow, displayReset
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    use :: HDF5_Access             , only : hdf5Access
    use :: IO_HDF5                 , only : hdf5Object
    use :: String_Handling         , only : operator(//)
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: mass         , time
    double precision                                       , intent(  out) :: rootVariance , rootVarianceLogarithmicGradient
    type            (hdf5Object                           ), save          :: errorFile
    type            (varying_string                       ), save          :: errorFileName
    !$omp threadprivate(errorFile,errorFileName)
    character       (len=1                                )                :: label
    double precision                                                       :: wavenumber   , rootVarianceGradient           , &
         &                                                                    interpolant  , h                              , &
         &                                                                    linearGrowth
    integer                                                                :: i            , j
    
    call self%retabulate(mass,time)
    if (self%growthIsMassDependent_) then
       call self%interpolantsTime(time,i,h)
    else
       i=1
       h=0.0d0
    end if
    rootVariance        =0.0d0
    rootVarianceGradient=0.0d0
    do j=i,i+1
       if (j == i) then
          interpolant=1.0d0-h
       else
          interpolant=      h
       end if
       if (interpolant == 0.0d0) cycle
       if (self%powerSpectrumWindowFunction_%amplitudeIsMassIndependent()) then
          ! For the case of a constant window function amplitude the logarithmic gradient can be found analytically.
          wavenumber          =+self%powerSpectrumWindowFunction_       %wavenumberMaximum(                mass         )
          rootVarianceGradient=+rootVarianceGradient                                                                         &
               &               -self%powerSpectrumPrimordialTransferred_%power            (wavenumber,self%times(j)     )    &
               &               *self%sigmaNormalization                                                                  **2 &
               &               *self%powerSpectrumWindowFunction_       %value            (wavenumber,     mass    ,time)**2 &
               &               *                                                           wavenumber                    **3 &
               &               /12.0d0                                                                                       &
               &               /Pi                                                                                       **2 &
               &               *interpolant
       else
          ! Compute the gradient by interpolation in the tabulated relation.
          rootVarianceGradient=+rootVarianceGradient                                &
               &               +self%rootVarianceTable(j)%interpolateGradient(mass) &
               &               *                                              mass  &
               &               *interpolant
       end if
       rootVariance=+rootVariance                                &
            &       +self%rootVarianceTable(j)%interpolate(mass) &
            &       *interpolant
    end do
    if (rootVariance /= 0.0d0) then
       if (self%powerSpectrumWindowFunction_%amplitudeIsMassIndependent()) then
          rootVarianceLogarithmicGradient=+rootVarianceGradient    &
               &                          /rootVariance        **2
       else
          rootVarianceLogarithmicGradient=+rootVarianceGradient    &
               &                          /rootVariance
       end if
    else
       if (rootVarianceGradient == 0.0d0) then
          rootVarianceLogarithmicGradient=0.0d0
       else
          call Error_Report('dσ/dM ≠ 0 but σ = 0 - can not compute dlogσ/dlogM'//{introspection:location})
       end if
    end if
    ! Scale by the linear growth factor if growth is not mass-dependent.
    if (.not.self%growthIsMassDependent_) rootVariance=+rootVariance                   &
         &                                             *self%linearGrowth_%value(time)
    ! Validate the logarithmic gradient.
    if (rootVarianceLogarithmicGradient > 0.0d0) then
       if (rootVarianceLogarithmicGradient < self%rootVarianceLogarithmicGradientTolerance) then
          ! Ignore small positive gradients which can occur due to rounding errors.
          rootVarianceLogarithmicGradient=0.0d0
          return
       end if
       ! Logarithmic gradient is positive, which should not happen.
       if (self%monotonicInterpolation) then
          ! Monotonic interpolation is being used - a positive logarithmic gradient should be impossible.
          ! Write a file with the interpolating data.
          if (.not.self%powerSpectrumWindowFunction_%amplitudeIsMassIndependent()) then
             errorFileName=self%fileName//".error."//GetPID()
             !$ call hdf5Access%set()
             call    errorFile%openFile(char(errorFileName),overWrite=.true.,objectsOverwritable=.true.)
             call    errorFile%writeAttribute(mass                           ,'mass'                           )
             call    errorFile%writeAttribute(rootVariance                   ,'rootVariance'                   )
             call    errorFile%writeAttribute(rootVarianceGradient           ,'rootVarianceGradient'           )
             call    errorFile%writeAttribute(rootVarianceLogarithmicGradient,'rootVarianceLogarithmicGradient')
             if (.not.self%growthIsMassDependent_) then
                linearGrowth=self%linearGrowth_%value(time)
                call errorFile%writeAttribute(linearGrowth                   ,'linearGrowth'                   )
             end if
             do j=i,i+1
                if (j == i) then
                   interpolant=1.0d0-h
                else
                   interpolant=      h
                end if
                if (interpolant == 0.0d0) cycle
                write (label,'(i1)') j-i+1
                rootVarianceGradient=self%rootVarianceTable(j)%interpolateGradient(mass)
                call errorFile%writeAttribute(     interpolant                 ,'interpolant'         //trim(adjustl(label)))
                call errorFile%writeAttribute(     rootVarianceGradient        ,'rootVarianceGradient'//trim(adjustl(label)))
                call errorFile%writeDataset  (self%rootVarianceTable   (j)%xs(),'massTable'           //trim(adjustl(label)))
                call errorFile%writeDataset  (self%rootVarianceTable   (j)%ys(),'rootVarianceTable'   //trim(adjustl(label)))
             end do
          end if
          call errorFile%close()
          !$ call hdf5Access%unset()
          call Error_Report('dlogσ/dlogM > 0 detected, but monotonic interpolation was used - this should not happen'//char(10)//'  table data written to:'//char(10)//'    '//char(errorFileName)//{introspection:location})
       else
          ! Recommend that monotonic interpolation be used.
          call Error_Report(                                                                                                                                                                                                                  &
               &            'dlogσ/dlogM > 0 detected'//char(10)//                                                                                                                                                                            &
               &            displayGreen()//'HELP:'//displayReset()//' set <'//displayBlue()//'monotonicInterpolation'//displayReset()//' '//displayYellow()//'value'//displayReset()//'='//displayGreen()//'"true"'//displayReset()//'/> '// &
               &            'to ensure monotonic interpolation of tabulated σ(M) function'//{introspection:location}                                                                                                                          &
               &           )
       end if
    end if
    return
  end subroutine filteredPowerRootVarianceAndLogarithmicGradient

  double precision function filteredPowerMass(self,rootVariance,time)
    !!{
    Return the mass corresponding to the given {\normalfont \ttfamily } root-variance of the cosmological density field.
    !!}
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: rootVariance      , time
    double precision                                                       :: rootVarianceActual
    double precision                                                       :: h                 , hTime     , &
         &                                                                    interpolantTime
    integer                                                                :: i                 , iBoundLeft, &
         &                                                                    iBoundRight       , k         , &
         &                                                                    j

    if (self%growthIsMassDependent_) then
       rootVarianceActual=rootVariance
    else
       rootVarianceActual=rootVariance/self%linearGrowth_%value(time)
    end if
    ! If the requested root-variance is below the lowest value tabulated, attempt to tabulate to higher mass (lower
    ! root-variance).
    call self%retabulate(time=time)
    do while (rootVarianceActual < self%rootVarianceTable(1)%y(-1))
       call self%retabulate(self%rootVarianceTable(1)%x(-1)*2.0d0,time)
    end do
    ! Get interpolants in time.
    if (self%growthIsMassDependent_) then
       call self%interpolantsTime(time,k,hTime)
    else
       k    =1
       hTime=0.0d0
    end if
    ! If σ exceeds the highest value tabulated, simply return the lowest tabulated mass.
    if (rootVarianceActual > self%rootVarianceTable(k)%y(1)) then
       filteredPowerMass=self%rootVarianceTable(k)%x(1)
    else
       ! Iterate over times.
       filteredPowerMass=0.0d0
       do j=k,k+1
          ! Compute interpolating factor.
          if (j == k) then
             interpolantTime=1.0d0-hTime
          else
             interpolantTime=      hTime
          end if
          if (interpolantTime == 0.0d0) cycle
          ! Find the largest mass corresponding to this σ.
          iBoundLeft =1
          iBoundRight=size(self%rootVarianceUniqueTable(j)%rootVariance)
          do while (iBoundLeft+1 < iBoundRight)
             i=int((iBoundLeft+iBoundRight)/2)
             if (self%rootVarianceUniqueTable(j)%rootVariance(i) < rootVarianceActual) then
                iBoundRight=i
             else
                iBoundLeft =i
             end if
          end do
          i                =self%rootVarianceUniqueTable(j)%index(iBoundRight)
          h                =+(     rootVarianceActual         -self%rootVarianceTable(j)%y(i)) &
               &            /(self%rootVarianceTable(j)%y(i-1)-self%rootVarianceTable(j)%y(i))
          filteredPowerMass=+filteredPowerMass                                    &
               &            +exp(                                                 &
               &                 +log(self%rootVarianceTable(j)%x(i  ))*(1.0d0-h) &
               &                 +log(self%rootVarianceTable(j)%x(i-1))*       h  &
               &                )                                                 &
               &            *interpolantTime
       end do
    end if
    return
  end function filteredPowerMass

  subroutine filteredPowerRetabulate(self,mass,time)
    !!{
    Tabulate the cosmological mass variance.
    !!}
    use :: Cosmology_Parameters    , only : hubbleUnitsLittleH
    use :: Display                 , only : displayIndent            , displayMagenta                   , displayMessage, displayReset, &
          &                                 displayUnindent          , verbosityLevelWorking
    use :: File_Utilities          , only : File_Lock                , File_Unlock                      , lockDescriptor
    use :: Error                   , only : Error_Report             , Warn
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Ranges        , only : Make_Range               , rangeTypeLogarithmic
    use :: Tables                  , only : table1DLogarithmicCSpline, table1DLogarithmicMonotoneCSpline
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout)               :: self
    double precision                                       , intent(in   ), optional     :: mass                      , time
    ! Radius for σ(M) normalization in Mpc/h.
    double precision                                       , parameter                   :: radiusNormalization =8.0d0
    integer                                                                              :: i                         , rootVarianceTableCount , &
         &                                                                                  j                         , rootVarianceUniqueCount, &
         &                                                                                  rootVarianceTimeCount     , k                      , &
         &                                                                                  countNewLower             , countNewUpper          , &
         &                                                                                  iMinimum
    double precision                                                                     :: sigma                     , smoothingMass          , &
         &                                                                                  massMinimum               , massMaximum            , &
         &                                                                                  timeMinimum               , timeMaximum            , &
         &                                                                                  sigmaMinimum
    logical                                                , allocatable  , dimension(:) :: rootVarianceIsUnique
    type            (varying_string                       ), save                        :: message
    character       (len=12                               )                              :: label                     , labelLow               , &
         &                                                                                  labelHigh                 , labelTarget
    ! The variable "message" is saved (and made threadprivate) as its destructor is expensive, and this function gets called a
    ! lot.
    !$omp threadprivate(message)

    if (self%remakeTable(mass,time)) then
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(self%fileName),self%fileLock,lockIsShared=.true.)
       call self%fileRead()
       call File_Unlock(self%fileLock,sync=.false.)
    end if
    if (self%remakeTable(mass,time)) then
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(self%fileName),self%fileLock,lockIsShared=.false.)
       ! Try again to read the file - another process/thread may have already created the file in which case we may not need to do so again.
       call self%fileRead()
       if (self%remakeTable(mass,time)) then
          ! Compute the mass at which the mass variance is normalized.
          smoothingMass=+(                                                               &
               &          +4.0d0                                                         &
               &          /3.0d0                                                         &
               &          *Pi                                                            &
               &         )                                                               &
               &        *  self%cosmologyParameters_%OmegaMatter    (                  ) &
               &        *  self%cosmologyParameters_%densityCritical(                  ) &
               &        *(                                                               &
               &          +radiusNormalization                                           &
               &          /self%cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH) &
               &         )**3
          ! Determine the normalization of the power spectrum.
          if (self%sigma8Value > 0.0d0) then
             ! σ₈ was given, so normalize to that.
             self%sigmaNormalization =+self%sigma8Value                                                                &
                  &                   /rootVariance(time_=self%cosmologyFunctions_%cosmicTime(1.0d0),useTopHat=.true.)
          else
             ! We must match power to a reference power spectrum at a specified wavenumber.
             !! Compute the normalization factor assuming linear scaling.
             self%sigmaNormalization=+sqrt(                                                                                                                                          &
                  &                        +self%cosmologicalMassVarianceReference          %powerNormalization(                                                                   ) &
                  &                        *self%powerSpectrumPrimordialTransferredReference%power             (self%wavenumberReference,self%cosmologyFunctions_%cosmicTime(1.0d0)) &
                  &                        /self%powerSpectrumPrimordialTransferred_        %power             (self%wavenumberReference,self%cosmologyFunctions_%cosmicTime(1.0d0)) &
                  &                       )
             !! Compute the value of our σ₈.
             self%sigma8Value=rootVariance(time_=self%cosmologyFunctions_%cosmicTime(1.0d0),useTopHat=.true.)*self%sigmaNormalization
          end if
          ! Find suitable range of masses to tabulate.
          if (present(mass)) then
             countNewLower=0
             countNewUpper=0
             if (self%initialized) then
                massMinimum     =min(mass/10.0d0,self%massMinimum)
                massMaximum     =max(mass*10.0d0,self%massMaximum)
             else
                self%massMinimum=    mass
                self%massMaximum=    mass
                massMinimum     =    mass/10.0d0
                massMaximum     =    mass*10.0d0
             end if
             ! Determine how many points the table must be extended by in each direction to span the new required range.
             if (self%massMinimum > massMinimum) countNewLower=int(+log10(self%massMinimum/massMinimum)*dble(pointsPerDecade)+1.0d0)
             if (self%massMaximum < massMaximum) countNewUpper=int(-log10(self%massMaximum/massMaximum)*dble(pointsPerDecade)+1.0d0)
             ! Adjust the limits of the table by an integer number of steps.
             self%massMinimum=self%massMinimum/10.0d0**(dble(countNewLower)/dble(pointsPerDecade))
             self%massMaximum=self%massMaximum*10.0d0**(dble(countNewUpper)/dble(pointsPerDecade))
          else if (.not.self%initialized) then
             ! No mass was given, but the tables are not initialized. Must provide some mass range.
             self%massMinimum=1.0d10
             self%massMaximum=1.0d15
          end if
          rootVarianceTableCount=int(                         &
               &                     +log10(                  &
               &                            +self%massMaximum &
               &                            /self%massMinimum &
               &                           )                  &
               &                     *dble(pointsPerDecade)   &
               &                    )
          ! Find suitable range of times to tabulate.
          if (self%growthIsMassDependent_) then
             if (present(time)) then
                countNewLower=0
                countNewUpper=0
                if (self%initialized) then
                   timeMinimum     =min(time/2.0d0,self%timeMinimum)
                   timeMaximum     =max(time*2.0d0,self%timeMaximum)
                else
                   self%timeMinimum=    time
                   self%timeMaximum=    time
                   timeMinimum     =    time/2.0d0
                   timeMaximum     =    time*2.0d0
                end if
                ! Determine how many points the table must be extended by in each direction to span the new required range.
                if (self%timeMinimum > timeMinimum) countNewLower=int(+log10(self%timeMinimum/timeMinimum)*dble(timePointsPerDecade)+1.0d0)
                if (self%timeMaximum < timeMaximum) countNewUpper=int(-log10(self%timeMaximum/timeMaximum)*dble(timePointsPerDecade)+1.0d0)
                ! Adjust the limits of the table by an integer number of steps.
                self%timeMinimum=self%timeMinimum/10.0d0**(dble(countNewLower)/dble(timePointsPerDecade))
                self%timeMaximum=self%timeMaximum*10.0d0**(dble(countNewUpper)/dble(timePointsPerDecade))
             else if (.not.self%initialized) then
                ! No time was given, but the tables are not initialized. Must provide some time range.
                self%timeMinimum=self%cosmologyFunctions_%cosmicTime(0.5d0)
                self%timeMaximum=self%cosmologyFunctions_%cosmicTime(1.0d0)
             end if
             rootVarianceTimeCount =int(                           &
                  &                     +log10(                    &
                  &                            +self%timeMaximum   &
                  &                            /self%timeMinimum   &
                  &                           )                    &
                  &                     *dble(timePointsPerDecade) &
                  &                    )
             self%timeMinimumLogarithmic     =                              log(                 self%timeMinimum)
             self%timeLogarithmicDeltaInverse=dble(rootVarianceTimeCount-1)/log(self%timeMaximum/self%timeMinimum)
          else
             ! Growth of the transferred power spectrum is independent of mass - we can therefore tabulate σ(M) at a single epoch
             ! and use the linear growth factor to transform it to other epochs.
             self%timeMinimum                =self%cosmologyFunctions_%cosmicTime(1.0d0)
             self%timeMaximum                =self%cosmologyFunctions_%cosmicTime(1.0d0)
             rootVarianceTimeCount           =1
             self%timeMinimumLogarithmic     =0.0d0
             self%timeLogarithmicDeltaInverse=0.0d0
          end if
          if (allocated(self%times                  )) deallocate(self%times                  )
          if (allocated(self%rootVarianceUniqueTable)) deallocate(self%rootVarianceUniqueTable)
          allocate(self%rootVarianceUniqueTable(rootVarianceTimeCount))
          allocate(self%times                  (rootVarianceTimeCount))
          if (self%growthIsMassDependent_) then
             self%times=Make_Range(self%timeMinimum,self%timeMaximum,rootVarianceTimeCount,rangeTypeLogarithmic)
          else
             self%times=self%timeMinimum
          end if
          ! Allocate table grid.
          if (allocated(self%rootVarianceTable)) then
             do i=1,size(self%rootVarianceTable)
                call self%rootVarianceTable(i)%destroy()
             end do
             deallocate(self%rootVarianceTable)
          end if
          if (self%monotonicInterpolation) then
             allocate(table1DLogarithmicMonotoneCSpline :: self%rootVarianceTable(rootVarianceTimeCount))
          else
             allocate(table1DLogarithmicCSpline         :: self%rootVarianceTable(rootVarianceTimeCount))
          end if
          call displayIndent("retabulating σ(M)",verbosityLevelWorking)
          write    (labelLow   ,'(e9.2)') self%massMinimum
          write    (labelHigh  ,'(e9.2)') self%massMaximum
          if (present(mass)) then
             write (labelTarget,'(e9.2)')      mass
          else
             labelTarget="unspecified"
          end if
          call displayMessage("mass range: "//labelLow//" < "//labelTarget//" < "//labelHigh//" M☉" ,verbosityLevelWorking)
          write    (labelLow   ,'(f9.4)') self%timeMinimum
          write    (labelHigh  ,'(f9.4)') self%timeMaximum
          if (present(time)) then
             write (labelTarget,'(f9.4)')      time
          else
             labelTarget="unspecified"
          end if
          call displayMessage("time range: "//labelLow//" < "//labelTarget//" < "//labelHigh//" Gyr",verbosityLevelWorking)
          do k=1,rootVarianceTimeCount
             call self%rootVarianceTable(k)%create(self%massMinimum,self%massMaximum,rootVarianceTableCount)
             allocate(rootVarianceIsUnique(rootVarianceTableCount))
             rootVarianceIsUnique=.true.
             ! Compute σ(M) at each tabulated point.
             massMinimum =-1.0d0
             sigmaMinimum=-1.0d0
             iMinimum    =-1
             do i=1,rootVarianceTableCount
                smoothingMass=+self        %rootVarianceTable(k)%x(                                    i)
                sigma        =+rootVariance                       (time_=self%times(k),useTopHat=.false.) &
                     &        *self%sigmaNormalization
                ! Enforce monotonicity.
                if (i > 1) then
                   if (sigma >= self%rootVarianceTable(k)%y(i-1)) then
                      iMinimum              =i
                      massMinimum            =smoothingMass
                      rootVarianceIsUnique(i)=.false.                      
                   end if
                   sigma=min(sigma,self%rootVarianceTable(k)%y(i-1))
                end if
                ! Store the value.
                call self%rootVarianceTable(k)%populate(sigma,i,computeSpline=(i == rootVarianceTableCount))
             end do
             ! Find unique values in the variance table.
             rootVarianceUniqueCount=count(rootVarianceIsUnique)
             allocate(self%rootVarianceUniqueTable(k)%rootVariance(rootVarianceUniqueCount))
             allocate(self%rootVarianceUniqueTable(k)%index       (rootVarianceUniqueCount))
             j=1
             do i=1,rootVarianceTableCount
                if (rootVarianceIsUnique(i)) then
                   self%rootVarianceUniqueTable(k)%rootVariance(j)=self%rootVarianceTable(k)%y(i)
                   self%rootVarianceUniqueTable(k)%index       (j)=i
                   j                                              =j+1
                end if
             end do
             deallocate(rootVarianceIsUnique)
             ! Abort or warn if σ(M) has no increase below some mass scale.
             if (massMinimum > 0.0d0) then
                if (self%nonMonotonicIsFatal) then
                   message=         ""
                else
                   message=         displayMagenta()//"WARNING: "//displayReset()
                end if
                write (label,'(e12.6)') massMinimum
                message=         "σ(M) is non-increasing below mass M="//label//"M☉"
                write (label,'(e12.6)') self%times(k)
                message=message//" at time t="//label//"Gyr."//char(10)
                message=message//"                 M  ="
                do i=max(1,iMinimum-2),min(rootVarianceTableCount,iMinimum+2)
                   write (label,'(e12.6)') self%rootVarianceTable(k)%x(i)
                   message=message//" "//label
                end do
                message=message//char(10)
                message=message//"   monotonized σ(M) ="
                do i=max(1,iMinimum-2),min(rootVarianceTableCount,iMinimum+2)
                   write (label,'(e12.6)') self%rootVarianceTable(k)%y(i)
                   message=message//" "//label
                end do
                message=message//char(10)
                message=message//"      original σ(M) ="
                do i=max(1,iMinimum-2),min(rootVarianceTableCount,iMinimum+2)
                   if (i == iMinimum) then
                      write (label,'(e12.6)') sigmaMinimum
                   else
                      write (label,'(a    )') repeat(" ",12)
                   end if
                   message=message//" "//label
                end do
                if (self%nonMonotonicIsFatal) then
                   call Error_Report(message//{introspection:location})
                else
                   message=message//char(10)//"         If problems occur consider not attempting to model structure below this mass scale."
                   if (.not.self%warnedNonIncreasing) then
                      call Warn     (message                          )
                      self%warnedNonIncreasing=.true.
                   end if
                end if
             end if
          end do
          call displayUnindent("done",verbosityLevelWorking)
          ! Table is now initialized.
          self%initialized=.true.
          ! Store file.
          call self%fileWrite()
       end if
       call File_Unlock(self%fileLock)
    end if
    return

  contains

    double precision function rootVariance(time_,useTopHat)
      !!{
      Compute the root-variance of mass in spheres enclosing the given {\normalfont \ttfamily mass} from the power spectrum.
      !!}
      use, intrinsic :: ISO_C_Binding           , only : c_size_t
      use            :: Display                 , only : displayIndent     , displayUnindent              , displayMessage, verbosityLevelStandard, &
           &                                             verbosityLevelWarn, enumerationVerbosityLevelType, displayReset  , displayGreen          , &
           &                                             displayBlue       , displayYellow
      use            :: Error                   , only : Error_Report      , Warn
      use            :: Interface_GSL           , only : GSL_EBadTol       , GSL_ETol                     , GSL_ERound    , GSL_Success           , &
           &                                             GSL_EMaxIter      , GSL_ESing                    , gslErrorDecode
      use            :: Numerical_Constants_Math, only : Pi
      use            :: Numerical_Integration   , only : GSL_Integ_Gauss15 , integrator
      use            :: Sorting                 , only : sort
      use            :: String_Handling         , only : operator(//)
      use            :: ISO_Varying_String      , only : var_str
      implicit none
      double precision                         , intent(in   ) :: time_
      logical                                  , intent(in   ) :: useTopHat
      double precision                         , parameter     :: wavenumberBAO                         =5.0d0 ! The wavenumber above which baryon acoustic oscillations are small - used to split the integral allowing the oscillating part to be handled robustly.
      integer         (c_size_t  )             , parameter     :: intervalsMaximum                      =10000_c_size_t
      double precision            , allocatable, dimension(:)  :: wavenumbers                                          , wavenumbersRestricted , &
           &                                                      wavenumbersLocalMinimaTransferFunction
      double precision                                         :: wavenumberMinimum                                    , wavenumberMaximum     , &
           &                                                      integrand                                            , integrandInterval     , &
           &                                                      wavenumberLower                                      , wavenumberUpper       , &
           &                                                      topHatRadius                                         , toleranceRelative
      integer                                                  :: i                                                    , status
      logical                                                  :: computeLogarithmically
      type            (integrator)                             :: integrator_                                          , integratorLogarithmic_

      time__      =time_
      topHatRadius=(                                             &
           &        +(                                           &
           &          +3.0d0                                     &
           &          /4.0d0                                     &
           &          /Pi                                        &
           &         )                                           &
           &        *smoothingMass                               &
           &        /self%cosmologyParameters_%OmegaMatter    () &
           &        /self%cosmologyParameters_%densityCritical() &
           &       )**(1.0d0/3.0d0)
      ! Find the minimum wavenumber for the integral.
      if (self%truncateAtParticleHorizon) then
         wavenumberMinimum=+1.0d0                                                                                                &
              &            /self%cosmologyFunctions_%distanceParticleHorizonComoving(self%cosmologyFunctions_%cosmicTime(1.0d0))
      else
         wavenumberMinimum=+0.0d0
      end if
      ! Find the maximum wavenumber for the integral (and establish the integrator).
      if (useTopHat) then
         wavenumberMaximum     =min(1.0d3/topHatRadius,self%powerSpectrumWindowFunctionTopHat_%wavenumberMaximum(smoothingMass))
         toleranceRelative     =                                                                +self%toleranceTopHat
         integrator_           =integrator(varianceIntegrandTopHat           ,toleranceRelative=+self%toleranceTopHat,integrationRule=GSL_Integ_Gauss15,intervalsMaximum=intervalsMaximum)
         integratorLogarithmic_=integrator(varianceIntegrandTopHatLogarithmic,toleranceRelative=+self%toleranceTopHat,integrationRule=GSL_Integ_Gauss15,intervalsMaximum=intervalsMaximum)
      else
         wavenumberMaximum     =min(1.0d3/topHatRadius,self%powerSpectrumWindowFunction_      %wavenumberMaximum(smoothingMass))
         toleranceRelative     =                                                                +self%tolerance
         integrator_           =integrator(varianceIntegrand                 ,toleranceRelative=+self%tolerance      ,integrationRule=GSL_Integ_Gauss15,intervalsMaximum=intervalsMaximum)
         integratorLogarithmic_=integrator(varianceIntegrandLogarithmic      ,toleranceRelative=+self%tolerance      ,integrationRule=GSL_Integ_Gauss15,intervalsMaximum=intervalsMaximum)
      end if
      ! Split the integral into subsets of the wavenumber range:
      !
      !!  1. The integral over the power spectrum is split at a wavenumber corresponding to the smallest scale at which BAO
      !!  features are significant (unless the upper limit of the integral is already below that wavenumber). This allows the
      !!  oscillatory part of the integral to be computed more accurately, without affecting the non-oscillatory part at larger
      !!  wavenumbers, and leads to an overall more accurate and robust determination of σ(M).
      !!
      !!  2. The integral is split at a modest multiple of the inverse top hat radius. We would expect the window function to cut
      !!  off above this scale, so it is useful to have a split just above this cut off.
      !!
      !!  3. For non-cold dark matter models with suppressed power spectrum on small scales, another splitting is done around the
      !!  half-mode wavenumber so that σ(M) is better resolved near the saturation value. This is particularly useful if the
      !!  transfer function oscillates on small scales.
      !!
      !!  4. Split at any local minima of the transfer function. This is useful for oscillating transfer functions to allow
      !!  integration to accurately capture each oscillation.
      call self%transferFunction_%wavenumbersLocalMinima(wavenumbersLocalMinimaTransferFunction)
      allocate(wavenumbers(4+size(wavenumbersLocalMinimaTransferFunction)))
      wavenumbers(1:4                )=[                                &
           &                                        wavenumberMinimum , &
           &                            +3.0d0*self%wavenumberHalfMode, &
           &                            +3.0d0/     topHatRadius      , &
           &                                        wavenumberBAO       &
           &                           ]
      wavenumbers(5:size(wavenumbers))=wavenumbersLocalMinimaTransferFunction
      ! Restrict the intervals to those between the minimum and maximum wavenumbers, and then sort them.
      wavenumbersRestricted=pack(wavenumbers,wavenumbers >= wavenumberMinimum .and. wavenumbers < wavenumberMaximum)
      call sort(wavenumbersRestricted)
      ! Iterate over intervals, accumulation the integral.
      integrand=0.0d0
      do i=1,size(wavenumbersRestricted)
         wavenumberLower=wavenumbersRestricted(i)
         if (i < size(wavenumbersRestricted)) then
            wavenumberUpper=wavenumbersRestricted(i+1)
         else
            wavenumberUpper=wavenumberMaximum
         end if
         integrandInterval=integrator_%integrate(wavenumberLower,wavenumberUpper,status=status)
         ! Decide if we need to try a logarithmic integral instead.
         if     (                         &
              &   (                       &
              &     status == GSL_EBadTol &
              &    .or.                   &
              &     status == GSL_ETol    &
              &    .or.                   &
              &     status == GSL_ERound  &
              &   )                       &
              &  .and.                    &
              &   wavenumberLower > 0.0d0 &
              & ) then
            ! The integration failed numerically, and the lower limit is non-zero, so try a logarithmic integral instead.
            computeLogarithmically=.true.
         else if (status /= GSL_Success) then
            ! Integration failed for some other reason, report an error.
            computeLogarithmically=.false.
            if (self%integrationFailureIsFatal) then
            block
              type     (varying_string) :: message
              character(len=12        ) :: label
              message=var_str('integration over interval failed [error number: ')//status//'; "'//gslErrorDecode(status)//'"]'//char(10)// &
                   &  displayGreen()//'HELP:'//displayReset()//' try increasing the value of <'//displayBlue()//'tolerance'
              if (useTopHat) then
                 message=message//"TopHat"
                 write (label,'(e12.6)') self%toleranceTopHat
              else
                 write (label,'(e12.6)') self%tolerance
              end if
              message=message//displayReset()//' '//displayYellow()//'value'//displayReset()//'='//displayGreen()//'"'//trim(adjustl(label))//'"'//displayReset()//'/> in the <'//displayBlue()//'cosmologicalMassVariance'//displayReset()//' '//displayYellow()//'value'//displayReset()//'='//displayGreen()//'"filteredPower"'//displayReset()//'> parameter'
              call Error_Report(message//{introspection:location})
            end block
            end if
         else if (integrandInterval <= 0.0d0 .and. wavenumberLower > 0.0d0) then
            ! Integration gave a zero result, and the lower limit is non-zero. This may be because the upper limit is large and
            ! the power is confined to small wavenumbers near the lower limit. This can happen, for example, if attempting to
            ! compute σ(M) for mass scales far below the cut off for power spectra with a small-scale cut off. In such cases
            ! attempt to evaluate the integral again, but integrating over log(wavenumber) such that more points in the integrand
            ! are placed at small wavenumber.
            computeLogarithmically=.true.
         else
            ! The integration was successful, no need to try a logarithmic integral.
            computeLogarithmically=.false.
         end if
         ! Recompute the integral using a logarithmic integral if necessary.
         if (computeLogarithmically) then
            integrandInterval=integratorLogarithmic_%integrate(log(wavenumberLower),log(wavenumberUpper),status=status)
            if ((status == GSL_EMaxIter .or. status == GSL_ERound) .and. integrandInterval < toleranceRelative*integrand) then
               ! Maximum iterations were exceeded, but the integrand is tiny - so proceed.
            else if (status /= GSL_Success) then
               if (self%integrationFailureIsFatal)                                                   &
                    & call Error_Report(var_str       ('variance integral failed [error number: ')// &
                    &                                  status                                     // &
                    &                           '; "'                                             // &
                    &                   gslErrorDecode(status                                    )// &
                    &                           '"]'                                              // &
                    &                   {introspection:location}                                     &
                    &                  )
            else if (integrandInterval <= 0.0d0) then
               ! The integrated power is still non-positive - check the integrand at the lower limit. If it is positive, we have a
               ! problem. Otherwise, the power in this interval seems to be actually zero, so proceed.
               if (varianceIntegrand(wavenumberLower) > 0.0d0) then
                  if (varianceIntegrand(wavenumberUpper) > 0.0d0) then
                     ! The integrand is non-zero at the upper limit also. A zero integral should therefore not be possible. First,
                     ! check if the integrand is so small that, when multiplied by the typical step size, it underflows to
                     ! zero. If it is, this is acceptable, and no error needs to be emitted. If such underflow did not occur, then
                     ! the integral should, definitely be positive. As it is not, this is an error situation.
                     if     (                                                                                             &
                          &    +max(                                                                                      &
                          &         +varianceIntegrandLogarithmic(log(wavenumberUpper)),                                  &
                          &         +varianceIntegrandLogarithmic(log(wavenumberLower))                                   &
                          &         )                                                                                     &
                          &    *    (                                                                                     &
                          &          +                            log(wavenumberUpper)                                    &
                          &          -                            log(wavenumberLower)                                    &
                          &         )                                                                                     &
                          &    /10.0d0                                                                                    &
                          &   >                                                                                           &
                          &    + 0.0d0                                                                                    &
                          &  .and.                                                                                        &
                          &   self%integrationFailureIsFatal                                                              &
                          & ) call Error_Report('no power in interval integrand - unexpected'//{introspection:location})
                  else
                     ! The integrand is zero at the upper limit - try reducing the upper limit until we get a non-zero result.
                     do while (wavenumberUpper > wavenumberLower .and. varianceIntegrand(wavenumberUpper) <= 0.0d0)
                        integrandInterval=integratorLogarithmic_%integrate(log(wavenumberLower),log(wavenumberUpper),status=status)
                        if (status == GSL_Success .and. integrandInterval > 0.0d0) exit                       
                        wavenumberUpper=sqrt(wavenumberUpper*wavenumberLower)
                     end do
                     ! If we still failed to get a non-zero power (or if the integral failed) - there's nothing more we can do.
                     if (integrandInterval <= 0.0d0 .or. status /= GSL_Success) then
                        block
                          type     (varying_string               ) :: message
                          character(len=13                       ) :: label
                          integer                                  :: j
                          type     (enumerationVerbosityLevelType) :: verbosityLevel
                          if (self%integrationFailureIsFatal) then
                             verbosityLevel=verbosityLevelStandard
                          else
                             verbosityLevel=verbosityLevelWarn
                          end if
                          call displayIndent ('Integration report:'  ,verbosityLevel)
                          call displayIndent ('Wavenumber intervals:',verbosityLevel)
                          call displayMessage('i    kₘᵢₙ,ᵢ  kₘₐₓ,ᵢ'    ,verbosityLevel)
                          do j=1,size(wavenumbersRestricted)
                             wavenumberLower=wavenumbersRestricted(j)
                             if (j < size(wavenumbersRestricted)) then
                                wavenumberUpper=wavenumbersRestricted(j+1)
                             else
                                wavenumberUpper=wavenumberMaximum
                             end if
                             write (label,'(i4)') j
                             message=trim(label)
                             write (label,'(e13.6)') wavenumberLower
                             message=message//" "//label
                             write (label,'(e13.6)') wavenumberUpper
                             message=message//" "//label
                             if (j == i) message=message//" ← failed"
                             call displayMessage(message,verbosityLevel)
                          end do
                          call displayUnindent('done' ,verbosityLevel)
                          write (label,'(e13.6)') topHatRadius
                          message="Top-hat radius integral = "//label
                          call displayMessage (message,verbosityLevel)
                          message="Current integral = "//label
                          call displayMessage (message,verbosityLevel)
                          call displayUnindent('done' ,verbosityLevel)
                        end block
                        if (self%integrationFailureIsFatal) then
                           call Error_Report('variance integration failed - unexpected'//{introspection:location})
                        else
                           call Warn        ('variance integration failed - unexpected'                          )
                        end if
                     end if
                  end if
               end if
            end if
         end if
         ! Accumulate the integral from this region.
         integrand=+integrand         &
              &    +integrandInterval
      end do
      ! Apply factors of 2 and π to compute the root variance.
      rootVariance=+integrand          &
           &       /2.0d0              &
           &       /Pi**2
      rootVariance=+sqrt(rootVariance)
      return
    end function rootVariance

    double precision function varianceIntegrand(wavenumber)
      !!{
      Integrand function used in computing the variance in (real space) top-hat spheres from the power spectrum.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber

      ! Return power spectrum multiplied by window function and volume element in k-space. Factors of 2 and π are included
      ! elsewhere.
      varianceIntegrand=+  self%powerSpectrumPrimordialTransferred_%power(wavenumber              ,time__) &
           &            *(                                                                                 &
           &              +self%powerSpectrumWindowFunction_       %value(wavenumber,smoothingMass,time__) &
           &              *                                               wavenumber                       &
           &             )**2
      return
    end function varianceIntegrand

    double precision function varianceIntegrandLogarithmic(wavenumberLogarithmic)
      !!{
      Integrand function used in computing the variance in (real space) top-hat spheres from the power spectrum. This version
      integrates with respect to $\log(k)$.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumberLogarithmic
      double precision                :: wavenumber

      ! Return power spectrum multiplied by window function and volume element in k-space. Factors of 2 and π are included
      ! elsewhere.
      wavenumber                  =+exp(wavenumberLogarithmic)
      varianceIntegrandLogarithmic=+  self%powerSpectrumPrimordialTransferred_%power(wavenumber              ,time__) &
           &                       *(                                                                                 &
           &                         +self%powerSpectrumWindowFunction_       %value(wavenumber,smoothingMass,time__) &
           &                         *                                               wavenumber                       &
           &                        )**2                                                                              &
           &                       *                                                 wavenumber
      return
    end function varianceIntegrandLogarithmic

    double precision function varianceIntegrandTopHat(wavenumber)
      !!{
      Integrand function used in computing the variance in (real space) top-hat spheres from the power spectrum.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber

      ! Return power spectrum multiplied by window function and volume element in k-space. Factors of 2 and π are included
      ! elsewhere.
      varianceIntegrandTopHat=+  self%powerSpectrumPrimordialTransferred_%power(wavenumber              ,time__) &
           &                  *(                                                                                 &
           &                    +self%powerSpectrumWindowFunctionTopHat_ %value(wavenumber,smoothingMass,time__) &
           &                    *                                               wavenumber                       &
           &                   )**2
      return
    end function varianceIntegrandTopHat

    double precision function varianceIntegrandTopHatLogarithmic(wavenumberLogarithmic)
      !!{
      Integrand function used in computing the variance in (real space) top-hat spheres from the power spectrum.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumberLogarithmic
      double precision                :: wavenumber

      ! Return power spectrum multiplied by window function and volume element in k-space. Factors of 2 and π are included
      ! elsewhere.
      wavenumber                        =+exp(wavenumberLogarithmic)
      varianceIntegrandTopHatLogarithmic=+  self%powerSpectrumPrimordialTransferred_%power(wavenumber              ,time__) &
           &                             *(                                                                                 &
           &                               +self%powerSpectrumWindowFunctionTopHat_ %value(wavenumber,smoothingMass,time__) &
           &                               *                                               wavenumber                       &
           &                              )**2                                                                              &
           &                             *                                                 wavenumber
      return
    end function varianceIntegrandTopHatLogarithmic

  end subroutine filteredPowerRetabulate

  subroutine filteredPowerInterpolantsTime(self,time,i,h)
    !!{
    Compute interpolants in time.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    integer                                                , intent(  out) :: i
    double precision                                       , intent(  out) :: h

    h=(log(time)-self%timeMinimumLogarithmic)*self%timeLogarithmicDeltaInverse+1.0d0
    i=  int (h)
    h=h-dble(i)
    if (i == size(self%times)) then
       ! Requested time must exactly equal the maximum tabulated time.
       i=size(self%times)-1
       h=1.0d0
    else if (i < 1) then
       call Error_Report('interpolant out of range'//{introspection:location})
    end if
    return
  end subroutine filteredPowerInterpolantsTime

  logical function filteredPowerGrowthIsMassDependent(self)
    !!{
    Return true if the growth rate of the variance is mass-dependent.
    !!}
    implicit none
    class(cosmologicalMassVarianceFilteredPower), intent(inout) :: self

    filteredPowerGrowthIsMassDependent=self%growthIsMassDependent_
    return
  end function filteredPowerGrowthIsMassDependent

  subroutine filteredPowerFileRead(self)
    !!{
    Read tabulated data on mass variance from file.
    !!}
    use :: Display       , only : displayMessage           , verbosityLevelWorking
    use :: File_Utilities, only : File_Exists
    use :: HDF5_Access   , only : hdf5Access
    use :: IO_HDF5       , only : hdf5Object
    use :: Tables        , only : table1DLogarithmicCSpline, table1DLogarithmicMonotoneCSpline
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout)               :: self
    double precision                                       , dimension(:  ), allocatable :: massTmp        , timesTmp
    double precision                                       , dimension(:,:), allocatable :: rootVarianceTmp, rootVarianceUniqueTmp
    integer                                                , dimension(:  ), allocatable :: uniqueSizeTmp
    integer                                                , dimension(:,:), allocatable :: indexTmp
    type            (hdf5Object                           )                              :: dataFile
    integer                                                                              :: i              , useCache

    !$omp critical(cosmologicalMassVarianceFilteredPowerCache)
    useCache=0
    if (countCache > 0) then
       do i=1,countCache
          if (cachedVariances(i)%fileName == self%fileName) then
             useCache=i
             exit
          end if
       end do
    end if    
    if (useCache /= 0) then
       timesTmp                        =cachedVariances(useCache)%timesTmp
       massTmp                         =cachedVariances(useCache)%massTmp
       rootVarianceTmp                 =cachedVariances(useCache)%rootVarianceTmp
       rootVarianceUniqueTmp           =cachedVariances(useCache)%rootVarianceUniqueTmp
       indexTmp                        =cachedVariances(useCache)%indexTmp
       uniqueSizeTmp                   =cachedVariances(useCache)%uniqueSizeTmp
       self%sigma8Value                =cachedVariances(useCache)%sigma8Value
       self%sigmaNormalization         =cachedVariances(useCache)%sigmaNormalization
       self%massMinimum                =cachedVariances(useCache)%massMinimum
       self%massMaximum                =cachedVariances(useCache)%massMaximum
       self%timeMinimum                =cachedVariances(useCache)%timeMinimum
       self%timeMaximum                =cachedVariances(useCache)%timeMaximum
       self%timeMinimumLogarithmic     =cachedVariances(useCache)%timeMinimumLogarithmic
       self%timeLogarithmicDeltaInverse=cachedVariances(useCache)%timeLogarithmicDeltaInverse
       self%initialized                =.true.
    end if
    !$omp end critical(cosmologicalMassVarianceFilteredPowerCache)
    if (useCache == 0) then
       ! Return if we are not using stored solutions.
       if (.not.self%storeTabulations) return
       ! Return if the file does not exist.
       if (.not.File_Exists(char(self%fileName))) return
       call displayMessage('reading σ(M) data from: '//self%fileName,verbosityLevelWorking)
       !$ call hdf5Access%set()
       call dataFile%openFile     (char(self%fileName)          ,overWrite                       =.false.,readOnly=.true.)
       call dataFile%readDataset  ('times'                      ,     timesTmp                                           )
       call dataFile%readDataset  ('mass'                       ,     massTmp                                            )
       call dataFile%readDataset  ('rootVariance'               ,     rootVarianceTmp                                    )
       call dataFile%readDataset  ('rootVarianceUnique'         ,     rootVarianceUniqueTmp                              )
       call dataFile%readDataset  ('indexUnique'                ,     indexTmp                                           )
       call dataFile%readDataset  ('uniqueSize'                 ,     uniqueSizeTmp                                      )
       call dataFile%readAttribute('sigma8'                     ,self%sigma8Value                                        )
       call dataFile%readAttribute('sigmaNormalization'         ,self%sigmaNormalization                                 )
       call dataFile%readAttribute('massMinimum'                ,self%massMinimum                                        )
       call dataFile%readAttribute('massMaximum'                ,self%massMaximum                                        )
       call dataFile%readAttribute('timeMinimum'                ,self%timeMinimum                                        )
       call dataFile%readAttribute('timeMaximum'                ,self%timeMaximum                                        )
       call dataFile%readAttribute('timeMinimumLogarithmic'     ,self%timeMinimumLogarithmic                             )
       call dataFile%readAttribute('timeLogarithmicDeltaInverse',self%timeLogarithmicDeltaInverse                        )
       call dataFile%close        (                                                                                      )
       !$ call hdf5Access%unset()
       ! Cache this variance for possible later reuse.
       !$omp critical(cosmologicalMassVarianceFilteredPowerCache)
       lastCache=lastCache+1
       if (lastCache > sizeCache) lastCache=1
       countCache=max(countCache,lastCache)
       cachedVariances(lastCache)%fileName                   =self%fileName
       cachedVariances(lastCache)%timesTmp                   =     timesTmp
       cachedVariances(lastCache)%massTmp                    =     massTmp
       cachedVariances(lastCache)%rootVarianceTmp            =     rootVarianceTmp
       cachedVariances(lastCache)%rootVarianceUniqueTmp      =     rootVarianceUniqueTmp
       cachedVariances(lastCache)%indexTmp                   =     indexTmp
       cachedVariances(lastCache)%uniqueSizeTmp              =     uniqueSizeTmp
       cachedVariances(lastCache)%sigma8Value                =self%sigma8Value
       cachedVariances(lastCache)%sigmaNormalization         =self%sigmaNormalization
       cachedVariances(lastCache)%massMinimum                =self%massMinimum
       cachedVariances(lastCache)%massMaximum                =self%massMaximum
       cachedVariances(lastCache)%timeMinimum                =self%timeMinimum
       cachedVariances(lastCache)%timeMaximum                =self%timeMaximum
       cachedVariances(lastCache)%timeMinimumLogarithmic     =self%timeMinimumLogarithmic
       cachedVariances(lastCache)%timeLogarithmicDeltaInverse=self%timeLogarithmicDeltaInverse
       !$omp end critical(cosmologicalMassVarianceFilteredPowerCache)
    end if
    if (allocated(self%times                  )) deallocate(self%times                  )
    if (allocated(self%rootVarianceTable      )) deallocate(self%rootVarianceTable      )
    if (allocated(self%rootVarianceUniqueTable)) deallocate(self%rootVarianceUniqueTable)
    allocate(self%times                  (size(timesTmp)))
    allocate(self%rootVarianceUniqueTable(size(timesTmp)))
    if (self%monotonicInterpolation) then
       allocate(table1DLogarithmicMonotoneCSpline :: self%rootVarianceTable(size(timesTmp)))
    else
       allocate(table1DLogarithmicCSpline         :: self%rootVarianceTable(size(timesTmp)))
    end if
    self%times=timesTmp
    do i=1,size(self%times)
       allocate(self%rootVarianceUniqueTable(i)%rootVariance(uniqueSizeTmp(i)))
       allocate(self%rootVarianceUniqueTable(i)%index       (uniqueSizeTmp(i)))
       call self%rootVarianceTable(i)%create  (self%massMinimum,self%massMaximum,size(massTmp))
       call self%rootVarianceTable(i)%populate(rootVarianceTmp(:,i))
       self%rootVarianceUniqueTable(i)%rootVariance=rootVarianceUniqueTmp(1:uniqueSizeTmp(i),i)
       self%rootVarianceUniqueTable(i)%index       =indexTmp             (1:uniqueSizeTmp(i),i)
    end do
    deallocate(rootVarianceTmp      )
    deallocate(rootVarianceUniqueTmp)
    deallocate(indexTmp             )
    deallocate(massTmp              )
    deallocate(timesTmp             )
    deallocate(uniqueSizeTmp        )
    self%initialized=.true.
    return
  end subroutine filteredPowerFileRead

  subroutine filteredPowerFileWrite(self)
    !!{
    Write tabulated data on mass variance to file.
    !!}
    use :: Display    , only : displayMessage, verbosityLevelWorking
    use :: HDF5       , only : hsize_t
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout)               :: self
    double precision                                       , dimension(:  ), allocatable :: massTmp
    double precision                                       , dimension(:,:), allocatable :: rootVarianceTmp, rootVarianceUniqueTmp
    integer                                                , dimension(:  ), allocatable :: uniqueSizeTmp
    integer                                                , dimension(:,:), allocatable :: indexTmp
    type            (hdf5Object                           )                              :: dataFile
    integer                                                                              :: i              , useCache

    ! Prepare data.
    allocate(massTmp              (self%rootVarianceTable(1)%size()                 ))
    allocate(rootVarianceTmp      (self%rootVarianceTable(1)%size(),size(self%times)))
    allocate(rootVarianceUniqueTmp(self%rootVarianceTable(1)%size(),size(self%times)))
    allocate(indexTmp             (self%rootVarianceTable(1)%size(),size(self%times)))
    allocate(uniqueSizeTmp        (                                 size(self%times)))
    rootVarianceUniqueTmp   ( :                , : )=     0.0d0
    indexTmp                ( :                , : )=     0
    massTmp                 ( :                    )=     self%rootVarianceTable      (1)%xs          ()
    do i=1,size(self%times)
       rootVarianceTmp      ( :                ,i:i)=     self%rootVarianceTable      (i)%ys          ()
       uniqueSizeTmp        (                   i  )=size(self%rootVarianceUniqueTable(i)%index         )
       rootVarianceUniqueTmp(1:uniqueSizeTmp(i),i  )=     self%rootVarianceUniqueTable(i)%rootVariance
       indexTmp             (1:uniqueSizeTmp(i),i  )=     self%rootVarianceUniqueTable(i)%index
    end do
    ! Cache this variance for possible later reuse.
    !$omp critical(cosmologicalMassVarianceFilteredPowerCache)
    useCache=0
    if (countCache > 0) then
       do i=1,countCache
          if (cachedVariances(i)%fileName == self%fileName) then
             useCache=i
             exit
          end if
       end do
    end if
    if (useCache == 0) then
       lastCache=lastCache+1
       if (lastCache > sizeCache) lastCache=1
       countCache=max(countCache,lastCache)
       useCache  =lastCache
    else
       deallocate(cachedVariances(useCache)%massTmp              )
       deallocate(cachedVariances(useCache)%timesTmp             )
       deallocate(cachedVariances(useCache)%rootVarianceTmp      )
       deallocate(cachedVariances(useCache)%rootVarianceUniqueTmp)
       deallocate(cachedVariances(useCache)%uniqueSizeTmp        )
       deallocate(cachedVariances(useCache)%indexTmp             )
    end if
    cachedVariances(useCache)%fileName                   =self%fileName
    cachedVariances(useCache)%timesTmp                   =self%times
    cachedVariances(useCache)%massTmp                    =     massTmp
    cachedVariances(useCache)%rootVarianceTmp            =     rootVarianceTmp
    cachedVariances(useCache)%rootVarianceUniqueTmp      =     rootVarianceUniqueTmp
    cachedVariances(useCache)%indexTmp                   =     indexTmp
    cachedVariances(useCache)%uniqueSizeTmp              =     uniqueSizeTmp
    cachedVariances(useCache)%sigma8Value                =self%sigma8Value
    cachedVariances(useCache)%sigmaNormalization         =self%sigmaNormalization
    cachedVariances(useCache)%massMinimum                =self%massMinimum
    cachedVariances(useCache)%massMaximum                =self%massMaximum
    cachedVariances(useCache)%timeMinimum                =self%timeMinimum
    cachedVariances(useCache)%timeMaximum                =self%timeMaximum
    cachedVariances(useCache)%timeMinimumLogarithmic     =self%timeMinimumLogarithmic
    cachedVariances(useCache)%timeLogarithmicDeltaInverse=self%timeLogarithmicDeltaInverse
    !$omp end critical(cosmologicalMassVarianceFilteredPowerCache)
    ! Store to file if requested.
    if (self%storeTabulations) then
       call displayMessage('writing σ(M) data to: '//self%fileName,verbosityLevelWorking)
       ! Open the data file.
       !$ call hdf5Access%set()
       call dataFile%openFile      (char(self%fileName)             ,overWrite                    =.true.,objectsOverwritable=.true.,chunkSize=100_hsize_t,compressionLevel=9)
       call dataFile%writeDataset  (self%times                      ,'times'                                                                                                 )
       call dataFile%writeDataset  (     massTmp                    ,'mass'                                                                                                  )
       call dataFile%writeDataset  (     rootVarianceTmp            ,'rootVariance'                                                                                          )
       call dataFile%writeDataset  (     rootVarianceUniqueTmp      ,'rootVarianceUnique'                                                                                    )
       call dataFile%writeDataset  (     indexTmp                   ,'indexUnique'                                                                                           )
       call dataFile%writeDataset  (     uniqueSizeTmp              ,'uniqueSize'                                                                                            )
       call dataFile%writeAttribute(self%sigma8Value                ,'sigma8'                                                                                                )
       call dataFile%writeAttribute(self%sigmaNormalization         ,'sigmaNormalization'                                                                                    )
       call dataFile%writeAttribute(self%massMinimum                ,'massMinimum'                                                                                           )
       call dataFile%writeAttribute(self%massMaximum                ,'massMaximum'                                                                                           )
       call dataFile%writeAttribute(self%timeMinimum                ,'timeMinimum'                                                                                           )
       call dataFile%writeAttribute(self%timeMaximum                ,'timeMaximum'                                                                                           )
       call dataFile%writeAttribute(self%timeMinimumLogarithmic     ,'timeMinimumLogarithmic'                                                                                )
       call dataFile%writeAttribute(self%timeLogarithmicDeltaInverse,'timeLogarithmicDeltaInverse'                                                                           )
       call dataFile%close         (                                                                                                                                         )
       !$ call hdf5Access%unset()
    end if
    deallocate(rootVarianceTmp      )
    deallocate(rootVarianceUniqueTmp)
    deallocate(indexTmp             )
    deallocate(massTmp              )
    deallocate(uniqueSizeTmp        )
    return
  end subroutine filteredPowerFileWrite

  logical function filteredPowerRemakeTable(self,mass,time)
    !!{
    Determine if the table should be remade.
    !!}
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout)           :: self
    double precision                                       , intent(in   ), optional :: mass, time

    if (self%initialized) then
       filteredPowerRemakeTable=.false.
       if (present(mass)) then
          filteredPowerRemakeTable=(                          &
               &                     mass < self%massMinimum  &
               &                    .or.                      &
               &                     mass > self%massMaximum  &
               &                   )
       end if
       if (present(time).and.self%growthIsMassDependent_) then
          filteredPowerRemakeTable=(                          &
               &                     filteredPowerRemakeTable &
               &                    .or.                      &
               &                     time < self%timeMinimum  &
               &                    .or.                      &
               &                     time > self%timeMaximum  &
               &                   )
       end if
    else
       filteredPowerRemakeTable=.true.
    end if
    return
  end function filteredPowerRemakeTable

  subroutine filteredPowerDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (cosmologicalMassVarianceFilteredPower), intent(inout)           :: self
    type     (inputParameters                      ), intent(inout)           :: descriptor
    logical                                         , intent(in   ), optional :: includeClass, includeFileModificationTimes
    type     (inputParameters                      )                          :: parameters

    call self%descriptorNormalizationOnly(descriptor,includeClass,includeFileModificationTimes)
    parameters=descriptor%subparameters('cosmologicalMassVariance')
    call self%powerSpectrumWindowFunction_%descriptor(parameters,includeClass,includeFileModificationTimes)
    return
  end subroutine filteredPowerDescriptor

  subroutine filteredPowerDescriptorNormalizationOnly(self,descriptor,includeClass,includeFileModificationTimes)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object, for power spectrum normalization usage
    only (i.e. we exclude the window function).
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (cosmologicalMassVarianceFilteredPower), intent(inout)           :: self
    type     (inputParameters                      ), intent(inout)           :: descriptor
    logical                                         , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=18                               )                          :: parameterLabel
    type     (inputParameters                      )                          :: parameters    , referenceParameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('cosmologicalMassVariance','filteredPower')
    parameters=descriptor%subparameters('cosmologicalMassVariance')
    if (self%normalizationSigma8) then
       write (parameterLabel,'(e17.10)') self%sigma8Value
       call parameters%addParameter('sigma_8'                ,trim(adjustl(parameterLabel)))
    else
       write (parameterLabel,'(e17.10)') self%wavenumberReference
       call parameters%addParameter('wavenumberReference'    ,trim(adjustl(parameterLabel)))
       call parameters%addParameter('reference','')
       referenceParameters=parameters%subparameters('reference')
       call self%cosmologicalMassVarianceReference          %descriptorNormalizationOnly(referenceParameters)
       call self%powerSpectrumPrimordialTransferredReference%descriptor                 (referenceParameters)
    end if
    write    (parameterLabel,'(e17.10)') self%toleranceTopHat
    call    parameters%addParameter('toleranceTopHat'       ,trim(adjustl(parameterLabel)))
    write    (parameterLabel,'(e17.10)') self%tolerance
    call    parameters%addParameter('tolerance'             ,trim(adjustl(parameterLabel)))
    write    (parameterLabel,'(l1)'    ) self%nonMonotonicIsFatal
    call    parameters%addParameter('nonMonotonicIsFatal'   ,trim(adjustl(parameterLabel)))
    write    (parameterLabel,'(l1)'    ) self%monotonicInterpolation
    call    parameters%addParameter('monotonicInterpolation',trim(adjustl(parameterLabel)))
    call    self%cosmologyParameters_                       %descriptor(parameters,includeClass,includeFileModificationTimes)
    call    self%cosmologyFunctions_                        %descriptor(parameters,includeClass,includeFileModificationTimes)
    call    self%powerSpectrumPrimordialTransferred_        %descriptor(parameters,includeClass,includeFileModificationTimes)
    call    self%linearGrowth_                              %descriptor(parameters,includeClass,includeFileModificationTimes)
    call    self%transferFunction_                          %descriptor(parameters,includeClass,includeFileModificationTimes)
    return
  end subroutine filteredPowerDescriptorNormalizationOnly
