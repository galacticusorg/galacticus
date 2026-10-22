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

  !+    Contributions to this file made by: Andrew Benson, Christoph Behrens, Xiaolong Du, Claude.

  !!{RST
  An implementation of cosmological density field mass variance computed using a filtered power spectrum.
  !!}

  !![
  <cosmologicalMassVariance name="cosmologicalMassVarianceFilteredPower" docformat="rst">
   <description>
   Mass variance of cosmological density fields computed from a filtered power spectrum:

   .. math::

      \sigma^2(M) = {1 \over 2 \pi^2} \int_0^\infty P(k) T^2(k) W^2(k) k^2 \mathrm{d}k

   where :math:`P(k) T^2(k)` is the (normalized) transferred primordial power spectrum (see :galacticus-class:`powerSpectrumPrimordialTransferred`), and :math:`W(k)` is the power spectrum variance window function (see :galacticus-class:`powerSpectrumWindowFunction`).

   The normalization of the mass variance is that of the transferred primordial power spectrum, which is determined by its :galacticus-class:`powerSpectrumNormalizer` (for example, to give a specified value of :math:`\sigma_8`).

   The mass variance, :math:`\sigma(M)`, is found by integration over the linear theory power spectrum, with the specified power spectrum window function. The fractional tolerance for this integration can be set via the ``[tolerance]`` parameter. This is tabulated across the required range.

   The value of :math:`\sigma_8` (the linear theory root-variance of the density field at the present day in spheres of radius :math:`8h^{-1}`\ Mpc) is computed from the power spectrum using a real-space top-hat window function (as per the definition of :math:`\sigma_8`), unless a different window function is explicitly defined via the ``[powerSpectrumWindowFunctionTopHat]`` parameter. For this integration the tolerance can be set via the ``[toleranceTopHat]`` parameter.

   Cubic spline interpolation is then used to interpolate in this table to give :math:`\sigma(M)` at any required value of :math:`M`. The tabulation is always forced to be monotonically decreasing with :math:`M`. However, the interpolation is not necessarily monotonic---for example in cases where :math:`\sigma(M)` becomes constant or close to constant as a function of :math:`M` the interpolation can become non-monotonic over some ranges of :math:`M`. If strict monotonicity is required set ``[monotonicInterpolation]``\ =\ ``true``. This causes a monotonic spline interpolator to be used instead which guarantees monotonicity.
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
  use :: Numerical_Ranges                    , only : rangeLattice

  !![
  <stateStorable class="uniqueTable"/>
  !!]

  type :: uniqueTable
     !!{RST
     Type used to store unique values of the mass variance.
     !!}
     private
     double precision, allocatable, dimension(:) :: rootVariance
     integer         , allocatable, dimension(:) :: index
  end type uniqueTable

  type, extends(cosmologicalMassVarianceClass) :: cosmologicalMassVarianceFilteredPower
     !!{RST
     A cosmological mass variance class computing variance from a filtered power spectrum.
     !!}
     private
     class           (cosmologyParametersClass               ), pointer                   :: cosmologyParameters_                => null()
     class           (cosmologyFunctionsClass                ), pointer                   :: cosmologyFunctions_                 => null()
     class           (powerSpectrumPrimordialTransferredClass), pointer                   :: powerSpectrumPrimordialTransferred_ => null()
     class           (linearGrowthClass                      ), pointer                   :: linearGrowth_                       => null()
     class           (transferFunctionClass                  ), pointer                   :: transferFunction_                   => null()
     class           (powerSpectrumWindowFunctionClass       ), pointer                   :: powerSpectrumWindowFunction_        => null() , powerSpectrumWindowFunctionTopHat_          => null()
     logical                                                                              :: initialized                         =  .false., nonMonotonicIsFatal                                  , &
          &                                                                                  integrationFailureIsFatal
     double precision                                                                     :: tolerance                                     , toleranceTopHat                                      , &
          &                                                                                  sigma8Value                                   , massMinimum                                          , &
          &                                                                                  massMaximum                                   , timeMinimum                                          , &
          &                                                                                  timeMaximum                                   , timeLogarithmicDeltaInverse                          , &
          &                                                                                  wavenumberHalfMode                            , rootVarianceLogarithmicGradientTolerance
     double precision                                         , allocatable, dimension(:) :: times
     class           (table1DLinearCSpline                   ), allocatable, dimension(:) :: rootVarianceTable
     ! Lattices to which the two axes of the tabulation are pinned. These are the source of truth for its extent: the limits
     ! above are derived from them (in fact from the tables built on them - see `filteredPowerRetabulate`), since they are read
     ! in many places. The time lattice is left undefined where growth is mass-independent, since σ(M) is then tabulated at the
     ! present day alone and there is no time axis to pin.
     type            (rangeLattice                           )                            :: latticeMass                                   , latticeTime
     type            (varying_string                         )                            :: fileName
     type            (lockDescriptor                         )                            :: fileLock
     ! Unique values in the variance table and their corresponding indices.
     type            (uniqueTable                            ), allocatable, dimension(:) :: rootVarianceUniqueTable
     logical                                                                              :: monotonicInterpolation                        , growthIsMassDependent_                               , &
          &                                                                                  truncateAtParticleHorizon                     , storeTabulations                                     , &
          &                                                                                  warnedNonIncreasing
   contains
     !![
     <methods docformat="rst">
       <method description="Tabulate cosmological mass variance."                                                                                                                                     method="retabulate"          />
       <method description="Compute the interpolating factors in time."                                                                                                                               method="interpolantsTime"    />
       <method description="Write the tabulated mass variance to file."                                                                                                                               method="fileWrite"           />
       <method description="Read the tabulated mass variance from file."                                                                                                                              method="fileRead"            />
       <method description="Return true if the table must be remade."                                                                                                                                 method="remakeTable"         />
       <method description="Compute the root-variance of the density field in a region containing the given mass at the given time, using either the window function or the top-hat window function." method="rootVarianceIntegral"/>
     </methods>
     !!]
     final     ::                                        filteredPowerDestructor
     procedure :: sigma8                              => filteredPowerSigma8
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
     procedure :: rootVarianceIntegral                => filteredPowerRootVarianceIntegral
     procedure :: descriptor                          => filteredPowerDescriptor
  end type cosmologicalMassVarianceFilteredPower

  interface cosmologicalMassVarianceFilteredPower
     !!{RST
     Constructors for the :galacticus-class:`cosmologicalMassVarianceFilteredPower` cosmological mass variance class.
     !!}
     module procedure filteredPowerConstructorParameters
     module procedure filteredPowerConstructorInternal
  end interface cosmologicalMassVarianceFilteredPower

  ! Number of points per decade to use in tabulation of σ(M).
  integer                         , parameter :: pointsPerDecade=10             , timePointsPerDecade=100

  ! Intervals, in lattice steps, to which the bounds of the two axes are pinned. The mass axis is pinned to whole decades: a
  ! decade of extra tabulation costs `pointsPerDecade` integrals per epoch, and the ×10 safety margin which this replaces was
  ! already of that order. The time axis is pinned to the lattice points themselves, because *every* extra epoch costs a
  ! complete σ(M) tabulation, and a whole decade of cosmic time - `timePointsPerDecade` epochs - would be extravagant beyond
  ! any measure of the accuracy it buys. Both choices are equally deterministic; they differ only in the granularity of the
  ! discrete set of ranges which can result.
  integer                         , parameter :: anchorEveryMass=pointsPerDecade, anchorEveryTime    =  1

  ! Format version of stored σ(M) tabulations. Version 2 marks tabulations made after transfer functions were required to be
  ! normalized to unity on large scales. The stored normalization of the power spectrum is relative to the normalization of the
  ! transfer function, so a tabulation made before that change would misnormalize the power spectrum. Such files can not be
  ! detected through the file name (whose source digest does not include that of the transfer function), so are instead
  ! rejected on this version. Version 3 marks tabulations made after the curvature parameter passed to CAMB and CLASS was
  ! corrected, and the curvature factor 1-3K/k² removed from tabulated transfer functions, both of which change the transfer
  ! function in non-flat cosmologies. Version 4 marks tabulations made after the normalization of the power spectrum was moved to
  ! the transferred primordial power spectrum (so that σ(M) is tabulated from the normalized power spectrum, and the
  ! normalization is no longer stored).
  integer                         , parameter :: filteredPowerFileFormatCurrent=4

  ! Seed range of masses to tabulate. Every tabulation contains this range, so any two which are never asked for a mass outside
  ! it span *exactly* the same range - which matters here beyond the usual guarantee that two ranges overlap and so can be
  ! merged. σ(M) is interpolated with a cubic spline, every coefficient of which depends on every tabulated value, so σ at a
  ! given mass is not a function of that mass alone: it depends on the range over which the tabulation was built. Since
  ! tabulations are cached to file and merged across runs and across models sharing the same parameters, a model which extended
  ! the range would otherwise obtain a different σ(M) from one which ran before it and did not - making results depend on what
  ! else had run, which is the dependence this pinning exists to remove. The range spans the halo masses which models
  ! plausibly reach; a mass outside it still extends the tabulation, and the dependence returns for models which straddle it.
  double precision                , parameter :: massTableSeedMinimum=1.0d6     , massTableSeedMaximum=1.0d18

  ! Cached copies of tabulated solutions. These are used to avoid re-reading from file if the same variance is requested multiple times.
  type :: cachedVariance
     type            (varying_string)                              :: fileName
     double precision                , dimension(:  ), allocatable :: massTmp        , timesTmp
     double precision                , dimension(:,:), allocatable :: rootVarianceTmp, rootVarianceUniqueTmp
     integer                         , dimension(:  ), allocatable :: uniqueSizeTmp
     integer                         , dimension(:,:), allocatable :: indexTmp
     type            (rangeLattice  )                              :: latticeMass    , latticeTime
  end type cachedVariance
  
  integer                , parameter            :: sizeCache      =25
  integer                                       :: countCache     = 0, lastCache=0
  type   (cachedVariance), dimension(sizeCache) :: cachedVariances
  
contains

  function filteredPowerConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`cosmologicalMassVarianceFilteredPower` cosmological mass variance class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (cosmologicalMassVarianceFilteredPower  )                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (cosmologyParametersClass               ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    class           (powerSpectrumPrimordialTransferredClass), pointer       :: powerSpectrumPrimordialTransferred_
    class           (powerSpectrumWindowFunctionClass       ), pointer       :: powerSpectrumWindowFunction_            , powerSpectrumWindowFunctionTopHat_
    class           (linearGrowthClass                      ), pointer       :: linearGrowth_
    class           (transferFunctionClass                  ), pointer       :: transferFunction_
    double precision                                                         :: tolerance                               , toleranceTopHat                   , &
         &                                                                      rootVarianceLogarithmicGradientTolerance
    logical                                                                  :: monotonicInterpolation                  , nonMonotonicIsFatal               , &
         &                                                                      truncateAtParticleHorizon               , storeTabulations                  , &
         &                                                                      integrationFailureIsFatal
    
    !![
    <objectBuilder    class="cosmologyParameters"                name="cosmologyParameters_"                source="parameters"                                                  />
    <objectBuilder    class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"                                                  />
    <objectBuilder    class="powerSpectrumPrimordialTransferred" name="powerSpectrumPrimordialTransferred_" source="parameters"                                                  />
    <objectBuilder    class="powerSpectrumWindowFunction"        name="powerSpectrumWindowFunction_"        source="parameters"                                                  />
    <objectBuilder    class="linearGrowth"                       name="linearGrowth_"                       source="parameters"                                                  />
    <objectBuilder    class="transferFunction"                   name="transferFunction_"                   source="parameters"                                                  />
    !!]
    if (parameters%isPresent('powerSpectrumWindowFunctionTopHat')) then
       !![
       <objectBuilder class="powerSpectrumWindowFunction"        name="powerSpectrumWindowFunctionTopHat_"  source="parameters" parameterName="powerSpectrumWindowFunctionTopHat"/>
       !!]
    else
       nullify(powerSpectrumWindowFunctionTopHat_)
    end if
    !![
    <inputParameter docformat="rst">
      <name>toleranceTopHat</name>
      <source>parameters</source>
      <defaultValue>1.0d-6</defaultValue>
      <description>
      The relative tolerance to use in integrating over the linear power spectrum using a top-hat (real space) window function to compute :math:`\sigma_8`.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>tolerance</name>
      <source>parameters</source>
      <defaultValue>4.0d-6</defaultValue>
      <description>
      The relative tolerance to use in integrating over the linear power spectrum to compute the cosmological mass variance.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>nonMonotonicIsFatal</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true any non-monotonicity in the tabulated :math:`\sigma(M)` is treated as a fatal error. Otherwise a only a warning is issued.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>integrationFailureIsFatal</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true any failed integrals when evaluating :math:`\sigma(M)` are treated as fatal errors. Otherwise a only a warning is issued.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>monotonicInterpolation</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>
      If true use a monotonic cubic spline interpolator to interpolate in the :math:`\sigma(M)` table. Otherwise use a standard cubic spline interpolator. Use of the monotonic interpolator can be helpful is :math:`\sigma(M)` must be strictly monotonic but becomes a very weak function of :math:`M` at low masses.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>rootVarianceLogarithmicGradientTolerance</name>
      <source>parameters</source>
      <defaultValue>1.0d-12</defaultValue>
      <description>
      The tolerance in :math:`\mathrm{d}\log\sigma/\mathrm{d}\log M` to allow before reporting errors when monotonic interpolation is used..
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>truncateAtParticleHorizon</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>
      If true then integration over the power spectrum is truncated at a wavenumber :math:`k=1/H_\mathrm{p}(t_0)`, where :math:`H_\mathrm{p}(t_0)` is the comoving distance to the particle horizon at the present epoch. Otherwise, integration continues to :math:`k=0`.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>storeTabulations</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true then tabulated :math:`\sigma(M)` results are stored to file for future re-use.
      </description>
    </inputParameter>
    !!]
    if (associated(powerSpectrumWindowFunctionTopHat_)) then
       self=filteredPowerConstructorInternal(                                                                                   &
            &                                tolerance                               =tolerance                               , &
            &                                toleranceTopHat                         =toleranceTopHat                         , &
            &                                nonMonotonicIsFatal                     =nonMonotonicIsFatal                     , &
            &                                integrationFailureIsFatal               =integrationFailureIsFatal               , &
            &                                monotonicInterpolation                  =monotonicInterpolation                  , &
            &                                rootVarianceLogarithmicGradientTolerance=rootVarianceLogarithmicGradientTolerance, &
            &                                truncateAtParticleHorizon               =truncateAtParticleHorizon               , &
            &                                storeTabulations                        =storeTabulations                        , &
            &                                cosmologyParameters_                    =cosmologyParameters_                    , &
            &                                cosmologyFunctions_                     =cosmologyFunctions_                     , &
            &                                linearGrowth_                           =linearGrowth_                           , &
            &                                transferFunction_                       =transferFunction_                       , &
            &                                powerSpectrumPrimordialTransferred_     =powerSpectrumPrimordialTransferred_     , &
            &                                powerSpectrumWindowFunction_            =powerSpectrumWindowFunction_            , &
            &                                powerSpectrumWindowFunctionTopHat_      =powerSpectrumWindowFunctionTopHat_        &
            &                               )
    else
       self=filteredPowerConstructorInternal(                                                                                   &
            &                                tolerance                               =tolerance                               , &
            &                                toleranceTopHat                         =toleranceTopHat                         , &
            &                                nonMonotonicIsFatal                     =nonMonotonicIsFatal                     , &
            &                                integrationFailureIsFatal               =integrationFailureIsFatal               , &
            &                                monotonicInterpolation                  =monotonicInterpolation                  , &
            &                                rootVarianceLogarithmicGradientTolerance=rootVarianceLogarithmicGradientTolerance, &
            &                                truncateAtParticleHorizon               =truncateAtParticleHorizon               , &
            &                                storeTabulations                        =storeTabulations                        , &
            &                                cosmologyParameters_                    =cosmologyParameters_                    , &
            &                                cosmologyFunctions_                     =cosmologyFunctions_                     , &
            &                                linearGrowth_                           =linearGrowth_                           , &
            &                                transferFunction_                       =transferFunction_                       , &
            &                                powerSpectrumPrimordialTransferred_     =powerSpectrumPrimordialTransferred_     , &
            &                                powerSpectrumWindowFunction_            =powerSpectrumWindowFunction_              &
            &                               )
    end if
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor    name="cosmologyParameters_"               />
    <objectDestructor    name="cosmologyFunctions_"                />
    <objectDestructor    name="linearGrowth_"                      />
    <objectDestructor    name="powerSpectrumPrimordialTransferred_"/>
    <objectDestructor    name="powerSpectrumWindowFunction_"       />
    <objectDestructor    name="transferFunction_"                  />
    !!]
    if (parameters%isPresent('powerSpectrumWindowFunctionTopHat')) then
       !![
       <objectDestructor name="powerSpectrumWindowFunctionTopHat_" />
       !!]
    end if
    return
  end function filteredPowerConstructorParameters

  function filteredPowerConstructorInternal(tolerance,toleranceTopHat,nonMonotonicIsFatal,integrationFailureIsFatal,monotonicInterpolation,rootVarianceLogarithmicGradientTolerance,truncateAtParticleHorizon,storeTabulations,cosmologyParameters_,cosmologyFunctions_,linearGrowth_,transferFunction_,powerSpectrumPrimordialTransferred_,powerSpectrumWindowFunction_,powerSpectrumWindowFunctionTopHat_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`cosmologicalMassVarianceFilteredPower` cosmological mass variance class.
    !!}
    use :: File_Utilities                   , only : Directory_Make                              , File_Path
    use :: Input_Paths                      , only : inputPath                                   , pathTypeDataDynamic
    use :: Power_Spectrum_Variance_Integrals, only : Power_Spectrum_Variance_Wavenumber_Half_Mode
    use :: Power_Spectrum_Window_Functions  , only : powerSpectrumWindowFunctionTopHat
    implicit none
    type            (cosmologicalMassVarianceFilteredPower  )                                  :: self
    double precision                                         , intent(in   )                   :: tolerance                               , toleranceTopHat       , &
         &                                                                                        rootVarianceLogarithmicGradientTolerance
    logical                                                  , intent(in   )                   :: nonMonotonicIsFatal                     , monotonicInterpolation, &
         &                                                                                        truncateAtParticleHorizon               , storeTabulations      , &
         &                                                                                        integrationFailureIsFatal
    class           (cosmologyParametersClass               ), intent(in   ), target           :: cosmologyParameters_
    class           (cosmologyFunctionsClass                ), intent(in   ), target           :: cosmologyFunctions_
    class           (powerSpectrumPrimordialTransferredClass), intent(in   ), target           :: powerSpectrumPrimordialTransferred_
    class           (powerSpectrumWindowFunctionClass       ), intent(in   ), target           :: powerSpectrumWindowFunction_
    class           (powerSpectrumWindowFunctionClass       ), intent(in   ), target, optional :: powerSpectrumWindowFunctionTopHat_
    class           (linearGrowthClass                      ), intent(in   ), target           :: linearGrowth_
    class           (transferFunctionClass                  ), intent(in   ), target, optional :: transferFunction_ 
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
    self%initialized           =.false.
    self%warnedNonIncreasing   =.false.
    self%sigma8Value           =-1.0d0
    self%growthIsMassDependent_=self%powerSpectrumPrimordialTransferred_%growthIsWavenumberDependent()
    self%fileName              =inputPath(pathTypeDataDynamic)                                                       // &
         &                      'largeScaleStructure/'                                                               // &
         &                      self%objectType      (                                                              )// &
         &                      '_'                                                                                  // &
         &                      self%hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)// &
         &                      '.hdf5'
    call Directory_Make(File_Path(self%fileName))
    if (present(transferFunction_)) then
       self%wavenumberHalfMode=Power_Spectrum_Variance_Wavenumber_Half_Mode(self%transferFunction_,self%cosmologyParameters_)
    else
       self%wavenumberHalfMode=-1.0d0
    end if
    return
  end function filteredPowerConstructorInternal

  subroutine filteredPowerDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`cosmologicalMassVarianceFilteredPower` cosmological mass variance class.
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
    if (self%initialized) then
       do i=1,size(self%rootVarianceTable)
          call self%rootVarianceTable(i)%destroy()
       end do
    end if
    return
  end subroutine filteredPowerDestructor

  double precision function filteredPowerSigma8(self)
    !!{RST
    Return the value of :math:`\sigma_8`, the linear theory root-variance of the density field at the present day in spheres of
    radius :math:`8h^{-1}` Mpc, computed (on the first call) using a top-hat window function.
    !!}
    use :: Cosmology_Parameters    , only : hubbleUnitsLittleH
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    ! Radius for σ₈ in Mpc/h.
    double precision                                       , parameter     :: radiusSigma8=8.0d0
    double precision                                                       :: massSigma8

    if (self%sigma8Value < 0.0d0) then
       massSigma8      =+(                                                               &
            &             +4.0d0                                                         &
            &             /3.0d0                                                         &
            &             *Pi                                                            &
            &            )                                                               &
            &           *  self%cosmologyParameters_%OmegaMatter    (                  ) &
            &           *  self%cosmologyParameters_%densityCritical(                  ) &
            &           *(                                                               &
            &             +radiusSigma8                                                  &
            &             /self%cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH) &
            &            )**3
       self%sigma8Value=self%rootVarianceIntegral(massSigma8,self%cosmologyFunctions_%cosmicTime(1.0d0),useTopHat=.true.)
    end if
    filteredPowerSigma8=self%sigma8Value
    return
  end function filteredPowerSigma8

  double precision function filteredPowerRootVariance(self,mass,time)
    !!{RST
    Return the root-variance of the cosmological density field in a spherical region containing the given ``mass`` on average.
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
    !!{RST
    Return the logarithmic gradient with respect to mass of the root-variance of the cosmological density field in a spherical region containing the given ``mass`` on average.
    !!}
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: mass        , time
    double precision                                                       :: rootVariance

    call self%rootVarianceAndLogarithmicGradient(mass,time,rootVariance,filteredPowerRootVarianceLogarithmicGradient)
    return
  end function filteredPowerRootVarianceLogarithmicGradient

  double precision function filteredPowerRootVarianceLogarithmicGradientTime(self,mass,time)
    !!{RST
    Return the logarithmic gradient with respect to time of the root-variance of the cosmological density field in a spherical region containing the given ``mass`` on average.
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
    !!{RST
    Return the value and logarithmic gradient with respect to mass of the root-variance of the cosmological density field in a spherical region containing the given ``mass`` on average.
    !!}
    use :: Display                 , only : displayGreen, displayReset
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    use :: HDF5_Access             , only : hdf5Access
    use :: IO_HDF5                 , only : hdf5File
    use :: String_Handling         , only : operator(//), stringXMLFormat
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: mass         , time
    double precision                                       , intent(  out) :: rootVariance , rootVarianceLogarithmicGradient
    type            (varying_string                       ), save          :: errorFileName
    !$omp threadprivate(errorFileName)
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
             !$ call hdf5Access%set()
             hdf5ErrorScope: block
               type(hdf5File  ) :: errorFile
               errorFileName=self%fileName//".error."//GetPID()
               errorFile=hdf5File(errorFileName,overWrite=.true.,objectsOverwritable=.true.)
               call errorFile%writeAttribute(mass                           ,'mass'                           )
               call errorFile%writeAttribute(rootVariance                   ,'rootVariance'                   )
               call errorFile%writeAttribute(rootVarianceGradient           ,'rootVarianceGradient'           )
               call errorFile%writeAttribute(rootVarianceLogarithmicGradient,'rootVarianceLogarithmicGradient')
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
             end block hdf5ErrorScope
             !$ call hdf5Access%unset()
          end if
          call Error_Report('dlogσ/dlogM > 0 detected, but monotonic interpolation was used - this should not happen'//char(10)//'  table data written to:'//char(10)//'    '//char(errorFileName)//{introspection:location})
       else
          ! Recommend that monotonic interpolation be used.
          call Error_Report(                                                                                                                                                                                                       &
               &            'dlogσ/dlogM > 0 detected'//char(10)//                                                                                                                                                                 &
               &            displayGreen()//'HELP:'//displayReset()//' set the highlighted option in your input parameter file as shown below: '                                                            //char(10)//char(10)// &
               &            stringXMLFormat('<cosmologicalMassVariance value="'//self%objectType(short=.true.)//'">**B<monotonicInterpolation value="true"/>**C</cosmologicalMassVariance>',indentInitial=3)//char(10)//char(10)// &
               &            'to ensure monotonic interpolation of tabulated σ(M) function'//{introspection:location}                                                                                                               &
               &           )
       end if
    end if
    return
  end subroutine filteredPowerRootVarianceAndLogarithmicGradient

  double precision function filteredPowerMass(self,rootVariance,time)
    !!{RST
    Return the mass corresponding to the given root-variance of the cosmological density field.
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
    !!{RST
    Tabulate the cosmological mass variance.
    !!}
    use :: Display                 , only : displayIndent            , displayMagenta                   , displayMessage     , displayReset, &
          &                                 displayUnindent          , verbosityLevelWorking
    use :: File_Utilities          , only : File_Lock                , File_Unlock                      , lockDescriptor
    use :: Error                   , only : Error_Report             , Warn
    use :: Numerical_Ranges        , only : Range_Pinned             , Range_Lattice_Offset             , gridSchemePerDecade
    use :: Tables                  , only : table1DLogarithmicCSpline, table1DLogarithmicMonotoneCSpline
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout)               :: self
    double precision                                       , intent(in   ), optional     :: mass                           , time
    integer                                                                              :: i                              , rootVarianceTableCount , &
         &                                                                                  j                              , rootVarianceUniqueCount, &
         &                                                                                  rootVarianceTimeCount          , k                      , &
         &                                                                                  offsetTime                     , countTimePrevious      , &
         &                                                                                  offsetMass                     , countMassPrevious      , &
         &                                                                                  kPrevious                      , iMinimum
    double precision                                                                     :: sigma                          , smoothingMass          , &
         &                                                                                  massMinimum                    , sigmaMinimum
    double precision                                                     , dimension(2) :: timeSeed
    type            (rangeLattice                         )                              :: latticeMass                    , latticeTime
    logical                                                                              :: carryOver
    logical                                                , allocatable  , dimension(:) :: rootVarianceIsUnique           , rootVarianceIsComputed
    class           (table1DLinearCSpline                 ), allocatable  , dimension(:) :: rootVarianceTablePrevious
    type            (varying_string                       ), save                        :: message
    character       (len=12                               )                              :: label                          , labelLow               , &
         &                                                                                  labelHigh                      , labelTarget
    ! The variable "message" is saved (and made threadprivate) as its destructor is expensive, and this function gets called a
    ! lot.
    !$omp threadprivate(message)

    if (self%remakeTable(mass,time)) then
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       if (self%storeTabulations) call File_Lock(self%fileName,self%fileLock,lockIsShared=.true.)
       call self%fileRead()
       if (self%storeTabulations) call File_Unlock(self%fileLock,sync=.false.)
    end if
    if (self%remakeTable(mass,time)) then
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       if (self%storeTabulations) call File_Lock(self%fileName,self%fileLock,lockIsShared=.false.)
       ! Try again to read the file - another process/thread may have already created the file in which case we may not need to do so again.
       call self%fileRead()
       if (self%remakeTable(mass,time)) then
          ! Find suitable range of masses to tabulate. The range is pinned to an absolute lattice, so that the abscissae of the
          ! tabulation - and therefore every value interpolated from it - depend only on which lattice points are spanned, and
          ! not on the sequence of masses which happened to be requested. Taking the union with the lattice already in use
          ! (`latticeCurrent`) guarantees that the range only ever grows, so that previously computed values can always be
          ! carried over; the request itself is passed as the target, never the current range, since folding the current range
          ! into the target would apply the safety margin to an already-margined bound and ratchet the range upward on every
          ! retabulation.
          if      (present(mass)               ) then
             latticeMass=Range_Pinned(                                                             &
                  &                                   [mass]                                     , &
                  &                                    pointsPerDecade                           , &
                  &                                    gridSchemePerDecade                       , &
                  &                   anchorEvery   =  anchorEveryMass                           , &
                  &                   rangeCurrent  = [massTableSeedMinimum,massTableSeedMaximum], &
                  &                   latticeCurrent=self%latticeMass                              &
                  &                  )
          else if (self%latticeMass%isDefined()) then
             latticeMass=self%latticeMass
          else
             ! No mass was given, and there is no tabulation to retain. The seed range serves as the target; it is already a whole
             ! number of decades, and so is taken with no safety margin.
             latticeMass=Range_Pinned(                                                             &
                  &                                   [massTableSeedMinimum,massTableSeedMaximum], &
                  &                                    pointsPerDecade                           , &
                  &                                    gridSchemePerDecade                       , &
                  &                   marginFactor  =  1.0d0                                     , &
                  &                   anchorEvery   =  anchorEveryMass                             &
                  &                  )
          end if
          rootVarianceTableCount=latticeMass%count
          ! Find suitable range of times to tabulate, pinned to an absolute lattice exactly as the mass range is.
          if (self%growthIsMassDependent_) then
             if      (present(time)               ) then
                latticeTime=Range_Pinned(                                      &
                     &                                   [time]              , &
                     &                                    timePointsPerDecade, &
                     &                                    gridSchemePerDecade, &
                     &                   anchorEvery   =  anchorEveryTime    , &
                     &                   latticeCurrent=self%latticeTime       &
                     &                  )
             else if (self%latticeTime%isDefined()) then
                latticeTime=self%latticeTime
             else
                ! No time was given, and there is no tabulation to retain. Must provide some time range.
                timeSeed   =[                                            &
                     &       self%cosmologyFunctions_%cosmicTime(0.5d0), &
                     &       self%cosmologyFunctions_%cosmicTime(1.0d0)  &
                     &      ]
                latticeTime=Range_Pinned(                                      &
                     &                                   timeSeed            , &
                     &                                    timePointsPerDecade, &
                     &                                    gridSchemePerDecade, &
                     &                   marginFactor  =  1.0d0              , &
                     &                   anchorEvery   =  anchorEveryTime      &
                     &                  )
             end if
             rootVarianceTimeCount=latticeTime%count
          else
             ! Growth of the transferred power spectrum is independent of mass - we can therefore tabulate σ(M) at a single epoch
             ! and use the linear growth factor to transform it to other epochs. There is then no time axis to pin, so the time
             ! lattice is left undefined.
             latticeTime                     =rangeLattice()
             rootVarianceTimeCount           =1
             self%timeMinimum                =self%cosmologyFunctions_%cosmicTime(1.0d0)
             self%timeMaximum                =self%timeMinimum
             self%timeLogarithmicDeltaInverse=0.0d0
          end if
          ! Establish the epochs at which σ(M) is to be tabulated. The abscissae are taken from the lattice rather than by
          ! subdividing the range, so that they are bit-identical to those of any other tabulation built on the same lattice,
          ! and the interpolating factor in time likewise comes from the lattice spacing - which is a pure function of the
          ! density of points - rather than from the extent of this particular range.
          if (allocated(self%times)) deallocate(self%times)
          allocate(self%times(rootVarianceTimeCount))
          if (self%growthIsMassDependent_) then
             self%times                      =latticeTime%values()
             self%timeMinimum                =self%times(                    1)
             self%timeMaximum                =self%times(rootVarianceTimeCount)
             self%timeLogarithmicDeltaInverse=1.0d0/latticeTime%stepLogarithmic()
          else
             self%times                      =self%timeMinimum
          end if
          if (allocated(self%rootVarianceUniqueTable)) deallocate(self%rootVarianceUniqueTable)
          allocate(self%rootVarianceUniqueTable(rootVarianceTimeCount))
          ! Determine which of the tables already built are to be carried over, and where they sit in the extended tabulation.
          ! Along the time axis a table is carried over whole; along the mass axis it is carried over as the block of points
          ! which its lattice occupies in the new one. Both offsets are found in exact integer arithmetic from the lattice
          ! indices, so no abscissa is ever compared or searched for.
          carryOver        =       allocated  (self%rootVarianceTable)                                       &
               &            .and.  self%latticeMass%isDefined      ()                                        &
               &            .and. (self%latticeTime%isDefined      () .or. .not.self%growthIsMassDependent_)
          offsetTime       =0
          offsetMass       =0
          countTimePrevious=0
          countMassPrevious=0
          if (carryOver) then
             countMassPrevious=self%latticeMass%count
             offsetMass       =Range_Lattice_Offset(self%latticeMass,latticeMass)
             if (self%growthIsMassDependent_) then
                countTimePrevious=self%latticeTime%count
                offsetTime       =Range_Lattice_Offset(self%latticeTime,latticeTime)
             else
                ! There is a single, fixed epoch, so the sole table is always carried over in place.
                countTimePrevious=1
             end if
             if (countTimePrevious /= size(self%rootVarianceTable)) call Error_Report('tabulation does not match the lattice on which it was built'//{introspection:location})
             call Move_Alloc(self%rootVarianceTable,rootVarianceTablePrevious)
          else if (allocated(self%rootVarianceTable)) then
             do i=1,size(self%rootVarianceTable)
                call self%rootVarianceTable(i)%destroy()
             end do
             deallocate(self%rootVarianceTable)
          end if
          self%latticeMass=latticeMass
          self%latticeTime=latticeTime
          ! Allocate table grid.
          if (self%monotonicInterpolation) then
             allocate(table1DLogarithmicMonotoneCSpline :: self%rootVarianceTable(rootVarianceTimeCount))
          else
             allocate(table1DLogarithmicCSpline         :: self%rootVarianceTable(rootVarianceTimeCount))
          end if
          call displayIndent("retabulating σ(M)",verbosityLevelWorking)
          write    (labelLow   ,'(e9.2)') latticeMass%minimum()
          write    (labelHigh  ,'(e9.2)') latticeMass%maximum()
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
             call self%rootVarianceTable(k)%extend(latticeMass,rootVarianceIsComputed)
             kPrevious=k-offsetTime
             if (kPrevious >= 1 .and. kPrevious <= countTimePrevious) then
                self%rootVarianceTable(k)%yv(offsetMass+1:offsetMass+countMassPrevious,1)=rootVarianceTablePrevious(kPrevious)%yv(:,1)
                rootVarianceIsComputed   (offsetMass+1:offsetMass+countMassPrevious    )=.true.
             end if
             allocate(rootVarianceIsUnique(rootVarianceTableCount))
             rootVarianceIsUnique=.true.
             ! Compute σ(M) at each tabulated point.
             massMinimum =-1.0d0
             sigmaMinimum=-1.0d0
             iMinimum    =-1
             do i=1,rootVarianceTableCount
                smoothingMass=+self%rootVarianceTable(k)%x(i)
                if (rootVarianceIsComputed(i)) then
                   ! This point was tabulated on an earlier pass. Its σ(M) integral - which is the whole cost of the tabulation -
                   ! is precisely what carrying values over exists to avoid repeating.
                   sigma     =+self%rootVarianceTable(k)%y(i)
                else
                   sigma     =+self%rootVarianceIntegral(smoothingMass,self%times(k),useTopHat=.false.)
                end if
                ! Enforce monotonicity. The clamp is re-applied across the whole axis, carried-over points included: a point
                ! which acquires a new predecessor when the range is extended downward must be clamped against it. Re-clamping
                ! can only lower a value, so the result is exactly that of a single-pass tabulation over the whole range - and
                ! it costs nothing, being a comparison rather than an integral.
                if (i > 1) then
                   if (sigma >= self%rootVarianceTable(k)%y(i-1)) then
                      iMinimum               =i
                      massMinimum            =smoothingMass
                      sigmaMinimum           =sigma
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
          if (allocated(rootVarianceTablePrevious)) then
             do k=1,size(rootVarianceTablePrevious)
                call rootVarianceTablePrevious(k)%destroy()
             end do
             deallocate(rootVarianceTablePrevious)
          end if
          call displayUnindent("done",verbosityLevelWorking)
          ! Record the limits of the tabulation. These are taken from the end points of the tables themselves rather than from
          ! the lattice: the abscissae of a logarithmically-spaced table are the exponentials of its internal, logarithmic
          ! abscissae, and 10**(k/N) is not bit-identical to exp((k/N)ln10). Taking the limits from the lattice could therefore
          ! place them a fraction of an ulp outside the table, so that a mass which passes the test in `remakeTable` still lies
          ! beyond the range which the table actually spans.
          self%massMinimum=self%rootVarianceTable(1)%x(+1)
          self%massMaximum=self%rootVarianceTable(1)%x(-1)
          ! Table is now initialized.
          self%initialized=.true.
          ! Store file.
          call self%fileWrite()
       end if
       if (self%storeTabulations) call File_Unlock(self%fileLock)
    end if
    return
  end subroutine filteredPowerRetabulate

  double precision function filteredPowerRootVarianceIntegral(self,smoothingMass,time,useTopHat) result(rootVariance)
    !!{RST
    Compute the root-variance of the density field in a region containing the given ``smoothingMass`` on average, at the given
    ``time``, by integration over the power spectrum. If ``useTopHat`` is true the top-hat window function (used to compute
    :math:`\sigma_8`) is used, otherwise the window function used to compute :math:`\sigma(M)` is used.
    !!}
    use :: Numerical_Constants_Math         , only : Pi
    use :: Power_Spectrum_Variance_Integrals, only : Power_Spectrum_Variance_Integrate
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout)               :: self
    double precision                                       , intent(in   )               :: smoothingMass               , time
    logical                                                , intent(in   )               :: useTopHat
    double precision                                       , allocatable  , dimension(:) :: wavenumbersLocalMinima
    class           (powerSpectrumWindowFunctionClass     ), pointer                     :: powerSpectrumWindowFunction_
    double precision                                                                     :: wavenumberMinimum           , wavenumberMaximum, &
         &                                                                                  topHatRadius                , toleranceRelative
    character       (len=15                               )                              :: toleranceName

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
    ! Select the window function, and find the maximum wavenumber for the integral.
    if (useTopHat) then
       powerSpectrumWindowFunction_ => self%powerSpectrumWindowFunctionTopHat_
       toleranceRelative            =  self%toleranceTopHat
       toleranceName                =  'toleranceTopHat'
    else
       powerSpectrumWindowFunction_ => self%powerSpectrumWindowFunction_
       toleranceRelative            =  self%tolerance
       toleranceName                =  'tolerance'
    end if
    wavenumberMaximum=min(1.0d3/topHatRadius,powerSpectrumWindowFunction_%wavenumberMaximum(smoothingMass))
    ! Split at any local minima of the transfer function. This is useful for oscillating transfer functions to allow integration
    ! to accurately capture each oscillation.
    call self%transferFunction_%wavenumbersLocalMinima(wavenumbersLocalMinima)
    rootVariance=Power_Spectrum_Variance_Integrate(                                                               &
         &                                         integrand                =integrand                          , &
         &                                         integrandLogarithmic     =integrandLogarithmic               , &
         &                                         wavenumberMinimum        =wavenumberMinimum                  , &
         &                                         wavenumberMaximum        =wavenumberMaximum                  , &
         &                                         radiusTopHat             =topHatRadius                       , &
         &                                         wavenumberHalfMode       =self%wavenumberHalfMode            , &
         &                                         wavenumbersLocalMinima   =wavenumbersLocalMinima             , &
         &                                         toleranceRelative        =toleranceRelative                  , &
         &                                         integrationFailureIsFatal=self%integrationFailureIsFatal     , &
         &                                         className                ='cosmologicalMassVariance'         , &
         &                                         classValue               =char(self%objectType(short=.true.)), &
         &                                         toleranceName            =trim(toleranceName)                  &
         &                                        )
    return

  contains

    double precision function integrand(wavenumber)
      !!{RST
      Integrand function used in computing the variance in the window function from the power spectrum.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber

      ! Return power spectrum multiplied by window function and volume element in k-space. Factors of 2 and π are included
      ! elsewhere.
      integrand=+  self%powerSpectrumPrimordialTransferred_%power(wavenumber              ,time) &
           &    *(                                                                               &
           &      +     powerSpectrumWindowFunction_       %value(wavenumber,smoothingMass,time) &
           &      *                                               wavenumber                     &
           &     )**2
      return
    end function integrand

    double precision function integrandLogarithmic(wavenumberLogarithmic)
      !!{RST
      Integrand function used in computing the variance in the window function from the power spectrum. This version integrates
      with respect to :math:`\log(k)`.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumberLogarithmic
      double precision                :: wavenumber

      ! Return power spectrum multiplied by window function and volume element in k-space. Factors of 2 and π are included
      ! elsewhere.
      wavenumber          =+exp(wavenumberLogarithmic)
      integrandLogarithmic=+  self%powerSpectrumPrimordialTransferred_%power(wavenumber              ,time) &
           &               *(                                                                               &
           &                 +     powerSpectrumWindowFunction_       %value(wavenumber,smoothingMass,time) &
           &                 *                                               wavenumber                     &
           &                )**2                                                                            &
           &               *                                                 wavenumber
      return
    end function integrandLogarithmic

  end function filteredPowerRootVarianceIntegral

  subroutine filteredPowerInterpolantsTime(self,time,i,h)
    !!{RST
    Compute interpolants in time.

    The position of the epoch along the axis is found as its coordinate on the *absolute* lattice---a quantity which depends
    only on the epoch and on the density of lattice points, never on which part of the lattice this particular tabulation
    spans---and is split there into the index of the lattice point below it and the fraction of the interval above that. Only
    then is the index of the first tabulated epoch subtracted, in exact integer arithmetic.

    The order matters. Forming the position relative to the first tabulated epoch *first* is exact in the subtraction, but the
    fractional part is then extracted from a number whose magnitude is the index within the tabulation, and so is rounded on a
    grid which coarsens as the tabulation grows: extending it would perturb every interpolated value in its last bits, which is
    exactly the dependence on the sequence of requests that pinning the range exists to remove.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    integer                                                , intent(  out) :: i
    double precision                                       , intent(  out) :: h
    double precision                                                       :: coordinate

    coordinate=log10(time)*dble(self%latticeTime%pointsPer)
    i         =floor(coordinate)
    h         =coordinate-dble(i)
    i         =i-self%latticeTime%indexMinimum+1
    if      (i == size(self%times)                        ) then
       ! Requested time must exactly equal the maximum tabulated time.
       i=size(self%times)-1
       h=1.0d0
    else if (i == 0                .and. h > 1.0d0-1.0d-6) then
       ! Requested time is the minimum tabulated time. The lattice coordinate of a lattice point need not evaluate to exactly
       ! its integer index, so the point which defines the lower bound of the tabulation can place itself a fraction of a step
       ! below it.
       i=1
       h=0.0d0
    else if (i <  1                                       ) then
       call Error_Report('interpolant out of range'//{introspection:location})
    end if
    return
  end subroutine filteredPowerInterpolantsTime

  logical function filteredPowerGrowthIsMassDependent(self)
    !!{RST
    Return true if the growth rate of the variance is mass-dependent.
    !!}
    implicit none
    class(cosmologicalMassVarianceFilteredPower), intent(inout) :: self

    filteredPowerGrowthIsMassDependent=self%growthIsMassDependent_
    return
  end function filteredPowerGrowthIsMassDependent

  subroutine filteredPowerFileRead(self)
    !!{RST
    Read tabulated data on mass variance from file.
    !!}
    use :: Display         , only : displayMessage           , verbosityLevelWorking
    use :: File_Utilities  , only : File_Exists
    use :: HDF5_Access     , only : hdf5Access
    use :: IO_HDF5         , only : hdf5File
    use :: Numerical_Ranges, only : gridSchemePerDecade
    use :: Table_Caches    , only : Table_Cache_Lattice_Read
    use :: Tables          , only : table1DLogarithmicCSpline, table1DLogarithmicMonotoneCSpline
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout)               :: self
    double precision                                       , dimension(:  ), allocatable :: massTmp        , timesTmp
    double precision                                       , dimension(:,:), allocatable :: rootVarianceTmp, rootVarianceUniqueTmp
    integer                                                , dimension(:  ), allocatable :: uniqueSizeTmp
    integer                                                , dimension(:,:), allocatable :: indexTmp
    type            (rangeLattice                         )                              :: latticeMass    , latticeTime
    logical                                                , dimension(:  ), allocatable :: isComputed
    integer                                                                              :: i              , useCache       , &
         &                                                                                  fileFormat

    latticeMass=rangeLattice()
    latticeTime=rangeLattice()
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
       timesTmp               =cachedVariances(useCache)%timesTmp
       massTmp                =cachedVariances(useCache)%massTmp
       rootVarianceTmp        =cachedVariances(useCache)%rootVarianceTmp
       rootVarianceUniqueTmp  =cachedVariances(useCache)%rootVarianceUniqueTmp
       indexTmp               =cachedVariances(useCache)%indexTmp
       uniqueSizeTmp          =cachedVariances(useCache)%uniqueSizeTmp
       latticeMass            =cachedVariances(useCache)%latticeMass
       latticeTime            =cachedVariances(useCache)%latticeTime
    end if
    !$omp end critical(cosmologicalMassVarianceFilteredPowerCache)
    if (useCache == 0) then
       ! Return if we are not using stored solutions.
       if (.not.self%storeTabulations) return
       ! Return if the file does not exist.
       if (.not.File_Exists(self%fileName)) return
       call displayMessage('reading σ(M) data from: '//self%fileName,verbosityLevelWorking)
       !$ call hdf5Access%set()
       hdf5ReadScope: block
         type(hdf5File  ) :: dataFile
         dataFile=hdf5File(self%fileName,overWrite=.false.,readOnly=.true.)
         ! Recover the lattices on which the stored tabulation was built. A file which does not record them, or which records
         ! lattices incommensurate with those this object would use, is simply ignored - and, since the file name carries a
         ! digest of this source file, such a file can in any case only be one written by a different build.
         ! Note that an undefined time lattice is a legitimate record - it is what is written where growth is mass-independent
         ! and sigma(M) is tabulated at a single epoch - so a missing record is told apart from that one through the mass axis,
         ! which is never undefined in a usable file.
         ! A file of any other format version is ignored by leaving the mass lattice undefined.
         fileFormat=0
         if (dataFile%hasAttribute('fileFormat')) call dataFile%readAttribute('fileFormat',fileFormat)
         if (fileFormat == filteredPowerFileFormatCurrent) then
            call Table_Cache_Lattice_Read(dataFile,'mass',gridSchemePerDecade,pointsPerDecade    ,latticeMass)
            call Table_Cache_Lattice_Read(dataFile,'time',gridSchemePerDecade,timePointsPerDecade,latticeTime)
         end if
         if (latticeMass%isDefined()) then
            call dataFile%readDataset  ('times'             ,     timesTmp             )
            call dataFile%readDataset  ('mass'              ,     massTmp              )
            call dataFile%readDataset  ('rootVariance'      ,     rootVarianceTmp      )
            call dataFile%readDataset  ('rootVarianceUnique',     rootVarianceUniqueTmp)
            call dataFile%readDataset  ('indexUnique'       ,     indexTmp             )
            call dataFile%readDataset  ('uniqueSize'        ,     uniqueSizeTmp        )
         end if
       end block hdf5ReadScope
       !$ call hdf5Access%unset()
       if (.not.latticeMass%isDefined()) return
       ! Reject a stored tabulation whose datasets do not match the lattices recorded alongside them, or whose time axis does not
       ! match the treatment of growth which this object requires.
       if     (                                                                                    &
            &   size(rootVarianceTmp,dim=1)  /= latticeMass%count                                  &
            &  .or.                                                                                &
            &   size(massTmp             )   /= latticeMass%count                                  &
            &  .or.                                                                                &
            &   size(rootVarianceTmp,dim=2)  /= size(timesTmp)                                     &
            &  .or.                                                                                &
            &   (     self%growthIsMassDependent_ .and. .not.latticeTime%isDefined()            )  &
            &  .or.                                                                                &
            &   (.not.self%growthIsMassDependent_ .and.      latticeTime%isDefined()            )  &
            &  .or.                                                                                &
            &   (     self%growthIsMassDependent_ .and.      latticeTime%count  /= size(timesTmp)) &
            &  .or.                                                                                &
            &   (.not.self%growthIsMassDependent_ .and.      size(timesTmp)     /= 1             ) &
            & ) return
       ! Cache this variance for possible later reuse.
       !$omp critical(cosmologicalMassVarianceFilteredPowerCache)
       lastCache=lastCache+1
       if (lastCache > sizeCache) lastCache=1
       countCache=max(countCache,lastCache)
       cachedVariances(lastCache)%fileName             =self%fileName
       cachedVariances(lastCache)%timesTmp             =     timesTmp
       cachedVariances(lastCache)%massTmp              =     massTmp
       cachedVariances(lastCache)%rootVarianceTmp      =     rootVarianceTmp
       cachedVariances(lastCache)%rootVarianceUniqueTmp=     rootVarianceUniqueTmp
       cachedVariances(lastCache)%indexTmp             =     indexTmp
       cachedVariances(lastCache)%uniqueSizeTmp        =     uniqueSizeTmp
       cachedVariances(lastCache)%latticeMass          =     latticeMass
       cachedVariances(lastCache)%latticeTime          =     latticeTime
       !$omp end critical(cosmologicalMassVarianceFilteredPowerCache)
    end if
    if (allocated(self%times                  )) deallocate(self%times                  )
    if (allocated(self%rootVarianceUniqueTable)) deallocate(self%rootVarianceUniqueTable)
    if (allocated(self%rootVarianceTable      )) then
       do i=1,size(self%rootVarianceTable)
          call self%rootVarianceTable(i)%destroy()
       end do
       deallocate(self%rootVarianceTable)
    end if
    allocate(self%times                  (size(timesTmp)))
    allocate(self%rootVarianceUniqueTable(size(timesTmp)))
    if (self%monotonicInterpolation) then
       allocate(table1DLogarithmicMonotoneCSpline :: self%rootVarianceTable(size(timesTmp)))
    else
       allocate(table1DLogarithmicCSpline         :: self%rootVarianceTable(size(timesTmp)))
    end if
    self%times      =timesTmp
    self%latticeMass=latticeMass
    self%latticeTime=latticeTime
    do i=1,size(self%times)
       allocate(self%rootVarianceUniqueTable(i)%rootVariance(uniqueSizeTmp(i)))
       allocate(self%rootVarianceUniqueTable(i)%index       (uniqueSizeTmp(i)))
       ! Build the table on the lattice recorded in the file, rather than by subdividing the range which that lattice spans, so
       ! that its abscissae are bit-identical to those of any other tabulation built on the same lattice - which is what allows
       ! the two to be spliced together exactly when the range is later extended.
       call self%rootVarianceTable(i)%extend  (latticeMass,isComputed)
       call self%rootVarianceTable(i)%populate(rootVarianceTmp(:,i)  )
       self%rootVarianceUniqueTable(i)%rootVariance=rootVarianceUniqueTmp(1:uniqueSizeTmp(i),i)
       self%rootVarianceUniqueTable(i)%index       =indexTmp             (1:uniqueSizeTmp(i),i)
    end do
    ! Recover the derived descriptions of the two axes. None of these is read from the file: each is a function of the lattices
    ! and of the epochs, so recomputing them here - exactly as `filteredPowerRetabulate` does - guarantees that a restored
    ! tabulation is described identically to a freshly built one.
    self%massMinimum=self%rootVarianceTable(1)%x( 1)
    self%massMaximum=self%rootVarianceTable(1)%x(-1)
    if (self%growthIsMassDependent_) then
       self%timeMinimum                =self%times(                 1)
       self%timeMaximum                =self%times(size(self%times)  )
       self%timeLogarithmicDeltaInverse=1.0d0/latticeTime%stepLogarithmic()
    else
       self%timeMinimum                =self%times(1)
       self%timeMaximum                =self%times(1)
       self%timeLogarithmicDeltaInverse=0.0d0
    end if
    self%initialized=.true.
    return
  end subroutine filteredPowerFileRead

  subroutine filteredPowerFileWrite(self)
    !!{RST
    Write tabulated data on mass variance to file.
    !!}
    use :: Display     , only : displayMessage           , verbosityLevelWorking
    use :: HDF5        , only : hsize_t
    use :: HDF5_Access , only : hdf5Access
    use :: IO_HDF5     , only : hdf5File
    use :: Table_Caches, only : Table_Cache_Lattice_Write
    implicit none
    class           (cosmologicalMassVarianceFilteredPower), intent(inout)               :: self
    double precision                                       , dimension(:  ), allocatable :: massTmp
    double precision                                       , dimension(:,:), allocatable :: rootVarianceTmp, rootVarianceUniqueTmp
    integer                                                , dimension(:  ), allocatable :: uniqueSizeTmp
    integer                                                , dimension(:,:), allocatable :: indexTmp
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
    cachedVariances(useCache)%fileName             =self%fileName
    cachedVariances(useCache)%timesTmp             =self%times
    cachedVariances(useCache)%massTmp              =     massTmp
    cachedVariances(useCache)%rootVarianceTmp      =     rootVarianceTmp
    cachedVariances(useCache)%rootVarianceUniqueTmp=     rootVarianceUniqueTmp
    cachedVariances(useCache)%indexTmp             =     indexTmp
    cachedVariances(useCache)%uniqueSizeTmp        =     uniqueSizeTmp
    cachedVariances(useCache)%latticeMass          =self%latticeMass
    cachedVariances(useCache)%latticeTime          =self%latticeTime
    !$omp end critical(cosmologicalMassVarianceFilteredPowerCache)
    ! Store to file if requested.
    if (self%storeTabulations) then
       call displayMessage('writing σ(M) data to: '//self%fileName,verbosityLevelWorking)
       ! Open the data file.
       !$ call hdf5Access%set()
       hdf5WriteScope: block
         type(hdf5File  ) :: dataFile
         dataFile=hdf5File(self%fileName,overWrite=.true.,objectsOverwritable=.true.,chunkSize=100_hsize_t,compressionLevel=9)
         call dataFile%writeDataset  (self%times                         ,'times'             )
         call dataFile%writeDataset  (     massTmp                       ,'mass'              )
         call dataFile%writeDataset  (     rootVarianceTmp               ,'rootVariance'      )
         call dataFile%writeDataset  (     rootVarianceUniqueTmp         ,'rootVarianceUnique')
         call dataFile%writeDataset  (     indexTmp                      ,'indexUnique'       )
         call dataFile%writeDataset  (     uniqueSizeTmp                 ,'uniqueSize'        )
         call dataFile%writeAttribute(     filteredPowerFileFormatCurrent,'fileFormat'        )
         ! Record the lattices on which the two axes are built. The bounds and interpolating factors which were formerly stored
         ! alongside them are not: each is a function of the lattices and of the epochs, and is recomputed when the file is read,
         ! so that a restored tabulation cannot come to be described differently from a freshly built one.
         call Table_Cache_Lattice_Write(dataFile,'mass',self%latticeMass)
         call Table_Cache_Lattice_Write(dataFile,'time',self%latticeTime)
       end block hdf5WriteScope
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
    !!{RST
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

  subroutine filteredPowerDescriptor(self,descriptor,includeClass,includeFileModificationTimes,parameterName)
    !!{RST
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters, descriptorParameterName
    implicit none
    class    (cosmologicalMassVarianceFilteredPower), intent(inout)           :: self
    type     (inputParameters                      ), intent(inout)           :: descriptor
    logical                                         , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=*                                ), intent(in   ), optional :: parameterName
    character(len=18                               )                          :: parameterLabel
    type     (inputParameters                      )                          :: parameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter(descriptorParameterName('cosmologicalMassVariance',parameterName),'filteredPower')
    parameters=descriptor%subparameters(descriptorParameterName('cosmologicalMassVariance',parameterName))
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
    call    self%powerSpectrumWindowFunction_               %descriptor(parameters,includeClass,includeFileModificationTimes)
    return
  end subroutine filteredPowerDescriptor
