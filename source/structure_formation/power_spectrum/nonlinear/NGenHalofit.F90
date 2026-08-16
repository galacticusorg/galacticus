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
Implements a nonlinear power spectrum class in which the nonlinear power spectrum is computed using the :term:`NGenHalofit` model of :cite:t:`smith_precision_2019`.
!!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: File_Utilities            , only : lockDescriptor
  use :: Numerical_Interpolation   , only : interpolator2D
  use :: Power_Spectra             , only : powerSpectrumClass
  use :: Power_Spectra_Primordial  , only : powerSpectrumPrimordialClass

  !![
  <powerSpectrumNonlinear name="powerSpectrumNonlinearNGenHalofit" docformat="rst">
   <description>
   Provides a nonlinear power spectrum class in which the nonlinear power spectrum is computed as :math:`P_\mathrm{NL}(k,t) =
   R(k,z) P_\mathrm{L}(k,t)`, where :math:`R(k,z)` is the ratio of nonlinear to linear power computed by the :term:`NGenHalofit`
   model of :cite:t:`smith_precision_2019`, and :math:`P_\mathrm{L}(k,t)` is the linear theory power spectrum. NGenHalofit
   recalibrates the halofit model against the Daemmerung suite of :math:`N`-body simulations, and supplements it with a 2-loop
   Multi-point Propagator Theory calculation on large scales. The NGenHalofit code, and the Cuba integration library on which it
   depends, will be downloaded, compiled and run as necessary if this option is utilized.

   The linear theory power spectrum used by NGenHalofit is that computed by Galacticus, and so is fully consistent with the
   linear power spectrum to which the ratio is applied. The dependence on cosmological parameters is captured by a second-order
   Taylor expansion about the Planck-2013-like fiducial cosmology of the Daemmerung suite---the range over which that expansion
   is calibrated is enforced (see the class source for the precise limits). The model is quoted as accurate to
   :math:`\lesssim3`\ % for :math:`k \le 9 h`/Mpc and :math:`z \le 3`. It is tabulated here to :math:`z=30`, but note that
   NGenHalofit itself transitions smoothly back to the halofit model of :cite:t:`takahashi_revising_2012` for :math:`a &lt; 0.25`
   (i.e. :math:`z > 3`), and for :math:`k` outside the range :math:`0.05 \le k \le 10 h`/Mpc. Below the smallest tabulated
   wavenumber the ratio is held fixed at its value there, while above the largest tabulated wavenumber it is extrapolated as a
   power law in :math:`k`.
   </description>
  </powerSpectrumNonlinear>
  !!]
  type, extends(powerSpectrumNonlinearClass) :: powerSpectrumNonlinearNGenHalofit
     !!{RST
     A nonlinear power spectrum class in which the nonlinear power spectrum is computed using the :term:`NGenHalofit` model of :cite:t:`smith_precision_2019`.
     !!}
     private
     double precision                               , allocatable, dimension(:  ) :: wavenumberTable                    , redshiftTable
     double precision                               , allocatable, dimension(:,:) :: ratioTable
     type            (interpolator2D               ), allocatable                 :: interpolator_
     double precision                                                             :: amplitudeScalar
     logical                                                                      :: amplitudeScalarIsSet               , tabulated
     type            (lockDescriptor               )                              :: fileLock                           , installLock
     class           (cosmologyFunctionsClass      ), pointer                     :: cosmologyFunctions_       => null()
     class           (cosmologyParametersClass     ), pointer                     :: cosmologyParameters_      => null()
     class           (powerSpectrumClass           ), pointer                     :: powerSpectrum_            => null()
     class           (powerSpectrumPrimordialClass ), pointer                     :: powerSpectrumPrimordial_  => null()
     class           (cosmologicalMassVarianceClass), pointer                     :: cosmologicalMassVariance_ => null()
   contains
     !![
     <methods docformat="rst">
       <method description="Tabulate the ratio of nonlinear to linear power, :math:`R(k,z)`, by running the NGenHalofit code." method="tabulate"/>
     </methods>
     !!]
     final     ::             nGenHalofitDestructor
     procedure :: value    => nGenHalofitValue
     procedure :: tabulate => nGenHalofitTabulate
  end type powerSpectrumNonlinearNGenHalofit

  interface powerSpectrumNonlinearNGenHalofit
     !!{RST
     Constructors for the :galacticus-class:`powerSpectrumNonlinearNGenHalofit` nonlinear power spectrum class.
     !!}
     module procedure nGenHalofitConstructorParameters
     module procedure nGenHalofitConstructorInternal
  end interface powerSpectrumNonlinearNGenHalofit

  ! Range of wavenumbers (in h/Mpc), and number of points spaced uniformly in ln(k), at which the ratio of nonlinear to linear
  ! power is tabulated.
  double precision                 , parameter :: wavenumberMinimum          =1.0d-3, wavenumberMaximum=1.0d+2
  integer                          , parameter :: countWavenumber            =400
  ! Maximum redshift at which the ratio is tabulated, along with the number of redshifts spaced uniformly in ln(1+z). The
  ! tabulation is limited to the range spanned by the Daemmerung simulation snapshots (which extend to a=0.0321). NGenHalofit
  ! emulates all redshifts in a single invocation, so a fine grid is inexpensive. With this number of redshifts the error
  ! introduced by interpolating in redshift is at most 0.21% (attained close to z=2.6), which is small compared to the ~1%
  ! precision of the model itself. Note that the interpolation error does not fall as the square of the grid spacing---the
  ! residual arises from small-scale structure in the redshift dependence of the ratio, not from smooth curvature---so little is
  ! gained by refining this grid further.
  double precision                 , parameter :: redshiftMaximum            =30.0d0
  integer                          , parameter :: countRedshift              =91
  ! Range of wavenumbers (in h/Mpc), and number of points spaced uniformly in ln(k), at which the linear theory power spectrum is
  ! written for input to NGenHalofit. This must not exceed NPOWTARGET=1000 in the NGenHalofit sources.
  double precision                 , parameter :: wavenumberLinearMinimum    =1.0d-4, wavenumberLinearMaximum=1.0d+3
  integer                          , parameter :: countWavenumberLinear      =701
  ! Tolerance within which a slightly-negative redshift arising from round-off in the conversion from cosmic time at the present
  ! epoch is treated as zero.
  double precision                 , parameter :: toleranceRedshift          =1.0d-6
  ! Tolerance within which the dark energy equation of state must match the CPL form, w(a)=w₀+w_a(1-a), assumed by NGenHalofit.
  double precision                 , parameter :: toleranceEquationState     =1.0d-3
  ! The expansion factor at which the CPL form of the dark energy equation of state is tested.
  double precision                 , parameter :: expansionFactorTest        =0.5d+0
  ! Tolerance within which the cosmology must be flat.
  double precision                 , parameter :: toleranceFlatness          =1.0d-3
  ! The pivot wavenumber (in Mpc¯¹) at which the primordial power spectrum index and its running are evaluated. This matches the
  ! `pivot_scalar` used for the CAMB linear power spectra against which NGenHalofit was calibrated.
  double precision                 , parameter :: wavenumberPivot            =0.05d0
  ! Logarithmic offset in wavenumber used to evaluate the running of the primordial power spectrum index by finite difference.
  double precision                 , parameter :: wavenumberOffsetLogarithmic=1.0d-2
  ! The amplitude of the primordial scalar power spectrum in the fiducial Daemmerung cosmology. NGenHalofit expects the amplitude
  ! to be given in units of this value.
  double precision                 , parameter :: amplitudeScalarFiducial    =2.14818d-9
  ! The fiducial cosmological parameters of the Daemmerung simulation suite, the step sizes by which each was varied, and the
  ! multiple of the step size beyond which NGenHalofit switches off the Taylor expansion in that parameter (see
  ! `setFidParAndDeltaPar()` and the `KILLTAYLORSTRICT` block of `NGenHalofit()` in the NGenHalofit sources). Beyond these limits
  ! NGenHalofit silently discards the *entire* cosmological parameter dependence---note that its own limits are one-sided for the
  ! parameters which are positive-definite, whereas the symmetric limits enforced here bound the region over which the Taylor
  ! expansion is actually calibrated.
  double precision                 , parameter, dimension(8) :: parameterFiducial  =[-1.00000d0,0.0d0,0.6928849d0,0.11889000d0,0.02216100d0,0.961100d0,1.0d0,0.00d0]
  double precision                 , parameter, dimension(8) :: parameterDelta     =[+0.10000d0,0.2d0,0.0500000d0,0.00594450d0,0.00110805d0,0.048055d0,0.1d0,0.01d0]
  double precision                 , parameter, dimension(8) :: parameterKillTaylor=[+0.50000d0,2.0d0,2.0000000d0,2.00000000d0,2.00000000d0,2.000000d0,2.0d0,2.00d0]
  character       (len=19         ), parameter, dimension(8) :: parameterName      =[                                                                   &
       &                                                                             "w_0                ","w_a                ","OmegaDarkEnergy    ", &
       &                                                                             "omegaColdDarkMatter","omegaBaryon        ","index              ", &
       &                                                                             "amplitudeScalar    ","running            "                        &
       &                                                                            ]

  ! Generate a source digest.
  !![
  <sourceDigest name="nGenHalofitSourceDigest"/>
  !!]

contains

  function nGenHalofitConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`powerSpectrumNonlinearNGenHalofit` nonlinear power spectrum class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumNonlinearNGenHalofit)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass         ), pointer       :: cosmologyParameters_
    class           (powerSpectrumClass               ), pointer       :: powerSpectrum_
    class           (powerSpectrumPrimordialClass     ), pointer       :: powerSpectrumPrimordial_
    class           (cosmologicalMassVarianceClass    ), pointer       :: cosmologicalMassVariance_
    double precision                                                   :: amplitudeScalar

    ! Check and read parameters.
    if (parameters%isPresent('amplitudeScalar')) then
       !![
       <inputParameter docformat="rst">
         <name>amplitudeScalar</name>
         <source>parameters</source>
         <description>
         The amplitude of the primordial scalar power spectrum, :math:`A_\mathrm{s}`, such that :math:`P_\chi(k) = A_\mathrm{s}
         (k/k_\mathrm{s0})^{n_\mathrm{s}-1}` with :math:`k_\mathrm{s0}=0.05` Mpc\ :math:`^{-1}`, to be passed to NGenHalofit. If
         not specified, :math:`A_\mathrm{s}` is computed from :math:`\sigma_8` using the ratio :math:`\sigma_8^2/A_\mathrm{s}`
         evaluated by :term:`CLASS`.
         </description>
       </inputParameter>
       !!]
    end if
    ! Construct required objects.
    !![
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="powerSpectrum"            name="powerSpectrum_"            source="parameters"/>
    <objectBuilder class="powerSpectrumPrimordial"  name="powerSpectrumPrimordial_"  source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    ! Call the internal constructor.
    !![
    <conditionalCall>
     <call>
      self=powerSpectrumNonlinearNGenHalofit(                                                     &amp;
       &amp;                                 cosmologyFunctions_      =cosmologyFunctions_      , &amp;
       &amp;                                 cosmologyParameters_     =cosmologyParameters_     , &amp;
       &amp;                                 powerSpectrum_           =powerSpectrum_           , &amp;
       &amp;                                 powerSpectrumPrimordial_ =powerSpectrumPrimordial_ , &amp;
       &amp;                                 cosmologicalMassVariance_=cosmologicalMassVariance_  &amp;
       &amp;                                 {conditions}                                         &amp;
       &amp;                                )
     </call>
     <argument name="amplitudeScalar" value="amplitudeScalar" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="powerSpectrum_"           />
    <objectDestructor name="powerSpectrumPrimordial_" />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function nGenHalofitConstructorParameters

  function nGenHalofitConstructorInternal(cosmologyFunctions_,cosmologyParameters_,powerSpectrum_,powerSpectrumPrimordial_,cosmologicalMassVariance_,amplitudeScalar) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`powerSpectrumNonlinearNGenHalofit` nonlinear power spectrum class.
    !!}
    use :: Error               , only : Error_Report
    use :: Numerical_Comparison, only : Values_Differ
    implicit none
    type            (powerSpectrumNonlinearNGenHalofit)                          :: self
    class           (cosmologyFunctionsClass          ), intent(in   ), target   :: cosmologyFunctions_
    class           (cosmologyParametersClass         ), intent(in   ), target   :: cosmologyParameters_
    class           (powerSpectrumClass               ), intent(in   ), target   :: powerSpectrum_
    class           (powerSpectrumPrimordialClass     ), intent(in   ), target   :: powerSpectrumPrimordial_
    class           (cosmologicalMassVarianceClass    ), intent(in   ), target   :: cosmologicalMassVariance_
    double precision                                   , intent(in   ), optional :: amplitudeScalar
    double precision                                                             :: equationOfState0        , equationOfStateA, &
         &                                                                          equationOfStateTest
    !![
    <constructorAssign variables="*cosmologyFunctions_, *cosmologyParameters_, *powerSpectrum_, *powerSpectrumPrimordial_, *cosmologicalMassVariance_"/>
    !!]

    ! Initialize state.
    self%tabulated=.false.
    if (present(amplitudeScalar)) then
       self%amplitudeScalar     = amplitudeScalar
       self%amplitudeScalarIsSet=.true.
    else
       self%amplitudeScalar     =-1.0d0
       self%amplitudeScalarIsSet=.false.
    end if
    ! Check that this is a flat cosmology. NGenHalofit assumes flatness, taking Ω_m = 1 - Ω_DE and deriving the Hubble parameter
    ! as h = √[(ω_c+ω_b)/Ω_m] - if the cosmology is not flat that derived value will not match the Hubble parameter used by
    ! Galacticus to compute the linear power spectrum supplied to NGenHalofit.
    if     (                                                                                      &
         &  Values_Differ(                                                                        &
         &                +self%cosmologyParameters_%OmegaMatter    ()                            &
         &                +self%cosmologyParameters_%OmegaDarkEnergy()                         ,  &
         &                       1.0d+0                                                        ,  &
         &                absTol=toleranceFlatness                                                &
         &               )                                                                        &
         & )                                                                                      &
         & call Error_Report(                                                                     &
         &                   'this method is applicable only to flat matter+dark energy models'// &
         &                    {introspection:location}                                            &
         &                  )
    ! Check that the dark energy equation of state follows the CPL form, w(a)=w₀+w_a(1-a), assumed by NGenHalofit. The coefficient
    ! w_a is identified by matching the derivative of w(a) at the present epoch, so the CPL form is exact for the
    ! `cosmologyFunctionsMatterDarkEnergyCPL` class (and for a constant equation of state), but not, in general, for the
    ! w(a)=w₀+w₁a(1-a) form of `cosmologyFunctionsMatterDarkEnergy`.
    equationOfState0   =self%cosmologyFunctions_%equationOfStateDarkEnergy(expansionFactor=1.0d0              )
    equationOfStateA   =nGenHalofitEquationOfStateA(self)
    equationOfStateTest=self%cosmologyFunctions_%equationOfStateDarkEnergy(expansionFactor=expansionFactorTest)
    if     (                                                                                                                   &
         &  Values_Differ(                                                                                                     &
         &                +equationOfStateTest                                                                               , &
         &                +equationOfState0                                                                                    &
         &                +equationOfStateA                                                                                    &
         &                *(1.0d0-expansionFactorTest)                                                                       , &
         &                absTol=toleranceEquationState                                                                        &
         &               )                                                                                                     &
         & )                                                                                                                   &
         & call Error_Report(                                                                                                  &
         &                   'this method is applicable only to dark energy equations of state of the form w(a)=w₀+w_a(1-a)'// &
         &                    {introspection:location}                                                                         &
         &                  )
    return
  end function nGenHalofitConstructorInternal

  subroutine nGenHalofitDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`powerSpectrumNonlinearNGenHalofit` nonlinear power spectrum class.
    !!}
    implicit none
    type(powerSpectrumNonlinearNGenHalofit), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%powerSpectrum_"           />
    <objectDestructor name="self%powerSpectrumPrimordial_" />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine nGenHalofitDestructor

  double precision function nGenHalofitEquationOfStateA(self) result(equationOfStateA)
    !!{RST
    Return the coefficient :math:`w_\mathrm{a}` of the CPL dark energy equation of state, :math:`w(a)=w_0+w_\mathrm{a}(1-a)`,
    found by matching the derivative :math:`\mathrm{d}w/\mathrm{d}a` at the present epoch.
    !!}
    implicit none
    class           (powerSpectrumNonlinearNGenHalofit), intent(inout) :: self
    double precision                                   , parameter     :: expansionFactorOffset=1.0d-3

    equationOfStateA=-(                                                                                                 &
         &             +self%cosmologyFunctions_%equationOfStateDarkEnergy(expansionFactor=1.0d0                      ) &
         &             -self%cosmologyFunctions_%equationOfStateDarkEnergy(expansionFactor=1.0d0-expansionFactorOffset) &
         &            )                                                                                                 &
         &           /expansionFactorOffset
    return
  end function nGenHalofitEquationOfStateA

  double precision function nGenHalofitRedshift(iRedshift) result(redshift)
    !!{RST
    Return the redshift of the ``iRedshift``\ :sup:`th` point in the grid at which the ratio of nonlinear to linear power is
    tabulated. The grid is uniform in :math:`\ln(1+z)`.
    !!}
    implicit none
    integer, intent(in   ) :: iRedshift

    redshift=+exp(                            &
         &        +log(1.0d0+redshiftMaximum) &
         &        *dble(iRedshift    -1)      &
         &        /dble(countRedshift-1)      &
         &       )                            &
         &   -1.0d0
    return
  end function nGenHalofitRedshift

  subroutine nGenHalofitTabulate(self)
    !!{RST
    Tabulate the ratio of nonlinear to linear power, :math:`R(k,z)`, by downloading, building, and running the
    :term:`NGenHalofit` code as necessary.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsLittleH
    use :: Dependencies        , only : dependencyVersion            , dependencyVersionLabel
    use :: Display             , only : displayMessage               , verbosityLevelWorking
    use :: Error               , only : Error_Report
    use :: File_Utilities      , only : Directory_Make               , File_Exists           , File_Lock        , File_Unlock , &
         &                              File_Name_Temporary
    use :: Hashes_Cryptographic, only : Hash_MD5
    use :: Input_Paths         , only : inputPath                    , pathTypeDataDynamic   , pathTypeTools
    use :: Interfaces_CLASS    , only : Interface_CLASS_Normalization
    use :: ISO_Varying_String  , only : varying_string               , assignment(=)         , char             , operator(//)
    use :: String_Handling     , only : String_C_To_Fortran
    use :: System_Command      , only : System_Command_Do            , shellEscape
    use :: System_Compilers    , only : compiler                     , compilerValidate      , languageC
    use :: System_Download     , only : download
    implicit none
    class           (powerSpectrumNonlinearNGenHalofit), intent(inout)               :: self
    double precision                                   , dimension(8)                :: parameterValue
    double precision                                   , allocatable, dimension(:  ) :: ratioRow            , wavenumberRead
    double precision                                   , allocatable, dimension(:,:) :: ratioRead
    type            (varying_string                   )                              :: ratioFile           , uniqueLabel          , &
         &                                                                              nGenHalofitVersion  , nGenHalofitPath      , &
         &                                                                              cubaVersion         , cubaPath             , &
         &                                                                              escapedArchiveFile  , escapedPath          , &
         &                                                                              escapedCompiler     , escapedCubaPath      , &
         &                                                                              escapedCubaLibs     , cubaLibraries        , &
         &                                                                              escapedParameters   , escapedExpansion     , &
         &                                                                              outputPath          , command              , &
         &                                                                              fileNameParameters  , fileNameExpansion    , &
         &                                                                              escapedOutputPath   , fileNamePowerSpectrum
    character       (len=32                           )                              :: parameterLabel      , fileLabel
    integer                                                                          :: ratioUnit           , iParameter           , &
         &                                                                              iRedshift           , iWavenumber          , &
         &                                                                              status              , parametersUnit       , &
         &                                                                              expansionUnit       , powerSpectrumUnit    , &
         &                                                                              outputUnit          , countWavenumberRead  , &
         &                                                                              countRedshiftRead
    double precision                                                                 :: hubbleConstant      , timePresent          , &
         &                                                                              wavenumber          , powerSpectrumLinear  , &
         &                                                                              powerLinear         , powerNonlinear       , &
         &                                                                              powerHalofit

    ! Assemble the cosmological parameters required by NGenHalofit. Note that NGenHalofit assumes a flat cosmology with massless
    ! neutrinos and no radiation, so the cold dark matter density is simply the total matter density less the baryon density.
    hubbleConstant   =self%cosmologyParameters_%HubbleConstant(units=hubbleUnitsLittleH)
    parameterValue(1)=self%cosmologyFunctions_ %equationOfStateDarkEnergy(expansionFactor=1.0d0)
    parameterValue(2)=nGenHalofitEquationOfStateA(self)
    parameterValue(3)=   self%cosmologyParameters_%OmegaDarkEnergy()
    parameterValue(4)=+(                                             &
         &              +self%cosmologyParameters_%OmegaMatter    () &
         &              -self%cosmologyParameters_%OmegaBaryon    () &
         &             )                                             &
         &            *hubbleConstant**2
    parameterValue(5)=+  self%cosmologyParameters_%OmegaBaryon    () &
         &            *hubbleConstant**2
    parameterValue(6)=self%powerSpectrumPrimordial_%logarithmicDerivative(wavenumberPivot)
    ! The running of the primordial power spectrum index, dn_s/dln k, evaluated by finite difference at the pivot wavenumber.
    parameterValue(8)=+(                                                                                                        &
         &              +self%powerSpectrumPrimordial_%logarithmicDerivative(wavenumberPivot*exp(+wavenumberOffsetLogarithmic)) &
         &              -self%powerSpectrumPrimordial_%logarithmicDerivative(wavenumberPivot*exp(-wavenumberOffsetLogarithmic)) &
         &             )                                                                                                        &
         &            /2.0d0                                                                                                    &
         &            /wavenumberOffsetLogarithmic
    ! Determine the amplitude of the primordial scalar power spectrum. Unless it was specified explicitly, invert the ratio
    ! σ₈²/A_s computed by CLASS - this is exact, as the linear power spectrum is directly proportional to A_s. NGenHalofit expects
    ! the amplitude in units of that of the fiducial Daemmerung cosmology.
    if (self%amplitudeScalarIsSet) then
       parameterValue(7)=+self%amplitudeScalar         &
            &            /     amplitudeScalarFiducial
    else
       parameterValue(7)=+self%cosmologicalMassVariance_%sigma8(                         )**2 &
            &            /Interface_CLASS_Normalization        (self%cosmologyParameters_)    &
            &            /amplitudeScalarFiducial
    end if
    ! Validate that the cosmological parameters fall within the range over which the Taylor expansion used by NGenHalofit to
    ! capture the dependence on cosmology is calibrated. Beyond these limits NGenHalofit switches that expansion off entirely -
    ! and does so for all parameters as soon as any one of them is out of range - silently returning a result which retains no
    ! dependence on cosmology beyond that of the linear power spectrum itself.
    do iParameter=1,8
       if     (                                                                                                     &
            &   abs(parameterValue(iParameter)-parameterFiducial(iParameter))                                       &
            &  >                                                                                                    &
            &   parameterKillTaylor(iParameter)*parameterDelta(iParameter)                                          &
            & ) then
          write (parameterLabel,'(e12.6)') parameterValue(iParameter)
          call Error_Report(                                                                                        &
               &            'parameter "'//trim(parameterName(iParameter))//'" = '//trim(adjustl(parameterLabel))// &
               &            ' is outside the range over which NGenHalofit is calibrated'                         // &
               &            {introspection:location}                                                                &
               &           )
       end if
    end do
    ! Build the name of the file to which the tabulated ratio is cached. The cache key includes the versions of NGenHalofit and
    ! Cuba, and a digest of this source file, so that a ratio computed by an earlier version of either code, or using a different
    ! tabulation grid, is not silently reused.
    nGenHalofitVersion=dependencyVersion("ngenhalofit")
    cubaVersion       =dependencyVersion("cuba"       )
    nGenHalofitPath   =inputPath(pathTypeTools)//"NGenHalofit-"//nGenHalofitVersion
    cubaPath          =inputPath(pathTypeTools)//"Cuba-"       //cubaVersion
    uniqueLabel       ="_sourceDigest:"//String_C_To_Fortran(nGenHalofitSourceDigest)//dependencyVersionLabel("ngenhalofit")//dependencyVersionLabel("cuba")
    do iParameter=1,8
       write (parameterLabel,'(e17.10)') parameterValue(iParameter)
       uniqueLabel=uniqueLabel//"_"//trim(parameterName(iParameter))//":"//trim(adjustl(parameterLabel))
    end do
    call Directory_Make(inputPath(pathTypeDataDynamic)//"largeScaleStructure")
    ratioFile=inputPath(pathTypeDataDynamic)//"largeScaleStructure/powerSpectrumRatioNGenHalofit_"//Hash_MD5(uniqueLabel)//".txt"
    ! Check for existence of the tabulated ratio, building it if necessary.
    call File_Lock(ratioFile,self%fileLock,lockIsShared=.true.)
    if (.not.File_Exists(ratioFile)) then
       call File_Unlock(self%fileLock)
       call File_Lock(ratioFile,self%fileLock,lockIsShared=.false.)
       ! NGenHalofit locates its simulation data, and both reads and writes its tabulations of the halofit effective spectral
       ! parameters, using paths relative to the current working directory. It must therefore be run from its own directory, and
       ! only one process may do so at a time - otherwise concurrent runs for different cosmologies would overwrite each other's
       ! `DATAEffectiveSpectra/GenerateEffectiveSpec.Target.dat`. Take an exclusive lock on the installation directory for the
       ! duration of the build and the run.
       call Directory_Make(inputPath(pathTypeTools))
       call File_Lock(nGenHalofitPath,self%installLock,lockIsShared=.false.)
       ! Check for presence of the executable.
       if (.not.File_Exists(nGenHalofitPath//"/NGenHalofit.exe")) then
          call compilerValidate(languageC,'NGenHalofit')
          ! Build the Cuba library, on which the Multi-point Propagator Theory module of NGenHalofit depends.
          if (.not.File_Exists(cubaPath//"/libcuba.a")) then
             if (.not.File_Exists(cubaPath//"/configure")) then
                if (.not.File_Exists(cubaPath//".tar.gz")) then
                   call displayMessage("downloading Cuba library....",verbosityLevelWorking)
                   call download("https://feynarts.de/cuba/Cuba-"//char(cubaVersion)//".tar.gz",char(cubaPath)//".tar.gz",status=status,retries=5,retryWait=60)
                   if (status /= 0 .or. .not.File_Exists(cubaPath//".tar.gz")) &
                        & call Error_Report("failed to download Cuba library"//{introspection:location})
                end if
                call displayMessage("unpacking Cuba library....",verbosityLevelWorking)
                call Directory_Make(cubaPath)
                escapedArchiveFile=cubaPath//".tar.gz"
                escapedArchiveFile=shellEscape(escapedArchiveFile)
                escapedCubaPath   =shellEscape(cubaPath          )
                call System_Command_Do("tar -xzf "//escapedArchiveFile//" -C "//escapedCubaPath//" --strip-components=1",status)
                if (status /= 0 .or. .not.File_Exists(cubaPath//"/configure")) &
                     & call Error_Report("failed to unpack Cuba library"//{introspection:location})
             end if
             ! Build the library. The `gnu11` standard is imposed as Cuba declares its own `bool` type, which conflicts with the
             ! `bool` keyword introduced in C23 - the default standard for GCC 15 and later.
             call displayMessage("compiling Cuba library....",verbosityLevelWorking)
             escapedCubaPath=shellEscape(cubaPath     )
             escapedCompiler=compiler   (languageC    )
             call System_Command_Do(                                                  &
                  &                 "cd "//escapedCubaPath//"; ./configure"        // &
                  &                 " CC='"//escapedCompiler//"'"                  // &
                  &                 " CFLAGS='-O2 -fomit-frame-pointer -std=gnu11'"// &
                  &                 " && make lib"                                  , &
                  &                 status                                            &
                  &                )
             if (status /= 0 .or. .not.File_Exists(cubaPath//"/libcuba.a")) &
                  & call Error_Report("failed to build Cuba library"//{introspection:location})
          end if
          ! Check for presence of the NGenHalofit source code.
          if (.not.File_Exists(nGenHalofitPath//"/main.c")) then
             if (.not.File_Exists(nGenHalofitPath//".tar.gz")) then
                call displayMessage("downloading NGenHalofit code....",verbosityLevelWorking)
                call download("https://bitbucket.org/ngenhalofitteam/ngenhalofitpublic/get/"//char(nGenHalofitVersion)//".tar.gz",char(nGenHalofitPath)//".tar.gz",status=status,retries=5,retryWait=60)
                if (status /= 0 .or. .not.File_Exists(nGenHalofitPath//".tar.gz")) &
                     & call Error_Report("failed to download NGenHalofit code"//{introspection:location})
             end if
             call displayMessage("unpacking NGenHalofit code....",verbosityLevelWorking)
             call Directory_Make(nGenHalofitPath)
             escapedArchiveFile=nGenHalofitPath//".tar.gz"
             escapedArchiveFile=shellEscape(escapedArchiveFile)
             escapedPath       =shellEscape(nGenHalofitPath   )
             call System_Command_Do("tar -xzf "//escapedArchiveFile//" -C "//escapedPath//" --strip-components=1",status)
             if (status /= 0 .or. .not.File_Exists(nGenHalofitPath//"/main.c")) &
                  & call Error_Report("failed to unpack NGenHalofit code"//{introspection:location})
          end if
          ! Build the code. The compiler, optimization flags, and the GSL and Cuba include and library paths are all supplied on
          ! the `make` command line - the NGenHalofit makefile hard-codes paths which do not, in general, exist. Note that the
          ! GSL paths are expanded by the shell, so must be enclosed in double (not single) quotes.
          call displayMessage("compiling NGenHalofit code....",verbosityLevelWorking)
          escapedPath    =shellEscape(nGenHalofitPath)
          cubaLibraries  ="-L"//cubaPath//" -lcuba"
          escapedCubaPath=shellEscape(cubaPath       )
          escapedCubaLibs=shellEscape(cubaLibraries  )
          escapedCompiler=compiler   (languageC      )
          call System_Command_Do(                                                                &
               &                 "cd "//escapedPath//"; gslPrefix=`gsl-config --prefix`; make"// &
               &                 " CC='"//escapedCompiler//"'"                                // &
               &                 " OPTIMIZE=-O2"                                              // &
               &                 ' GSL_INCL="-I${gslPrefix}/include"'                         // &
               &                 ' GSL_LIBS="-L${gslPrefix}/lib -lgsl -lgslcblas"'            // &
               &                 " CUBA_INCL=-I"//escapedCubaPath                             // &
               &                 " CUBA_LIBS="//escapedCubaLibs                               ,  &
               &                 status                                                          &
               &                )
          if (status /= 0 .or. .not.File_Exists(nGenHalofitPath//"/NGenHalofit.exe")) &
               & call Error_Report("failed to build NGenHalofit code"//{introspection:location})
       end if
       ! Generate the tabulations of the halofit effective spectral parameters for the Daemmerung simulations. These depend only
       ! on the simulation cosmologies, not on the target cosmology, and so are computed once and cached in the installation
       ! directory - exactly as done by the `run.sh` script supplied with NGenHalofit.
       if (.not.File_Exists(nGenHalofitPath//"/DATAEffectiveSpectra/GenerateEffectiveSpec.iVAR_16.dat")) then
          call displayMessage("initializing NGenHalofit effective spectra....",verbosityLevelWorking)
          escapedPath=shellEscape(nGenHalofitPath)
          call System_Command_Do("cd "//escapedPath//"; ./NGenHalofit.exe parameterfiles/setup.dat parameterfiles/setup_expansion.dat",status)
          if (status /= 0 .or. .not.File_Exists(nGenHalofitPath//"/DATAEffectiveSpectra/GenerateEffectiveSpec.iVAR_16.dat")) &
               & call Error_Report("failed to initialize NGenHalofit effective spectra"//{introspection:location})
       end if
       ! Write the input files for this cosmology to a temporary directory.
       outputPath=File_Name_Temporary("nGenHalofitOutput",char(inputPath(pathTypeDataDynamic)//"largeScaleStructure"))
       call Directory_Make(outputPath)
       fileNameParameters   =outputPath//"/parameters.txt"
       fileNameExpansion    =outputPath//"/expansionFactors.txt"
       fileNamePowerSpectrum=outputPath//"/powerSpectrumLinear.txt"
       ! Write the linear theory power spectrum at the present day, converting from the Mpc¯¹ and Mpc³ units used by Galacticus to
       ! the h/Mpc and (Mpc/h)³ units expected by NGenHalofit.
       timePresent=self%cosmologyFunctions_%cosmicTime(1.0d0)
       open(newUnit=powerSpectrumUnit,file=char(fileNamePowerSpectrum),status='unknown',form='formatted')
       do iWavenumber=1,countWavenumberLinear
          wavenumber=+wavenumberLinearMinimum            &
               &     *(                                  &
               &       +wavenumberLinearMaximum          &
               &       /wavenumberLinearMinimum          &
               &      )**(                               &
               &          +dble(iWavenumber          -1) &
               &          /dble(countWavenumberLinear-1) &
               &         )
          powerSpectrumLinear=self%powerSpectrum_%power(wavenumber*hubbleConstant,timePresent)
          write (powerSpectrumUnit,'(2e20.12)') wavenumber,powerSpectrumLinear*hubbleConstant**3
       end do
       close(powerSpectrumUnit)
       ! Write the list of expansion factors at which the power spectrum is required.
       open(newUnit=expansionUnit,file=char(fileNameExpansion),status='unknown',form='formatted')
       do iRedshift=1,countRedshift
          write (expansionUnit,'(e20.12)') 1.0d0/(1.0d0+nGenHalofitRedshift(iRedshift))
       end do
       close(expansionUnit)
       ! Write the NGenHalofit parameter file.
       open(newUnit=parametersUnit,file=char(fileNameParameters),status='unknown',form='formatted')
       write (parametersUnit,'(a,e20.12)') "aexp              ",1.0d0
       write (parametersUnit,'(a,e20.12)') "w0                ",parameterValue(1)
       write (parametersUnit,'(a,e20.12)') "w1                ",parameterValue(2)
       write (parametersUnit,'(a,e20.12)') "om_DE0            ",parameterValue(3)
       write (parametersUnit,'(a,e20.12)') "om_ch20           ",parameterValue(4)
       write (parametersUnit,'(a,e20.12)') "om_bh20           ",parameterValue(5)
       write (parametersUnit,'(a,e20.12)') "pindex            ",parameterValue(6)
       write (parametersUnit,'(a,e20.12)') "As                ",parameterValue(7)
       write (parametersUnit,'(a,e20.12)') "running           ",parameterValue(8)
       write (parametersUnit,'(a,i8    )') "nPkOut            ",countWavenumber
       write (parametersUnit,'(a,e20.12)') "rkOutMIN          ",wavenumberMinimum
       write (parametersUnit,'(a,e20.12)') "rkOutMAX          ",wavenumberMaximum
       write (parametersUnit,'(a,i8    )') "iLogOrLin         ",1
       write (parametersUnit,'(a,i8    )') "iGenEffSpecTarget ",1
       write (parametersUnit,'(a,i8    )') "iGenEffSpecVar    ",0
       write (parametersUnit,'(a,a     )') "OutputDir         ",char(outputPath)//"/"
       write (parametersUnit,'(a,a     )') "OutputFileBase    ","powerSpectrum"
       write (parametersUnit,'(a,a     )') "OutputMPTFileBase ","powerSpectrumMPT"
       write (parametersUnit,'(a,a     )') "PowDirTarget      ",char(outputPath)
       write (parametersUnit,'(a,a     )') "PowFileTarget     ","powerSpectrumLinear.txt"
       write (parametersUnit,'(a,a     )') "PowDirVar         ","DATAPlinCAMB/"
       write (parametersUnit,'(a,a     )') "PowFileBaseVar    ","Planck2013.Step_ByHand.HighAcc_matterpower"
       close(parametersUnit)
       ! Run NGenHalofit from its own directory.
       escapedPath      =shellEscape(nGenHalofitPath   )
       escapedParameters=shellEscape(fileNameParameters)
       escapedExpansion =shellEscape(fileNameExpansion )
       command="cd "//escapedPath//"; ./NGenHalofit.exe "//escapedParameters//" "//escapedExpansion
       call System_Command_Do(char(command),status)
       fileLabel=nGenHalofitFileIndex(countRedshift-1)
       if (status /= 0 .or. .not.File_Exists(outputPath//"/powerSpectrum."//trim(fileLabel)//".dat")) &
            & call Error_Report("failed to run NGenHalofit code"//{introspection:location})
       ! Release the lock on the installation directory - all further work uses only the temporary output directory.
       call File_Unlock(self%installLock)
       ! Read the results and store the ratio of nonlinear to linear power to the cache file. NGenHalofit writes one file per
       ! expansion factor, each containing the wavenumber, the linear power spectrum, the halofit power spectrum, and the
       ! NGenHalofit power spectrum. Taking the ratio of the last to the second of these removes the linear growth factor, so that
       ! the result is independent of any small difference between the linear growth computed by NGenHalofit and by Galacticus.
       allocate(wavenumberRead(countWavenumber              ))
       allocate(ratioRead     (countWavenumber,countRedshift))
       do iRedshift=1,countRedshift
          fileLabel=nGenHalofitFileIndex(iRedshift-1)
          open(newUnit=outputUnit,file=char(outputPath)//"/powerSpectrum."//trim(fileLabel)//".dat",status='old',form='formatted')
          do iWavenumber=1,countWavenumber
             read (outputUnit,*) wavenumber,powerLinear,powerHalofit,powerNonlinear
             wavenumberRead(iWavenumber          )=wavenumber
             ratioRead     (iWavenumber,iRedshift)=powerNonlinear/powerLinear
          end do
          close(outputUnit)
       end do
       open(newUnit=ratioUnit,file=char(ratioFile),status='unknown',form='formatted')
       write (ratioUnit,'(2i8)') countWavenumber,countRedshift
       do iWavenumber=1,countWavenumber
          write (ratioUnit,'(1000e24.16)') wavenumberRead(iWavenumber),ratioRead(iWavenumber,:)
       end do
       close(ratioUnit)
       deallocate(wavenumberRead)
       deallocate(ratioRead     )
       ! Remove the temporary directory. NGenHalofit writes a large number of files into this directory (two per expansion factor,
       ! plus records of the parameter values used), so the directory is removed recursively rather than by removing each file
       ! individually - this avoids the run failing at this final stage if a future version of NGenHalofit should write some
       ! additional file.
       escapedOutputPath=shellEscape(outputPath)
       call System_Command_Do("rm -rf "//escapedOutputPath,status)
       if (status /= 0) call Error_Report("failed to remove temporary directory '"//outputPath//"'"//{introspection:location})
    end if
    ! Read the tabulated ratio.
    if (allocated(self%wavenumberTable)) deallocate(self%wavenumberTable)
    if (allocated(self%redshiftTable  )) deallocate(self%redshiftTable  )
    if (allocated(self%ratioTable     )) deallocate(self%ratioTable     )
    allocate(self%wavenumberTable(countWavenumber              ))
    allocate(self%redshiftTable  (                countRedshift))
    allocate(self%ratioTable     (countWavenumber,countRedshift))
    allocate(     ratioRow       (                countRedshift))
    open(newUnit=ratioUnit,file=char(ratioFile),status='old',form='formatted')
    read (ratioUnit,*) countWavenumberRead,countRedshiftRead
    if (countWavenumberRead /= countWavenumber .or. countRedshiftRead /= countRedshift) &
         & call Error_Report("cached NGenHalofit tabulation has unexpected dimensions"//{introspection:location})
    do iWavenumber=1,countWavenumber
       read (ratioUnit,*) self%wavenumberTable(iWavenumber),ratioRow
       self%ratioTable(iWavenumber,:)=ratioRow
    end do
    close(ratioUnit)
    deallocate(ratioRow)
    call File_Unlock(self%fileLock)
    ! Convert wavenumbers from the h/Mpc units used by NGenHalofit to the Mpc¯¹ used by Galacticus, and tabulate in terms of ln(k)
    ! and ln(1+z), interpolating in ln(R).
    do iRedshift=1,countRedshift
       self%redshiftTable(iRedshift)=log(1.0d0+nGenHalofitRedshift(iRedshift))
    end do
    self%wavenumberTable=log(self%wavenumberTable*hubbleConstant)
    self%ratioTable     =log(self%ratioTable                    )
    ! Build the interpolator.
    if (allocated(self%interpolator_)) deallocate(self%interpolator_)
    allocate(self%interpolator_)
    self%interpolator_=interpolator2D(self%wavenumberTable,self%redshiftTable,self%ratioTable)
    self%tabulated    =.true.
    return
  end subroutine nGenHalofitTabulate

  function nGenHalofitFileIndex(indexFile) result(label)
    !!{RST
    Return the zero-based index used by NGenHalofit to label the output file for each expansion factor, formatted as a string.
    !!}
    implicit none
    character(len=32)                :: label
    integer          , intent(in   ) :: indexFile

    write (label,'(i0)') indexFile
    return
  end function nGenHalofitFileIndex

  double precision function nGenHalofitValue(self,waveNumber,time) result(powerSpectrum)
    !!{RST
    Return a nonlinear power spectrum computed by applying the ratio of nonlinear to linear power computed by the
    :term:`NGenHalofit` model of :cite:t:`smith_precision_2019` to the linear theory power spectrum.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (powerSpectrumNonlinearNGenHalofit), intent(inout) :: self
    double precision                                   , intent(in   ) :: time                   , waveNumber
    double precision                                                   :: redshift               , ratio              , &
         &                                                                wavenumberLogarithmic  , redshiftLogarithmic, &
         &                                                                ratioLogarithmic       , gradient           , &
         &                                                                ratioLogarithmicMaximum

    ! Tabulate the ratio of nonlinear to linear power if necessary.
    if (.not.self%tabulated) then
       !$omp critical(powerSpectrumNonlinearNGenHalofitTabulate)
       if (.not.self%tabulated) call self%tabulate()
       !$omp end critical(powerSpectrumNonlinearNGenHalofitTabulate)
    end if
    ! Find the redshift corresponding to this time. Slightly-negative redshifts arising from round-off in the conversion from
    ! cosmic time at the present epoch are treated as zero.
    redshift=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    if (redshift < 0.0d0 .and. redshift > -toleranceRedshift) redshift=0.0d0
    if (redshift < 0.0d0 .or.  redshift >  redshiftMaximum  )                                         &
         & call Error_Report(                                                                         &
         &                   'redshift is outside the range tabulated for NGenHalofit (0 ≤ z ≤ 30)'// &
         &                    {introspection:location}                                                &
         &                  )
    ! Clamp to the tabulated range - the end points of the tabulation are computed from the same expression as the redshifts
    ! passed to NGenHalofit, so can differ from the value computed here at the level of round-off.
    redshiftLogarithmic  =min(                                       &
         &                    max(                                   &
         &                        log(1.0d0+redshift)              , &
         &                        self%redshiftTable(            1)  &
         &                       )                                 , &
         &                        self%redshiftTable(countRedshift)  &
         &                   )
    wavenumberLogarithmic=log(wavenumber)
    ! Evaluate the ratio of nonlinear to linear power. Below the smallest tabulated wavenumber the ratio is held at its value at
    ! that wavenumber - it is very close to unity there in any case - while above the largest tabulated wavenumber it is
    ! extrapolated using the logarithmic gradient at that wavenumber.
    if      (wavenumberLogarithmic <  self%wavenumberTable(                          1)) then
       ratioLogarithmic       =self%interpolator_%interpolate(self%wavenumberTable(                          1),redshiftLogarithmic)
    else if (wavenumberLogarithmic >  self%wavenumberTable(size(self%wavenumberTable) )) then
       ratioLogarithmicMaximum=self%interpolator_%interpolate(self%wavenumberTable(size(self%wavenumberTable) ),redshiftLogarithmic)
       gradient               =+(                                                                                                        &
            &                    +ratioLogarithmicMaximum                                                                                &
            &                    -self%interpolator_%interpolate(self%wavenumberTable(size(self%wavenumberTable)-1),redshiftLogarithmic) &
            &                   )                                                                                                        &
            &                  /(                                                                                                        &
            &                    +self%wavenumberTable(size(self%wavenumberTable)  )                                                     &
            &                    -self%wavenumberTable(size(self%wavenumberTable)-1)                                                     &
            &                   )
       ratioLogarithmic       =+ratioLogarithmicMaximum                                                                                  &
            &                  +gradient                                                                                                 &
            &                  *(                                                                                                        &
            &                    +wavenumberLogarithmic                                                                                  &
            &                    -self%wavenumberTable(size(self%wavenumberTable))                                                       &
            &                   )
    else
       ratioLogarithmic       =self%interpolator_%interpolate(wavenumberLogarithmic                                ,redshiftLogarithmic)
    end if
    ratio        =exp(ratioLogarithmic)
    powerSpectrum=+ratio                                    &
         &        *self%powerSpectrum_%power(wavenumber,time)
    return
  end function nGenHalofitValue
