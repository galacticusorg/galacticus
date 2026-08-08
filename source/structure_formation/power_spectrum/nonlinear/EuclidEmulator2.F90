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

!+    Contributions to this file made by: Andrew Benson. The interface to EuclidEmulator2 requested in issue #59 was
!+    drafted with assistance from Claude, and reviewed and verified by Andrew Benson.

!!{RST
Implements a nonlinear power spectrum class in which the nonlinear power spectrum is computed by applying the nonlinear correction factor emulated by the :term:`EuclidEmulator2` code of :cite:t:`euclid_collaboration_euclid_2021` to the linear theory power spectrum.
!!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: File_Utilities            , only : lockDescriptor
  use :: Numerical_Interpolation   , only : interpolator2D
  use :: Power_Spectra             , only : powerSpectrumClass
  use :: Power_Spectra_Primordial  , only : powerSpectrumPrimordialClass

  !![
  <powerSpectrumNonlinear name="powerSpectrumNonlinearEuclidEmulator2" docformat="rst">
   <description>
   Provides a nonlinear power spectrum class in which the nonlinear power spectrum is computed as :math:`P_\mathrm{NL}(k,t) =
   B(k,z) P_\mathrm{L}(k,t)`, where :math:`B(k,z)` is the nonlinear correction factor (or "boost factor") emulated by the
   :term:`EuclidEmulator2` code of :cite:t:`euclid_collaboration_euclid_2021`, and :math:`P_\mathrm{L}(k,t)` is the linear theory
   power spectrum. The EuclidEmulator2 code will be downloaded, compiled and run as necessary if this option is utilized. The
   emulation is applicable only for :math:`0 \le z \le 10`; below the smallest emulated wavenumber, :math:`k = 8.73 \times 10^{-3}
   h`/Mpc, the correction factor is held fixed at its value there, while above the largest emulated wavenumber, :math:`k = 9.41
   h`/Mpc, it is extrapolated as a power law in :math:`k`.
   </description>
  </powerSpectrumNonlinear>
  !!]
  type, extends(powerSpectrumNonlinearClass) :: powerSpectrumNonlinearEuclidEmulator2
     !!{RST
     A nonlinear power spectrum class in which the nonlinear power spectrum is computed using the :term:`EuclidEmulator2` code of :cite:t:`euclid_collaboration_euclid_2021`.
     !!}
     private
     double precision                               , allocatable, dimension(:  ) :: wavenumberTable                   , redshiftTable
     double precision                               , allocatable, dimension(:,:) :: boostTable
     type            (interpolator2D               ), allocatable                 :: interpolator_
     double precision                                                             :: amplitudeScalar
     logical                                                                      :: amplitudeScalarIsSet              , tabulated
     type            (lockDescriptor               )                              :: fileLock
     class           (cosmologyFunctionsClass      ), pointer                     :: cosmologyFunctions_       => null()
     class           (cosmologyParametersClass     ), pointer                     :: cosmologyParameters_      => null()
     class           (powerSpectrumClass           ), pointer                     :: powerSpectrum_            => null()
     class           (powerSpectrumPrimordialClass ), pointer                     :: powerSpectrumPrimordial_  => null()
     class           (cosmologicalMassVarianceClass), pointer                     :: cosmologicalMassVariance_ => null()
   contains
     !![
     <methods docformat="rst">
       <method description="Tabulate the nonlinear correction factor, :math:`B(k,z)`, by running the EuclidEmulator2 code." method="tabulate"/>
     </methods>
     !!]
     final     ::             euclidEmulator2Destructor
     procedure :: value    => euclidEmulator2Value
     procedure :: tabulate => euclidEmulator2Tabulate
  end type powerSpectrumNonlinearEuclidEmulator2

  interface powerSpectrumNonlinearEuclidEmulator2
     !!{RST
     Constructors for the :galacticus-class:`powerSpectrumNonlinearEuclidEmulator2` nonlinear power spectrum class.
     !!}
     module procedure euclidEmulator2ConstructorParameters
     module procedure euclidEmulator2ConstructorInternal
  end interface powerSpectrumNonlinearEuclidEmulator2

  ! Wavenumber range used for testing the shape of the primordial power spectrum.
  double precision                , parameter :: wavenumberLong        =0.01d0, wavenumberShort=1.0d0
  ! Maximum redshift emulated by EuclidEmulator2, along with the number of redshifts, spaced uniformly in ln(1+z), at which the
  ! nonlinear correction factor is tabulated. EuclidEmulator2 emulates all redshifts in a single invocation, so a fine grid costs
  ! essentially nothing.
  double precision                , parameter :: redshiftMaximum       =10.0d0
  integer                         , parameter :: countRedshift         =101
  ! Tolerance within which a slightly-negative redshift arising from round-off in the conversion from cosmic time at the present
  ! epoch is treated as zero.
  double precision                , parameter :: toleranceRedshift     =1.0d-6
  ! Tolerance within which the dark energy equation of state must match the CPL form, w(a)=w₀+w_a(1-a), assumed by
  ! EuclidEmulator2.
  double precision                , parameter :: toleranceEquationState=1.0d-3
  ! The expansion factor at which the CPL form of the dark energy equation of state is tested.
  double precision                , parameter :: expansionFactorTest   =0.5d+0
  ! Ranges of the cosmological parameters over which EuclidEmulator2 is applicable (see Table 1 of Euclid Collaboration (2021;
  ! https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.2840E). Note that, while w_a is drawn from the range [-0.7,+0.7],
  ! EuclidEmulator2 itself rejects values exceeding +0.5.
  double precision                , parameter, dimension(8) :: parameterMinimum=[0.04d+0,0.24d0,0.00d0,0.92d0,0.61d0,-1.3d0,-0.7d0,1.7d-9]
  double precision                , parameter, dimension(8) :: parameterMaximum=[0.06d+0,0.40d0,0.15d0,1.00d0,0.73d0,-0.7d0,+0.5d0,2.5d-9]
  character       (len=18        ), parameter, dimension(8) :: parameterName   =[                                                                &
       &                                                                         "OmegaBaryon       ","OmegaMatter       ","neutrinoMassSummed", &
       &                                                                         "index             ","HubbleConstant    ","w_0               ", &
       &                                                                         "w_a               ","amplitudeScalar   "                       &
       &                                                                        ]

  ! Generate a source digest.
  !![
  <sourceDigest name="euclidEmulator2SourceDigest"/>
  !!]

contains

  function euclidEmulator2ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`powerSpectrumNonlinearEuclidEmulator2` nonlinear power spectrum class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumNonlinearEuclidEmulator2)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass              ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass             ), pointer       :: cosmologyParameters_
    class           (powerSpectrumClass                   ), pointer       :: powerSpectrum_
    class           (powerSpectrumPrimordialClass         ), pointer       :: powerSpectrumPrimordial_
    class           (cosmologicalMassVarianceClass        ), pointer       :: cosmologicalMassVariance_
    double precision                                                       :: amplitudeScalar

    ! Check and read parameters.
    if (parameters%isPresent('amplitudeScalar')) then
       !![
       <inputParameter docformat="rst">
         <name>amplitudeScalar</name>
         <source>parameters</source>
         <description>
         The amplitude of the primordial scalar power spectrum, :math:`A_\mathrm{s}`, such that :math:`P_\chi(k) = A_\mathrm{s}
         (k/k_\mathrm{s0})^{n_\mathrm{s}-1}` with :math:`k_\mathrm{s0}=0.05` Mpc\ :math:`^{-1}`, to be passed to
         EuclidEmulator2. If not specified, :math:`A_\mathrm{s}` is computed from :math:`\sigma_8` using the ratio
         :math:`\sigma_8^2/A_\mathrm{s}` evaluated by :term:`CLASS`.
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
      self=powerSpectrumNonlinearEuclidEmulator2(                                                     &amp;
       &amp;                                     cosmologyFunctions_      =cosmologyFunctions_      , &amp;
       &amp;                                     cosmologyParameters_     =cosmologyParameters_     , &amp;
       &amp;                                     powerSpectrum_           =powerSpectrum_           , &amp;
       &amp;                                     powerSpectrumPrimordial_ =powerSpectrumPrimordial_ , &amp;
       &amp;                                     cosmologicalMassVariance_=cosmologicalMassVariance_  &amp;
       &amp;                                     {conditions}                                         &amp;
       &amp;                                    )
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
  end function euclidEmulator2ConstructorParameters

  function euclidEmulator2ConstructorInternal(cosmologyFunctions_,cosmologyParameters_,powerSpectrum_,powerSpectrumPrimordial_,cosmologicalMassVariance_,amplitudeScalar) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`powerSpectrumNonlinearEuclidEmulator2` nonlinear power spectrum class.
    !!}
    use :: Error               , only : Error_Report
    use :: Numerical_Comparison, only : Values_Differ
    implicit none
    type            (powerSpectrumNonlinearEuclidEmulator2)                          :: self
    class           (cosmologyFunctionsClass              ), intent(in   ), target   :: cosmologyFunctions_
    class           (cosmologyParametersClass             ), intent(in   ), target   :: cosmologyParameters_
    class           (powerSpectrumClass                   ), intent(in   ), target   :: powerSpectrum_
    class           (powerSpectrumPrimordialClass         ), intent(in   ), target   :: powerSpectrumPrimordial_
    class           (cosmologicalMassVarianceClass        ), intent(in   ), target   :: cosmologicalMassVariance_
    double precision                                       , intent(in   ), optional :: amplitudeScalar
    double precision                                                                 :: equationOfState0         , equationOfStateA, &
         &                                                                              equationOfStateTest
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
    ! Check that this is a flat cosmology - EuclidEmulator2 assumes that the dark energy density parameter is given by
    ! Ω_DE = 1 - Ω_m - Ω_γ - Ω_ν.
    if     (                                                                                      &
         &  Values_Differ(                                                                        &
         &                +self%cosmologyParameters_%OmegaMatter    ()                            &
         &                +self%cosmologyParameters_%OmegaDarkEnergy()                         ,  &
         &                       1.0d+0                                                        ,  &
         &                absTol=1.0d-3                                                           &
         &               )                                                                        &
         & )                                                                                      &
         & call Error_Report(                                                                     &
         &                   'this method is applicable only to flat matter+dark energy models'// &
         &                    {introspection:location}                                            &
         &                  )
    ! Check that the primordial power spectrum has no running of the spectral index.
    if     (                                                                                                    &
         &  Values_Differ(                                                                                      &
         &                self%powerSpectrumPrimordial_%logarithmicDerivative(wavenumberShort),                 &
         &                self%powerSpectrumPrimordial_%logarithmicDerivative(wavenumberLong ),                 &
         &                relTol=1.0d-3                                                                         &
         &               )                                                                                      &
         & )                                                                                                    &
         & call Error_Report(                                                                                   &
         &                   'this method is applicable only to models with no running of the spectral index'// &
         &                    {introspection:location}                                                          &
         &                  )
    ! Check that the dark energy equation of state follows the CPL form, w(a)=w₀+w_a(1-a), assumed by EuclidEmulator2. The
    ! coefficient w_a is identified by matching the derivative of w(a) at the present epoch, so the CPL form is exact for a
    ! constant equation of state, but not, in general, for the w(a)=w₀+w₁a(1-a) form of `cosmologyFunctionsMatterDarkEnergy`.
    equationOfState0   =self%cosmologyFunctions_%equationOfStateDarkEnergy(expansionFactor=1.0d0              )
    equationOfStateA   =euclidEmulator2EquationOfStateA(self)
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
  end function euclidEmulator2ConstructorInternal

  subroutine euclidEmulator2Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`powerSpectrumNonlinearEuclidEmulator2` nonlinear power spectrum class.
    !!}
    implicit none
    type(powerSpectrumNonlinearEuclidEmulator2), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%powerSpectrum_"           />
    <objectDestructor name="self%powerSpectrumPrimordial_" />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine euclidEmulator2Destructor

  double precision function euclidEmulator2EquationOfStateA(self) result(equationOfStateA)
    !!{RST
    Return the coefficient :math:`w_\mathrm{a}` of the CPL dark energy equation of state, :math:`w(a)=w_0+w_\mathrm{a}(1-a)`,
    found by matching the derivative :math:`\mathrm{d}w/\mathrm{d}a` at the present epoch.
    !!}
    implicit none
    class           (powerSpectrumNonlinearEuclidEmulator2), intent(inout) :: self
    double precision                                       , parameter     :: expansionFactorOffset=1.0d-3

    equationOfStateA=-(                                                                                                 &
         &             +self%cosmologyFunctions_%equationOfStateDarkEnergy(expansionFactor=1.0d0                      ) &
         &             -self%cosmologyFunctions_%equationOfStateDarkEnergy(expansionFactor=1.0d0-expansionFactorOffset) &
         &            )                                                                                                 &
         &           /expansionFactorOffset
    return
  end function euclidEmulator2EquationOfStateA

  subroutine euclidEmulator2Tabulate(self)
    !!{RST
    Tabulate the nonlinear correction factor, :math:`B(k,z)`, by downloading, building, and running the :term:`EuclidEmulator2`
    code as necessary.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsLittleH
    use :: Dependencies        , only : dependencyVersion            , dependencyVersionLabel
    use :: Display             , only : displayMessage               , verbosityLevelWorking
    use :: Error               , only : Error_Report
    use :: File_Utilities      , only : Count_Lines_In_File          , Directory_Make        , Directory_Remove , File_Exists , &
         &                              File_Lock                    , File_Name_Temporary   , File_Unlock
    use :: Hashes_Cryptographic, only : Hash_MD5
    use :: Input_Paths         , only : inputPath                    , pathTypeDataDynamic   , pathTypeTools
    use :: Interfaces_CLASS    , only : Interface_CLASS_Normalization
    use :: ISO_Varying_String  , only : varying_string               , assignment(=)         , char             , operator(//)
    use :: String_Handling     , only : String_C_To_Fortran
    use :: System_Command      , only : System_Command_Do            , shellEscape
    use :: System_Compilers    , only : compiler                     , compilerValidate      , languageCPlusPlus
    use :: System_Download     , only : download
    implicit none
    class           (powerSpectrumNonlinearEuclidEmulator2), intent(inout)                :: self
    double precision                                       , dimension(8)                 :: parameterValue
    double precision                                       , allocatable  , dimension(:)  :: boostRow
    type            (varying_string                       )                               :: boostFile             , uniqueLabel        , &
         &                                                                                   euclidEmulator2Version, euclidEmulator2Path, &
         &                                                                                   escapedZipFile        , escapedToolsPath   , &
         &                                                                                   escapedPath           , escapedOutputPath  , &
         &                                                                                   escapedBoostFile      , escapedCompiler    , &
         &                                                                                   outputPath            , command            , &
         &                                                                                   arguments
    character       (len=32                               )                               :: parameterLabel
    integer                                                                               :: boostUnit             , iParameter         , &
         &                                                                                   iRedshift             , iWavenumber        , &
         &                                                                                   status
    double precision                                                                      :: hubbleConstant

    ! Assemble the cosmological parameters required by EuclidEmulator2. Note that Galacticus does not currently model massive
    ! neutrinos, so the summed neutrino mass is set to zero, and the matter density parameter is passed through unchanged (in
    ! EuclidEmulator2 the matter density parameter is the total, including neutrinos).
    hubbleConstant    =self%cosmologyParameters_%HubbleConstant(units=hubbleUnitsLittleH)
    parameterValue(1)=self%cosmologyParameters_ %OmegaBaryon              (                     )
    parameterValue(2)=self%cosmologyParameters_ %OmegaMatter              (                     )
    parameterValue(3)=0.0d0
    parameterValue(4)=self%powerSpectrumPrimordial_%logarithmicDerivative (wavenumberShort      )
    parameterValue(5)=hubbleConstant
    parameterValue(6)=self%cosmologyFunctions_  %equationOfStateDarkEnergy(expansionFactor=1.0d0)
    parameterValue(7)=euclidEmulator2EquationOfStateA(self)
    ! Determine the amplitude of the primordial scalar power spectrum. Unless it was specified explicitly, invert the ratio
    ! σ₈²/A_s computed by CLASS - this is exact, as the linear power spectrum is directly proportional to A_s.
    if (self%amplitudeScalarIsSet) then
       parameterValue(8)=self%amplitudeScalar
    else
       parameterValue(8)=+self%cosmologicalMassVariance_%sigma8(                         )**2 &
            &            /Interface_CLASS_Normalization        (self%cosmologyParameters_)
    end if
    ! Validate that the cosmological parameters fall within the range over which EuclidEmulator2 is applicable. EuclidEmulator2
    ! itself simply aborts if they do not, so check here to be able to report a useful error.
    do iParameter=1,8
       if     (                                                                                                     &
            &   parameterValue(iParameter) < parameterMinimum(iParameter)                                           &
            &  .or.                                                                                                 &
            &   parameterValue(iParameter) > parameterMaximum(iParameter)                                           &
            & ) then
          write (parameterLabel,'(e12.6)') parameterValue(iParameter)
          call Error_Report(                                                                                        &
               &            'parameter "'//trim(parameterName(iParameter))//'" = '//trim(adjustl(parameterLabel))// &
               &            ' is outside the range over which EuclidEmulator2 is applicable'                     // &
               &            {introspection:location}                                                                &
               &           )
       end if
    end do
    ! Build the name of the file to which the tabulated nonlinear correction factor is cached. The cache key includes the version
    ! of EuclidEmulator2 and a digest of this source file, so that a correction factor emulated by an earlier version of
    ! EuclidEmulator2, or using a different tabulation grid, is not silently reused.
    euclidEmulator2Version=dependencyVersion("euclidemulator2")
    euclidEmulator2Path   =inputPath(pathTypeTools)//"EuclidEmulator2-"//euclidEmulator2Version
    uniqueLabel           ="_sourceDigest:"//String_C_To_Fortran(euclidEmulator2SourceDigest)//dependencyVersionLabel("euclidemulator2")
    do iParameter=1,8
       write (parameterLabel,'(e17.10)') parameterValue(iParameter)
       uniqueLabel=uniqueLabel//"_"//trim(parameterName(iParameter))//":"//trim(adjustl(parameterLabel))
    end do
    call Directory_Make(inputPath(pathTypeDataDynamic)//"largeScaleStructure")
    boostFile=inputPath(pathTypeDataDynamic)//"largeScaleStructure/powerSpectrumBoostEuclidEmulator2_"//Hash_MD5(uniqueLabel)//".txt"
    ! Check for existence of the tabulated correction factor, building it if necessary.
    call File_Lock(boostFile,self%fileLock,lockIsShared=.true.)
    if (.not.File_Exists(boostFile)) then
       call File_Unlock(self%fileLock)
       call File_Lock(boostFile,self%fileLock,lockIsShared=.false.)
       ! Check for presence of the executable.
       if (.not.File_Exists(euclidEmulator2Path//"/ee2.exe")) then
          call compilerValidate(languageCPlusPlus,'EuclidEmulator2')
          ! Check for presence of the source code.
          if (.not.File_Exists(euclidEmulator2Path//"/src/main.cxx")) then
             ! Download the code.
             if (.not.File_Exists(euclidEmulator2Path//".zip")) then
                call displayMessage("downloading EuclidEmulator2 code....",verbosityLevelWorking)
                call download("https://github.com/miknab/EuclidEmulator2/archive/v"//char(euclidEmulator2Version)//".zip",char(euclidEmulator2Path)//".zip",status=status,retries=5,retryWait=60)
                if (status /= 0 .or. .not.File_Exists(euclidEmulator2Path//".zip")) &
                     & call Error_Report("failed to download EuclidEmulator2 code"//{introspection:location})
             end if
             ! Unpack the code.
             call displayMessage("unpacking EuclidEmulator2 code....",verbosityLevelWorking)
             escapedZipFile  =euclidEmulator2Path//".zip"
             escapedZipFile  =shellEscape(escapedZipFile   )
             escapedToolsPath=inputPath  (pathTypeTools    )
             escapedToolsPath=shellEscape(escapedToolsPath )
             call System_Command_Do("unzip -o "//escapedZipFile//" -d "//escapedToolsPath,status)
             if (status /= 0 .or. .not.File_Exists(euclidEmulator2Path//"/src/main.cxx")) &
                  & call Error_Report("failed to unpack EuclidEmulator2 code"//{introspection:location})
          end if
          ! Build the code. The GSL include and library paths are supplied on the `make` command line - the EuclidEmulator2
          ! makefile hard-codes paths that do not, in general, exist. The `cstdint` and `cmath` headers are force-included as
          ! the EuclidEmulator2 sources (and the version of `cxxopts` that they bundle) rely on those headers being included
          ! transitively, which is no longer the case for GCC 13 and later.
          call displayMessage("compiling EuclidEmulator2 code....",verbosityLevelWorking)
          escapedPath    =shellEscape(euclidEmulator2Path)
          escapedCompiler=compiler   (languageCPlusPlus  )
          call System_Command_Do(                                                                &
               &                 "cd "//escapedPath//"; make"                                 // &
               &                 " CC='"//escapedCompiler//" -include cstdint -include cmath'"// &
               &                 " I_GSL=-I`gsl-config --prefix`/include"                     // &
               &                 " LIBPATH=-L`gsl-config --prefix`/lib"                        , &
               &                 status                                                          &
               &                )
          if (status /= 0 .or. .not.File_Exists(euclidEmulator2Path//"/ee2.exe")) &
               & call Error_Report("failed to build EuclidEmulator2 code"//{introspection:location})
       end if
       ! Build the list of command line arguments. The redshifts are spaced uniformly in ln(1+z) over the full range emulated by
       ! EuclidEmulator2, and are all evaluated in a single invocation.
       arguments=''
       write (parameterLabel,'(f9.6 )') parameterValue(1)
       arguments=arguments//" -b "//trim(adjustl(parameterLabel))
       write (parameterLabel,'(f9.6 )') parameterValue(2)
       arguments=arguments//" -m "//trim(adjustl(parameterLabel))
       write (parameterLabel,'(f9.6 )') parameterValue(3)
       arguments=arguments//" -s "//trim(adjustl(parameterLabel))
       write (parameterLabel,'(f9.6 )') parameterValue(4)
       arguments=arguments//" -n "//trim(adjustl(parameterLabel))
       write (parameterLabel,'(f9.6 )') parameterValue(5)
       arguments=arguments//" -H "//trim(adjustl(parameterLabel))
       write (parameterLabel,'(f9.6 )') parameterValue(6)
       arguments=arguments//" -W "//trim(adjustl(parameterLabel))
       write (parameterLabel,'(f9.6 )') parameterValue(7)
       arguments=arguments//" -w "//trim(adjustl(parameterLabel))
       write (parameterLabel,'(e15.8)') parameterValue(8)
       arguments=arguments//" -A "//trim(adjustl(parameterLabel))
       do iRedshift=1,countRedshift
          write (parameterLabel,'(f12.8)') euclidEmulator2Redshift(iRedshift)
          arguments=arguments//" -z "//trim(adjustl(parameterLabel))
       end do
       ! Run EuclidEmulator2. It must be run from its own directory as it opens its data file using a hard-coded relative path -
       ! if that file can not be opened it will segfault rather than reporting an error. Output is written to a temporary
       ! directory, from which the resulting file is moved to the cache file.
       outputPath=File_Name_Temporary("euclidEmulator2Output",char(inputPath(pathTypeDataDynamic)//"largeScaleStructure"))
       call Directory_Make(outputPath)
       escapedPath      =shellEscape(euclidEmulator2Path)
       escapedOutputPath=shellEscape(outputPath         )
       escapedBoostFile =shellEscape(boostFile          )
       command="cd "//escapedPath//"; ./ee2.exe"//arguments//" -d "//escapedOutputPath//"/ -o boost"
       call System_Command_Do(char(command),status)
       if (status /= 0 .or. .not.File_Exists(outputPath//"/boost0.dat")) &
            & call Error_Report("failed to run EuclidEmulator2 code"//{introspection:location})
       escapedOutputPath=shellEscape(outputPath//"/boost0.dat")
       call System_Command_Do("mv "//escapedOutputPath//" "//escapedBoostFile,status)
       if (status /= 0 .or. .not.File_Exists(boostFile)) &
            & call Error_Report("failed to store EuclidEmulator2 output"//{introspection:location})
       call Directory_Remove(outputPath)
    end if
    ! Read the tabulated nonlinear correction factor.
    if (allocated(self%wavenumberTable)) deallocate(self%wavenumberTable)
    if (allocated(self%redshiftTable  )) deallocate(self%redshiftTable  )
    if (allocated(self%boostTable     )) deallocate(self%boostTable     )
    allocate(self%wavenumberTable(Count_Lines_In_File(boostFile,comment_char='#')              ))
    allocate(self%redshiftTable  (                                                countRedshift))
    allocate(self%boostTable     (size(self%wavenumberTable)                     ,countRedshift))
    allocate(boostRow            (                                                countRedshift))
    open(newUnit=boostUnit,file=char(boostFile),status='old',form='formatted')
    read (boostUnit,*)
    do iWavenumber=1,size(self%wavenumberTable)
       read (boostUnit,*) self%wavenumberTable(iWavenumber),boostRow
       self%boostTable(iWavenumber,:)=boostRow
    end do
    close(boostUnit)
    deallocate(boostRow)
    call File_Unlock(self%fileLock)
    ! Convert wavenumbers from the h/Mpc units used by EuclidEmulator2 to the Mpc¯¹ used by Galacticus, and tabulate in terms of
    ! ln(k) and ln(1+z), interpolating in ln(B).
    do iRedshift=1,countRedshift
       self%redshiftTable(iRedshift)=log(1.0d0+euclidEmulator2Redshift(iRedshift))
    end do
    self%wavenumberTable=log(self%wavenumberTable*hubbleConstant)
    self%boostTable     =log(self%boostTable                    )
    ! Build the interpolator.
    if (allocated(self%interpolator_)) deallocate(self%interpolator_)
    allocate(self%interpolator_)
    self%interpolator_=interpolator2D(self%wavenumberTable,self%redshiftTable,self%boostTable)
    self%tabulated    =.true.
    return
  end subroutine euclidEmulator2Tabulate

  double precision function euclidEmulator2Redshift(iRedshift) result(redshift)
    !!{RST
    Return the redshift of the ``iRedshift``\ :sup:`th` point in the grid at which the nonlinear correction factor is tabulated.
    The grid is uniform in :math:`\ln(1+z)` over the range emulated by :term:`EuclidEmulator2`.
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
  end function euclidEmulator2Redshift

  double precision function euclidEmulator2Value(self,waveNumber,time) result(powerSpectrum)
    !!{RST
    Return a nonlinear power spectrum computed by applying the nonlinear correction factor emulated by the :term:`EuclidEmulator2`
    code of :cite:t:`euclid_collaboration_euclid_2021` to the linear theory power spectrum.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (powerSpectrumNonlinearEuclidEmulator2), intent(inout) :: self
    double precision                                       , intent(in   ) :: time                   , waveNumber
    double precision                                                       :: redshift               , boost              , &
         &                                                                    wavenumberLogarithmic  , redshiftLogarithmic, &
         &                                                                    boostLogarithmic       , gradient           , &
         &                                                                    boostLogarithmicMaximum

    ! Tabulate the nonlinear correction factor if necessary.
    if (.not.self%tabulated) then
       !$omp critical(powerSpectrumNonlinearEuclidEmulator2Tabulate)
       if (.not.self%tabulated) call self%tabulate()
       !$omp end critical(powerSpectrumNonlinearEuclidEmulator2Tabulate)
    end if
    ! Find the redshift corresponding to this time. Slightly-negative redshifts arising from round-off in the conversion from
    ! cosmic time at the present epoch are treated as zero.
    redshift=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    if (redshift < 0.0d0 .and. redshift > -toleranceRedshift) redshift=0.0d0
    if (redshift < 0.0d0 .or.  redshift >  redshiftMaximum  )                                           &
         & call Error_Report(                                                                           &
         &                   'redshift is outside the range emulated by EuclidEmulator2 (0 ≤ z ≤ 10)'// &
         &                    {introspection:location}                                                  &
         &                  )
    ! Clamp to the tabulated range - the end points of the tabulation are computed from the same expression as the redshifts
    ! passed to EuclidEmulator2, so can differ from the value computed here at the level of round-off.
    redshiftLogarithmic  =min(                                       &
         &                    max(                                   &
         &                        log(1.0d0+redshift)              , &
         &                        self%redshiftTable(            1)  &
         &                       )                                 , &
         &                        self%redshiftTable(countRedshift)  &
         &                   )
    wavenumberLogarithmic=log(wavenumber)
    ! Evaluate the nonlinear correction factor. Below the smallest emulated wavenumber the correction factor is held at its value
    ! at that wavenumber, while above the largest emulated wavenumber it is extrapolated using the logarithmic gradient at that
    ! wavenumber.
    if      (wavenumberLogarithmic <  self%wavenumberTable(                          1)) then
       boostLogarithmic       =self%interpolator_%interpolate(self%wavenumberTable(                          1),redshiftLogarithmic)
    else if (wavenumberLogarithmic >  self%wavenumberTable(size(self%wavenumberTable) )) then
       boostLogarithmicMaximum=self%interpolator_%interpolate(self%wavenumberTable(size(self%wavenumberTable) ),redshiftLogarithmic)
       gradient               =+(                                                                                                        &
            &                    +boostLogarithmicMaximum                                                                                &
            &                    -self%interpolator_%interpolate(self%wavenumberTable(size(self%wavenumberTable)-1),redshiftLogarithmic) &
            &                   )                                                                                                        &
            &                  /(                                                                                                        &
            &                    +self%wavenumberTable(size(self%wavenumberTable)  )                                                     &
            &                    -self%wavenumberTable(size(self%wavenumberTable)-1)                                                     &
            &                   )
       boostLogarithmic       =+boostLogarithmicMaximum                                                                                  &
            &                  +gradient                                                                                                 &
            &                  *(                                                                                                        &
            &                    +wavenumberLogarithmic                                                                                  &
            &                    -self%wavenumberTable(size(self%wavenumberTable))                                                       &
            &                   )
    else
       boostLogarithmic       =self%interpolator_%interpolate(wavenumberLogarithmic                                ,redshiftLogarithmic)
    end if
    boost        =exp(boostLogarithmic)
    powerSpectrum=+boost                                            &
         &        *self%powerSpectrum_%power(wavenumber,time)
    return
  end function euclidEmulator2Value
