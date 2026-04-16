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

  !+    Contributions to this file made by: Andrew Benson, Copilot.
  
  !     Implementation planning was assisted by Copilot, which then made the initial implementation. All code checked, physics
  !     assumptions verified, documented, and cited by Andrew Benson.
  
  !!{
  An intrinsic three-component thin-disk AGN SED implementation with strict
  energy-closure scaling, based on the model of \cite{done_intrinsic_2012}.
  !!}

  use :: Black_Hole_Accretion_Rates, only : blackHoleAccretionRateClass
  use :: Accretion_Disks           , only : accretionDisksClass
  use :: Numerical_Interpolation   , only : interpolator
  use :: Kind_Numbers              , only : kind_int8

  !![
  <accretionDiskSpectra name="accretionDiskSpectraThinDisk">
   <description>
    Accretion disk spectra computed intrinsically for a three-component thin-disk AGN SED based
    on the model of \cite{done_intrinsic_2012} and comprising (1) a multitemperature blackbody
    disk following the \cite{shakura_black_1973} temperature profile, (2) a soft Comptonized warm corona,
    and (3) a hard Comptonized hot corona extending to hard X-rays.  Strict energy closure is
    enforced: fractions \mono{[fractionHot]} and \mono{[fractionWarm]} of the bolometric
    luminosity are allocated to the hot and warm coronae respectively, with the remainder $(1 -
    f_\mathrm{hot} - f_\mathrm{warm})$ going to the disk.

    The disk temperature normalization is derived self-consistently from energy conservation so
    that the SED integrates exactly to the bolometric luminosity. Specifically, the disk
    temperature profile follows that for a \cite{shakura_black_1973} disk (their equation~2.6):
    \begin{equation}
    T^4(x) = T_0^4 x^{-3} (1-x^{-1/2}),
    \end{equation}
    where $x=R/R_\mathrm{ISCO}$. Integrating the thermal emission over the disk from the ISCO to
    infinity and setting this equal to the disk luminosity gives
    \begin{equation}
    2 \times 2 \pi R_\mathrm{ISCO}^2 \sigma T_0^4 \int_1^\infty \mathrm{d}x x^{-2} (1-x^{-1/2}) = L_\mathrm{disk},
    \end{equation}
    where $\sigma$ is the Stefan-Boltzmann constant, and the first factor of 2 on the left
    arises from that fact that the disk has upper and lower surfaces. This results in
    \begin{equation}
    T_0^4 = \frac{3 L_\mathrm{disk}}{4 \pi \sigma R_\mathrm{ISCO}^2}.
    \end{equation}
    
    The warm and hot corona spectra are power laws with exponential cutoffs, normalized by the
    upper incomplete gamma function to achieve the same strict closure. This functional form is
    commonly used to model hot coronae \citep[e.g.][]{fabian_properties_2015}, and is used here
    as a phenomenological model for the warm bump also. The model uses only the thin-disk
    luminosity, and returns zero outside the thin-disk regime.
   </description>
  </accretionDiskSpectra>
  !!]
  type, extends(accretionDiskSpectraClass) :: accretionDiskSpectraThinDisk
     !!{
     Three-component thin-disk AGN SED with strict energy-closure scaling.
     !!}
     private
     class           (blackHoleAccretionRateClass), pointer                     :: blackHoleAccretionRate_ => null()
     class           (accretionDisksClass        ), pointer                     :: accretionDisks_         => null()
     ! Model parameters.
     double precision                                                           :: fractionHot                      , fractionWarm         , &
          &                                                                        temperatureHot                   , temperatureWarm      , &
          &                                                                        spectralIndexHot                 , spectralIndexWarm    , &
          &                                                                        energyMinimumHot                 , energyMinimumWarm    , &
          &                                                                        massBlackHoleFiducial            , spinBlackHoleFiducial
     ! Wavelength grid (logarithmic, set at construction time).
     integer                                                                    :: wavelengthCount
     double precision                             , allocatable, dimension(:  ) :: wavelengthTable
     type            (interpolator               )                              :: interpolatorWavelength
     ! Pre-computed corona normalization integrals (fixed for given parameters).
     double precision                                                           :: normalizationHot                 , normalizationWarm
     ! Cache: SED and metadata for the most recently computed node.
     integer         (kind_int8                  )                              :: lastUniqueID
     logical                                                                    :: lastSEDComputed
     double precision                             , allocatable, dimension(:  ) :: lastSED
   contains
     !![
     <methods>
       <method method="calculationReset" description="Reset the cached SED."                                          />
       <method method="computeSED"       description="Compute the SED for given BH properties and populate the cache."/>
     </methods>
     !!]
     final     ::                     thinDiskDestructor
     procedure :: autoHook         => thinDiskAutoHook
     procedure :: calculationReset => thinDiskCalculationReset
     procedure :: spectrumNode     => thinDiskSpectrumNode
     procedure :: spectrumMassRate => thinDiskSpectrumMassRate
     procedure :: wavelengths      => thinDiskWavelengths
     procedure :: computeSED       => thinDiskComputeSED
  end type accretionDiskSpectraThinDisk

  interface accretionDiskSpectraThinDisk
     !!{
     Constructors for the \refClass{accretionDiskSpectraThinDisk} accretion disk spectra class.
     !!}
     module procedure thinDiskConstructorParameters
     module procedure thinDiskConstructorInternal
  end interface accretionDiskSpectraThinDisk

contains

  function thinDiskConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{accretionDiskSpectraThinDisk} accretion disk spectra
    class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (accretionDiskSpectraThinDisk)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (blackHoleAccretionRateClass ), pointer       :: blackHoleAccretionRate_
    class           (accretionDisksClass         ), pointer       :: accretionDisks_
    double precision                                              :: fractionHot            , fractionWarm         , &
         &                                                           temperatureHot         , temperatureWarm      , &
         &                                                           spectralIndexHot       , spectralIndexWarm    , &
         &                                                           energyMinimumHot       , energyMinimumWarm    , &
         &                                                           massBlackHoleFiducial  , spinBlackHoleFiducial
    !![
    <inputParameter>
      <name>fractionHot</name>
      <source>parameters</source>
      <defaultValue>0.02d0</defaultValue>
      <defaultSource>chosen to produce a plausible AGN SED---in reality will depend on viewing angle, geometry, and energy partition</defaultSource>
      <description>The fraction of the bolometric luminosity allocated to the hot Comptonized corona.</description>
    </inputParameter>
    <inputParameter>
      <name>fractionWarm</name>
      <source>parameters</source>
      <defaultValue>0.05d0</defaultValue>
      <defaultSource>chosen to produce a plausible AGN SED---in reality will depend on viewing angle, geometry, and energy partition</defaultSource>
      <description>The fraction of the bolometric luminosity allocated to the warm Comptonized corona.</description>
    </inputParameter>
    <inputParameter>
      <name>temperatureHot</name>
      <source>parameters</source>
      <defaultValue>200.0d0</defaultValue>
      <defaultSource>\citep[][typical of the fits in their Table~1]{fabian_properties_2015}</defaultSource>
      <description>The exponential cutoff energy (in keV) of the hot corona power-law spectrum.</description>
    </inputParameter>
    <inputParameter>
      <name>temperatureWarm</name>
      <source>parameters</source>
      <defaultValue>0.2d0</defaultValue>
      <defaultSource>\citep[][their fits show a range of around 0.2--0.3~keV]{done_intrinsic_2012}</defaultSource>
      <description>The exponential cutoff energy (in keV) of the warm corona power-law spectrum.</description>
    </inputParameter>
    <inputParameter>
      <name>spectralIndexHot</name>
      <source>parameters</source>
      <defaultValue>1.9d0</defaultValue>
      <defaultSource>\citep[][typical value from their Table~2]{tortosa_nustar_2018}</defaultSource>
      <description>The photon spectral index $\Gamma$ of the hot corona power-law spectrum, so that $S_\nu \propto \nu^{-\Gamma}$.</description>
    </inputParameter>
    <inputParameter>
      <name>spectralIndexWarm</name>
      <source>parameters</source>
      <defaultValue>2.6d0</defaultValue>
      <defaultSource>\citep[][\S4.3.2; reported a range of 2.6--2.7]{petrucci_testing_2018}</defaultSource>
      <description>The photon spectral index $\Gamma$ of the warm corona power-law spectrum.</description>
    </inputParameter>
    <inputParameter>
      <name>energyMinimumHot</name>
      <source>parameters</source>
      <defaultValue>0.1d0</defaultValue>
      <defaultSource>chosen to restrict the hot component to the X-ray regime where it is typically measured</defaultSource>
      <description>The minimum seed-photon energy (in keV) below which the hot corona spectrum is set to zero.</description>
    </inputParameter>
    <inputParameter>
      <name>energyMinimumWarm</name>
      <source>parameters</source>
      <defaultValue>1.0d-2</defaultValue>
      <defaultSource>typical of the UV/EUV disk photons which seed the warm Comptonization, e.g. \cite{done_intrinsic_2012}</defaultSource>
      <description>The minimum seed-photon energy (in keV) below which the warm corona spectrum is set to zero.</description>
    </inputParameter>
    <inputParameter>
      <name>massBlackHoleFiducial</name>
      <source>parameters</source>
      <defaultValue>1.0d8</defaultValue>
      <description>The fiducial black hole mass (in $\mathrm{M}_\odot$) used when \mono{spectrumMassRate} is called without a specific node (i.e.\ when only the accretion rate and radiative efficiency are available).</description>
    </inputParameter>
    <inputParameter>
      <name>spinBlackHoleFiducial</name>
      <source>parameters</source>
      <defaultValue>0.5d0</defaultValue>
      <description>The fiducial dimensionless black hole spin used when \mono{spectrumMassRate} is called without a specific node.</description>
    </inputParameter>
    <objectBuilder class="blackHoleAccretionRate" name="blackHoleAccretionRate_" source="parameters"/>
    <objectBuilder class="accretionDisks"         name="accretionDisks_"         source="parameters"/>
    !!]
    self=accretionDiskSpectraThinDisk(fractionHot,fractionWarm,temperatureHot,temperatureWarm,spectralIndexHot,spectralIndexWarm,energyMinimumHot,energyMinimumWarm,massBlackHoleFiducial,spinBlackHoleFiducial,blackHoleAccretionRate_,accretionDisks_) 
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleAccretionRate_"/>
    <objectDestructor name="accretionDisks_"        />
    !!]
    return
  end function thinDiskConstructorParameters

  function thinDiskConstructorInternal(fractionHot,fractionWarm,temperatureHot,temperatureWarm,spectralIndexHot,spectralIndexWarm,energyMinimumHot,energyMinimumWarm,massBlackHoleFiducial,spinBlackHoleFiducial,blackHoleAccretionRate_,accretionDisks_) result(self)
    !!{
    Internal constructor for the \refClass{accretionDiskSpectraThinDisk} accretion disk
    spectra class.
    !!}
    use :: Gamma_Functions             , only : Gamma_Function_Incomplete_Unnormalized
    use :: Numerical_Constants_Physical, only : plancksConstant                       , speedLight
    use :: Numerical_Constants_Units   , only : electronVolt                          , metersToAngstroms
    use :: Numerical_Constants_Prefixes, only : kilo
    use :: Numerical_Ranges            , only : Make_Range                            , rangeTypeLogarithmic
    use :: Table_Labels                , only : extrapolationTypeZero
    implicit none
    type            (accretionDiskSpectraThinDisk)                        :: self
    class           (blackHoleAccretionRateClass ), intent(in   ), target :: blackHoleAccretionRate_
    class           (accretionDisksClass         ), intent(in   ), target :: accretionDisks_
    double precision                              , intent(in   )         :: fractionHot                     , fractionWarm           , &
         &                                                                   temperatureHot                  , temperatureWarm        , &
         &                                                                   spectralIndexHot                , spectralIndexWarm      , &
         &                                                                   energyMinimumHot                , energyMinimumWarm      , &
         &                                                                   massBlackHoleFiducial           , spinBlackHoleFiducial
    integer                                       , parameter             :: wavelengthCountPerDecade  =100
    double precision                              , parameter             :: wavelengthMaximumAngstroms=1.0d7
    double precision                                                      :: energyHotJoules                 , energyWarmJoules       , &
         &                                                                   energyMinimumHotJoules          , energyMinimumWarmJoules, &
         &                                                                   energyRatioHot                  , energyRatioWarm        , &
         &                                                                   wavelengthMinimumAngstroms
    !![
    <constructorAssign variables="fractionHot, fractionWarm, temperatureHot, temperatureWarm, spectralIndexHot, spectralIndexWarm, energyMinimumHot, energyMinimumWarm, massBlackHoleFiducial, spinBlackHoleFiducial, *blackHoleAccretionRate_, *accretionDisks_"/>
    !!]

    ! Initialize the per-node SED cache.
    self%lastUniqueID   =-1_kind_int8
    self%lastSEDComputed=.false.
    ! Build the wavelength table.
    !! Minimum wavelength: hc/5kT_hot – extends well into the hard X-ray regime, beyond the hot-corona exponential cutoff energy.
    !! Maximum wavelength: 10⁷ Å (far-infrared, well past all disk emission).
    energyHotJoules           =+temperatureHot    &
         &                     *kilo              &
         &                     *electronVolt
    wavelengthMinimumAngstroms=+plancksConstant   &
         &                     *speedLight        &
         &                     *metersToAngstroms &
         &                     /5.0d0             &
         &                     /energyHotJoules
    self%wavelengthCount     = int(dble(wavelengthCountPerDecade)*log10(wavelengthMaximumAngstroms/wavelengthMinimumAngstroms))+1
    allocate(self%wavelengthTable(self%wavelengthCount))
    allocate(self%lastSED        (self%wavelengthCount))
    self%wavelengthTable       =Make_Range  (                                                   &
         &                                                          wavelengthMinimumAngstroms, &
         &                                                          wavelengthMaximumAngstroms, &
         &                                                     self%wavelengthCount           , &
         &                                                          rangeTypeLogarithmic        &
         &                                  )
    self%interpolatorWavelength=interpolator(                                                   &
         &                                                     self%wavelengthTable           , &
         &                                   extrapolationType=     extrapolationTypeZero       &
         &                                  )
    ! Pre-compute corona normalization integrals.
    !! For S(ν) = A ν^{-Γ} exp(−h ν / E_c):
    !!   norm = ∫_{ν_min}^∞ ν^{-Γ} exp(−h ν / E_c) dν
    !!        = (E_c/h)^{1-Γ}  Γ_inc(1-Γ, h ν_min / E_c)
    !! where Γ_inc(a,x) = ∫_x^∞ t^{a-1} e^{-t} dt (upper incomplete gamma function).
    energyHotJoules        =+temperatureHot   *kilo*electronVolt
    energyWarmJoules       =+temperatureWarm  *kilo*electronVolt
    energyMinimumHotJoules =+energyMinimumHot *kilo*electronVolt
    energyMinimumWarmJoules=+energyMinimumWarm*kilo*electronVolt
    energyRatioHot         =+energyMinimumHotJoules /energyHotJoules
    energyRatioWarm        =+energyMinimumWarmJoules/energyWarmJoules
    self%normalizationHot  =+(energyHotJoules /plancksConstant)**  (1.0d0-spectralIndexHot                 ) &
         &                  *Gamma_Function_Incomplete_Unnormalized(1.0d0-spectralIndexHot ,energyRatioHot )
    self%normalizationWarm =+(energyWarmJoules/plancksConstant)**  (1.0d0-spectralIndexWarm                ) &
         &                  *Gamma_Function_Incomplete_Unnormalized(1.0d0-spectralIndexWarm,energyRatioWarm)
    ! Initialize the SED cache to zero.
    self%lastSED           =-huge(0.0d0)
    return
  end function thinDiskConstructorInternal

  subroutine thinDiskDestructor(self)
    !!{
    Destructor for the \refClass{accretionDiskSpectraThinDisk} accretion disk spectra
    class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(accretionDiskSpectraThinDisk), intent(inout) :: self

    !![
    <objectDestructor name="self%blackHoleAccretionRate_"/>
    <objectDestructor name="self%accretionDisks_"        />
    !!]
    if (calculationResetEvent%isAttached(self,thinDiskCalculationReset)) &
         & call calculationResetEvent%detach(self,thinDiskCalculationReset)
    return
  end subroutine thinDiskDestructor

  subroutine thinDiskAutoHook(self)
    !!{
    Attach to the calculation reset event so that the SED cache is cleared when a new
    node is processed.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(accretionDiskSpectraThinDisk), intent(inout) :: self

    call calculationResetEvent%attach(self,thinDiskCalculationReset,openMPThreadBindingAllLevels,label='accretionDiskSpectraThinDisk')
    return
  end subroutine thinDiskAutoHook

  subroutine thinDiskCalculationReset(self,node,uniqueID)
    !!{
    Reset the cached SED (triggered by the calculation-reset event).
    !!}
    use :: Kind_Numbers    , only : kind_int8
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    class  (accretionDiskSpectraThinDisk), intent(inout) :: self
    type   (treeNode                    ), intent(inout) :: node
    integer(kind_int8                   ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%lastUniqueID   =uniqueID
    self%lastSEDComputed=.false.
    self%lastSED        =-huge(0.0d0)
    return
  end subroutine thinDiskCalculationReset

  double precision function thinDiskSpectrumNode(self,node,wavelength) result(spectrum)
    !!{
    Return the accretion disk SED (in $L_\odot\,\mathrm{Hz}^{-1}$) at \mono{wavelength}
    (in \AA) for the black hole in \mono{node}, or zero if outside the thin-disk regime.
    !!}
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Accretion_Disks                 , only : accretionDiskTypeThin
    use            :: Black_Hole_Fundamentals         , only : Black_Hole_Eddington_Accretion_Rate, Black_Hole_ISCO_Radius
    use            :: Galacticus_Nodes                , only : nodeComponentBlackHole
    use            :: Numerical_Constants_Physical    , only : gravitationalConstant              , speedLight
    use            :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class           (accretionDiskSpectraThinDisk), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: wavelength
    class           (nodeComponentBlackHole      ), pointer       :: blackHole
    double precision                                              :: rateAccretionSpheroid          , rateAccretionHotHalo , &
         &                                                           rateAccretionNuclearStarCluster, rateAccretion        , &
         &                                                           efficiencyRadiative            , luminosityBolometric , &
         &                                                           massBlackHole                  , spinBlackHole        , &
         &                                                           radiusISCO
    integer         (c_size_t                    )                :: iWavelength
    double precision                             , dimension(0:1) :: hWavelength

    spectrum  =  0.0d0
    blackHole => node%blackHole()
    ! Get the total accretion rate onto this black hole.
    call self%blackHoleAccretionRate_%rateAccretion(blackHole,rateAccretionSpheroid,rateAccretionHotHalo,rateAccretionNuclearStarCluster)
    rateAccretion=+rateAccretionSpheroid           &
         &         +rateAccretionHotHalo            &
         &         +rateAccretionNuclearStarCluster
    ! Return zero for non-positive accretion rates.
    if (rateAccretion <= 0.0d0) return
    ! Force a SED reset if the node has changed.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! (Re)compute the full SED for this node if needed.
    if (.not.self%lastSEDComputed) then
       massBlackHole      =blackHole                %mass               (                                             )
       spinBlackHole      =blackHole                %spin               (                                             )
       efficiencyRadiative=self     %accretionDisks_%efficiencyRadiative(blackHole,rateAccretion,accretionDiskTypeThin)
       ! ISCO radius in meters: r_isco = x_isco GM•/c², where x_isco is the ISCO radius in gravitational units.
       radiusISCO         =+Black_Hole_ISCO_Radius(spinBlackHole) &
            &              *gravitationalConstant                 &
            &              *massBlackHole                         &
            &              *massSolar                             &
            &              /speedLight**2
       ! Bolometric luminosity in Galacticus-internal units (M☉/Gyr (m/s)²);
       ! the function thinDiskComputeSED() applies the massSolar/gigaYear conversion internally.
       luminosityBolometric =+efficiencyRadiative    &
            &                *rateAccretion          &
            &                *speedLight         **2
       call self%computeSED(luminosityBolometric,radiusISCO)
       self%lastSEDComputed=.true.
       self%lastUniqueID   =node%uniqueID()
    end if
    ! Interpolate the cached SED at the requested wavelength.
    call self%interpolatorWavelength%linearFactors(wavelength,iWavelength,hWavelength)
    spectrum=+self%lastSED(iWavelength  )*hWavelength(0) &
         &   +self%lastSED(iWavelength+1)*hWavelength(1)
    spectrum=max(spectrum,0.0d0)
    return
  end function thinDiskSpectrumNode

  double precision function thinDiskSpectrumMassRate(self,accretionRate,efficiencyRadiative,wavelength) result(spectrum)
    !!{
    Return the accretion disk SED (in $L_\odot\,\mathrm{Hz}^{-1}$) at \mono{wavelength}
    (in \AA) for the given \mono{accretionRate} ($\mathrm{M}_\odot\,\mathrm{Gyr}^{-1}$)
    and \mono{efficiencyRadiative}, using the fiducial black hole mass and spin.
    Returns zero if outside the thin-disk regime.
    !!}
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Black_Hole_Fundamentals         , only : Black_Hole_ISCO_Radius
    use            :: Numerical_Constants_Astronomical, only : massSolar             , gigaYear
    use            :: Numerical_Constants_Atomic      , only : massHydrogenAtom
    use            :: Numerical_Constants_Math        , only : Pi
    use            :: Numerical_Constants_Physical    , only : gravitationalConstant , speedLight, &
         &                                                     thomsonCrossSection
    implicit none
    class           (accretionDiskSpectraThinDisk), intent(inout)              :: self
    double precision                              , intent(in   )              :: accretionRate, efficiencyRadiative , &
         &                                                                        wavelength
    double precision                                                           :: radiusISCO   , luminosityBolometric
    integer         (c_size_t                    )                             :: iWavelength
    double precision                             , dimension(0:1)              :: hWavelength
    double precision                             , dimension( : ), allocatable :: tmpSED

    spectrum=0.0d0
    if (accretionRate       <= 0.0d0) return
    if (efficiencyRadiative <= 0.0d0) return
    ! ISCO radius in meters using the fiducial spin.
    radiusISCO=+Black_Hole_ISCO_Radius(self%spinBlackHoleFiducial) &
         &     *gravitationalConstant                              &
         &     *self%massBlackHoleFiducial                         &
         &     *massSolar                                          &
         &     /speedLight**2
    ! Bolometric luminosity (in Galacticus-internal units: M☉/Gyr × (m/s)²).
    luminosityBolometric=+efficiencyRadiative &
         &               *accretionRate       &
         &               *speedLight**2
    ! Compute SED into a temporary array.
    allocate(tmpSED(self%wavelengthCount))
    tmpSED=0.0d0
    call self%computeSED(luminosityBolometric,radiusISCO,sedOut=tmpSED)
    ! Interpolate at the requested wavelength.
    call self%interpolatorWavelength%linearFactors(wavelength,iWavelength,hWavelength)
    spectrum=+tmpSED(iWavelength  )*hWavelength(0) &
         &   +tmpSED(iWavelength+1)*hWavelength(1)
    spectrum=max(spectrum,0.0d0)
    deallocate(tmpSED)
    return
  end function thinDiskSpectrumMassRate

  subroutine thinDiskWavelengths(self,wavelengthsCount,wavelengths)
    !!{
    Return the wavelength grid at which the thin-disk SED is tabulated.
    !!}
    implicit none
    class           (accretionDiskSpectraThinDisk)                           , intent(inout) :: self
    integer                                                                  , intent(  out) :: wavelengthsCount
    double precision                              , allocatable, dimension(:), intent(  out) :: wavelengths

    wavelengthsCount=self%wavelengthCount
    allocate(wavelengths(wavelengthsCount))
    wavelengths=self%wavelengthTable
    return
  end subroutine thinDiskWavelengths

  subroutine thinDiskComputeSED(self,luminosityBolometric,radiusISCO,sedOut)
    !!{
    Compute the three-component thin-disk SED (in $L_\odot\,\mathrm{Hz}^{-1}$) and
    store the result. If \mono{sedOut} is present the result is written there;
    otherwise it is stored in \mono{self\%lastSED}.

    Strict energy closure is enforced for all three components:
    \begin{itemize}
      \item Disk: multitemperature blackbody with $T_0$ derived from the
            relation $L_\mathrm{disk} = 4\pi\sigma r_\mathrm{isco}^2 T_0^4 / 3$,
            where $L_\mathrm{disk} = (1 - f_\mathrm{hot} - f_\mathrm{warm})\,L_\mathrm{bol}$.
            The spectral integral is performed numerically over disk radii
            with the logarithmic change of variable $t = \ln(r/r_\mathrm{isco})$.
      \item Warm/hot corona: power law $\propto\nu^{-\Gamma}\mathrm{e}^{-h\nu/E_c}$
            normalized so the integral equals $f_\mathrm{warm/hot}\,L_\mathrm{bol}$,
            using the upper incomplete gamma function computed at construction time.
    \end{itemize}
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar         , gigaYear             , luminositySolar
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Physical    , only : plancksConstant   , speedLight           , stefanBoltzmannConstant
    use :: Numerical_Constants_Units       , only : electronVolt      , metersToAngstroms
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Integration           , only : integrator
    use :: Thermodynamics_Radiation        , only : Blackbody_Emission, radianceTypeFrequency
    implicit none
    class           (accretionDiskSpectraThinDisk)                                           , intent(inout) :: self
    double precision                                                                         , intent(in   ) :: luminosityBolometric        , radiusISCO
    double precision                              , optional, dimension(self%wavelengthCount), intent(  out) :: sedOut
    double precision                              ,           dimension(self%wavelengthCount)                :: sed
    double precision                              , parameter                                                :: radiusMaximum         =1.0d6 ! Maximum radius of the disk in r_ISCO units.
    double precision                              , parameter                                                :: argumentMaximum       =3.0d2 ! Maximum argument for exponential terms.
    double precision                                                                                         :: luminosityDisk              , temperatureNormalization, &
         &                                                                                                      luminosityWarm              , luminosityHot           , &
         &                                                                                                      normalizationDisk           , frequency               , &
         &                                                                                                      energyHotJoules             , energyWarmJoules        , &
         &                                                                                                      energyMinimumHotJoules      , energyMinimumWarmJoules , &
         &                                                                                                      uArgument                   , fractionDisk            , &
         &                                                                                                      normalizationWarm           , normalizationHot        , &
         &                                                                                                      radiusISCO2                 , wavelengthCurrent
    integer                                                                                                  :: i
    type            (integrator                  )                                                           :: integrator_

    sed=0.0d0
    ! Guard against degenerate geometry.
    if (radiusISCO <= 0.0d0) then
       if (present(sedOut)) then
          sedOut      =sed
       else
          self%lastSED=sed
       end if
       return
    end if
    ! Luminosity fractions (all converted to Watts).
    fractionDisk  =+1.0d0             &
         &         -self%fractionHot  &
         &         -self%fractionWarm
    luminosityDisk=     fractionDisk*luminosityBolometric*massSolar/gigaYear
    luminosityWarm=self%fractionWarm*luminosityBolometric*massSolar/gigaYear
    luminosityHot =self%fractionHot *luminosityBolometric*massSolar/gigaYear
    ! Disk: T₀ from energy conservation assuming a Shakura-Sunyaev disk (see class description):
    !       T₀ = ( 3 L_disk / (4π σ r_isco²) )^{1/4}
    radiusISCO2=radiusISCO**2
    if (luminosityDisk > 0.0d0) then
       temperatureNormalization=+(                         &
            &                     +3.0d0                   &
            &                     *luminosityDisk          &
            &                     /+4.0d0                  &
            &                     /Pi                      &
            &                     /stefanBoltzmannConstant &
            &                     /radiusISCO2             &
            &                    )**(1.0d0/4.0d0)
    else
       temperatureNormalization=+0.0d0
    end if
    ! Corona spectral normalization factors.
    ! S_corona(ν) = normalization * ν^{-Γ} * exp(-hν / E_c)  [L☉ / Hz]
    ! normalization = L_corona / (luminositySolar * normIntegral)
    normalizationWarm=+0.0d0
    normalizationHot =+0.0d0
    if (self%normalizationWarm > 0.0d0)                                            &
         & normalizationWarm=luminosityWarm/luminositySolar/self%normalizationWarm
    if (self%normalizationHot  > 0.0d0)                                            &
         & normalizationHot =luminosityHot /luminositySolar/self%normalizationHot
    ! Disk normalization: L_ν_disk = normalizationDisk * integral [L☉ / Hz]. The factor of 4π² arises from a geometric factor of
    ! 2π (surface area element of the disk), a factor of 2 accounting for the upper and lower surfaces of the disk, and a factor
    ! of π representing the total solid angle into which each surface is radiating.
    normalizationDisk=+4.0d0           &
         &            *Pi**2           &
         &            *radiusISCO2     &
         &            /luminositySolar
    ! Corona cutoff energies (Joules).
    energyHotJoules        =self%temperatureHot   *kilo*electronVolt
    energyWarmJoules       =self%temperatureWarm  *kilo*electronVolt
    energyMinimumHotJoules =self%energyMinimumHot *kilo*electronVolt
    energyMinimumWarmJoules=self%energyMinimumWarm*kilo*electronVolt
    ! Iterate over the wavelength table.
    integrator_=integrator(diskIntegrand,toleranceRelative=1.0d-3)
    do i=1,self%wavelengthCount
       wavelengthCurrent=self%wavelengthTable(i)
       ! (1) Disk: multitemperature blackbody
       if (temperatureNormalization > 0.0d0) then
          ! L_ν_disk = normalizationDisk * int_0^{ln(x_max)} B_ν(T(exp(t))) exp(2t) dt
          ! where x = r/r_isco, T(x) = T₀ x^{-3/4} (1-x^{-1/2})^{1/4}.
          sed(i)=normalizationDisk*integrator_%integrate(0.0d0,log(radiusMaximum))
       end if
       ! (2) Warm corona:
       frequency=+speedLight        &
            &    *metersToAngstroms &
            &    /wavelengthCurrent
       if (plancksConstant*frequency >= energyMinimumWarmJoules) then
          uArgument=plancksConstant*frequency/energyWarmJoules
          if (uArgument < argumentMaximum)                    &
               & sed(i)=+sed(i)                               &
               &        +normalizationWarm                    &
               &        *frequency**(-self%spectralIndexWarm) &
               &        *exp(-uArgument)
       end if
       ! (3) Hot corona:
       if (plancksConstant*frequency >= energyMinimumHotJoules) then
          uArgument=plancksConstant*frequency/energyHotJoules
          if (uArgument < argumentMaximum)                   &
               & sed(i)=+sed(i)                              &
               &        +normalizationHot                    &
               &        *frequency**(-self%spectralIndexHot) &
               &        *exp(-uArgument)
       end if
    end do
    ! Store or return result.
    if (present(sedOut)) then
       sedOut      =sed
    else
       self%lastSED=sed
    end if
    return

  contains

    double precision function diskIntegrand(t)
      !!{
      Integrand for the disk radial quadrature with $t = \ln(r/r_\mathrm{isco})$.
      The Shakura--Sunyaev temperature profile is
      $T(x) = T_\mathrm{max}\,x^{-3/4}\,(1 - x^{-1/2})^{1/4}$ where $x = e^t$ (Shakura \& Sunyaev; 1983; A\&A; 24; 337; eq. 3.5).
      The factor $x^2 = e^{2t}$ arises from the change of variables $\mathrm{d}r = r\,\mathrm{d}t$.
      !!}
      implicit none
      double precision, intent(in   ) :: t
      double precision                :: x, temperature

      x=exp(t)
      ! Temperature profile is zero at the ISCO (x = 1) and peaks near x = 49/36.
      if (x <= 1.0d0) then
         diskIntegrand = 0.0d0
         return
      end if
      temperature  =+temperatureNormalization     &
           &        *       x**(-0.75d0)          &
           &        *(1.0d0-x**(-0.50d0))**0.25d0
      diskIntegrand=+Blackbody_Emission(wavelengthCurrent,temperature,radianceTypeFrequency)&
           &        *x**2
      return
    end function diskIntegrand

  end subroutine thinDiskComputeSED

