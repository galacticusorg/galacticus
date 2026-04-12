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
  An intrinsic three-component thin-disk AGN SED implementation with strict
  energy-closure scaling.  The model comprises (1) a multitemperature blackbody
  accretion disk, (2) a soft Comptonized warm corona, and (3) a hard Comptonized hot
  corona.  The bolometric luminosity fractions allocated to the warm and hot coronae are
  set by \mono{[fractionWarm]} and \mono{[fractionHot]}, with the remainder going to
  the disk.  Strict energy closure is imposed: the integral of the SED over all
  frequencies equals the disk bolometric luminosity exactly.  The model is active only
  in the thin-disk regime, i.e.\ when $\dot{M}/\dot{M}_\mathrm{Edd} \geq$
  \mono{[mdotThinMinimum]}.
  !!}

  use :: Black_Hole_Accretion_Rates, only : blackHoleAccretionRateClass
  use :: Accretion_Disks           , only : accretionDisksClass
  use :: Numerical_Interpolation   , only : interpolator
  use :: Kind_Numbers              , only : kind_int8

  !![
  <accretionDiskSpectra name="accretionDiskSpectraThinDisk">
   <description>
    Accretion disk spectra computed intrinsically for a three-component thin-disk AGN SED
    comprising (1) a multitemperature blackbody disk following the Novikov--Thorne
    temperature profile, (2) a soft Comptonized warm corona, and (3) a hard Comptonized
    hot corona extending to hard X-rays.  Strict energy closure is enforced: fractions
    \mono{[fractionHot]} and \mono{[fractionWarm]} of the bolometric luminosity are
    allocated to the hot and warm coronae respectively, with the remainder $(1 -
    f_\mathrm{hot} - f_\mathrm{warm})$ going to the disk.  The disk peak temperature is
    derived self-consistently from energy conservation so that the SED integrates exactly
    to the bolometric luminosity.  The warm and hot corona spectra are power laws with
    exponential cutoffs, normalised by the upper incomplete gamma function to achieve the
    same strict closure.  The model returns zero outside the thin-disk regime, i.e.\ when
    $\dot{M}/\dot{M}_\mathrm{Edd} &lt;$ \mono{[mdotThinMinimum]}.
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
     double precision                                                            :: mdotThinMinimum                   , &
          &                                                                        fractionHot                       , &
          &                                                                        fractionWarm                      , &
          &                                                                        temperatureHot                    , &
          &                                                                        temperatureWarm                   , &
          &                                                                        spectralIndexHot                  , &
          &                                                                        spectralIndexWarm                 , &
          &                                                                        energyMinimumHot                  , &
          &                                                                        energyMinimumWarm                 , &
          &                                                                        massBlackHoleFiducial             , &
          &                                                                        spinBlackHoleFiducial
     ! Wavelength grid (logarithmic, set at construction time).
     integer                                                                     :: wavelengthCount
     double precision                             , allocatable, dimension(:  ) :: wavelengthTable
     type            (interpolator               )                              :: interpolatorWavelength
     ! Pre-computed corona normalization integrals (fixed for given parameters).
     double precision                                                            :: normHot                           , &
          &                                                                        normWarm
     ! Cache: SED and metadata for the most recently computed node.
     integer         (kind_int8                  )                              :: lastUniqueID
     logical                                                                    :: lastSEDComputed
     double precision                             , allocatable, dimension(:  ) :: lastSED
   contains
     !![
     <methods>
       <method description="Reset the cached SED." method="calculationReset" />
       <method description="Compute the SED for given BH properties and populate the cache." method="computeSED" />
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
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (accretionDiskSpectraThinDisk)                :: self
    type (inputParameters             ), intent(inout) :: parameters
    class(blackHoleAccretionRateClass ), pointer       :: blackHoleAccretionRate_
    class(accretionDisksClass         ), pointer       :: accretionDisks_
    double precision                                   :: mdotThinMinimum      , fractionHot          , &
         &                                                fractionWarm          , temperatureHot       , &
         &                                                temperatureWarm       , spectralIndexHot     , &
         &                                                spectralIndexWarm     , energyMinimumHot     , &
         &                                                energyMinimumWarm     , massBlackHoleFiducial, &
         &                                                spinBlackHoleFiducial

    !![
    <inputParameter>
      <name>mdotThinMinimum</name>
      <source>parameters</source>
      <defaultValue>0.01d0</defaultValue>
      <description>The minimum Eddington-scaled accretion rate, $\dot{M}/\dot{M}_\mathrm{Edd}$, below which the thin-disk SED returns zero.</description>
    </inputParameter>
    <inputParameter>
      <name>fractionHot</name>
      <source>parameters</source>
      <defaultValue>0.02d0</defaultValue>
      <description>The fraction of the bolometric luminosity allocated to the hot Comptonized corona.</description>
    </inputParameter>
    <inputParameter>
      <name>fractionWarm</name>
      <source>parameters</source>
      <defaultValue>0.05d0</defaultValue>
      <description>The fraction of the bolometric luminosity allocated to the warm Comptonized corona.</description>
    </inputParameter>
    <inputParameter>
      <name>temperatureHot</name>
      <source>parameters</source>
      <defaultValue>100.0d0</defaultValue>
      <description>The exponential cutoff energy (in keV) of the hot corona power-law spectrum.</description>
    </inputParameter>
    <inputParameter>
      <name>temperatureWarm</name>
      <source>parameters</source>
      <defaultValue>0.2d0</defaultValue>
      <description>The exponential cutoff energy (in keV) of the warm corona power-law spectrum.</description>
    </inputParameter>
    <inputParameter>
      <name>spectralIndexHot</name>
      <source>parameters</source>
      <defaultValue>1.7d0</defaultValue>
      <description>The photon spectral index $\Gamma$ of the hot corona power-law spectrum, so that $S_\nu \propto \nu^{-\Gamma}$.</description>
    </inputParameter>
    <inputParameter>
      <name>spectralIndexWarm</name>
      <source>parameters</source>
      <defaultValue>2.5d0</defaultValue>
      <description>The photon spectral index $\Gamma$ of the warm corona power-law spectrum.</description>
    </inputParameter>
    <inputParameter>
      <name>energyMinimumHot</name>
      <source>parameters</source>
      <defaultValue>0.1d0</defaultValue>
      <description>The minimum seed-photon energy (in keV) below which the hot corona spectrum is set to zero.</description>
    </inputParameter>
    <inputParameter>
      <name>energyMinimumWarm</name>
      <source>parameters</source>
      <defaultValue>1.0d-3</defaultValue>
      <description>The minimum seed-photon energy (in keV) below which the warm corona spectrum is set to zero.</description>
    </inputParameter>
    <inputParameter>
      <name>massBlackHoleFiducial</name>
      <source>parameters</source>
      <defaultValue>1.0d8</defaultValue>
      <description>The fiducial black hole mass (in $\mathrm{M}_\odot$) used when \refMethod{spectrumMassRate} is called without a specific node (i.e.\ when only the accretion rate and radiative efficiency are available).</description>
    </inputParameter>
    <inputParameter>
      <name>spinBlackHoleFiducial</name>
      <source>parameters</source>
      <defaultValue>0.5d0</defaultValue>
      <description>The fiducial dimensionless black hole spin used when \refMethod{spectrumMassRate} is called without a specific node.</description>
    </inputParameter>
    <objectBuilder class="blackHoleAccretionRate" name="blackHoleAccretionRate_" source="parameters"/>
    <objectBuilder class="accretionDisks"         name="accretionDisks_"         source="parameters"/>
    !!]
    self=accretionDiskSpectraThinDisk(                      &
         &                            mdotThinMinimum      , &
         &                            fractionHot          , &
         &                            fractionWarm         , &
         &                            temperatureHot       , &
         &                            temperatureWarm      , &
         &                            spectralIndexHot     , &
         &                            spectralIndexWarm    , &
         &                            energyMinimumHot     , &
         &                            energyMinimumWarm    , &
         &                            massBlackHoleFiducial, &
         &                            spinBlackHoleFiducial, &
         &                            blackHoleAccretionRate_, &
         &                            accretionDisks_         &
         &                           )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleAccretionRate_"/>
    <objectDestructor name="accretionDisks_"        />
    !!]
    return
  end function thinDiskConstructorParameters

  function thinDiskConstructorInternal(                      &
       &                               mdotThinMinimum      , &
       &                               fractionHot          , &
       &                               fractionWarm         , &
       &                               temperatureHot       , &
       &                               temperatureWarm      , &
       &                               spectralIndexHot     , &
       &                               spectralIndexWarm    , &
       &                               energyMinimumHot     , &
       &                               energyMinimumWarm    , &
       &                               massBlackHoleFiducial, &
       &                               spinBlackHoleFiducial, &
       &                               blackHoleAccretionRate_, &
       &                               accretionDisks_         &
       &                              ) result(self)
    !!{
    Internal constructor for the \refClass{accretionDiskSpectraThinDisk} accretion disk
    spectra class.
    !!}
    use :: Gamma_Functions             , only : Gamma_Function_Incomplete_Unnormalized
    use :: Numerical_Constants_Physical, only : plancksConstant, speedLight
    use :: Numerical_Constants_Units   , only : electronVolt   , metersToAngstroms
    use :: Numerical_Constants_Prefixes, only : kilo
    use :: Numerical_Ranges            , only : Make_Range                             , rangeTypeLogarithmic
    use :: Table_Labels                , only : extrapolationTypeZero
    implicit none
    type (accretionDiskSpectraThinDisk)                        :: self
    class(blackHoleAccretionRateClass ), target, intent(in   ) :: blackHoleAccretionRate_
    class(accretionDisksClass         ), target, intent(in   ) :: accretionDisks_
    double precision                           , intent(in   ) :: mdotThinMinimum      , fractionHot          , &
         &                                                        fractionWarm          , temperatureHot       , &
         &                                                        temperatureWarm       , spectralIndexHot     , &
         &                                                        spectralIndexWarm     , energyMinimumHot     , &
         &                                                        energyMinimumWarm     , massBlackHoleFiducial, &
         &                                                        spinBlackHoleFiducial
    double precision                                           :: energyHotJoules       , energyWarmJoules     , &
         &                                                        energyMinHotJoules    , energyMinWarmJoules  , &
         &                                                        uMinHot               , uMinWarm             , &
         &                                                        lambdaMinAngstroms
    integer                                                    :: wavelengthCountPerDecade
    !![
    <constructorAssign variables="mdotThinMinimum, fractionHot, fractionWarm, temperatureHot, temperatureWarm, spectralIndexHot, spectralIndexWarm, energyMinimumHot, energyMinimumWarm, massBlackHoleFiducial, spinBlackHoleFiducial, *blackHoleAccretionRate_, *accretionDisks_"/>
    !!]

    ! Initialize the per-node SED cache.
    self%lastUniqueID   = -1_kind_int8
    self%lastSEDComputed = .false.

    ! -----------------------------------------------------------------------
    ! Build the wavelength table.
    ! -----------------------------------------------------------------------
    ! Minimum wavelength: hc / (5 * kT_hot) – extends well into the hard X-ray
    ! regime, beyond the hot-corona exponential cutoff energy.
    ! Maximum wavelength: 1e7 Å (far-infrared, well past all disk emission).
    energyHotJoules      = temperatureHot * kilo * electronVolt
    lambdaMinAngstroms   = plancksConstant                                    &
         &               * speedLight * metersToAngstroms                     &
         &               / (5.0d0 * energyHotJoules)
    wavelengthCountPerDecade = 40
    self%wavelengthCount     = int(                                           &
         &   dble(wavelengthCountPerDecade)                                   &
         & * log10(1.0d7 / lambdaMinAngstroms)                               &
         & ) + 1
    allocate(self%wavelengthTable(self%wavelengthCount))
    allocate(self%lastSED        (self%wavelengthCount))
    self%wavelengthTable = Make_Range(                       &
         &                            lambdaMinAngstroms   , &
         &                            1.0d7                , &
         &                            self%wavelengthCount , &
         &                            rangeTypeLogarithmic   &
         &                           )
    self%interpolatorWavelength = interpolator(                      &
         &                                      self%wavelengthTable , &
         &                                      extrapolationType=extrapolationTypeZero &
         &                                     )

    ! -----------------------------------------------------------------------
    ! Pre-compute corona normalization integrals.
    ! -----------------------------------------------------------------------
    ! For S(ν) = A ν^{-Γ} exp(−h ν / E_c):
    !   norm = ∫_{ν_min}^∞ ν^{-Γ} exp(−h ν / E_c) dν
    !        = (E_c/h)^{1-Γ}  Γ_inc(1-Γ, h ν_min / E_c)
    ! where Γ_inc(a,x) = ∫_x^∞ t^{a-1} e^{-t} dt  (upper incomplete gamma).
    energyHotJoules    = temperatureHot  * kilo * electronVolt
    energyWarmJoules   = temperatureWarm * kilo * electronVolt
    energyMinHotJoules = energyMinimumHot  * kilo * electronVolt
    energyMinWarmJoules= energyMinimumWarm * kilo * electronVolt
    uMinHot  = energyMinHotJoules  / energyHotJoules
    uMinWarm = energyMinWarmJoules / energyWarmJoules
    self%normHot  = (energyHotJoules  / plancksConstant)**(1.0d0 - spectralIndexHot ) &
         &        * Gamma_Function_Incomplete_Unnormalized(1.0d0 - spectralIndexHot , uMinHot )
    self%normWarm = (energyWarmJoules / plancksConstant)**(1.0d0 - spectralIndexWarm) &
         &        * Gamma_Function_Incomplete_Unnormalized(1.0d0 - spectralIndexWarm, uMinWarm)

    ! Initialize the SED cache to zero.
    self%lastSED = 0.0d0
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

    call calculationResetEvent%attach(self,thinDiskCalculationReset,openMPThreadBindingAllLevels, &
         &                            label='accretionDiskSpectraThinDisk')
    return
  end subroutine thinDiskAutoHook

  subroutine thinDiskCalculationReset(self,node,uniqueID)
    !!{
    Reset the cached SED (triggered by the calculation-reset event).
    !!}
    use :: Kind_Numbers   , only : kind_int8
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    class  (accretionDiskSpectraThinDisk), intent(inout) :: self
    type   (treeNode                    ), intent(inout) :: node
    integer(kind_int8                   ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%lastUniqueID    = uniqueID
    self%lastSEDComputed = .false.
    self%lastSED         = 0.0d0
    return
  end subroutine thinDiskCalculationReset

  double precision function thinDiskSpectrumNode(self,node,wavelength) result(spectrum)
    !!{
    Return the accretion disk SED (in $L_\odot\,\mathrm{Hz}^{-1}$) at \mono{wavelength}
    (in \AA) for the black hole in \mono{node}, or zero if outside the thin-disk regime.
    !!}
    use, intrinsic :: ISO_C_Binding      , only : c_size_t
    use :: Black_Hole_Fundamentals      , only : Black_Hole_Eddington_Accretion_Rate , &
         &                                       Black_Hole_ISCO_Radius
    use :: Galacticus_Nodes             , only : nodeComponentBlackHole
    use :: Numerical_Constants_Physical , only : gravitationalConstant              , &
         &                                       speedLight
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class           (accretionDiskSpectraThinDisk), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: wavelength
    class           (nodeComponentBlackHole      ), pointer       :: blackHole
    double precision                                              :: rateAccretionSpheroid            , rateAccretionHotHalo             , &
         &                                                           rateAccretionNuclearStarCluster  , rateAccretion                    , &
         &                                                           efficiencyRadiative              , rateAccretionEddington           , &
         &                                                           massBlackHole                    , spinBlackHole                    , &
         &                                                           radiusISCO                       , luminosityBolometric
    integer         (c_size_t                    )                :: iWavelength
    double precision                             , dimension(0:1) :: hWavelength

    spectrum  =  0.0d0
    blackHole => node%blackHole()
    ! Get the total accretion rate onto this black hole.
    call self%blackHoleAccretionRate_%rateAccretion(blackHole,rateAccretionSpheroid,rateAccretionHotHalo,rateAccretionNuclearStarCluster)
    rateAccretion = +rateAccretionSpheroid           &
         &          +rateAccretionHotHalo            &
         &          +rateAccretionNuclearStarCluster
    ! Return zero for non-positive accretion rates.
    if (rateAccretion <= 0.0d0) return
    ! Check thin-disk regime.
    rateAccretionEddington = Black_Hole_Eddington_Accretion_Rate(blackHole)
    if (rateAccretionEddington <= 0.0d0) return
    if (rateAccretion / rateAccretionEddington < self%mdotThinMinimum) return
    ! Force a SED reset if the node has changed.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! (Re)compute the full SED for this node if needed.
    if (.not. self%lastSEDComputed) then
       massBlackHole       = blackHole%mass()
       spinBlackHole       = blackHole%spin()
       efficiencyRadiative = self%accretionDisks_%efficiencyRadiative(blackHole,rateAccretion)
       ! ISCO radius in metres: r_isco = iscoFactor * G * M / c^2.
       radiusISCO          = Black_Hole_ISCO_Radius(spinBlackHole)             &
            &              * gravitationalConstant * massBlackHole * massSolar &
            &              / speedLight**2
       ! Bolometric luminosity in Galacticus-internal units (M_sun/Gyr * (m/s)^2);
       ! thinDiskComputeSED applies the massSolar/gigaYear conversion internally.
       luminosityBolometric = efficiencyRadiative * rateAccretion * speedLight**2
       call self%computeSED(luminosityBolometric,radiusISCO)
       self%lastSEDComputed = .true.
       self%lastUniqueID    = node%uniqueID()
    end if
    ! Interpolate the cached SED at the requested wavelength.
    call self%interpolatorWavelength%linearFactors(wavelength,iWavelength,hWavelength)
    spectrum = +self%lastSED(iWavelength  ) * hWavelength(0) &
         &     +self%lastSED(iWavelength+1) * hWavelength(1)
    spectrum = max(spectrum, 0.0d0)
    return
  end function thinDiskSpectrumNode

  double precision function thinDiskSpectrumMassRate(self,accretionRate,efficiencyRadiative,wavelength) result(spectrum)
    !!{
    Return the accretion disk SED (in $L_\odot\,\mathrm{Hz}^{-1}$) at \mono{wavelength}
    (in \AA) for the given \mono{accretionRate} ($\mathrm{M}_\odot\,\mathrm{Gyr}^{-1}$)
    and \mono{efficiencyRadiative}, using the fiducial black hole mass and spin.
    Returns zero if outside the thin-disk regime.
    !!}
    use, intrinsic :: ISO_C_Binding         , only : c_size_t
    use :: Black_Hole_Fundamentals         , only : Black_Hole_ISCO_Radius
    use :: Numerical_Constants_Astronomical, only : massSolar            , gigaYear
    use :: Numerical_Constants_Atomic      , only : massHydrogenAtom
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Physical    , only : gravitationalConstant, speedLight, &
         &                                          thomsonCrossSection
    implicit none
    class           (accretionDiskSpectraThinDisk), intent(inout) :: self
    double precision                              , intent(in   ) :: accretionRate      , efficiencyRadiative, &
         &                                                           wavelength
    double precision                                              :: radiusISCO         , luminosityBolometric , &
         &                                                           rateAccretionEddington
    integer         (c_size_t                    )                :: iWavelength
    double precision                             , dimension(0:1) :: hWavelength
    double precision                             , allocatable    , dimension(:) :: tmpSED

    spectrum = 0.0d0
    if (accretionRate       <= 0.0d0) return
    if (efficiencyRadiative <= 0.0d0) return
    ! Eddington accretion rate for the fiducial black hole mass (M_sun/Gyr).
    ! Uses the same formula as Black_Hole_Eddington_Accretion_Rate: no massSolar
    ! factor because blackHole%mass() in M_sun and the massSolar / massSolar cancel.
    rateAccretionEddington = 4.0d0 * Pi * gravitationalConstant             &
         &                 * self%massBlackHoleFiducial                      &
         &                 * massHydrogenAtom * gigaYear                     &
         &                 / thomsonCrossSection / speedLight
    if (rateAccretionEddington <= 0.0d0) return
    if (accretionRate / rateAccretionEddington < self%mdotThinMinimum) return
    ! ISCO radius in metres using the fiducial spin.
    radiusISCO = Black_Hole_ISCO_Radius(self%spinBlackHoleFiducial)               &
         &     * gravitationalConstant * self%massBlackHoleFiducial * massSolar   &
         &     / speedLight**2
    ! Bolometric luminosity (Galacticus-internal: M_sun/Gyr × (m/s)^2).
    luminosityBolometric = efficiencyRadiative * accretionRate * speedLight**2
    ! Compute SED into a temporary array.
    allocate(tmpSED(self%wavelengthCount))
    tmpSED = 0.0d0
    call self%computeSED(luminosityBolometric,radiusISCO,sedOut=tmpSED)
    ! Interpolate at the requested wavelength.
    call self%interpolatorWavelength%linearFactors(wavelength,iWavelength,hWavelength)
    spectrum = +tmpSED(iWavelength  ) * hWavelength(0) &
         &     +tmpSED(iWavelength+1) * hWavelength(1)
    spectrum = max(spectrum,0.0d0)
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

    wavelengthsCount = self%wavelengthCount
    allocate(wavelengths(wavelengthsCount))
    wavelengths = self%wavelengthTable
    return
  end subroutine thinDiskWavelengths

  subroutine thinDiskComputeSED(self,luminosityBolometric,radiusISCO,sedOut)
    !!{
    Compute the three-component thin-disk SED (in $L_\odot\,\mathrm{Hz}^{-1}$) and
    store the result.  If \mono{sedOut} is present the result is written there;
    otherwise it is stored in \mono{self\%lastSED}.

    Strict energy closure is enforced for all three components:
    \begin{itemize}
      \item Disk: multitemperature blackbody with $T_\mathrm{max}$ derived from the
            relation $L_\mathrm{disk} = 4\pi\sigma r_\mathrm{isco}^2 T_\mathrm{max}^4 / 3$,
            where $L_\mathrm{disk} = (1 - f_\mathrm{hot} - f_\mathrm{warm})\,L_\mathrm{bol}$.
            The spectral integral is performed by numerical quadrature over disk radii
            with the logarithmic change of variable $t = \ln(r/r_\mathrm{isco})$.
      \item Warm/hot corona: power law $\propto\nu^{-\Gamma}\mathrm{e}^{-h\nu/E_c}$
            normalised so the integral equals $f_\mathrm{warm/hot}\,L_\mathrm{bol}$,
            using the upper incomplete gamma function computed at construction time.
    \end{itemize}
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar              , gigaYear              , luminositySolar
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Physical    , only : plancksConstant        , speedLight            , &
         &                                          stefanBoltzmannConstant
    use :: Numerical_Constants_Units       , only : electronVolt           , metersToAngstroms
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Integration           , only : integrator
    use :: Thermodynamics_Radiation        , only : Blackbody_Emission     , radianceTypeFrequency
    implicit none
    class           (accretionDiskSpectraThinDisk)                               , intent(inout) :: self
    double precision                                                             , intent(in   ) :: luminosityBolometric , radiusISCO
    double precision                              , optional, dimension(self%wavelengthCount), intent(  out) :: sedOut
    double precision                              ,           dimension(self%wavelengthCount)                :: sed
    double precision                                                                            :: luminosityDisk       , temperatureMax , &
         &                                                                                         luminosityWarm       , luminosityHot  , &
         &                                                                                         normFactorDisk       , frequencyNu    , &
         &                                                                                         energyHotJoules      , energyWarmJoules, &
         &                                                                                         energyMinHotJoules   , energyMinWarmJoules, &
         &                                                                                         uArg                 , fractionDisk   , &
         &                                                                                         normFactorWarm       , normFactorHot  , &
         &                                                                                         radiusISCOSq
    integer                                                                                     :: i
    type            (integrator             )                                                   :: integrator_
    ! Thread-private scalars used inside the contained integrand function.
    double precision                                                                            :: wavelengthCurrent    , temperatureMaxLocal

    sed = 0.0d0

    ! Guard against degenerate geometry.
    if (radiusISCO <= 0.0d0) then
       if (present(sedOut)) then
          sedOut = sed
       else
          self%lastSED = sed
       end if
       return
    end if

    ! -----------------------------------------------------------------------
    ! Luminosity fractions (all converted to Watts).
    ! -----------------------------------------------------------------------
    fractionDisk   = 1.0d0 - self%fractionHot - self%fractionWarm
    luminosityDisk = fractionDisk          * luminosityBolometric * massSolar / gigaYear
    luminosityWarm = self%fractionWarm     * luminosityBolometric * massSolar / gigaYear
    luminosityHot  = self%fractionHot      * luminosityBolometric * massSolar / gigaYear

    ! -----------------------------------------------------------------------
    ! Disk: T_max from energy conservation.
    ! L_disk = 4 pi sigma r_isco^2 T_max^4 / 3
    ! => T_max = ( 3 L_disk / (4 pi sigma r_isco^2) )^{1/4}
    ! -----------------------------------------------------------------------
    radiusISCOSq = radiusISCO**2
    if (luminosityDisk > 0.0d0) then
       temperatureMax = ( 3.0d0 * luminosityDisk                      &
            &           / (4.0d0 * Pi * stefanBoltzmannConstant        &
            &              * radiusISCOSq)                             &
            &           )**0.25d0
    else
       temperatureMax = 0.0d0
    end if

    ! -----------------------------------------------------------------------
    ! Corona spectral normalization factors.
    ! S_corona(nu) = normFactor * nu^{-Gamma} * exp(-h nu / E_c)  [L_sun / Hz]
    ! normFactor = L_corona / (luminositySolar * normIntegral)
    ! -----------------------------------------------------------------------
    normFactorWarm = 0.0d0
    normFactorHot  = 0.0d0
    if (self%normWarm > 0.0d0) &
         & normFactorWarm = luminosityWarm / (luminositySolar * self%normWarm)
    if (self%normHot  > 0.0d0) &
         & normFactorHot  = luminosityHot  / (luminositySolar * self%normHot )
    ! Disk normalization: L_nu_disk = normFactorDisk * integral [L_sun / Hz].
    normFactorDisk = 4.0d0 * Pi**2 * radiusISCOSq / luminositySolar

    ! Corona cutoff energies (Joules).
    energyHotJoules     = self%temperatureHot    * kilo * electronVolt
    energyWarmJoules    = self%temperatureWarm   * kilo * electronVolt
    energyMinHotJoules  = self%energyMinimumHot  * kilo * electronVolt
    energyMinWarmJoules = self%energyMinimumWarm * kilo * electronVolt

    ! -----------------------------------------------------------------------
    ! Loop over the wavelength table.
    ! -----------------------------------------------------------------------
    temperatureMaxLocal = temperatureMax
    integrator_         = integrator(diskIntegrand, toleranceRelative=1.0d-3)

    do i = 1, self%wavelengthCount
       wavelengthCurrent = self%wavelengthTable(i)

       ! ---- (1) Disk: multitemperature blackbody ----
       if (temperatureMaxLocal > 0.0d0) then
          ! L_nu_disk = normFactorDisk * int_0^{ln(x_max)} B_nu(T(exp(t))) exp(2t) dt
          ! where x = r/r_isco, T(x) = T_max x^{-3/4} (1-x^{-1/2})^{1/4}.
          sed(i) = normFactorDisk * integrator_%integrate(0.0d0, log(1.0d6))
       end if

       ! ---- (2) Warm corona ----
       frequencyNu = speedLight * metersToAngstroms / wavelengthCurrent
       if (plancksConstant * frequencyNu >= energyMinWarmJoules) then
          uArg = plancksConstant * frequencyNu / energyWarmJoules
          if (uArg < 300.0d0) &
               & sed(i) = sed(i) + normFactorWarm &
               &        * frequencyNu**(-self%spectralIndexWarm) * exp(-uArg)
       end if

       ! ---- (3) Hot corona ----
       if (plancksConstant * frequencyNu >= energyMinHotJoules) then
          uArg = plancksConstant * frequencyNu / energyHotJoules
          if (uArg < 300.0d0) &
               & sed(i) = sed(i) + normFactorHot &
               &        * frequencyNu**(-self%spectralIndexHot) * exp(-uArg)
       end if

    end do

    ! -----------------------------------------------------------------------
    ! Store or return result.
    ! -----------------------------------------------------------------------
    if (present(sedOut)) then
       sedOut = sed
    else
       self%lastSED = sed
    end if
    return

  contains

    double precision function diskIntegrand(t)
      !!{
      Integrand for the disk radial quadrature with $t = \ln(r/r_\mathrm{isco})$.
      The Novikov--Thorne temperature profile is
      $T(x) = T_\mathrm{max}\,x^{-3/4}\,(1 - x^{-1/2})^{1/4}$ where $x = e^t$.
      The factor $x^2 = e^{2t}$ arises from the change of variables $\mathrm{d}r = r\,\mathrm{d}t$.
      !!}
      implicit none
      double precision, intent(in   ) :: t
      double precision                :: x, temperature

      x = exp(t)
      ! Temperature profile is zero at the ISCO (x = 1) and peaks near x = 49/36.
      if (x <= 1.0d0) then
         diskIntegrand = 0.0d0
         return
      end if
      temperature   = temperatureMaxLocal * x**(-0.75d0) * (1.0d0 - x**(-0.5d0))**0.25d0
      diskIntegrand = Blackbody_Emission(wavelengthCurrent,temperature,radianceTypeFrequency) * x**2
      return
    end function diskIntegrand

  end subroutine thinDiskComputeSED

