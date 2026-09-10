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
  Implementation of the :cite:t:`blitz_role_2006` star formation rate surface density law for galactic disks.
  !!}

  use :: Kind_Numbers       , only : kind_int8
  use :: Root_Finder        , only : rootFinder
  use :: Math_Exponentiation, only : fastExponentiator
  use :: Numerical_Ranges   , only : rangeLattice

  ! Floor on the disk gas mass, below which the disk is treated as gas-free.
  !
  ! The stellar pressure-boost coefficient scales as factorBoostStellarCoefficient proportional to
  ! 1/massGas, and is subsequently raised to the fourth power when locating the molecular/atomic
  ! transition radius. A vanishingly small but strictly positive gas mass therefore drives the
  ! coefficient (or its fourth power) past HUGE(0.0d0), raising an overflow floating point exception.
  ! Writing the finite remainder of the coefficient as C = σ_gas*2*π*R²*√(M_star/[2*π²*G*
  ! h_star*R³]) (C <~ 1e10 for any physical disk), the coefficient itself overflows for
  ! massGas <~ C/HUGE ~ 1e-298, and its fourth power for the far less extreme
  ! massGas <~ C/HUGE**(1/4) ~ 1e-68. A gas mass of ~1e-30 Msun is numerical residue rather than gas -
  ! its star formation rate is already zero - so flooring here removes the exception with tens of orders
  ! of margin above the binding (fourth-power) bound while remaining ~30 orders below any physical disk
  ! gas mass. Note: this floor must be applied identically in blitz2006ComputeFactors, blitz2006Rate,
  ! and blitz2006Intervals, since the coefficients are only *set* (in computeFactors) under
  ! massGas > massGasFloor and must only be *used* (in Rate/Intervals) under the complementary condition.
  double precision, parameter :: massGasFloor=1.0d-30

  ! Densities of tabulation points along the three axes of the tabulated solution, and the intervals - in lattice steps - to
  ! which their bounds are pinned. The anchoring is to *tenths* of a decade rather than to whole decades because this tabulation
  ! is three-dimensional and its cost is the product of its three extents: rounding each axis outward to a whole decade would,
  ! for the ranges these axes actually reach, inflate the number of points to be evaluated - each of which is a numerical
  ! integral - by about a half, whereas tenths of a decade cost a few percent. Note that the anchor interval must divide the
  ! density of points, so that whole numbers of the lattice coordinate are themselves anchor points and a bound lying on one -
  ! such as the hard lower limits below - can actually be reached.
  integer         , parameter :: pointsPerDecadeFactorBoost       =30, pointsPerDecadeFactorBoostStellar=30, &
       &                         pointsPerDecadeRadius            =30
  integer         , parameter :: anchorEveryFactorBoost           = 3, anchorEveryFactorBoostStellar    = 3, &
       &                         anchorEveryRadius                = 3

  ! Generate a source digest. The tabulated integral depends on only two of this class's parameters, so its file name is built
  ! from those alone rather than from a full hashed descriptor - see the constructor. But it depends also on the code which
  ! generates it: the integrand, the integration tolerance and rule, the floors applied to the two boost coefficients, and the
  ! densities of tabulation points above. None of those are visible to a descriptor of the parameters, so the digest is folded
  ! into the name to keep a tabulation built by earlier code from being silently reused.
  !![
  <sourceDigest name="blitz2006SourceDigest"/>
  !!]

  !![
  <starFormationRateSurfaceDensityDisks name="starFormationRateSurfaceDensityDisksBlitz2006" docformat="rst">
   <description>
   A star formation rate surface density class which assumes that the star formation rate is given by :cite:p:`blitz_role_2006`:

   .. math::

      \dot{\Sigma}_\star(R) = \nu_\mathrm{SF}(R) \Sigma_\mathrm{H_2, disk}(R),

   where :math:`\nu_\mathrm{SF}` is a frequency given by

   .. math::

      \nu_\mathrm{SF}(R) = \nu_\mathrm{SF,0} \left[ 1 + \left({\Sigma_\mathrm{HI}\over \Sigma_0}\right)^q \right],

   where :math:`q=`\ ``[surfaceDensityExponent]`` and :math:`\Sigma_0=`\ ``[surfaceDensityCritical]`` are parameters, the surface density of molecular gas :math:`\Sigma_\mathrm{H_2} = (P_\mathrm{ext}/P_0)^\alpha \Sigma_\mathrm{HI}`, where :math:`\alpha=`\ ``[pressureExponent]`` and :math:`P_0=`\ ``[pressureCharacteristic]`` are parameters, and the hydrostatic pressure in the disk plane assuming locally isothermal gas and stellar components is given by

   .. math::

      P_\mathrm{ext} \approx {\pi\over 2} \G \Sigma_\mathrm{gas} \left[ \Sigma_\mathrm{gas} + \left({\sigma_\mathrm{gas}\over
      \sigma_\star}\right)\Sigma_\star\right]

   where we assume that the velocity dispersion in the gas is fixed at :math:`\sigma_\mathrm{gas}=`\ ``[velocityDispersionDiskGas]`` and, assuming :math:`\Sigma_\star \gg \Sigma_\mathrm{gas}`, we can write the stellar velocity dispersion in terms of the disk scale height, :math:`h_\star`, as

   .. math::

      \sigma_\star = \sqrt{\pi \G h_\star \Sigma_\star}

   where we assume :math:`h_\star/R_\mathrm{disk}=`\ ``[heightToRadialScaleDiskBlitzRosolowsky]``.
   </description>
  </starFormationRateSurfaceDensityDisks>
  !!]
  type, extends(starFormationRateSurfaceDensityDisksClass) :: starFormationRateSurfaceDensityDisksBlitz2006
     !!{RST
     Implementation of the :cite:t:`blitz_role_2006` star formation rate surface density law for galactic disks.
     !!}
     private
     integer         (kind_int8             )                                :: lastUniqueID
     logical                                                                 :: factorsComputed                     , assumeMonotonicSurfaceDensity       , &
          &                                                                     isExponentialDisk                   , useTabulation
     double precision                                                        :: heightToRadialScaleDisk             , pressureCharacteristic              , &
          &                                                                     pressureExponent                    , starFormationFrequencyNormalization , &
          &                                                                     surfaceDensityCritical              , surfaceDensityExponent              , &
          &                                                                     velocityDispersionDiskGas           , radiusDisk                          , &
          &                                                                     massGas                             , hydrogenMassFraction                , &
          &                                                                     massStellar                         , massGasPrevious                     , &
          &                                                                     massStellarPrevious                 , hydrogenMassFractionPrevious        , &
          &                                                                     radiusDiskPrevious                  , radiusCritical                      , &
          &                                                                     radiusCriticalPrevious              , factorBoostStellarCoefficient       , &
          &                                                                     pressureRatioCoefficient
     type            (rootFinder            )                                :: finder
     type            (fastExponentiator     )                                :: pressureRatioExponentiator
     double precision                        , allocatable, dimension(:,:,:) :: integralPartiallyMolecularTable
     logical                                                                 :: tableInitialized
     double precision                                                        :: coefficientFactorBoostMinimum       , coefficientFactorBoostMaximum       , &
          &                                                                     coefficientFactorBoostStellarMinimum, coefficientFactorBoostStellarMaximum, &
          &                                                                     radiusScaleFreeMinimum              , radiusScaleFreeMaximum
     ! Lattices to which the three axes of the tabulation are pinned. These are the source of truth for its extent: the limits
     ! above are derived from them, and so are the offset and inverse step which the interpolation uses, which are therefore no
     ! longer held on the object.
     type            (rangeLattice          )                                :: latticeFactorBoost                  , latticeFactorBoostStellar           , &
          &                                                                     latticeRadiusScaleFree
     type            (varying_string        )                                :: filenameTable
   contains
     !![
     <methods docformat="rst">
       <method description="Reset memoized calculations."                                                                       method="calculationReset"          />
       <method description="Compute various factors."                                                                           method="computeFactors"            />
       <method description="Compute the pressure ratio."                                                                        method="pressureRatio"             />
       <method description="Compute the integral of the star formation rate surface density in the fully-molecular regime."     method="integralFullyMolecular"    />
       <method description="Compute the integral of the star formation rate surface density in the partially-molecular regime." method="integralPartiallyMolecular"/>
     </methods>
     !!]
     final     ::                               blitz2006Destructor
     procedure :: autoHook                   => blitz2006AutoHook
     procedure :: calculationReset           => blitz2006CalculationReset
     procedure :: rate                       => blitz2006Rate
     procedure :: computeFactors             => blitz2006ComputeFactors
     procedure :: unchanged                  => blitz2006Unchanged
     procedure :: intervals                  => blitz2006Intervals
     procedure :: pressureRatio              => blitz2006PressureRatio
     procedure :: integralFullyMolecular     => blitz2006IntegralFullyMolecular
     procedure :: integralPartiallyMolecular => blitz2006IntegralPartiallyMolecular
  end type starFormationRateSurfaceDensityDisksBlitz2006

  interface starFormationRateSurfaceDensityDisksBlitz2006
     !!{RST
     Constructors for the :galacticus-class:`starFormationRateSurfaceDensityDisksBlitz2006` star formation surface density rate in disks class.
     !!}
     module procedure blitz2006ConstructorParameters
     module procedure blitz2006ConstructorInternal
  end interface starFormationRateSurfaceDensityDisksBlitz2006

  ! Submodule-scope pointer to the active node.
  class           (starFormationRateSurfaceDensityDisksBlitz2006), pointer   :: self_
  type            (treeNode                                     ), pointer   :: node_
  !$omp threadprivate(self_,node_)

contains

  function blitz2006ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`starFormationRateSurfaceDensityDisksBlitz2006` star formation surface density rate in disks class which takes a parameter set as input.
    !!}
    implicit none
    type            (starFormationRateSurfaceDensityDisksBlitz2006)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    double precision                                                               :: velocityDispersionDiskGas          , heightToRadialScaleDisk, &
         &                                                                            surfaceDensityCritical             , surfaceDensityExponent , &
         &                                                                            starFormationFrequencyNormalization, pressureCharacteristic , &
         &                                                                            pressureExponent
    logical                                                                        :: assumeMonotonicSurfaceDensity      , useTabulation

    !![
    <inputParameter docformat="rst">
      <name>velocityDispersionDiskGas</name>
      <defaultSource>
      :cite:p:`leroy_star_2008`
      </defaultSource>
      <defaultValue>10.0d0</defaultValue>
      <description>
      The velocity dispersion of gas in galactic disks (in km/s), used to compute the hydrostatic midplane pressure that determines the molecular-to-atomic gas ratio in the :cite:t:`blitz_role_2006` star formation model.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>heightToRadialScaleDisk</name>
      <defaultSource>
      :cite:p:`kregel_flattening_2002`
      </defaultSource>
      <defaultValue>0.137d0</defaultValue>
      <description>
      The ratio of scale height to scale radius for disks in the :cite:t:`blitz_role_2006` star formation timescale calculation.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>surfaceDensityCritical</name>
      <defaultSource>
      :cite:p:`bigiel_star_2008`
      </defaultSource>
      <defaultValue>200.0d0</defaultValue>
      <description>
      The surface density (in units of :math:`\mathrm{M}_\odot` pc\ :math:`^{-2}`) in the :cite:t:`blitz_role_2006` star formation timescale calculation at which low-density truncation begins.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>surfaceDensityExponent</name>
      <defaultSource>
      :cite:p:`bigiel_star_2008`
      </defaultSource>
      <defaultValue>0.4d0</defaultValue>
      <description>
      The exponent for surface density in the :cite:t:`blitz_role_2006` star formation timescale calculation at in the high density regime.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>starFormationFrequencyNormalization</name>
      <defaultSource>
      :cite:p:`leroy_star_2008`
      </defaultSource>
      <defaultValue>5.25d-10</defaultValue>
      <description>
      The star formation frequency (in the low-density limit and in units of yr\ :math:`^{-1}`) in the :cite:t:`blitz_role_2006` star formation timescale calculation.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>pressureCharacteristic</name>
      <defaultSource>
      :cite:p:`blitz_role_2006`
      </defaultSource>
      <defaultValue>4.54d0</defaultValue>
      <description>
      The characteristic pressure (given as :math:`P_0/k_\mathrm{B}` in units of K cm\ :math:`^{-3}`) in the scaling relation of molecular hydrogen fraction with disk pressure in the :cite:t:`blitz_role_2006` star formation timescale calculation.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>pressureExponent</name>
      <defaultSource>
      :cite:p:`blitz_role_2006`
      </defaultSource>
      <defaultValue>0.92d0</defaultValue>
      <description>
      The exponent in the scaling relation of molecular hydrogen fraction with disk pressure in the :cite:t:`blitz_role_2006` star formation timescale calculation.
      </description>
      <source>parameters</source>
      <minimum>0.0</minimum>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>assumeMonotonicSurfaceDensity</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, assume that the surface density in disks is always monotonically decreasing.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>useTabulation</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, then use tabulated solutions to the integrated star formation rate.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=starFormationRateSurfaceDensityDisksBlitz2006(velocityDispersionDiskGas,heightToRadialScaleDisk,surfaceDensityCritical,surfaceDensityExponent,starFormationFrequencyNormalization,pressureCharacteristic,pressureExponent,assumeMonotonicSurfaceDensity,useTabulation)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function blitz2006ConstructorParameters

  function blitz2006ConstructorInternal(velocityDispersionDiskGas,heightToRadialScaleDisk,surfaceDensityCritical,surfaceDensityExponent,starFormationFrequencyNormalization,pressureCharacteristic,pressureExponent,assumeMonotonicSurfaceDensity,useTabulation) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`starFormationRateSurfaceDensityDisksBlitz2006` star formation surface density rate in disks class.
    !!}
    use :: Error                           , only : Error_Report
    use :: Input_Paths                     , only : inputPath                , pathTypeDataDynamic
    use :: Numerical_Constants_Astronomical, only : massSolar                , megaParsec
    use :: Numerical_Constants_Physical    , only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes    , only : giga                     , hecto                        , kilo                         , mega
    use :: Root_Finder                     , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    use :: Hashes_Cryptographic            , only : Hash_MD5
    use :: String_Handling                 , only : String_C_To_Fortran
    implicit none
    type            (starFormationRateSurfaceDensityDisksBlitz2006)                :: self
    double precision                                               , intent(in   ) :: velocityDispersionDiskGas          , heightToRadialScaleDisk, &
         &                                                                            surfaceDensityCritical             , surfaceDensityExponent , &
         &                                                                            starFormationFrequencyNormalization, pressureCharacteristic , &
         &                                                                            pressureExponent
    logical                                                        , intent(in   ) :: assumeMonotonicSurfaceDensity      , useTabulation
    type            (varying_string                               )                :: descriptorString
    character       (len=17                                       )                :: parameterLabel
    !![
    <constructorAssign variables="velocityDispersionDiskGas, heightToRadialScaleDisk, surfaceDensityCritical, surfaceDensityExponent, starFormationFrequencyNormalization, pressureCharacteristic, pressureExponent, assumeMonotonicSurfaceDensity, useTabulation"/>
    !!]

    self%lastUniqueID   =-1_kind_int8
    self%factorsComputed=.false.
    ! Validate
    if (pressureExponent < 0.0d0) call Error_Report('pressureExponent < 0 violates assumptions'//{introspection:location})
    ! Convert parameters to internal units.
    self%surfaceDensityCritical             =self%surfaceDensityCritical*(mega**2)                                                    ! Convert to M☉/Mpc².
    self%starFormationFrequencyNormalization=self%starFormationFrequencyNormalization*giga                                            ! Convert to Gyr⁻¹.
    self%pressureCharacteristic             =self%pressureCharacteristic*boltzmannsConstant*((hecto*megaParsec)**3)/massSolar/kilo**2 ! Convert to M☉(km/s)²/Mpc.
    ! Build fast exponentiator.
    self%pressureRatioExponentiator         =fastExponentiator(0.0d0,1.0d0,pressureExponent,1000.0d0,.false.)
    ! Build root finder.
    self%finder=rootFinder(                                                             &
         &                 rootFunction                 =blitz2006CriticalDensityRoot , &
         &                 toleranceAbsolute            =0.0d+0                       , &
         &                 toleranceRelative            =1.0d-4                       , &
         &                 rangeExpandUpward            =2.0d+0                       , &
         &                 rangeExpandDownward          =0.5d+0                       , &
         &                 rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
         &                 rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
         &                 rangeExpandType              =rangeExpandMultiplicative      &
         &                )
    ! Initialize memoized values.
    self%massGasPrevious                     =-huge(0.0d0)
    self%massStellarPrevious                 =-huge(0.0d0)
    self%radiusDiskPrevious                  =-huge(0.0d0)
    self%hydrogenMassFractionPrevious        =-huge(0.0d0)
    self%radiusCriticalPrevious              =-huge(0.0d0)
    ! Initialize table.
    self%tableInitialized                    =.false.
    self%coefficientFactorBoostMinimum       =+huge(0.0d0)
    self%coefficientFactorBoostMaximum       =-huge(0.0d0)
    self%coefficientFactorBoostStellarMinimum=+huge(0.0d0)
    self%coefficientFactorBoostStellarMaximum=-huge(0.0d0)
    self%radiusScaleFreeMinimum              =+huge(0.0d0)
    self%radiusScaleFreeMaximum              =-huge(0.0d0)
    ! Generate a file name for the table using the two parameters upon which it depends to create a hashed descriptor suffix.
    ! Deliberately only these two: the tabulated integral is scale free in every other parameter of this class, so including
    ! them would force the table - which is three-dimensional, and each of whose points is a numerical integral - to be rebuilt
    ! on changes which cannot affect it. The source digest is included because the tabulated values do depend on the code which
    ! generates them, which no descriptor of the parameters would capture.
    descriptorString=""
    write (parameterLabel,'(e17.10)') self%pressureExponent
    descriptorString=descriptorString//parameterLabel//" "
    write (parameterLabel,'(e17.10)') self%surfaceDensityExponent
    descriptorString=descriptorString//parameterLabel//" "
    descriptorString=descriptorString//String_C_To_Fortran(blitz2006SourceDigest)
    self%filenameTable                       =     inputPath  (pathTypeDataDynamic)// &
         &                                    'starFormation/'                     // &
         &                                    self%objectType (                   )// &
         &                                    '_'                                  // &
         &                                    Hash_MD5        (descriptorString   )// &
         &                                    '.hdf5'
    return
  end function blitz2006ConstructorInternal

  subroutine blitz2006AutoHook(self)
    !!{RST
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self

    call calculationResetEvent%attach(self,blitz2006CalculationReset,openMPThreadBindingAllLevels,label='starFormationRateSurfaceDensityDisksBlitz2006')
    return
  end subroutine blitz2006AutoHook

  subroutine blitz2006Destructor(self)
    !!{RST
    Destructor for the blitz2006 cooling radius class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self

    if (calculationResetEvent%isAttached(self,blitz2006CalculationReset)) call calculationResetEvent%detach(self,blitz2006CalculationReset)
    return
  end subroutine blitz2006Destructor

  subroutine blitz2006CalculationReset(self,node,uniqueID)
    !!{RST
    Reset the Kennicutt-Schmidt relation calculation.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self
    type   (treeNode                                     ), intent(inout) :: node
    integer(kind_int8                                    ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%factorsComputed       =.false.
    self%lastUniqueID          =uniqueID
    self%radiusCriticalPrevious=-huge(0.0d0)
    return
  end subroutine blitz2006CalculationReset

  double precision function blitz2006Rate(self,node,radius)
    !!{RST
    Returns the star formation rate surface density (in :math:`\mathrm{M}_\odot` Gyr\ :math:`^{-1}` Mpc\ :math:`^{-2}`) for star formation in the galactic disk of ``node``. The disk is assumed to obey the :cite:t:`blitz_role_2006` star formation rule.
    !!}
    implicit none
    class           (starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    double precision                                               , intent(in   ) :: radius
    double precision                                                               :: molecularFraction, pressureRatio, &
         &                                                                            surfaceDensityGas, factorBoost

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! Compute factors.
    call self%computeFactors(node)
    ! Return zero rate for non-positive radius or negligible mass.
    if (self%massGas <= massGasFloor .or. self%massStellar < 0.0d0 .or. self%radiusDisk <= 0.0d0) then
       blitz2006Rate=0.0d0
       return
    end if
    ! Compute the pressure ratio that Blitz & Rosolowsky (2006) use to compute the molecular fraction.    
    pressureRatio=self%pressureRatio(node,radius,surfaceDensityGas)
    ! Compute the molecular fraction, limited to 100% molecular.
    if (pressureRatio >= 1.0d0) then
       molecularFraction=                                                                1.0d0
    else
       molecularFraction=min(self%pressureRatioExponentiator%exponentiate(pressureRatio),1.0d0)
    end if
    ! Compute the star formation rate surface density.
    factorBoost  =+self%hydrogenMassFraction                  &
         &        *surfaceDensityGas                          &
         &        /self%surfaceDensityCritical
    blitz2006Rate=+surfaceDensityGas                          &
         &        *self%hydrogenMassFraction                  &
         &        *molecularFraction                          &
         &        *self%starFormationFrequencyNormalization   &
         &        *(                                          &
         &          +1.0d0                                    &
         &          +factorBoost**self%surfaceDensityExponent &
         &         )
    return
  end function blitz2006Rate

  logical function blitz2006Unchanged(self,node)
    !!{RST
    Determine if the surface rate density of star formation is unchanged.
    !!}
    implicit none
    class(starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self
    type (treeNode                                     ), intent(inout) :: node

    call self%computeFactors(node)
    blitz2006Unchanged= self%massGas              == self%massGasPrevious              &
         &             .and.                                                           &
         &              self%massStellar          == self%massStellarPrevious          &
         &             .and.                                                           &
         &              self%radiusDisk           == self%radiusDiskPrevious           &
         &             .and.                                                           &
         &              self%hydrogenMassFraction == self%hydrogenMassFractionPrevious
    if (.not.blitz2006Unchanged) then
       self%massGasPrevious             =self%massGas
       self%massStellarPrevious         =self%massStellar
       self%radiusDiskPrevious          =self%radiusDisk
       self%hydrogenMassFractionPrevious=self%hydrogenMassFraction
    end if
    return
  end function blitz2006Unchanged
  
  subroutine blitz2006ComputeFactors(self,node)
    !!{RST
    Compute various factors for the ``blitz2006`` star formation rate surface density calculation.
    !!}
    use :: Abundances_Structure            , only : abundances
    use :: Galacticus_Nodes                , only : nodeComponentDisk
    use :: Galactic_Structure_Options      , only : componentTypeDisk             , massTypeGaseous                  , massTypeStellar
    use :: Mass_Distributions              , only : massDistributionClass         , massDistributionCylindricalScaler, massDistributionExponentialDisk
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class(starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self
    type (treeNode                                     ), intent(inout) :: node
    class(nodeComponentDisk                            ), pointer       :: disk
    class(massDistributionClass                        ), pointer       :: massDistributionGaseous, massDistributionStellar, &
         &                                                                 massDistribution_
    type (abundances                                   ), save          :: abundancesFuel
    !$omp threadprivate(abundancesFuel)

    ! Check if factors have been precomputed.
    if (.not.self%factorsComputed) then
       ! Get the disk properties.
       disk         => node%disk   ()
       self%massGas =  disk%massGas()
       if (self%massGas > massGasFloor) then
          self%massStellar=disk%massStellar()
          self%radiusDisk =disk%radius     ()
          ! Find the hydrogen fraction in the disk gas of the fuel supply.
          abundancesFuel=disk%abundancesGas()
          call abundancesFuel%massToMassFraction(self%massGas)
          self%hydrogenMassFraction=abundancesFuel%hydrogenMassFraction()
          ! Determine if we have an exponential disk.
          massDistributionGaseous => node%massDistribution(componentType=componentTypeDisk,massType=massTypeGaseous)
          massDistributionStellar => node%massDistribution(componentType=componentTypeDisk,massType=massTypeStellar)
          self%isExponentialDisk  =  .true.
          select type (massDistributionGaseous)
          class is (massDistributionExponentialDisk)
             ! The disk is exponential - no change needed.
          class is (massDistributionCylindricalScaler     )
             ! Check the unscale distribution.
             massDistribution_ => massDistributionGaseous%unscaled()
             select type (massDistribution_)
             class is (massDistributionExponentialDisk)
                ! The disk is exponential - no change needed.
             class default
                self%isExponentialDisk=.false.
             end select
          class default
             ! Not an exponential distribution.
             self%isExponentialDisk=.false.
          end select
          select type (massDistributionStellar)
          class is (massDistributionExponentialDisk)
             ! The disk is exponential - no change needed.
          class is (massDistributionCylindricalScaler     )
             ! Check the unscale distribution.
             massDistribution_ => massDistributionStellar%unscaled()
             select type (massDistribution_)
             class is (massDistributionExponentialDisk)
                ! The disk is exponential - no change needed.
             class default
                self%isExponentialDisk=.false.
             end select
          class default
             ! Not an exponential distribution.
             self%isExponentialDisk=.false.
          end select
          !![
	  <objectDestructor name="massDistributionGaseous"/>
	  <objectDestructor name="massDistributionStellar"/>
	  !!]
          ! Properties required for exponential disks.
          if (self%isExponentialDisk .and. self%massStellar >= 0.0d0 .and. self%radiusDisk > 0.0d0) then
             self%pressureRatioCoefficient     =+gravitationalConstant_internal          &
                  &                             /8.0d0                                   &
                  &                             /Pi                                      &
                  &                             *self%massGas                        **2 &
                  &                             /self%pressureCharacteristic             &
                  &                             /self%radiusDisk                     **4
             self%factorBoostStellarCoefficient=+self%velocityDispersionDiskGas          &
                  &                             *2.0d0                                   &
                  &                             *Pi                                      &
                  &                             *self%radiusDisk                     **2 &
                  &                             /self%massGas                            &
                  &                             *sqrt(                                   &
                  &                                   +self%massStellar                  &
                  &                                   /2.0d0                             &
                  &                                   /Pi                            **2 &
                  &                                   /gravitationalConstant_internal    &
                  &                                   /self%heightToRadialScaleDisk      &
                  &                                   /self%radiusDisk               **3 &
                  &                                  ) 
          end if
       else
          ! No gas mass, so other factors are irrelevant.
          self%massStellar         =0.0d0
          self%radiusDisk          =0.0d0
          self%hydrogenMassFraction=0.0d0
       end if
       ! Record that factors have now been computed.
       self%factorsComputed=.true.
    end if
    return
  end subroutine blitz2006ComputeFactors

  function blitz2006Intervals(self,node,radiusInner,radiusOuter,intervalIsAnalytic,integralsAnalytic)
    !!{RST
    Returns intervals to use for integrating the :cite:t:`krumholz_star_2009` star formation rate over a galactic disk.
    !!}
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Galactic_Structure_Options      , only : componentTypeDisk             , massTypeGaseous, massTypeStellar
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (starFormationRateSurfaceDensityDisksBlitz2006), intent(inout), target                      :: self
    double precision                                                              , allocatable, dimension(:,:) :: blitz2006Intervals
    type            (treeNode                                     ), intent(inout), target                      :: node
    double precision                                               , intent(in   )                              :: radiusInner                             , radiusOuter
    logical                                                        , intent(inout), allocatable, dimension(  :) :: intervalIsAnalytic
    double precision                                               , intent(inout), allocatable, dimension(  :) :: integralsAnalytic
    class           (massDistributionClass                        ), pointer                                    :: massDistributionGaseous                 , massDistributionStellar
    double precision                                               , parameter                                  :: factorBoostStellarCoefficientTiny=1.0d-6
    double precision                                                                                            :: coefficientNormalization                , coefficientFactorBoost       , &
         &                                                                                                         coefficientFactorBoostStellar           , coefficientMolecular         , &
         &                                                                                                         rootValueInner                          , rootValueOuter               , &
         &                                                                                                         radiusAnalytic                          , sqrtTerm
    logical                                                                                                     :: thresholdCondition                      , assumeMonotonicSurfaceDensity

    ! Check if we can assume a monotonic surface density.
    massDistributionGaseous       =>  node                   %massDistribution                       (componentType=componentTypeDisk,massType=massTypeGaseous)
    massDistributionStellar       =>  node                   %massDistribution                       (componentType=componentTypeDisk,massType=massTypeStellar)
    assumeMonotonicSurfaceDensity =   massDistributionGaseous%assumeMonotonicDecreasingSurfaceDensity(                                                        ) &
         &                           .and.                                                                                                                      &
         &                            massDistributionStellar%assumeMonotonicDecreasingSurfaceDensity(                                                        )
    !![
    <objectDestructor name="massDistributionGaseous"/>
    <objectDestructor name="massDistributionStellar"/>
    !!]
    if (assumeMonotonicSurfaceDensity) then
       ! Set the critical radius to a very negative value so that pressure ratio is always computed.
       self%radiusCritical=-huge(0.0d0)
       ! Compute factors.
       call self%computeFactors(node)
       ! Set zero intervals for non-positive radius or negligible mass.
       if (self%massGas <= massGasFloor .or. self%massStellar < 0.0d0 .or. self%radiusDisk <= 0.0d0) then
          allocate(blitz2006Intervals(2,0))
          self%radiusCritical=-huge(0.0d0)          
       else
          self_ => self
          node_ => node
          ! Test if the inner radius is below the pressure threshold.
          if (self%isExponentialDisk) then
             ! For exponential disks this condition has a simple analytic form.
             rootValueInner       =-huge(0.0d0)
             if (self%pressureRatioCoefficient > 0.0d0 .and. -exponent(self%pressureRatioCoefficient) < maxExponent(0.0d0)) then
                thresholdCondition=1.0d0/self%pressureRatioCoefficient-self%factorBoostStellarCoefficient >= 1.0d0
             else
                thresholdCondition=.true.
             end if
           else
             ! For generic disks test this numerically.
             rootValueInner       =blitz2006CriticalDensityRoot(radiusInner)
             thresholdCondition   =rootValueInner                                                         <= 0.0d0
          end if          
          if (thresholdCondition) then
             ! The entire disk is below the pressure threshold so use a single interval.
             allocate(blitz2006Intervals(2,1))
             allocate(intervalIsAnalytic(  1))
             if (self%isExponentialDisk.and.self%useTabulation) then
                call computeCoefficients()
                allocate(integralsAnalytic(1))
                intervalIsAnalytic=.true.
                integralsAnalytic =self%integralPartiallyMolecular(coefficientNormalization*coefficientMolecular,coefficientFactorBoost,coefficientFactorBoostStellar,radiusInner,radiusOuter)
             else
                intervalIsAnalytic=.false.
             end if
             blitz2006Intervals =reshape([radiusInner,radiusOuter],[2,1])
             self%radiusCritical=-huge(0.0d0)
          else
             ! Compute coefficients needed for analytic solutions.
             if (self%isExponentialDisk.and.self%useTabulation) call computeCoefficients()
             ! Test the surface density at the outer radius.
             rootValueOuter=blitz2006CriticalDensityRoot(radiusOuter)
             if (rootValueOuter >= 0.0d0) then
                ! Entire disk is above the pressure threshold so use a single interval.
                allocate(blitz2006Intervals(2,1))
                allocate(intervalIsAnalytic(  1))
                if (self%isExponentialDisk.and.self%useTabulation) then
                   allocate(integralsAnalytic (  1))
                   intervalIsAnalytic =.true.
                   integralsAnalytic  =self%integralFullyMolecular(coefficientNormalization,coefficientFactorBoost,radiusInner,radiusOuter)
                else
                   intervalIsAnalytic=.false.
                end if
                blitz2006Intervals =reshape([radiusInner,radiusOuter],[2,1])
                self%radiusCritical=radiusOuter
             else
                ! The disk transitions the pressure threshold - attempt to locate the radius at which this happens and use two
                ! intervals split at this point.
                if (self%isExponentialDisk) then
                   ! For exponential disks we have an analytic solution for the transition radius.
                   if (self%factorBoostStellarCoefficient <= factorBoostStellarCoefficientTiny) then
                      radiusAnalytic=+0.5d0*log(self%pressureRatioCoefficient)
                   else
                      sqrtTerm      =+(                                                                                      &
                           &                 +  9.0d0*self%pressureRatioCoefficient**2*self%factorBoostStellarCoefficient**2 &
                           &           +sqrt(                                                                                &
                           &                 +  3.0d0                                                                        &
                           &                )                                                                                &
                           &           *sqrt(                                                                                &
                           &                 +256.0d0*self%pressureRatioCoefficient**3                                       &
                           &                 + 27.0d0*self%pressureRatioCoefficient**4*self%factorBoostStellarCoefficient**4 &
                           &                )                                                                                &
                           &          )**(1.0d0/3.0d0)
                      radiusAnalytic=+2.0d0                                                                                                                  &
                           &         *log(                                                                                                                   &
                           &              +0.5d0                                                                                                             &
                           &              *sqrt(                                                                                                             &
                           &                          -4.0d0* (2.0d0/3.0d0)**(1.0d0/3.0d0)                      *self%pressureRatioCoefficient     /sqrtTerm &
                           &                          +1.0d0/( 2.0d0       **(1.0d0/3.0d0)*3.0d0**(2.0d0/3.0d0))                                   *sqrtTerm &
                           &                   )                                                                                                             &
                           &              +0.5d0                                                                                                             &
                           &              *sqrt(                                                                                                             &
                           &                          +4.0d0* (2.0d0/3.0d0)**(1.0d0/3.0d0)                      *self%pressureRatioCoefficient     /sqrtTerm &
                           &                          -1.0d0/( 2.0d0       **(1.0d0/3.0d0)*3.0d0**(2.0d0/3.0d0))                                   *sqrtTerm &
                           &                          +2.0d0                                                                                                 &
                           &                                                                                    *self%pressureRatioCoefficient               &
                           &                                                                                    *self%factorBoostStellarCoefficient          &
                           &                    /sqrt(                                                                                                       &
                           &                          -4.0d0* (2.0d0/3.0d0)**(1.0d0/3.0d0)                      *self%pressureRatioCoefficient     /sqrtTerm &
                           &                          +1.0d0/( 2.0d0       **(1.0d0/3.0d0)*3.0d0**(2.0d0/3.0d0))                                   *sqrtTerm &
                           &                         )                                                                                                       &
                           &                   )                                                                                                             &
                           &             )
                   end if
                   self%radiusCritical=+     radiusAnalytic &
                        &              *self%radiusDisk
                else
                   ! For non-exponential disks, seek a solution numerically.
                   if (self%radiusCriticalPrevious > 0.0d0) then
                      self%radiusCritical=self%finder%find(rootGuess=self%radiusCriticalPrevious)
                   else
                      self%radiusCritical=self%finder%find(rootRange=[radiusInner,radiusOuter],rootRangeValues=[rootValueInner,rootValueOuter])
                   end if
                end if
                self%radiusCriticalPrevious=self%radiusCritical
                allocate(blitz2006Intervals(2,2))
                allocate(intervalIsAnalytic(  2))
                if (self%isExponentialDisk.and.self%useTabulation) then
                   allocate(integralsAnalytic (  2))
                   intervalIsAnalytic =.true.
                   integralsAnalytic  =[                                                                                                                                                                             &
                        &               self%integralFullyMolecular    (coefficientNormalization                     ,coefficientFactorBoost                              ,     radiusInner   ,self%radiusCritical), &
                        &               self%integralPartiallyMolecular(coefficientNormalization*coefficientMolecular,coefficientFactorBoost,coefficientFactorBoostStellar,self%radiusCritical,     radiusOuter   )  &
                        &              ]                   
                else
                   intervalIsAnalytic=.false.
                end if
                blitz2006Intervals=reshape([radiusInner,self%radiusCritical,self%radiusCritical,radiusOuter],[2,2])
             end if
          end if
       end if
    else
       ! Disk pressure can not be assumed to be monotonic - use a single interval.
       allocate(blitz2006Intervals(2,1))
       allocate(intervalIsAnalytic(  1))
       intervalIsAnalytic=.false.
       blitz2006Intervals=reshape([radiusInner,radiusOuter],[2,1])
       self%radiusCritical=radiusInner
    end if
    return

  contains

    subroutine computeCoefficients()
      !!{RST
      Compute coefficients needed in analytic and tabulated solutions.
      !!}
      implicit none
      
      coefficientNormalization     =    +self%massGas/2.0d0/Pi                       *self%hydrogenMassFraction*self%starFormationFrequencyNormalization
      coefficientFactorBoost       =    +self%massGas/2.0d0/Pi/self%radiusDisk**2    *self%hydrogenMassFraction/self%surfaceDensityCritical
      coefficientMolecular         =+(                                                &
           &                          +(+self%massGas/2.0d0/Pi/self%radiusDisk**2)**2 &
           &                          *gravitationalConstant_internal                 &
           &                          *Pi                                             &
           &                          /2.0d0                                          &
           &                          /self%pressureCharacteristic                    &
           &                         )**self%pressureExponent
      coefficientFactorBoostStellar=+self%velocityDispersionDiskGas       &
           &                        *2.0d0                                &
           &                        *Pi                                   &
           &                        *self%radiusDisk**2                   &
           &                        /self%massGas                         &
           &                        *sqrt(                                &
           &                              +self%massStellar               &
           &                              /2.0d0                          &
           &                              /Pi                             &
           &                              /self%radiusDisk**2             &
           &                              /Pi                             &
           &                              /gravitationalConstant_internal &
           &                              /self%heightToRadialScaleDisk   &
           &                              /self%radiusDisk                &
           &                             ) 
      return
    end subroutine computeCoefficients

  end function blitz2006Intervals

  double precision function blitz2006IntegralFullyMolecular(self,coefficientNormalization,coefficientFactorBoost,radiusInner,radiusOuter)
    !!{RST
    Evaluate the integral of the star formation rate surface density in the fully-molecular regime.
    !!}
    implicit none
    class           (starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self
    double precision                                               , intent(in   ) :: coefficientNormalization, coefficientFactorBoost, &
         &                                                                            radiusInner             , radiusOuter
    
    blitz2006IntegralFullyMolecular=+integralAnalyticFullyMolecular(radiusOuter) &
         &                          -integralAnalyticFullyMolecular(radiusInner)
    return
    
  contains
    
    double precision function integralAnalyticFullyMolecular(r)
      !!{RST
      Analytic solution to the improper integral of the star formation rate surface density over an exponential disk.
      !!}
      implicit none
      double precision, intent(in   ) :: r
      double precision                :: x

      x                             =+r                                                  &
           &                         /self%radiusDisk
      integralAnalyticFullyMolecular=+coefficientNormalization                           &
           &                         *exp(-x)                                            &
           &                         *(                                                  &
           &                           -1.0d0                                            &
           &                           -x                                                &
           &                           -(                                                &
           &                             +coefficientFactorBoost                         &
           &                             *exp(-x)                                        &
           &                            )**self%surfaceDensityExponent                   &
           &                           *(1.0d0+x*(1.0d0+self%surfaceDensityExponent)   ) &
           &                           /         (1.0d0+self%surfaceDensityExponent)**2  &
           &                          )
      return
    end function integralAnalyticFullyMolecular
    
  end function blitz2006IntegralFullyMolecular
  
  double precision function blitz2006IntegralPartiallyMolecular(self,coefficientNormalization,coefficientFactorBoost,coefficientFactorBoostStellar,radiusInner,radiusOuter)
    !!{RST
    Evaluate the integral of the star formation rate surface density in the fully-molecular regime.
    !!}
    implicit none
    class           (starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self
    double precision                                               , intent(in   ) :: radiusInner                              , radiusOuter                             , &
         &                                                                            coefficientNormalization                 , coefficientFactorBoost                  , &
         &                                                                            coefficientFactorBoostStellar
    double precision                                               , parameter     :: coefficientFactorBoostStellarLarge=1.0d+4, radiusScaleFreeTiny              =1.0d-3
    double precision                                               , parameter     :: coefficientFactorBoostTiny        =1.0d-6, coefficientFactorBoostStellarTiny=1.0d-6
    double precision                                               , save          :: coefficientFactorBoost_                  , coefficientFactorBoostStellar_
    !$omp threadprivate(coefficientFactorBoost_,coefficientFactorBoostStellar_)
    integer                                                                        :: i
    double precision                                                               :: multiplier                               , radiusScaleFree                         , &
         &                                                                            radiusScaleFreeInner                     , radiusScaleFreeOuter
    
    ! Compute scale free radii.
    radiusScaleFreeInner=radiusInner/self%radiusDisk
    radiusScaleFreeOuter=radiusOuter/self%radiusDisk
    ! Handle cases where we have an (approximate) analytic solution.
    if      (coefficientFactorBoostStellar <= 0.0d0                            ) then
       ! If the stellar boost is zero, we have an analytic solution.
       blitz2006IntegralPartiallyMolecular=+integralAnalyticPartiallyMolecularZeroStellarBoost (radiusScaleFreeOuter) &
            &                              -integralAnalyticPartiallyMolecularZeroStellarBoost (radiusScaleFreeInner)
       return
    else if (coefficientFactorBoostStellar > coefficientFactorBoostStellarLarge) then
       ! If the stellar boost is large, we have an approximate analytic solution.
       blitz2006IntegralPartiallyMolecular=+integralAnalyticPartiallyMolecularLargeStellarBoost(radiusScaleFreeOuter) &
            &                              -integralAnalyticPartiallyMolecularLargeStellarBoost(radiusScaleFreeInner)
       return
    else
       ! General case - no analytic solution available.
       blitz2006IntegralPartiallyMolecular=0.0d0
       do i=1,2
         select case (i)
         case (1)
            radiusScaleFree=radiusInner/self%radiusDisk
            multiplier     =-1.0d0
         case (2)
            radiusScaleFree=radiusOuter/self%radiusDisk
            multiplier     =+1.0d0
         end select
         ! Use series or tabulated solutions.
         if (radiusScaleFree < radiusScaleFreeTiny) then
            ! Small radius limit - use a series solution.
            blitz2006IntegralPartiallyMolecular=+blitz2006IntegralPartiallyMolecular                           &
                 &                              +multiplier                                                    &
                 &                              *integralAnalyticPartiallyMolecularSmallRadii(radiusScaleFree)
         else
            ! No approximation available - use the tabulated solution.
            blitz2006IntegralPartiallyMolecular=+blitz2006IntegralPartiallyMolecular                           &
                 &                              +multiplier                                                    &
                 &                              *integralAnalyticPartiallyMolecularGeneric   (radiusScaleFree)
         end if
      end do
   end if
   return

  contains

    double precision function integralAnalyticPartiallyMolecularZeroStellarBoost(x)
      !!{RST
      Analytic solution to the improper integral of the star formation rate surface density over an exponential disk for the case of zero stellar boost factor..
      !!}
      implicit none
      double precision, intent(in   ) :: x

      integralAnalyticPartiallyMolecularZeroStellarBoost=+coefficientNormalization                                                       &
           &                                             *exp(-x*(1.0d0+2.0d0*self%pressureExponent))                                    &
           &                                             *(                                                                              &
           &                                               -1.0d0/(1.0d0+2.0d0*self%pressureExponent)**2                                 &
           &                                               -x    /(1.0d0+2.0d0*self%pressureExponent)                                    &
           &                                               -(                                                                            &
           &                                                 +coefficientFactorBoost                                                     &
           &                                                 *exp(-x)                                                                    &
           &                                                )**self%surfaceDensityExponent                                               &
           &                                               *(1.0d0+x*(1.0d0+self%surfaceDensityExponent+2.0d0*self%pressureExponent)   ) &
           &                                               /         (1.0d0+self%surfaceDensityExponent+2.0d0*self%pressureExponent)**2  &
           &                                              )
      return
    end function integralAnalyticPartiallyMolecularZeroStellarBoost
    
    double precision function integralAnalyticPartiallyMolecularLargeStellarBoost(x)
      !!{RST
      Analytic solution to the improper integral of the star formation rate surface density over an exponential disk for the case of large stellar boost factor.
      !!}
      implicit none
      double precision, intent(in   ) :: x

      integralAnalyticPartiallyMolecularLargeStellarBoost=+coefficientNormalization                                                            &
           &                                             *exp(-x*(1.0d0+2.0d0*self%pressureExponent))                                          &
           &                                             *(+coefficientFactorBoostStellar*exp(x/2))**self%pressureExponent                     &
           &                                             *(                                                                                    &
           &                                               -4.0d0  /(2.0d0+3.0d0*self%pressureExponent)**2                                     & 
           &                                               -2.0d0*x/(2.0d0+3.0d0*self%pressureExponent)                                        &
           &                                               -2.0d0                                                                              &
           &                                               *(                                                                                  &
           &                                                 +coefficientFactorBoost                                                           &
           &                                                 *exp(-x)                                                                          &
           &                                                )**self%surfaceDensityExponent                                                     &
           &                                               *(2.0d0+x*(2.0d0+2.0d0*self%surfaceDensityExponent+3.0d0*self%pressureExponent)   ) &
           &                                               /         (2.0d0+2.0d0*self%surfaceDensityExponent+3.0d0*self%pressureExponent)**2  &
           &                                              )
      return
    end function integralAnalyticPartiallyMolecularLargeStellarBoost
   
    double precision function integralAnalyticPartiallyMolecularSmallRadii(x)
      !!{RST
      Analytic solution to the improper integral of the star formation rate surface density over an exponential disk for the case of small radii. Uses a series solution.
      !!}
      implicit none
      double precision, intent(in   ) :: x

      integralAnalyticPartiallyMolecularSmallRadii=-coefficientNormalization                                                                                                 &
           &                                       /6.0d0                                                                                                                    &
           &                                       *(1.0d0+coefficientFactorBoostStellar)**(-1.0d0+self%pressureExponent)                                                    &
           &                                       *x**2                                                                                                                     &
           &                                       *(                                                                                                                        &
           &                                         +(1.0d0+      coefficientFactorBoostStellar)*(-3.0d0+2.0d0*x)                                                           &
           &                                         +(4.0d0+3.0d0*coefficientFactorBoostStellar)*x*self%pressureExponent                                                    &
           &                                         +coefficientFactorBoost**self%surfaceDensityExponent                                                                    &
           &                                         *(                                                                                                                      &
           &                                           -3.0d0                                                                                                                &
           &                                           +                                      2.0d0*x*(1.0d0+2.0d0*self%pressureExponent+      self%surfaceDensityExponent)  &
           &                                           +      coefficientFactorBoostStellar*(-3.0d0+x*(2.0d0+3.0d0*self%pressureExponent+2.0d0*self%surfaceDensityExponent)) &
           &                                          )                                                                                                                      &
           &                                        )
      return
    end function integralAnalyticPartiallyMolecularSmallRadii
   
    double precision function integralAnalyticPartiallyMolecularGeneric(radiusScaleFree)
      !!{RST
      Analytic solution to the improper integral of the star formation rate surface density over an exponential disk for the general case.
      !!}
      use :: Display              , only : displayCounter          , displayCounterClear      , displayIndent      , displayMessage, &
           &                               displayUnindent         , verbosityLevelWorking
      use :: Numerical_Integration, only : integrator              , GSL_Integ_Gauss61
      use :: HDF5_Access          , only : hdf5Access
      use :: IO_HDF5              , only : hdf5File
      use :: File_Utilities       , only : File_Exists             , File_Lock                , File_Unlock        , lockDescriptor, &
           &                               Directory_Make          , File_Path
      use :: Numerical_Ranges     , only : Range_Pinned            , Range_Lattice_Offset     , gridSchemePerDecade
      use :: Table_Caches         , only : Table_Cache_Lattice_Read, Table_Cache_Lattice_Write
      implicit none
      double precision                , intent(in   )                 :: radiusScaleFree
      double precision                , parameter                     :: toleranceRelative                 =1.0d-6
      integer                                                         :: countFactorBoost                         , countFactorBoostStellar        , &
           &                                                             countRadii                               , i                              , &
           &                                                             j                                        , k                              , &
           &                                                             ii                                       , jj                             , &
           &                                                             kk                                       , loopCount                      , &
           &                                                             loopCountTotal                           , offsetFactorBoost              , &
           &                                                             offsetFactorBoostStellar                 , offsetRadius                   , &
           &                                                             countFactorBoostPrevious                 , countFactorBoostStellarPrevious, &
           &                                                             countRadiiPrevious
      double precision                                                :: integral                                                                  , &
           &                                                             hi                                       , hj                             , &
           &                                                             hk                                       , hhi                            , &
           &                                                             hhj                                      , hhk                            , &
           &                                                             radiusMinimum                            , radiusMaximum                  , &
           &                                                             coordinateFactorBoost                    , coordinateFactorBoostStellar   , &
           &                                                             coordinateRadiusScaleFree
      type            (rangeLattice  )                                :: latticeFactorBoost                       , latticeFactorBoostStellar      , &
           &                                                             latticeRadiusScaleFree
      double precision                , allocatable, dimension(:    ) :: valuesFactorBoost                        , valuesFactorBoostStellar       , &
           &                                                             valuesRadiusScaleFree
      double precision                , allocatable, dimension(:,:,:) :: integralPartiallyMolecularPrevious
      logical                                                         :: carryOver                                , isCarriedFactorBoost           , &
           &                                                             isCarriedFactorBoostStellar
      character       (len= 8        )                                :: tableSize
      character       (len= 8        )                                :: countSteps
      character       (len=12        )                                :: rangeLower                               , rangeUpper
      type            (integrator    ), allocatable                   :: integrator_
      type            (varying_string), save                          :: message
      type            (lockDescriptor), save                          :: fileLock
      logical                                                         :: haveLock
      !$omp threadprivate(message,fileLock)

      ! If our table is insufficient (or does not yet exist), attempt to read the table from file.
      haveLock=.false.
      if (tableIsInsufficient()) then
         call Directory_Make(File_Path(self%filenameTable))
         call File_Lock(self%filenameTable,fileLock,lockIsShared=.false.)
         haveLock=.true.
         if (File_Exists(self%filenameTable)) then
            message='Reading Blitz2006 star formation rate tabulation from file: '//self%filenameTable
            call displayMessage(message,verbosityLevelWorking)
            ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
            !$ call hdf5Access%set()
            hdf5FileScopeRead: block
              type(hdf5File  ) :: file
              file=hdf5File(self%filenameTable,readOnly=.true.)
              ! Recover the lattices on which the stored tabulation was built. A file which records none, or which records
              ! lattices this object would not use, is ignored rather than misread.
              call Table_Cache_Lattice_Read(file,'coefficientFactorBoost'       ,gridSchemePerDecade,pointsPerDecadeFactorBoost       ,self%latticeFactorBoost       )
              call Table_Cache_Lattice_Read(file,'coefficientFactorBoostStellar',gridSchemePerDecade,pointsPerDecadeFactorBoostStellar,self%latticeFactorBoostStellar)
              call Table_Cache_Lattice_Read(file,'radiusScaleFree'              ,gridSchemePerDecade,pointsPerDecadeRadius            ,self%latticeRadiusScaleFree   )
              if     (                                             &
                   &   self%latticeFactorBoost       %isDefined()  &
                   &  .and.                                        &
                   &   self%latticeFactorBoostStellar%isDefined()  &
                   &  .and.                                        &
                   &   self%latticeRadiusScaleFree   %isDefined()  &
                   & ) call file%readDataset('integral',self%integralPartiallyMolecularTable)
            end block hdf5FileScopeRead
            !$ call hdf5Access%unset()
            ! Adopt the stored tabulation only if its lattices were recovered and the array stored alongside them matches.
            if     (                                                 &
                 &   self%latticeFactorBoost       %isDefined()      &
                 &  .and.                                            &
                 &   self%latticeFactorBoostStellar%isDefined()      &
                 &  .and.                                            &
                 &   self%latticeRadiusScaleFree   %isDefined()      &
                 &  .and.                                            &
                 &   allocated(self%integralPartiallyMolecularTable) &
                 & ) then
               if     (                                                                                          &
                    &   size(self%integralPartiallyMolecularTable,dim=1) == self%latticeFactorBoost       %count &
                    &  .and.                                                                                     &
                    &   size(self%integralPartiallyMolecularTable,dim=2) == self%latticeFactorBoostStellar%count &
                    &  .and.                                                                                     &
                    &   size(self%integralPartiallyMolecularTable,dim=3) == self%latticeRadiusScaleFree   %count &
                    & ) then
                  ! Recover the limits from the lattices rather than reading them from the file, so that a restored tabulation
                  ! cannot come to be described differently from a freshly built one.
                  self%coefficientFactorBoostMinimum       =self%latticeFactorBoost       %minimum()
                  self%coefficientFactorBoostMaximum       =self%latticeFactorBoost       %maximum()
                  self%coefficientFactorBoostStellarMinimum=self%latticeFactorBoostStellar%minimum()
                  self%coefficientFactorBoostStellarMaximum=self%latticeFactorBoostStellar%maximum()
                  self%radiusScaleFreeMinimum              =self%latticeRadiusScaleFree   %minimum()
                  self%radiusScaleFreeMaximum              =self%latticeRadiusScaleFree   %maximum()
                  self%tableInitialized                    =.true.
               end if
            end if
            if (.not.self%tableInitialized) then
               ! The stored tabulation is unusable - discard whatever was read of it.
               if (allocated(self%integralPartiallyMolecularTable)) deallocate(self%integralPartiallyMolecularTable)
               self%latticeFactorBoost       =rangeLattice()
               self%latticeFactorBoostStellar=rangeLattice()
               self%latticeRadiusScaleFree   =rangeLattice()
            end if
         end if
      end if
      ! Having read the table from file (if it exists), check again to see if it is sufficient. If it is not, we must retabulate.
      if (tableIsInsufficient()) then
         ! Obtain a file lock if we don't already have one.
         if (.not.haveLock) then
            call File_Lock(self%filenameTable,fileLock,lockIsShared=.false.)
            haveLock=.true.
         end if
         ! Find the range to tabulate, pinning each axis to an absolute lattice so that the points evaluated - and therefore every
         ! value interpolated between them - depend only on which lattice points are spanned, and not on the sequence of values
         ! which happened to be requested. Each request is passed as the target and the range already tabulated is unioned in
         ! through `latticeCurrent`; folding the latter into the target instead would apply the factor-of-two margin to an
         ! already-margined bound and ratchet the range outward on every retabulation. The hard lower limits enter as
         ! `limitMinimum`. Note that each request is first clamped *up* to its hard limit - exactly as the interpolation below
         ! clamps it - rather than the limit being added to the target: adding it would make it the smallest target value, and so
         ! would drag the lower bound of every tabulation down to the limit instead of leaving it where the request puts it.
         latticeFactorBoost       =Range_Pinned(                                                                                       &
              &                                                [max(coefficientFactorBoost       ,coefficientFactorBoostTiny       )], &
              &                                                 pointsPerDecadeFactorBoost                                           , &
              &                                                 gridSchemePerDecade                                                  , &
              &                                 marginFactor  = 2.0d0                                                                , &
              &                                 anchorEvery   = anchorEveryFactorBoost                                               , &
              &                                 limitMinimum  = coefficientFactorBoostTiny                                           , &
              &                                 latticeCurrent=self%latticeFactorBoost                                                 &
              &                                )
         latticeFactorBoostStellar=Range_Pinned(                                                                                       &
              &                                                [max(coefficientFactorBoostStellar,coefficientFactorBoostStellarTiny)], &
              &                                                 pointsPerDecadeFactorBoostStellar                                    , &
              &                                                 gridSchemePerDecade                                                  , &
              &                                 marginFactor  = 2.0d0                                                                , &
              &                                 anchorEvery   = anchorEveryFactorBoostStellar                                        , &
              &                                 limitMinimum  = coefficientFactorBoostStellarTiny                                    , &
              &                                 latticeCurrent=self%latticeFactorBoostStellar                                          &
              &                                )
         latticeRadiusScaleFree   =Range_Pinned(                                                                                       &
              &                                                [radiusScaleFree]                                                     , &
              &                                                 pointsPerDecadeRadius                                                , &
              &                                                 gridSchemePerDecade                                                  , &
              &                                 marginFactor  = 2.0d0                                                                , &
              &                                 anchorEvery   = anchorEveryRadius                                                    , &
              &                                 latticeCurrent=self%latticeRadiusScaleFree                                             &
              &                                )
         countFactorBoost       =latticeFactorBoost       %count
         countFactorBoostStellar=latticeFactorBoostStellar%count
         countRadii             =latticeRadiusScaleFree   %count
         ! Record where the table already in hand sits within the extended one, so that the solutions it holds can be carried
         ! over. Every offset is found in exact integer arithmetic from the lattice indices, so no abscissa is compared.
         carryOver              =       self%tableInitialized                      &
              &                  .and.  self%latticeFactorBoost       %isDefined() &
              &                  .and.  self%latticeFactorBoostStellar%isDefined() &
              &                  .and.  self%latticeRadiusScaleFree   %isDefined()
         offsetFactorBoost       =0
         offsetFactorBoostStellar=0
         offsetRadius            =0
         countFactorBoostPrevious       =0
         countFactorBoostStellarPrevious=0
         countRadiiPrevious             =0
         if (carryOver) then
            offsetFactorBoost              =Range_Lattice_Offset(self%latticeFactorBoost       ,latticeFactorBoost       )
            offsetFactorBoostStellar       =Range_Lattice_Offset(self%latticeFactorBoostStellar,latticeFactorBoostStellar)
            offsetRadius                   =Range_Lattice_Offset(self%latticeRadiusScaleFree   ,latticeRadiusScaleFree   )
            countFactorBoostPrevious       =self%latticeFactorBoost       %count
            countFactorBoostStellarPrevious=self%latticeFactorBoostStellar%count
            countRadiiPrevious             =self%latticeRadiusScaleFree   %count
            call Move_Alloc(self%integralPartiallyMolecularTable,integralPartiallyMolecularPrevious)
         end if
         self%latticeFactorBoost       =latticeFactorBoost
         self%latticeFactorBoostStellar=latticeFactorBoostStellar
         self%latticeRadiusScaleFree   =latticeRadiusScaleFree
         ! Take the abscissae from the lattices. They must come from there, and never from an exponentiation open-coded here: the
         ! lattice evaluates them through a single, deliberately un-inlined path, so that a given lattice point is bit-identical
         ! between one tabulation and another regardless of how many points each spans.
         valuesFactorBoost       =latticeFactorBoost       %values()
         valuesFactorBoostStellar=latticeFactorBoostStellar%values()
         valuesRadiusScaleFree   =latticeRadiusScaleFree   %values()
         ! Record the limits. The offset and inverse step which the interpolation uses are not recorded: both are functions of
         ! the lattice, and are taken from it where they are used.
         self%coefficientFactorBoostMinimum       =latticeFactorBoost       %minimum()
         self%coefficientFactorBoostMaximum       =latticeFactorBoost       %maximum()
         self%coefficientFactorBoostStellarMinimum=latticeFactorBoostStellar%minimum()
         self%coefficientFactorBoostStellarMaximum=latticeFactorBoostStellar%maximum()
         self%radiusScaleFreeMinimum              =latticeRadiusScaleFree   %minimum()
         self%radiusScaleFreeMaximum              =latticeRadiusScaleFree   %maximum()
         !! Allocate the table, and carry over the block of solutions already found.
         allocate(self%integralPartiallyMolecularTable(countFactorBoost,countFactorBoostStellar,countRadii))
         self%integralPartiallyMolecularTable=0.0d0
         if (carryOver) then
            self%integralPartiallyMolecularTable(                                                                                     &
                 &                               offsetFactorBoost       +1:offsetFactorBoost       +countFactorBoostPrevious       , &
                 &                               offsetFactorBoostStellar+1:offsetFactorBoostStellar+countFactorBoostStellarPrevious, &
                 &                               offsetRadius            +1:offsetRadius            +countRadiiPrevious               &
                 &                              )=integralPartiallyMolecularPrevious
            deallocate(integralPartiallyMolecularPrevious)
         end if
         !! Populate the table.
         call displayIndent("tabulating solutions for Blitz2006 star formation rate in exponential disks",verbosityLevelWorking)
         call displayIndent("table ranges"                                                               ,verbosityLevelWorking)
         write (rangeLower,'(e12.6)') self%coefficientFactorBoostMinimum
         write (rangeUpper,'(e12.6)') self%coefficientFactorBoostMaximum
         write (countSteps,'(i6   )')      countFactorBoost
         message=rangeLower//" ≤ boost gas     ≤ "//rangeUpper//" ["//countSteps//" steps]"
         call displayMessage(message,verbosityLevelWorking)
         write (rangeLower,'(e12.6)') self%coefficientFactorBoostStellarMinimum
         write (rangeUpper,'(e12.6)') self%coefficientFactorBoostStellarMaximum
         write (countSteps,'(i6   )')      countFactorBoostStellar
         message=rangeLower//" ≤ boost stellar ≤ "//rangeUpper//" ["//countSteps//" steps]"
         call displayMessage(message,verbosityLevelWorking)
         write (rangeLower,'(e12.6)') self%radiusScaleFreeMinimum
         write (rangeUpper,'(e12.6)') self%radiusScaleFreeMaximum
         write (countSteps,'(i6   )')      countRadii
         message=rangeLower//" ≤ radius        ≤ "//rangeUpper//" ["//countSteps//" steps]"
         call displayMessage(message,verbosityLevelWorking)
         call displayUnindent("",verbosityLevelWorking)
         call displayIndent("table size",verbosityLevelWorking)
         write (tableSize,'(f8.4)') dble(sizeof(self%integralPartiallyMolecularTable))/1024.0d0**2
         message=trim(adjustl(tableSize))//" MiB"
         call displayMessage(message,verbosityLevelWorking)
         call displayUnindent("",verbosityLevelWorking)
         loopCount     =+0
         loopCountTotal=+countFactorBoost                                                            &
              &         *countFactorBoostStellar                                                     &
              &         *countRadii                                                                  &
              &         -countFactorBoostPrevious*countFactorBoostStellarPrevious*countRadiiPrevious
         !$omp parallel private(i,j,k,integrator_,radiusMinimum,radiusMaximum,isCarriedFactorBoost,isCarriedFactorBoostStellar)
         allocate(integrator_)
         integrator_   = integrator(integrand,toleranceRelative=toleranceRelative,integrationRule=GSL_Integ_Gauss61)
         do i=1,countFactorBoost
            coefficientFactorBoost_          =valuesFactorBoost       (i)
            isCarriedFactorBoost             =         carryOver                                                           &
                 &                            .and. i >           offsetFactorBoost                                        &
                 &                            .and. i <=          offsetFactorBoost+countFactorBoostPrevious
            do j=1,countFactorBoostStellar
               coefficientFactorBoostStellar_=valuesFactorBoostStellar(j)
               isCarriedFactorBoostStellar   =         isCarriedFactorBoost                                                &
                    &                         .and. j >           offsetFactorBoostStellar                                 &
                    &                         .and. j <=          offsetFactorBoostStellar+countFactorBoostStellarPrevious
               !$omp do
               do k=1,countRadii
                  ! Skip the block of solutions carried over from the tabulation already in hand - evaluating them again would
                  ! merely reproduce them, at the cost of a numerical integral apiece.
                  if     (                                           &
                       &        isCarriedFactorBoostStellar          &
                       &  .and. k >  offsetRadius                    &
                       &  .and. k <= offsetRadius+countRadiiPrevious &
                       & ) cycle
                  radiusMinimum                              =0.0d0
                  radiusMaximum                              =valuesRadiusScaleFree   (k)
                  self%integralPartiallyMolecularTable(i,j,k)=log(integrator_%integrate(radiusMinimum,radiusMaximum))
                  !$omp critical(blitz2006Tabulation)
                  call displayCounter(int(100.0d0*dble(loopCount)/dble(loopCountTotal)),loopCount==0,verbosityLevelWorking)
                  loopCount=loopCount+1
                  !$omp end critical(blitz2006Tabulation)
               end do
               !$omp end do
            end do
         end do
         deallocate(integrator_)
         !$omp end parallel
         call displayCounterClear(       verbosityLevelWorking)
         call displayUnindent    ("done",verbosityLevelWorking)
         self%tableInitialized=.true.
         ! Write the table to file.
         message='Writing Blitz2006 star formation rate tabulation to file: '//self%filenameTable
         call displayMessage(message,verbosityLevelWorking)
         call Directory_Make(File_Path(self%filenameTable))
         !$ call hdf5Access%set()
         hdf5FileScopeWrite: block
           type(hdf5File  ) :: file
           file=hdf5File(self%filenameTable,overWrite=.true.,readOnly=.false.)
           ! Record the lattices on which the three axes are built. The limits, offsets and inverse steps formerly stored
           ! alongside them are not: each is a function of the lattices, and is recovered from them when the file is read.
           call Table_Cache_Lattice_Write(file,'coefficientFactorBoost'       ,self%latticeFactorBoost       )
           call Table_Cache_Lattice_Write(file,'coefficientFactorBoostStellar',self%latticeFactorBoostStellar)
           call Table_Cache_Lattice_Write(file,'radiusScaleFree'              ,self%latticeRadiusScaleFree   )
           call file%writeDataset(self%integralPartiallyMolecularTable,'integral')
         end block hdf5FileScopeWrite
         !$ call hdf5Access%unset()
      end if
      if (haveLock) then
         call File_Unlock(fileLock)
         haveLock=.false.
      end if
      ! Interpolate in table.
      integralAnalyticPartiallyMolecularGeneric=0.0d0
      coefficientFactorBoost_            =max(coefficientFactorBoost       ,coefficientFactorBoostTiny       )
      coefficientFactorBoostStellar_     =max(coefficientFactorBoostStellar,coefficientFactorBoostStellarTiny)
      ! Find the position of each value along its axis. This is done as the coordinate of the value on the *absolute* lattice -
      ! a quantity which depends only on the value and on the density of lattice points, never on which part of the lattice this
      ! particular table happens to span - split there into the index of the lattice point below it and the fraction of the
      ! interval above that, with the index of the first tabulated point subtracted only afterwards, in exact integer
      ! arithmetic.
      !
      ! The order matters. Forming the position relative to the first tabulated point *first* is exact in the subtraction, but
      ! the fractional part is then extracted from a number whose magnitude is the index within the table, and so is rounded on
      ! a grid which coarsens as the table grows: extending the table would perturb every interpolated value in its last bits,
      ! which is exactly the dependence on the sequence of requests that pinning the ranges exists to remove.
      coordinateFactorBoost       =log10(coefficientFactorBoost_       )*dble(pointsPerDecadeFactorBoost       )
      coordinateFactorBoostStellar=log10(coefficientFactorBoostStellar_)*dble(pointsPerDecadeFactorBoostStellar)
      coordinateRadiusScaleFree   =log10(radiusScaleFree               )*dble(pointsPerDecadeRadius            )
      i                           =floor(coordinateFactorBoost       )
      j                           =floor(coordinateFactorBoostStellar)
      k                           =floor(coordinateRadiusScaleFree   )
      hi                          =coordinateFactorBoost       -dble(i)
      hj                          =coordinateFactorBoostStellar-dble(j)
      hk                          =coordinateRadiusScaleFree   -dble(k)
      i                           =i-self%latticeFactorBoost       %indexMinimum+1
      j                           =j-self%latticeFactorBoostStellar%indexMinimum+1
      k                           =k-self%latticeRadiusScaleFree   %indexMinimum+1
      ! Confine the indices to the table. Each is guaranteed to lie within it by the test which decided that the table need not
      ! be extended, but only up to rounding: a value which sits within an ulp of a bound can index one point beyond it, and the
      ! interpolation below reads the *next* point along every axis. Where an index is moved the interpolating factor is moved
      ! with it, to the end of the interval nearest the value, so that the interpolation returns the value at the edge of the
      ! table. Moving the index alone - as this formerly did - shifts the result by a whole interval at the bounds.
      if (i <                                                  1) hi=0.0d0
      if (i > size(self%integralPartiallyMolecularTable,dim=1)-1) hi=1.0d0
      if (j <                                                  1) hj=0.0d0
      if (j > size(self%integralPartiallyMolecularTable,dim=2)-1) hj=1.0d0
      if (k <                                                  1) hk=0.0d0
      if (k > size(self%integralPartiallyMolecularTable,dim=3)-1) hk=1.0d0
      i=min(max(i,1),size(self%integralPartiallyMolecularTable,dim=1)-1)
      j=min(max(j,1),size(self%integralPartiallyMolecularTable,dim=2)-1)
      k=min(max(k,1),size(self%integralPartiallyMolecularTable,dim=3)-1)
      integral=0.0d0
      do ii=0,1
         if (ii == 0) then
            hhi=+1.0d0-hi
         else
            hhi=      +hi
         end if
         do jj=0,1
            if (jj == 0) then
               hhj=+1.0d0-hj
            else
               hhj=      +hj
            end if
            do kk=0,1
               if (kk == 0) then
                  hhk=+1.0d0-hk
               else
                  hhk=      +hk
               end if
               integral=+integral                                             &
                    &   +self%integralPartiallyMolecularTable(i+ii,j+jj,k+kk) &
                    &   *                                     hhi             &
                    &   *                                          hhj        &
                    &   *                                               hhk
            end do
         end do
      end do
      integralAnalyticPartiallyMolecularGeneric=+integralAnalyticPartiallyMolecularGeneric &
           &                                    +exp(integral)                             &
           &                                    *coefficientNormalization
      return
    end function integralAnalyticPartiallyMolecularGeneric
    
    logical function tableIsInsufficient()
      !!{RST
      Determine if the current table is insufficient for our purposes.
      !!}
      implicit none
      
      tableIsInsufficient=                                    .not. self%tableInitialized                     &
           &              .or.                                                                                &
           &               (                                                                                  &
           &                 coefficientFactorBoost            <    self%coefficientFactorBoostMinimum        &
           &                .and.                                                                             &
           &                 coefficientFactorBoostTiny        <    self%coefficientFactorBoostMinimum        &
           &               )                                                                                  &
           &              .or.                                                                                &
           &                 coefficientFactorBoost            >    self%coefficientFactorBoostMaximum        &
           &              .or.                                                                                &
           &               (                                                                                  &
           &                 coefficientFactorBoostStellar     <    self%coefficientFactorBoostStellarMinimum &
           &                .and.                                                                             &
           &                 coefficientFactorBoostStellarTiny <    self%coefficientFactorBoostStellarMinimum &
           &               )                                                                                  &
           &              .or.                                                                                &
           &                 coefficientFactorBoostStellar     >    self%coefficientFactorBoostStellarMaximum &
           &              .or.                                                                                &
           &                 radiusScaleFree                   <    self%radiusScaleFreeMinimum               &
           &              .or.                                                                                &
           &                 radiusScaleFree                   >    self%radiusScaleFreeMaximum
      return
    end function tableIsInsufficient
    
    double precision function integrand(radiusScaleFree)
      !!{RST
      Integrand for the partially molecular case.
      !!}
      implicit none
      double precision, intent(in   ) :: radiusScaleFree

      integrand=+exp(-(1.0d0+2.0d0*self%pressureExponent)*radiusScaleFree) &
           &    *(                                                         &
           &      +1.0d0                                                   &
           &      +coefficientFactorBoostStellar_                          &
           &      *exp(0.5d0*radiusScaleFree)                              &
           &     )**self%pressureExponent                                  &
           &    *(                                                         &
           &      +1.0d0                                                   &
           &      +(                                                       &
           &        +coefficientFactorBoost_                               &
           &        *exp(-radiusScaleFree)                                 &
           &       )**self%surfaceDensityExponent                          &
           &     )                                                         &
           &    *radiusScaleFree
      return
    end function integrand
    
  end function blitz2006IntegralPartiallyMolecular
  
  double precision function blitz2006CriticalDensityRoot(radius)
    !!{RST
    Root function used in finding the radius in a disk where the pressure ratio exceeds the critical ratio.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    blitz2006CriticalDensityRoot=self_%pressureRatio(node_,radius)-1.0d0
    return
  end function blitz2006CriticalDensityRoot

  double precision function blitz2006PressureRatio(self,node,radius,surfaceDensityGas) result(pressureRatio)
    !!{RST
    Root function used in finding the radius in a disk where the pressure ratio exceeds the critical ratio.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Galactic_Structure_Options      , only : componentTypeDisk             , coordinateSystemCylindrical, massTypeGaseous, massTypeStellar
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Coordinates, only : coordinateCylindrical, assignment(=)
    implicit none
    class           (starFormationRateSurfaceDensityDisksBlitz2006), intent(inout)           :: self
    type            (treeNode                                     ), intent(inout)           :: node
    double precision                                               , intent(in   )           :: radius
    double precision                                               , intent(  out), optional :: surfaceDensityGas
    class           (massDistributionClass                        ), pointer                 :: massDistribution_
    type            (coordinateCylindrical                        )                          :: coordinates
    double precision                                                                         :: surfaceDensityGas_, surfaceDensityStellar, &
         &                                                                                      factorBoostStellar

    ! Get gas surface density.
    coordinates        =  [radius,0.0d0,0.0d0]
    massDistribution_  => node             %massDistribution(componentType=componentTypeDisk,massType=massTypeGaseous)
    surfaceDensityGas_ =  massDistribution_%surfaceDensity  (              coordinates                               )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    if (present(surfaceDensityGas)) surfaceDensityGas=surfaceDensityGas_
    ! If the radius is less than the critical radius the pressure radius is above 1 by definition, so simply pin it to that value.
    if (radius <= self%radiusCritical) then
       pressureRatio=1.0d0
    else
       ! Compute the pressure ratio that Blitz & Rosolowsky (2006) use to compute the molecular fraction.    
       !! We first compute the pressure ratio ignoring the boost from the stellar mass. If this already exceeds the characteristic
       !! limit then we will not need to compute the stellar contribution.
       pressureRatio=+0.5d0                             &
            &        *Pi                                &
            &        *gravitationalConstant_internal    &
            &        *surfaceDensityGas_            **2 &
            &        /self%pressureCharacteristic
       if (pressureRatio > 0.0d0 .and. pressureRatio < 1.0d0) then
          ! Compute the stellar boost factor.
          massDistribution_     =>  node             %massDistribution(componentType=componentTypeDisk,massType=massTypeStellar)
          surfaceDensityStellar =  +massDistribution_%surfaceDensity  (              coordinates                               )
          !![
	  <objectDestructor name="massDistribution_"/>
          !!]
          factorBoostStellar   =+1.0d0                                &
               &                +self%velocityDispersionDiskGas       &
               &                /surfaceDensityGas_                   &
               &                *sqrt(                                &
               &                      +surfaceDensityStellar          &
               &                      /Pi                             &
               &                      /gravitationalConstant_internal &
               &                      /self%heightToRadialScaleDisk   &
               &                      /self%radiusDisk                &
               &                     ) 
          pressureRatio        =+pressureRatio                        &
               &                *factorBoostStellar
       end if
    end if
    return
  end function blitz2006PressureRatio

