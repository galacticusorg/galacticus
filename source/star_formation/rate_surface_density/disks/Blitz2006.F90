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

  !![
  <starFormationRateSurfaceDensityDisks name="starFormationRateSurfaceDensityDisksBlitz2006" docformat="rst">
   <description>
   A star formation rate surface density class which assumes that the star formation rate is given by :cite:p:`blitz_role_2006`:

   .. math::

      \dot{\Sigma}_\star(R) = \nu_\mathrm{SF}(R) \Sigma_\mathrm{H_2, disk}(R),

   where :math:`\nu_\mathrm{SF}` is a frequency given by

   .. math::

      \nu_\mathrm{SF}(R) = \nu_\mathrm{SF,0} \left[ 1 + \left({\Sigma_\mathrm{H}\over \Sigma_0}\right)^q \right],

   where :math:`q=`\ ``[surfaceDensityExponent]`` and :math:`\Sigma_0=`\ ``[surfaceDensityCritical]`` are parameters, and :math:`\Sigma_\mathrm{H}` is the surface density of hydrogen, atomic and molecular together - the density at which :cite:t:`bigiel_star_2008`, the source of the default :math:`\Sigma_0`, find the star formation law to steepen. The ratio of molecular to atomic hydrogen is :math:`R_\mathrm{mol} = \Sigma_\mathrm{H_2}/\Sigma_\mathrm{HI} = (P_\mathrm{ext}/P_0)^\alpha`, where :math:`\alpha=`\ ``[pressureExponent]`` and :math:`P_0=`\ ``[pressureCharacteristic]`` are parameters, so that the surface density of molecular gas is :math:`\Sigma_\mathrm{H_2} = f_\mathrm{H_2} \Sigma_\mathrm{H}`, with :math:`f_\mathrm{H_2} = R_\mathrm{mol}/(1+R_\mathrm{mol})` the molecular fraction of eqn. (21) of :cite:t:`blitz_role_2006`. The hydrostatic pressure in the disk plane assuming locally isothermal gas and stellar components is given by

   .. math::

      P_\mathrm{ext} \approx {\pi\over 2} \G \Sigma_\mathrm{gas} \left[ \Sigma_\mathrm{gas} + \left({\sigma_\mathrm{gas}\over
      \sigma_\star}\right)\Sigma_\star\right]

   where we assume that the velocity dispersion in the gas is fixed at :math:`\sigma_\mathrm{gas}=`\ ``[velocityDispersionDiskGas]`` and, assuming :math:`\Sigma_\star \gg \Sigma_\mathrm{gas}`, we can write the stellar velocity dispersion in terms of the disk scale height, :math:`h_\star`, as

   .. math::

      \sigma_\star = \sqrt{\pi \G h_\star \Sigma_\star}

   where we assume :math:`h_\star/R_\mathrm{disk}=`\ ``[heightToRadialScaleDisk]``.
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
          &                                                                     isExponentialDisk
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
   contains
     !![
     <methods docformat="rst">
       <method description="Reset memoized calculations." method="calculationReset"/>
       <method description="Compute various factors."     method="computeFactors"  />
       <method description="Compute the pressure ratio."  method="pressureRatio"   />
     </methods>
     !!]
     final     ::                     blitz2006Destructor
     procedure :: autoHook         => blitz2006AutoHook
     procedure :: calculationReset => blitz2006CalculationReset
     procedure :: rate             => blitz2006Rate
     procedure :: computeFactors   => blitz2006ComputeFactors
     procedure :: unchanged        => blitz2006Unchanged
     procedure :: intervals        => blitz2006Intervals
     procedure :: pressureRatio    => blitz2006PressureRatio
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
    logical                                                                        :: assumeMonotonicSurfaceDensity

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
      <defaultValue>3.5d4</defaultValue>
      <description>
      The characteristic pressure (given as :math:`P_0/k_\mathrm{B}` in units of K cm\ :math:`^{-3}`) in the scaling relation of molecular hydrogen fraction with disk pressure in the :cite:t:`blitz_role_2006` star formation timescale calculation. The default is the mean value, :math:`P_0/k_\mathrm{B} = 3.5 \times 10^4` K cm\ :math:`^{-3}`, of Table 2 (and eqn. 12) of :cite:t:`blitz_role_2006`, which corresponds to the default :math:`\alpha=0.92` of ``[pressureExponent]``.
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
    !!]
    self=starFormationRateSurfaceDensityDisksBlitz2006(velocityDispersionDiskGas,heightToRadialScaleDisk,surfaceDensityCritical,surfaceDensityExponent,starFormationFrequencyNormalization,pressureCharacteristic,pressureExponent,assumeMonotonicSurfaceDensity)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function blitz2006ConstructorParameters

  function blitz2006ConstructorInternal(velocityDispersionDiskGas,heightToRadialScaleDisk,surfaceDensityCritical,surfaceDensityExponent,starFormationFrequencyNormalization,pressureCharacteristic,pressureExponent,assumeMonotonicSurfaceDensity) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`starFormationRateSurfaceDensityDisksBlitz2006` star formation surface density rate in disks class.
    !!}
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : massSolar                , megaParsec
    use :: Numerical_Constants_Physical    , only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes    , only : giga                     , hecto                        , kilo                         , mega
    use :: Root_Finder                     , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    implicit none
    type            (starFormationRateSurfaceDensityDisksBlitz2006)                :: self
    double precision                                               , intent(in   ) :: velocityDispersionDiskGas          , heightToRadialScaleDisk, &
         &                                                                            surfaceDensityCritical             , surfaceDensityExponent , &
         &                                                                            starFormationFrequencyNormalization, pressureCharacteristic , &
         &                                                                            pressureExponent
    logical                                                        , intent(in   ) :: assumeMonotonicSurfaceDensity
    !![
    <constructorAssign variables="velocityDispersionDiskGas, heightToRadialScaleDisk, surfaceDensityCritical, surfaceDensityExponent, starFormationFrequencyNormalization, pressureCharacteristic, pressureExponent, assumeMonotonicSurfaceDensity"/>
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
    self%massGasPrevious             =-huge(0.0d0)
    self%massStellarPrevious         =-huge(0.0d0)
    self%radiusDiskPrevious          =-huge(0.0d0)
    self%hydrogenMassFractionPrevious=-huge(0.0d0)
    self%radiusCriticalPrevious      =-huge(0.0d0)
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
         &                                                                            surfaceDensityGas, factorBoost  , &
         &                                                                            ratioMolecular

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
    ! Compute the ratio of molecular to atomic hydrogen, R_mol=Σ_H₂/Σ_HI=(P_ext/P₀)^α (eqn. 11 of Blitz & Rosolowsky 2006), and
    ! from it the molecular fraction, f_H₂=R_mol/(1+R_mol) (their eqn. 21). The fast exponentiator is tabulated over the range
    ! [0,1] only, so is used only where the pressure ratio falls within that range.
    if (pressureRatio >= 1.0d0) then
       ratioMolecular   =+                                             pressureRatio **self%pressureExponent
    else
       ratioMolecular   =+self%pressureRatioExponentiator%exponentiate(pressureRatio)
    end if
    molecularFraction   =+       ratioMolecular  &
         &               /(1.0d0+ratioMolecular)
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
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Galactic_Structure_Options, only : componentTypeDisk    , massTypeGaseous, massTypeStellar
    implicit none
    class           (starFormationRateSurfaceDensityDisksBlitz2006), intent(inout), target                      :: self
    double precision                                                              , allocatable, dimension(:,:) :: blitz2006Intervals
    type            (treeNode                                     ), intent(inout), target                      :: node
    double precision                                               , intent(in   )                              :: radiusInner                             , radiusOuter
    logical                                                        , intent(inout), allocatable, dimension(  :) :: intervalIsAnalytic
    double precision                                               , intent(inout), allocatable, dimension(  :) :: integralsAnalytic
    class           (massDistributionClass                        ), pointer                                    :: massDistributionGaseous                 , massDistributionStellar
    double precision                                               , parameter                                  :: factorBoostStellarCoefficientTiny=1.0d-6
    double precision                                                                                            :: rootValueInner                          , rootValueOuter               , &
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
             intervalIsAnalytic =.false.
             blitz2006Intervals =reshape([radiusInner,radiusOuter],[2,1])
             self%radiusCritical=-huge(0.0d0)
          else
             ! Test the surface density at the outer radius.
             rootValueOuter=blitz2006CriticalDensityRoot(radiusOuter)
             if (rootValueOuter >= 0.0d0) then
                ! Entire disk is above the pressure threshold so use a single interval.
                allocate(blitz2006Intervals(2,1))
                allocate(intervalIsAnalytic(  1))
                intervalIsAnalytic =.false.
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
                intervalIsAnalytic=.false.
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
  end function blitz2006Intervals

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
    use :: Coordinates                     , only : coordinateCylindrical         , assignment(=)
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
    ! Compute the pressure ratio that Blitz & Rosolowsky (2006) use to compute the molecular fraction. The molecular fraction,
    ! f_H₂=R_mol/(1+R_mol), is not capped at unity, so the ratio must be computed in full at every radius: neither pinning the
    ! ratio to unity inside the critical radius, nor omitting the stellar boost where the gas-only ratio already exceeds unity,
    ! is valid for that form.
    pressureRatio=+0.5d0                             &
         &        *Pi                                &
         &        *gravitationalConstant_internal    &
         &        *surfaceDensityGas_            **2 &
         &        /self%pressureCharacteristic
    if (pressureRatio > 0.0d0) then
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
    return
  end function blitz2006PressureRatio

