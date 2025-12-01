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

  !!{
  Implementation of the \cite{blitz_role_2006} star formation rate surface density law for galactic disks.
  !!}

  use :: Kind_Numbers       , only : kind_int8
  use :: Root_Finder        , only : rootFinder
  use :: Math_Exponentiation, only : fastExponentiator
  
  !![
  <starFormationRateSurfaceDensityDisks name="starFormationRateSurfaceDensityDisksBlitz2006">
   <description>
    A star formation rate surface density class which assumes that the star formation rate is given by \citep{blitz_role_2006}:
    \begin{equation}
     \dot{\Sigma}_\star(R) = \nu_\mathrm{SF}(R) \Sigma_\mathrm{H_2, disk}(R),
    \end{equation}
    where $\nu_\mathrm{SF}$ is a frequency given by
    \begin{equation}
     \nu_\mathrm{SF}(R) = \nu_\mathrm{SF,0} \left[ 1 + \left({\Sigma_\mathrm{HI}\over \Sigma_0}\right)^q \right],
    \end{equation}
    where $q=${\normalfont \ttfamily [surfaceDensityExponent]} and $\Sigma_0=${\normalfont \ttfamily [surfaceDensityCritical]}
    are parameters, the surface density of molecular gas $\Sigma_\mathrm{H_2} = (P_\mathrm{ext}/P_0)^\alpha
    \Sigma_\mathrm{HI}$, where $\alpha=${\normalfont \ttfamily [pressureExponent]} and $P_0=${\normalfont \ttfamily
    [pressureCharacteristic]} are parameters, and the hydrostatic pressure in the disk plane assuming locally isothermal gas
    and stellar components is given by
    \begin{equation}
     P_\mathrm{ext} \approx {\pi\over 2} \G \Sigma_\mathrm{gas} \left[ \Sigma_\mathrm{gas} + \left({\sigma_\mathrm{gas}\over
     \sigma_\star}\right)\Sigma_\star\right]
    \end{equation}
    where we assume that the velocity dispersion in the gas is fixed at $\sigma_\mathrm{gas}=${\normalfont \ttfamily
    [velocityDispersionDiskGas]} and, assuming $\Sigma_\star \gg \Sigma_\mathrm{gas}$, we can write the stellar velocity
    dispersion in terms of the disk scale height, $h_\star$, as
    \begin{equation}
     \sigma_\star = \sqrt{\pi \G h_\star \Sigma_\star}
    \end{equation}
    where we assume $h_\star/R_\mathrm{disk}=${\normalfont \ttfamily [heightToRadialScaleDiskBlitzRosolowsky]}.
   </description>
  </starFormationRateSurfaceDensityDisks>
  !!]
  type, extends(starFormationRateSurfaceDensityDisksClass) :: starFormationRateSurfaceDensityDisksBlitz2006
     !!{
     Implementation of the \cite{blitz_role_2006} star formation rate surface density law for galactic disks.
     !!}
     private
     integer         (kind_int8             )                                :: lastUniqueID
     logical                                                                 :: factorsComputed                                       , assumeMonotonicSurfaceDensity                      , &
          &                                                                     isExponentialDisk                                     , useTabulation
     double precision                                                        :: heightToRadialScaleDisk                               , pressureCharacteristic                             , &
          &                                                                     pressureExponent                                      , starFormationFrequencyNormalization                , &
          &                                                                     surfaceDensityCritical                                , surfaceDensityExponent                             , &
          &                                                                     velocityDispersionDiskGas                             , radiusDisk                                         , &
          &                                                                     massGas                                               , hydrogenMassFraction                               , &
          &                                                                     massStellar                                           , massGasPrevious                                    , &
          &                                                                     massStellarPrevious                                   , hydrogenMassFractionPrevious                       , &
          &                                                                     radiusDiskPrevious                                    , radiusCritical                                     , &
          &                                                                     radiusCriticalPrevious                                , factorBoostStellarCoefficient                      , &
          &                                                                     pressureRatioCoefficient
     type            (rootFinder            )                                :: finder
     type            (fastExponentiator     )                                :: pressureRatioExponentiator
     double precision                        , allocatable, dimension(:,:,:) :: integralPartiallyMolecularTable
     logical                                                                 :: tableInitialized
     double precision                                                        :: coefficientFactorBoostMinimum                         , coefficientFactorBoostMaximum                      , &
          &                                                                     coefficientFactorBoostStellarMinimum                  , coefficientFactorBoostStellarMaximum               , &
          &                                                                     radiusScaleFreeMinimum                                , radiusScaleFreeMaximum                             , &
          &                                                                     coefficientFactorBoostLogarithmicStepInverse          , coefficientFactorBoostStellarLogarithmicStepInverse, &
          &                                                                     radiusScaleFreeLogarithmicStepInverse                 , radiusScaleFreeLogarithmicOffset                   , &
          &                                                                     coefficientFactorBoostLogarithmicOffset               , coefficientFactorBoostStellarLogarithmicOffset
     type            (varying_string        )                                :: filenameTable
   contains
     !![
     <methods>
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
     !!{
     Constructors for the \refClass{starFormationRateSurfaceDensityDisksBlitz2006} star formation surface density rate in disks class.
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
    !!{
    Constructor for the \refClass{starFormationRateSurfaceDensityDisksBlitz2006} star formation surface density rate in disks class which takes a parameter set as input.
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
    <inputParameter>
      <name>velocityDispersionDiskGas</name>
      <defaultSource>\citep{leroy_star_2008}</defaultSource>
      <defaultValue>10.0d0</defaultValue>
      <description>The velocity dispersion of gas in disks.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>heightToRadialScaleDisk</name>
      <defaultSource>\citep{kregel_flattening_2002}</defaultSource>
      <defaultValue>0.137d0</defaultValue>
      <description>The ratio of scale height to scale radius for disks in the \cite{blitz_role_2006} star formation timescale calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>surfaceDensityCritical</name>
      <defaultSource>\citep{bigiel_star_2008}</defaultSource>
      <defaultValue>200.0d0</defaultValue>
      <description>The surface density (in units of $M_\odot$ pc$^{-2}$) in the \cite{blitz_role_2006} star formation timescale calculation at which low-density truncation begins.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>surfaceDensityExponent</name>
      <defaultSource>\citep{bigiel_star_2008}</defaultSource>
      <defaultValue>0.4d0</defaultValue>
      <description>The exponent for surface density in the \cite{blitz_role_2006} star formation timescale calculation at in the high density regime.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>starFormationFrequencyNormalization</name>
      <defaultSource>\citep{leroy_star_2008}</defaultSource>
      <defaultValue>5.25d-10</defaultValue>
      <description>The star formation frequency (in the low-density limit and in units of yr$^{-1}$) in the \cite{blitz_role_2006} star formation timescale calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>pressureCharacteristic</name>
      <defaultSource>\citep{blitz_role_2006}</defaultSource>
      <defaultValue>4.54d0</defaultValue>
      <description>The characteristic pressure (given as $P_0/k_\mathrm{B}$ in units of K cm$^{-3}$) in the scaling relation of molecular hydrogen fraction with disk pressure in the \cite{blitz_role_2006} star formation timescale calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>pressureExponent</name>
      <defaultSource>\citep{blitz_role_2006}</defaultSource>
      <defaultValue>0.92d0</defaultValue>
      <description>The exponent in the scaling relation of molecular hydrogen fraction with disk pressure in the \cite{blitz_role_2006} star formation timescale calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>assumeMonotonicSurfaceDensity</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, assume that the surface density in disks is always monotonically decreasing.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>useTabulation</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, then use tabulated solutions to the integrated star formation rate.</description>
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
    !!{
    Internal constructor for the \refClass{starFormationRateSurfaceDensityDisksBlitz2006} star formation surface density rate from disks class.
    !!}
    use :: Error                           , only : Error_Report
    use :: Input_Paths                     , only : inputPath                , pathTypeDataDynamic
    use :: Numerical_Constants_Astronomical, only : massSolar                , megaParsec
    use :: Numerical_Constants_Physical    , only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes    , only : giga                     , hecto                        , kilo                         , mega
    use :: Root_Finder                     , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    use :: Hashes_Cryptographic            , only : Hash_MD5
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
    descriptorString=""
    write (parameterLabel,'(e17.10)') self%pressureExponent
    descriptorString=descriptorString//parameterLabel//" "
    write (parameterLabel,'(e17.10)') self%surfaceDensityExponent
    descriptorString=descriptorString//parameterLabel//" "
    self%filenameTable                       =     inputPath  (pathTypeDataDynamic)// &
         &                                    'starFormation/'                     // &
         &                                    self%objectType (                   )// &
         &                                    '_'                                  // &
         &                                    Hash_MD5        (descriptorString   )// &
         &                                    '.hdf5'
    return
  end function blitz2006ConstructorInternal

  subroutine blitz2006AutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self

    call calculationResetEvent%attach(self,blitz2006CalculationReset,openMPThreadBindingAllLevels,label='starFormationRateSurfaceDensityDisksBlitz2006')
    return
  end subroutine blitz2006AutoHook

  subroutine blitz2006Destructor(self)
    !!{
    Destructor for the blitz2006 cooling radius class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self

    if (calculationResetEvent%isAttached(self,blitz2006CalculationReset)) call calculationResetEvent%detach(self,blitz2006CalculationReset)
    return
  end subroutine blitz2006Destructor

  subroutine blitz2006CalculationReset(self,node,uniqueID)
    !!{
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
    !!{
    Returns the star formation rate surface density (in $M_\odot$ Gyr$^{-1}$ Mpc$^{-2}$) for star formation
    in the galactic disk of {\normalfont \ttfamily node}. The disk is assumed to obey the
    \cite{blitz_role_2006} star formation rule.
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
    ! Return zero rate for non-positive radius or mass.
    if (self%massGas <= 0.0d0 .or. self%massStellar < 0.0d0 .or. self%radiusDisk <= 0.0d0) then
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
    !!{
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
    !!{
    Compute various factors for the {\normalfont \ttfamily blitz2006} star formation rate surface density calculation.
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
       if (self%massGas > 0.0d0) then
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
    !!{
    Returns intervals to use for integrating the \cite{krumholz_star_2009} star formation rate over a galactic disk.
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
       ! Set zero intervals for non-positive radius or mass.
       if (self%massGas <= 0.0d0 .or. self%massStellar < 0.0d0 .or. self%radiusDisk <= 0.0d0) then
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
      !!{
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
    !!{
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
      !!{
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
    !!{
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
      !!{
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
      !!{
      Analytic solution to the improper integral of the star formation rate surface density over an exponential disk for the case
      of large stellar boost factor.
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
      !!{
      Analytic solution to the improper integral of the star formation rate surface density over an exponential disk for the case
      of small radii. Uses a series solution.
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
      !!{
      Analytic solution to the improper integral of the star formation rate surface density over an exponential disk for the general case.
      !!}
      use :: Display              , only : displayCounter , displayCounterClear  , displayIndent, displayMessage, &
           &                               displayUnindent, verbosityLevelWorking
      use :: Numerical_Integration, only : integrator     , GSL_Integ_Gauss61
      use :: HDF5_Access          , only : hdf5Access
      use :: IO_HDF5              , only : hdf5Object
      use :: File_Utilities       , only : File_Exists    , File_Lock            , File_Unlock  , lockDescriptor, &
           &                               Directory_Make , File_Path
      implicit none
      double precision                , intent(in   ) :: radiusScaleFree
      double precision                , parameter     :: toleranceRelative                    =1.0d-6
      integer                         , parameter     :: pointsPerDecadeFactorBoost           =30    , pointsPerDecadeFactorBoostStellar           =30    , &
           &                                             pointsPerDecadeRadius                =30
      integer                                         :: countFactorBoost                            , countFactorBoostStellar                            , &
           &                                             countRadii                                  , i                                                  , &
           &                                             j                                           , k                                                  , &
           &                                             ii                                          , jj                                                 , &
           &                                             kk                                          , loopCount                                         ,  &
           &                                             loopCountTotal
      double precision                                :: radiusScaleFreeLogarithmicStep              , integral                                           , &
           &                                             coefficientFactorBoostLogarithmicStep       , coefficientFactorBoostStellarLogarithmicStep       , &
           &                                             hi                                          , hj                                                 , &
           &                                             hk                                          , hhi                                                , &
           &                                             hhj                                         , hhk                                                , &
           &                                             radiusMinimum                               , radiusMaximum
      character       (len= 8        )                :: tableSize
      character       (len= 8        )                :: countSteps
      character       (len=12        )                :: rangeLower                                  , rangeUpper
      type            (integrator    ), allocatable   :: integrator_
      type            (varying_string), save          :: message
      type            (hdf5Object    ), save          :: file
      type            (lockDescriptor), save          :: fileLock
      logical                                         :: haveLock
      !$omp threadprivate(message,file,fileLock)
      
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
            file=hdf5Object(char(self%filenameTable),readOnly=.true.)
            call file%readAttribute('coefficientFactorBoostMinimum'                      ,self%coefficientFactorBoostMinimum                      )
            call file%readAttribute('coefficientFactorBoostMaximum'                      ,self%coefficientFactorBoostMaximum                      )
            call file%readAttribute('coefficientFactorBoostStellarMinimum'               ,self%coefficientFactorBoostStellarMinimum               )
            call file%readAttribute('coefficientFactorBoostStellarMaximum'               ,self%coefficientFactorBoostStellarMaximum               )
            call file%readAttribute('radiusScaleFreeMinimum'                             ,self%radiusScaleFreeMinimum                             )
            call file%readAttribute('radiusScaleFreeMaximum'                             ,self%radiusScaleFreeMaximum                             )
            call file%readAttribute('coefficientFactorBoostLogarithmicOffset'            ,self%coefficientFactorBoostLogarithmicOffset            )
            call file%readAttribute('coefficientFactorBoostStellarLogarithmicOffset'     ,self%coefficientFactorBoostStellarLogarithmicOffset     )
            call file%readAttribute('radiusScaleFreeLogarithmicOffset'                   ,self%radiusScaleFreeLogarithmicOffset                   )
            call file%readAttribute('coefficientFactorBoostLogarithmicStepInverse'       ,self%coefficientFactorBoostLogarithmicStepInverse       )
            call file%readAttribute('coefficientFactorBoostStellarLogarithmicStepInverse',self%coefficientFactorBoostStellarLogarithmicStepInverse)
            call file%readAttribute('radiusScaleFreeLogarithmicStepInverse'              ,self%radiusScaleFreeLogarithmicStepInverse              )
            call file%readDataset  ('integral'                                           ,self%integralPartiallyMolecularTable                    )
            !$ call hdf5Access%unset()
            self%tableInitialized=.true.
         end if
      end if
      ! Having read the table from file (if it exists), check again to see if it is sufficient. If it is not, we must retabulate.
      if (tableIsInsufficient()) then
         ! Obtain a file lock if we don't already have one.
         if (.not.haveLock) then
            call File_Lock(char(self%filenameTable),fileLock,lockIsShared=.false.)
            haveLock=.true.
         end if
         ! Find range encompassing the existing table and the new point (with some buffer).
         self%coefficientFactorBoostMinimum       =max(min(self%coefficientFactorBoostMinimum       ,coefficientFactorBoost       /2.0d0),coefficientFactorBoostTiny             )
         self%coefficientFactorBoostMaximum       =    max(self%coefficientFactorBoostMaximum       ,coefficientFactorBoost       *2.0d0 ,coefficientFactorBoostTiny       *2.0d0)
         self%coefficientFactorBoostStellarMinimum=max(min(self%coefficientFactorBoostStellarMinimum,coefficientFactorBoostStellar/2.0d0),coefficientFactorBoostStellarTiny      )
         self%coefficientFactorBoostStellarMaximum=    max(self%coefficientFactorBoostStellarMaximum,coefficientFactorBoostStellar*2.0d0 ,coefficientFactorBoostStellarTiny*2.0d0)
         self%radiusScaleFreeMinimum              =    min(self%radiusScaleFreeMinimum              ,radiusScaleFree              /2.0d0                                         )
         self%radiusScaleFreeMaximum              =    max(self%radiusScaleFreeMaximum              ,radiusScaleFree              *2.0d0                                         )
         ! Create the table.
         !! Number of points in each axis.
         countFactorBoost       =int(log10(self%coefficientFactorBoostMaximum       /self%coefficientFactorBoostMinimum       )*dble(pointsPerDecadeFactorBoost       ))+1
         countFactorBoostStellar=int(log10(self%coefficientFactorBoostStellarMaximum/self%coefficientFactorBoostStellarMinimum)*dble(pointsPerDecadeFactorBoostStellar))+1
         countRadii             =int(log10(self%radiusScaleFreeMaximum              /self%radiusScaleFreeMinimum              )*dble(pointsPerDecadeRadius            ))+1
         !! Step sizes.
         coefficientFactorBoostLogarithmicStep       =log10(self%coefficientFactorBoostMaximum       /self%coefficientFactorBoostMinimum       )/dble(countFactorBoost       -1)
         coefficientFactorBoostStellarLogarithmicStep=log10(self%coefficientFactorBoostStellarMaximum/self%coefficientFactorBoostStellarMinimum)/dble(countFactorBoostStellar-1)
         radiusScaleFreeLogarithmicStep              =log10(self%radiusScaleFreeMaximum              /self%radiusScaleFreeMinimum              )/dble(countRadii             -1)
         !! Inverse step sizes.
         self%coefficientFactorBoostLogarithmicOffset            =log10(self%coefficientFactorBoostMinimum       )
         self%coefficientFactorBoostStellarLogarithmicOffset     =log10(self%coefficientFactorBoostStellarMinimum)
         self%radiusScaleFreeLogarithmicOffset                   =log10(self%radiusScaleFreeMinimum              )
         self%coefficientFactorBoostLogarithmicStepInverse       =1.0d0/coefficientFactorBoostLogarithmicStep
         self%coefficientFactorBoostStellarLogarithmicStepInverse=1.0d0/coefficientFactorBoostStellarLogarithmicStep
         self%radiusScaleFreeLogarithmicStepInverse              =1.0d0/radiusScaleFreeLogarithmicStep
         !! Allocate the table.
         if (self%tableInitialized) deallocate(self%integralPartiallyMolecularTable)
         allocate(self%integralPartiallyMolecularTable(countFactorBoost,countFactorBoostStellar,countRadii))
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
         loopCountTotal=+countFactorBoost        &
              &         *countFactorBoostStellar &
              &         *countRadii
         !$omp parallel private(i,j,k,integrator_,radiusMinimum,radiusMaximum)
         allocate(integrator_)
         integrator_   = integrator(integrand,toleranceRelative=toleranceRelative,integrationRule=GSL_Integ_Gauss61)
         do i=1,countFactorBoost
            coefficientFactorBoost_                          =10.0d0**(dble(i-1)*coefficientFactorBoostLogarithmicStep       +self%coefficientFactorBoostLogarithmicOffset       )
            do j=1,countFactorBoostStellar
               coefficientFactorBoostStellar_                =10.0d0**(dble(j-1)*coefficientFactorBoostStellarLogarithmicStep+self%coefficientFactorBoostStellarLogarithmicOffset)
               !$omp do
               do k=1,countRadii
                  radiusMinimum                              = 0.0d0
                  radiusMaximum                              =10.0d0**(dble(k-1)*radiusScaleFreeLogarithmicStep              +self%radiusScaleFreeLogarithmicOffset              )
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
         file=hdf5Object(char(self%filenameTable),overWrite=.true.,readOnly=.false.)
         call file%writeAttribute(self%coefficientFactorBoostMinimum                       ,'coefficientFactorBoostMinimum'                      )
         call file%writeAttribute(self%coefficientFactorBoostMaximum                       ,'coefficientFactorBoostMaximum'                      )
         call file%writeAttribute(self%coefficientFactorBoostStellarMinimum                ,'coefficientFactorBoostStellarMinimum'               )
         call file%writeAttribute(self%coefficientFactorBoostStellarMaximum                ,'coefficientFactorBoostStellarMaximum'               )
         call file%writeAttribute(self%radiusScaleFreeMinimum                              ,'radiusScaleFreeMinimum'                             )
         call file%writeAttribute(self%radiusScaleFreeMaximum                              ,'radiusScaleFreeMaximum'                             )
         call file%writeAttribute(self%coefficientFactorBoostLogarithmicOffset             ,'coefficientFactorBoostLogarithmicOffset'            )
         call file%writeAttribute(self%coefficientFactorBoostStellarLogarithmicOffset      ,'coefficientFactorBoostStellarLogarithmicOffset'     )
         call file%writeAttribute(self%radiusScaleFreeLogarithmicOffset                    ,'radiusScaleFreeLogarithmicOffset'                   )
         call file%writeAttribute(self%coefficientFactorBoostLogarithmicStepInverse        ,'coefficientFactorBoostLogarithmicStepInverse'       )
         call file%writeAttribute(self%coefficientFactorBoostStellarLogarithmicStepInverse ,'coefficientFactorBoostStellarLogarithmicStepInverse')
         call file%writeAttribute(self%radiusScaleFreeLogarithmicStepInverse               ,'radiusScaleFreeLogarithmicStepInverse'              )
         call file%writeDataset  (self%integralPartiallyMolecularTable                     ,'integral'                                           )
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
      i                                  =int((log10(coefficientFactorBoost_       )-self%coefficientFactorBoostLogarithmicOffset       )*self%coefficientFactorBoostLogarithmicStepInverse       )+1
      j                                  =int((log10(coefficientFactorBoostStellar_)-self%coefficientFactorBoostStellarLogarithmicOffset)*self%coefficientFactorBoostStellarLogarithmicStepInverse)+1
      hi                                 =    (log10(coefficientFactorBoost_       )-self%coefficientFactorBoostLogarithmicOffset       )*self%coefficientFactorBoostLogarithmicStepInverse        +1.0d0-dble(i)
      hj                                 =    (log10(coefficientFactorBoostStellar_)-self%coefficientFactorBoostStellarLogarithmicOffset)*self%coefficientFactorBoostStellarLogarithmicStepInverse +1.0d0-dble(j)
      k                                  =int((log10(radiusScaleFree               )-self%radiusScaleFreeLogarithmicOffset              )*self%radiusScaleFreeLogarithmicStepInverse              )+1
      hk                                 =    (log10(radiusScaleFree               )-self%radiusScaleFreeLogarithmicOffset              )*self%radiusScaleFreeLogarithmicStepInverse               +1.0d0-dble(k)
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
                    &   *                                      hhi            &
                    &   *                                           hhj       &
                    &   *                                                hhk
            end do
         end do
      end do
      integralAnalyticPartiallyMolecularGeneric=+integralAnalyticPartiallyMolecularGeneric &
           &                                    +exp(integral)                             &
           &                                    *coefficientNormalization
      return
    end function integralAnalyticPartiallyMolecularGeneric
    
    logical function tableIsInsufficient()
      !!{
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
      !!{
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
    !!{
    Root function used in finding the radius in a disk where the pressure ratio exceeds the critical ratio.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    blitz2006CriticalDensityRoot=self_%pressureRatio(node_,radius)-1.0d0
    return
  end function blitz2006CriticalDensityRoot

  double precision function blitz2006PressureRatio(self,node,radius,surfaceDensityGas) result(pressureRatio)
    !!{
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
