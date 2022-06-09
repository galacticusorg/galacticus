!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  use :: Galactic_Structure , only : galacticStructureClass
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
    are parameters and the surface density of molecular gas $\Sigma_\mathrm{H_2} = (P_\mathrm{ext}/P_0)^\alpha
    \Sigma_\mathrm{HI}$, where $\alpha=${\normalfont \ttfamily [pressureExponent]} and $P_0=${\normalfont \ttfamily
    [pressureCharacteristic]} are parameters and the hydrostatic pressure in the disk plane assuming location isothermal gas
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
     class           (galacticStructureClass), pointer :: galacticStructure_         => null()
     integer         (kind_int8             )          :: lastUniqueID
     logical                                           :: factorsComputed                     , assumeMonotonicSurfaceDensity
     double precision                                  :: heightToRadialScaleDisk             , pressureCharacteristic             , &
          &                                               pressureExponent                    , starFormationFrequencyNormalization, &
          &                                               surfaceDensityCritical              , surfaceDensityExponent             , &
          &                                               velocityDispersionDiskGas           , radiusDisk                         , &
          &                                               massGas                             , hydrogenMassFraction               , &
          &                                               massStellar                         , massGasPrevious                    , &
          &                                               massStellarPrevious                 , hydrogenMassFractionPrevious       , &
          &                                               radiusDiskPrevious                  , radiusCritical
     type            (rootFinder            )          :: finder
     type            (fastExponentiator     )          :: pressureRatioExponentiator
   contains
     !![
     <methods>
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
     !!{
     Constructors for the {\normalfont \ttfamily blitz2006} star formation surface density rate in disks class.
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
    Constructor for the {\normalfont \ttfamily blitz2006} star formation surface density rate in disks class which takes a parameter set as input.
    !!}
    implicit none
    type            (starFormationRateSurfaceDensityDisksBlitz2006)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    class           (galacticStructureClass                       ), pointer       :: galacticStructure_
    double precision                                                               :: velocityDispersionDiskGas          , heightToRadialScaleDisk, &
         &                                                                            surfaceDensityCritical             , surfaceDensityExponent , &
         &                                                                            starFormationFrequencyNormalization, pressureCharacteristic , &
         &                                                                            pressureExponent
    logical                                                                        :: assumeMonotonicSurfaceDensity

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
      <description>The ratio of scale height to scale radius for disks in the ``Blitz-Rosolowsky2006'' star formation timescale calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>surfaceDensityCritical</name>
      <defaultSource>\citep{bigiel_star_2008}</defaultSource>
      <defaultValue>200.0d0</defaultValue>
      <description>The surface density (in units of $M_\odot$ pc$^{-2}$) in the ``Blitz-Rosolowsky2006'' star formation timescale calculation at which low-density truncation begins.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>surfaceDensityExponent</name>
      <defaultSource>\citep{bigiel_star_2008}</defaultSource>
      <defaultValue>0.4d0</defaultValue>
      <description>The exponent for surface density in the ``Blitz-Rosolowsky2006'' star formation timescale calculation at in the high density regime.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>starFormationFrequencyNormalization</name>
      <defaultSource>\citep{leroy_star_2008}</defaultSource>
      <defaultValue>5.25d-10</defaultValue>
      <description>The star formation frequency (in the low-density limit and in units of yr$^{-1}$) in the ``Blitz-Rosolowsky2006'' star formation timescale calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>pressureCharacteristic</name>
      <defaultSource>\citep{blitz_role_2006}</defaultSource>
      <defaultValue>4.54d0</defaultValue>
      <description>The characteristic pressure (given as $P_0/k_\mathrm{B}$ in units of K cm$^{-3}$) in the scaling relation of molecular hydrogen fraction with disk pressure in the ``Blitz-Rosolowsky2006'' star formation timescale calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>pressureExponent</name>
      <defaultSource>\citep{blitz_role_2006}</defaultSource>
      <defaultValue>0.92d0</defaultValue>
      <description>The exponent in the scaling relation of molecular hydrogen fraction with disk pressure in the ``Blitz-Rosolowsky2006'' star formation timescale calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>assumeMonotonicSurfaceDensity</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, assume that the surface density in disks is always monotonically decreasing.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="galacticStructure" name="galacticStructure_" source="parameters"/>
    !!]
    self=starFormationRateSurfaceDensityDisksBlitz2006(velocityDispersionDiskGas,heightToRadialScaleDisk,surfaceDensityCritical,surfaceDensityExponent,starFormationFrequencyNormalization,pressureCharacteristic,pressureExponent,assumeMonotonicSurfaceDensity,galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticStructure_"/>
    !!]
    return
  end function blitz2006ConstructorParameters

  function blitz2006ConstructorInternal(velocityDispersionDiskGas,heightToRadialScaleDisk,surfaceDensityCritical,surfaceDensityExponent,starFormationFrequencyNormalization,pressureCharacteristic,pressureExponent,assumeMonotonicSurfaceDensity,galacticStructure_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily blitz2006} star formation surface density rate from disks class.
    !!}
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : massSolar                , megaParsec
    use :: Numerical_Constants_Physical    , only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes    , only : giga                     , hecto                        , kilo                         , mega
    use :: Root_Finder                     , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    implicit none
    type            (starFormationRateSurfaceDensityDisksBlitz2006)                        :: self
    class           (galacticStructureClass                       ), intent(in   ), target :: galacticStructure_
    double precision                                               , intent(in   )         :: velocityDispersionDiskGas          , heightToRadialScaleDisk, &
         &                                                                                    surfaceDensityCritical             , surfaceDensityExponent , &
         &                                                                                    starFormationFrequencyNormalization, pressureCharacteristic , &
         &                                                                                    pressureExponent
    logical                                                        , intent(in   )         :: assumeMonotonicSurfaceDensity
     !![
     <constructorAssign variables="velocityDispersionDiskGas, heightToRadialScaleDisk, surfaceDensityCritical, surfaceDensityExponent, starFormationFrequencyNormalization, pressureCharacteristic, pressureExponent, assumeMonotonicSurfaceDensity, *galacticStructure_"/>
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
         &                 rangeExpandUpward            =2.0d0                        , &
         &                 rangeExpandDownward          =0.5d0                        , &
         &                 rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
         &                 rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
         &                 rangeExpandType              =rangeExpandMultiplicative      &
         &                )
    ! Initialize memoized values.
    self%massGasPrevious             =-huge(0.0d0)
    self%massStellarPrevious         =-huge(0.0d0)
    self%radiusDiskPrevious          =-huge(0.0d0)
    self%hydrogenMassFractionPrevious=-huge(0.0d0)
    return
  end function blitz2006ConstructorInternal

  subroutine blitz2006AutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self

    call calculationResetEvent%attach(self,blitz2006CalculationReset,openMPThreadBindingAllLevels)
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
    !![
    <objectDestructor name="self%galacticStructure_"/>
    !!]
    return
  end subroutine blitz2006Destructor

  subroutine blitz2006CalculationReset(self,node)
    !!{
    Reset the Kennicutt-Schmidt relation calculation.
    !!}
    implicit none
    class(starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self
    type (treeNode                                     ), intent(inout) :: node

    self%factorsComputed=.false.
    self%lastUniqueID   =node%uniqueID()
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
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Compute factors.
    call self%computeFactors(node)
    ! Return zero rate for non-positive radius or mass.
    if (self%massGas <= 0.0d0 .or. self%massStellar < 0.0d0 .or. self%radiusDisk <= 0.0d0) then
       blitz2006Rate=0.0d0
       return
    end if
    ! Compute the pressure ratio that Blitz & Rosolowsky (2006) use to compute the molecular fraction.    
    !! We first compute the pressure ratio ignoring the boost from the stellar mass. If this already exceeds the characteristic
    !! limit then we will not need to compute the stellar contribution.
    pressureRatio=self%pressureRatio(node,radius,surfaceDensityGas)
    ! Compute the molecular fraction, limited to 100% molecular.
    if (pressureRatio >= 1.0d0) then
       molecularFraction=                                         1.0d0
    else
       molecularFraction=min(self%pressureRatioExponentiator%exponentiate(pressureRatio),1.0d0)
    end if
    ! Compute the star formation rate surface density.
    factorBoost  =+self%hydrogenMassFraction                &
         &        *surfaceDensityGas                        &
         &        /self%surfaceDensityCritical
    blitz2006Rate=+surfaceDensityGas                        &
         &        *self%hydrogenMassFraction                &
         &        *molecularFraction                        &
         &        *self%starFormationFrequencyNormalization &
         &        *(                                        &
         &          +1.0d0                                  &
         &          +(                                      &
         &            +self%hydrogenMassFraction            &
         &            *surfaceDensityGas                    &
         &            /self%surfaceDensityCritical          &
         &           )**self%surfaceDensityExponent         &
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
    use :: Abundances_Structure, only : abundances
    use :: Galacticus_Nodes    , only : nodeComponentDisk
    implicit none
    class(starFormationRateSurfaceDensityDisksBlitz2006), intent(inout) :: self
    type (treeNode                                     ), intent(inout) :: node
    class(nodeComponentDisk                            ), pointer       :: disk
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

  function blitz2006Intervals(self,node,radiusInner,radiusOuter)
    !!{
    Returns intervals to use for integrating the \cite{krumholz_star_2009} star formation rate over a galactic disk.
    !!}
    implicit none
    class           (starFormationRateSurfaceDensityDisksBlitz2006), intent(inout), target         :: self
    double precision                                               , allocatable  , dimension(:,:) :: blitz2006Intervals
    type            (treeNode                                     ), intent(inout), target         :: node
    double precision                                               , intent(in   )                 :: radiusInner       , radiusOuter

    ! Check if we can assume a monotonic surface density.
    if (self%assumeMonotonicSurfaceDensity) then
       ! Set the critical radius to a very negative value so that pressure ratio is always computed.
       self%radiusCritical=-huge(0.0d0)
       ! Compute factors.
       call self%computeFactors(node)
       ! Set single interval for non-positive radius or mass.
       if (self%massGas <= 0.0d0 .or. self%massStellar < 0.0d0 .or. self%radiusDisk <= 0.0d0) then
          allocate(blitz2006Intervals(2,1))
          blitz2006Intervals=reshape([radiusInner,radiusOuter],[2,1])
          self%radiusCritical=-huge(0.0d0)
       else
          ! Test if the inner radius is below the pressure threshold.
          if (self%pressureRatio(node,radiusInner) <= 1.0d0) then
             ! The entire disk is below the pressure threshold so use a single interval.
             allocate(blitz2006Intervals(2,1))
             blitz2006Intervals=reshape([radiusInner,radiusOuter],[2,1])
             self%radiusCritical=-huge(0.0d0)
          else
             ! Test the surface density at the outer radius.
             if (self%pressureRatio(node,radiusOuter) >= 1.0d0) then
                ! Entire disk is above the pressure threshold so use a single interval.
                allocate(blitz2006Intervals(2,1))
                blitz2006Intervals=reshape([radiusInner,radiusOuter],[2,1])
                self%radiusCritical=radiusOuter
             else
                ! The disk transitions the pressure ratio - attempt to locate the radius at which this happens and use two
                ! intervals split at this point.
                self_                => self
                node_                => node
                self %radiusCritical =  self%finder%find(rootRange=[radiusInner,radiusOuter])
                allocate(blitz2006Intervals(2,2))
                blitz2006Intervals=reshape([radiusInner,self%radiusCritical,self%radiusCritical,radiusOuter],[2,2])
             end if
          end if
       end if
    else
       ! Disk pressure can not be assumed to be monotonic - use a single interval.
       allocate(blitz2006Intervals(2,1))
       blitz2006Intervals=reshape([radiusInner,radiusOuter],[2,1])
       self%radiusCritical=radiusInner
    end if
    return
  end function blitz2006Intervals

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
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Galactic_Structure_Options      , only : componentTypeDisk              , coordinateSystemCylindrical, massTypeGaseous, massTypeStellar
    implicit none
    class           (starFormationRateSurfaceDensityDisksBlitz2006), intent(inout)           :: self
    type            (treeNode                                     ), intent(inout)           :: node
    double precision                                               , intent(in   )           :: radius
    double precision                                               , intent(  out), optional :: surfaceDensityGas
    double precision                                                                         :: surfaceDensityGas_, surfaceDensityStellar, &
         &                                                                                      factorBoostStellar

    ! Get gas surface density.
    surfaceDensityGas_=self%galacticStructure_%surfaceDensity(node,[radius,0.0d0,0.0d0],coordinateSystem=coordinateSystemCylindrical,componentType=componentTypeDisk,massType=massTypeGaseous)
    if (present(surfaceDensityGas)) surfaceDensityGas=surfaceDensityGas_
    ! If the radius is less than the critical radius the pressure radius is above 1 by definition, so simply pin it to that value.
    if (radius <= self%radiusCritical) then
       pressureRatio=1.0d0
    else
       ! Compute the pressure ratio that Blitz & Rosolowsky (2006) use to compute the molecular fraction.    
       !! We first compute the pressure ratio ignoring the boost from the stellar mass. If this already exceeds the characteristic
       !! limit then we will not need to compute the stellar contribution.
       pressureRatio=+0.5d0                              &
            &        *Pi                                 &
            &        *gravitationalConstantGalacticus    &
            &        *surfaceDensityGas_             **2 &
            &        /self%pressureCharacteristic
       if (pressureRatio > 0.0d0 .and. pressureRatio < 1.0d0) then
          ! Compute the stellar boost factor.
          surfaceDensityStellar=+self%galacticStructure_%surfaceDensity(node,[radius,0.0d0,0.0d0],coordinateSystem=coordinateSystemCylindrical,componentType=componentTypeDisk,massType=massTypeStellar) 
          factorBoostStellar   =+1.0d0                                 &
               &                +self%velocityDispersionDiskGas        &
               &                /surfaceDensityGas_                    &
               &                *sqrt(                                 &
               &                      +surfaceDensityStellar           &
               &                      /Pi                              &
               &                      /gravitationalConstantGalacticus &
               &                      /self%heightToRadialScaleDisk    &
               &                      /self%radiusDisk                 &
               &                     ) 
          pressureRatio        =+pressureRatio                         &
               &                *factorBoostStellar
       end if
    end if
    return
  end function blitz2006PressureRatio
  
