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
  Implementation of the \cite{krumholz_star_2009} star formation rate surface density law for galactic disks.
  !!}

  use :: Abundances_Structure, only : abundances
  use :: Kind_Numbers        , only : kind_int8
  use :: Math_Exponentiation , only : fastExponentiator
  use :: Tables              , only : table1DLinearLinear
  use :: Root_Finder         , only : rootFinder

  !![
  <starFormationRateSurfaceDensityDisks name="starFormationRateSurfaceDensityDisksKrumholz2009">
   <description>
    A star formation rate surface density class implementing the model of \citep{krumholz_star_2009}:
    \begin{equation}
     \dot{\Sigma}_\star(R) = \nu_\mathrm{SF} f_\mathrm{H_2}(R)\Sigma_\mathrm{HI, disk}(R) \left\{ \begin{array}{ll}
     (\Sigma_\mathrm{HI}/\Sigma_0)^{-1/3}, &amp; \hbox{ if } \Sigma_\mathrm{HI}/\Sigma_0 \le 1 \\
     (\Sigma_\mathrm{HI}/\Sigma_0)^{1/3}, &amp; \hbox{ if } \Sigma_\mathrm{HI}/\Sigma_0 &gt; 1 \end{array} \right. ,
    \end{equation}
    where $\nu_\mathrm{SF}=${\normalfont \ttfamily [frequencyStarFormation]} is a frequency and $\Sigma_0=85 M_\odot
    \hbox{pc}^{-2}$. The molecular fraction is given by
    \begin{equation}
     f_\mathrm{H_2} = 1 - \left( 1 + \left[ { 3 s \over 4 (1+\delta)} \right]^{-5} \right)^{-1/5},
    \end{equation}
    where
    \begin{equation}
     \delta = 0.0712 \left[ 0.1 s^{-1} + 0.675 \right]^{-2.8},
    \end{equation}
    and
    \begin{equation}
     s = {\ln(1+0.6\chi+0.01\chi^2) \over 0.04 \Sigma_\mathrm{comp,0} Z^\prime},
    \end{equation}
    with
    \begin{equation}
     \chi = 0.77 \left[ 1 + 3.1 Z^{\prime 0.365} \right],
    \end{equation}
    and $\Sigma_\mathrm{comp,0}=c \Sigma_\mathrm{HI}/M_\odot \hbox{pc}^{-2}$ where $c=${\normalfont \ttfamily
    [clumpingFactorMolecularComplex]} is a density enhancement factor relating the surface density of molecular complexes to
    the gas density on larger scales. Alternatively, if {\normalfont \ttfamily [molecularFractionFast]} is set to true, the
    molecular fraction will be computed using the faster (but less accurate at low molecular fraction) formula
    \begin{equation}
     f_\mathrm{H_2} = 1 - { 3s/4 \over (1 + s/4)}.
    \end{equation}
   </description>
  </starFormationRateSurfaceDensityDisks>
  !!]
  type, extends(starFormationRateSurfaceDensityDisksClass) :: starFormationRateSurfaceDensityDisksKrumholz2009
     !!{
     Implementation of the \cite{krumholz_star_2009} star formation rate surface density law for galactic disks.
     !!}
     private
     integer         (kind_int8          )                  :: lastUniqueID
     logical                                                :: factorsComputed
     double precision                                       :: massGasPrevious                   , radiusPrevious                , &
          &                                                    radiusCriticalPrevious            , radiusMaximumPrevious
     type            (abundances         )                  :: abundancesFuelPrevious
     double precision                                       :: chi                               , radiusDisk                    , &
          &                                                    massGas                           , hydrogenMassFraction          , &
          &                                                    metallicityRelativeToSolar        , sNormalization                , &
          &                                                    sigmaMolecularComplexNormalization, clumpingFactorMolecularComplex, &
          &                                                    frequencyStarFormation
     logical                                                :: assumeMonotonicSurfaceDensity     , molecularFractionFast
     type            (rootFinder         )                  :: finderCritical                    , finderMolecules
     type            (fastExponentiator  )                  :: surfaceDensityExponentiator
     type            (table1DLinearLinear)                  :: molecularFractionTable
     procedure       (double precision   ), nopass, pointer :: molecularFractionFunction
   contains
     !![
     <methods>
       <method description="Reset memoized calculations."              method="calculationReset"      />
       <method description="Compute constant factors required."        method="computeFactors"        />
       <method description="Compute surface density factors required." method="surfaceDensityFactors" />
       <method description="Compute the molecular fraction."           method="molecularFraction"     />
     </methods>
     !!]
     final     ::                          krumholz2009Destructor
     procedure :: autoHook              => krumholz2009AutoHook
     procedure :: computeFactors        => krumholz2009ComputeFactors
     procedure :: surfaceDensityFactors => krumholz2009SurfaceDensityFactors
     procedure :: calculationReset      => krumholz2009CalculationReset
     procedure :: rate                  => krumholz2009Rate
     procedure :: unchanged             => krumholz2009Unchanged
     procedure :: intervals             => krumholz2009Intervals
     procedure :: molecularFraction     => krumholz2009MolecularFraction 
  end type starFormationRateSurfaceDensityDisksKrumholz2009

  interface starFormationRateSurfaceDensityDisksKrumholz2009
     !!{
     Constructors for the \refClass{starFormationRateSurfaceDensityDisksKrumholz2009} star formation surface density rate in disks class.
     !!}
     module procedure krumholz2009ConstructorParameters
     module procedure krumholz2009ConstructorInternal
  end interface starFormationRateSurfaceDensityDisksKrumholz2009

  ! Module-scope pointer to the active node.
  class           (starFormationRateSurfaceDensityDisksKrumholz2009), pointer   :: self_
  type            (treeNode                                        ), pointer   :: node_
  !$omp threadprivate(self_,node_)

  ! Range of s-parameter to tabulate.
  double precision                                                  , parameter :: sMinimum=0.0d+0, sMaximum=10.0d0

contains

  function krumholz2009ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{starFormationRateSurfaceDensityDisksKrumholz2009} star formation surface density rate in disks class which takes a parameter set as input.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (starFormationRateSurfaceDensityDisksKrumholz2009)                :: self
    type            (inputParameters                                 ), intent(inout) :: parameters
    double precision                                                                  :: frequencyStarFormation, clumpingFactorMolecularComplex
    logical                                                                           :: molecularFractionFast , assumeMonotonicSurfaceDensity

    !![
    <inputParameter>
      <name>frequencyStarFormation</name>
      <defaultSource>\citep{krumholz_star_2009}</defaultSource>
      <defaultValue>0.385d0</defaultValue>
      <description>The star formation frequency (in units of Gyr$^{-1}$) in the ``Krumholz-McKee-Tumlinson'' star formation timescale calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>clumpingFactorMolecularComplex</name>
      <defaultValue>5.0d0</defaultValue>
      <description>The density enhancement (relative to mean disk density) for molecular complexes in the ``Krumholz-McKee-Tumlinson'' star formation timescale calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>molecularFractionFast</name>
      <defaultValue>.false.</defaultValue>
      <description>Selects whether the fast (but less accurate) fitting formula for molecular hydrogen should be used in the ``Krumholz-McKee-Tumlinson'' star formation timescale calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>assumeMonotonicSurfaceDensity</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, assume that the surface density in disks is always monotonically decreasing.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=starFormationRateSurfaceDensityDisksKrumholz2009(frequencyStarFormation,clumpingFactorMolecularComplex,molecularFractionFast,assumeMonotonicSurfaceDensity)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function krumholz2009ConstructorParameters

  function krumholz2009ConstructorInternal(frequencyStarFormation,clumpingFactorMolecularComplex,molecularFractionFast,assumeMonotonicSurfaceDensity) result(self)
    !!{
    Internal constructor for the \refClass{starFormationRateSurfaceDensityDisksKrumholz2009} star formation surface density rate from disks class.
    !!}
    use :: Abundances_Structure                      , only : unitAbundances
    use :: Table_Labels                              , only : extrapolationTypeFix
    use :: Root_Finder                               , only : rangeExpandMultiplicative        , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    use :: Star_Formation_Rate_Krumholz2009_Utilities, only : krumholz2009MolecularFractionSlow
    implicit none
    type            (starFormationRateSurfaceDensityDisksKrumholz2009)                :: self
    double precision                                                  , intent(in   ) :: frequencyStarFormation      , clumpingFactorMolecularComplex
    logical                                                           , intent(in   ) :: molecularFractionFast       , assumeMonotonicSurfaceDensity
    integer                                                           , parameter     :: sCount                =1000
    integer                                                                           :: i
    !![
    <constructorAssign variables="frequencyStarFormation, clumpingFactorMolecularComplex, molecularFractionFast, assumeMonotonicSurfaceDensity"/>
    !!]

    self%lastUniqueID          =-1_kind_int8
    self%massGasPrevious       =-1.0d0
    self%radiusPrevious        =-1.0d0
    self%radiusCriticalPrevious=-1.0d0
    self%radiusMaximumPrevious =-1.0d0
    self%abundancesFuelPrevious=-unitAbundances
    self%factorsComputed   =.false.
    ! Set a pointer to the molecular hydrogen fraction fitting function to be used.
    select case (molecularFractionFast)
    case (.true.)
       self%molecularFractionFunction => krumholz2009MolecularFractionFast
    case(.false.)
       self%molecularFractionFunction => krumholz2009MolecularFractionSlow
    end select
    ! Build a table of molecular fraction for fast look-up.
    call self%molecularFractionTable%create(sMinimum,sMaximum,sCount,extrapolationType=[extrapolationTypeFix,extrapolationTypeFix])
    do i=1,sCount
       call self%molecularFractionTable%populate(self%molecularFractionFunction(self%molecularFractionTable%x(i)),i)
    end do
    ! Initialize exponentiator.
    self%surfaceDensityExponentiator=fastExponentiator(1.0d0,1000.0d0,0.33d0,100.0d0,.false.)
    ! Build root finders.
    self%finderCritical =rootFinder(                                                                 &
         &                          rootFunction                 =krumholz2009CriticalDensityRoot  , &
         &                          toleranceAbsolute            =0.0d+0                           , &
         &                          toleranceRelative            =1.0d-4                           , &
         &                          rangeExpandUpward            =2.0d0                            , &
         &                          rangeExpandDownward          =0.5d0                            , &
         &                          rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative    , &
         &                          rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive    , &
         &                          rangeExpandType              =rangeExpandMultiplicative          &
         &                         )
    self%finderMolecules=rootFinder(                                                                 &
         &                          rootFunction                 =krumholz2009MolecularFractionRoot, &
         &                          toleranceAbsolute            =0.0d+0                           , &
         &                          toleranceRelative            =1.0d-4                           , &
         &                          rangeExpandUpward            =2.0d0                            , &
         &                          rangeExpandDownward          =0.5d0                            , &
         &                          rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative    , &
         &                          rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive    , &
         &                          rangeExpandType              =rangeExpandMultiplicative          &
         &                         )
    return
  end function krumholz2009ConstructorInternal

  subroutine krumholz2009AutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(starFormationRateSurfaceDensityDisksKrumholz2009), intent(inout) :: self

    call calculationResetEvent%attach(self,krumholz2009CalculationReset,openMPThreadBindingAllLevels,label='starFormationRateSurfaceDensityDisksKrumholz2009')
    return
  end subroutine krumholz2009AutoHook

  subroutine krumholz2009Destructor(self)
    !!{
    Destructor for the \refClass{starFormationRateSurfaceDensityDisksKrumholz2009} star formation surface density rate from disks class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(starFormationRateSurfaceDensityDisksKrumholz2009), intent(inout) :: self

    if (calculationResetEvent%isAttached(self,krumholz2009CalculationReset))                            &
         & call calculationResetEvent%detach                        (self,krumholz2009CalculationReset)
    call        self                 %molecularFractionTable%destroy(                                 )
    return
  end subroutine krumholz2009Destructor

  subroutine krumholz2009CalculationReset(self,node,uniqueID)
    !!{
    Reset the Kennicutt-Schmidt relation calculation.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (starFormationRateSurfaceDensityDisksKrumholz2009), intent(inout) :: self
    type   (treeNode                                        ), intent(inout) :: node
    integer(kind_int8                                       ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%factorsComputed       =.false.
    self%radiusCriticalPrevious=-1.0d0
    self%radiusMaximumPrevious =-1.0d0
    self%lastUniqueID          =uniqueID
    return
  end subroutine krumholz2009CalculationReset

  subroutine krumholz2009ComputeFactors(self,node)
    !!{
    Compute constant factors needed in the \cite{krumholz_star_2009} star formation rule.
    !!}
    use :: Abundances_Structure        , only : metallicityTypeLinearByMassSolar
    use :: Galacticus_Nodes            , only : nodeComponentDisk               , treeNode
    use :: Numerical_Constants_Prefixes, only : mega
    implicit none
    class(starFormationRateSurfaceDensityDisksKrumholz2009), intent(inout) :: self
    type (treeNode                                        ), intent(inout) :: node
    class(nodeComponentDisk                               ), pointer       :: disk
    type (abundances                                      ), save          :: abundancesFuel
    !$omp threadprivate(abundancesFuel)

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! Check if factors have been precomputed.
    if (.not.self%factorsComputed) then
       ! Get the disk properties.
       disk            => node%disk   ()
       self%massGas    =  disk%massGas()
       self%radiusDisk =  disk%radius ()
       ! Find the hydrogen fraction in the disk gas of the fuel supply.
       abundancesFuel=disk%abundancesGas()
       call abundancesFuel%massToMassFraction(self%massGas)
       self%hydrogenMassFraction=abundancesFuel%hydrogenMassFraction()
       ! Get the metallicity in Solar units, and related quantities.
       self%metallicityRelativeToSolar=abundancesFuel%metallicity(metallicityTypeLinearByMassSolar)
       if (self%metallicityRelativeToSolar > 0.0d0) then
          self%chi                               =0.77d0*(1.0d0+3.1d0*self%metallicityRelativeToSolar**0.365d0)
          self%sigmaMolecularComplexNormalization=self%hydrogenMassFraction*self%clumpingFactorMolecularComplex/mega**2
          self%sNormalization                    =log(1.0d0+0.6d0*self%chi+0.01d0*self%chi**2)/(0.04d0*self%metallicityRelativeToSolar)
       end if
       ! Record that factors have now been computed.
       self%factorsComputed=.true.
    end if
    return
  end subroutine krumholz2009ComputeFactors

  double precision function krumholz2009Rate(self,node,radius)
    !!{
    Returns the star formation rate surface density (in $M_\odot$ Gyr$^{-1}$ Mpc$^{-2}$) for star formation
    in the galactic disk of {\normalfont \ttfamily node}. The disk is assumed to obey the
    \cite{krumholz_star_2009} star formation rule.
    !!}
    implicit none
    class           (starFormationRateSurfaceDensityDisksKrumholz2009), intent(inout) :: self
    type            (treeNode                                        ), intent(inout) :: node
    double precision                                                  , intent(in   ) :: radius
    double precision                                                  , parameter     :: surfaceDensityGasTiny=1.0d-100
    double precision                                                                  :: surfaceDensityFactor          , molecularFraction             , &
         &                                                                               s                             , sigmaMolecularComplex         , &
         &                                                                               surfaceDensityGas             , surfaceDensityGasDimensionless
    
    ! Compute factors.
    call self%computeFactors(node)
    ! Check if the disk is physical.
    if     (                                                  &
         &   self%massGas                            <= 0.0d0 &
         &  .or.                                              &
         &   self%radiusDisk                         <= 0.0d0 &
         &  .or.                                              &
         &   self%metallicityRelativeToSolar         <= 0.0d0 &
         &  .or.                                              &
         &   self%sigmaMolecularComplexNormalization <= 0.0d0 &
         & ) then
       ! It is not, so return zero rate.
       krumholz2009Rate=0.0d0
    else
       ! Get surface density and related quantities.
       call self%surfaceDensityFactors(node,radius,surfaceDensityGas,surfaceDensityGasDimensionless)
       ! Check for non-positive gas mass.
       if (surfaceDensityGas <= surfaceDensityGasTiny) then
          krumholz2009Rate=0.0d0
       else
          ! Compute the molecular fraction.
          sigmaMolecularComplex=self%sigmaMolecularComplexNormalization*surfaceDensityGas
          s                    =self%sNormalization/sigmaMolecularComplex
          molecularFraction    =self%molecularFraction(s)
          ! Compute the cloud density factor.
          if      (surfaceDensityGasDimensionless <= 0.0d0) then
             surfaceDensityFactor=0.0d0
          else if (surfaceDensityGasDimensionless <  1.0d0) then
             surfaceDensityFactor=self%surfaceDensityExponentiator%exponentiate(1.0d0/surfaceDensityGasDimensionless)
          else
             surfaceDensityFactor=self%surfaceDensityExponentiator%exponentiate(      surfaceDensityGasDimensionless)
          end if
          ! Compute the star formation rate surface density.
          krumholz2009Rate=+self%frequencyStarFormation &
               &           *surfaceDensityGas           &
               &           *surfaceDensityFactor        &
               &           *molecularFraction
       end if
    end if
    return
  end function krumholz2009Rate

  double precision function krumholz2009MolecularFraction(self,s)
    !!{
    Compute the molecular fraction as a function of the $s$ parameter.
    !!}
    implicit none
    class           (starFormationRateSurfaceDensityDisksKrumholz2009), intent(inout) :: self
    double precision                                                  , intent(in   ) :: s
    
    if (s > sMaximum) then
       krumholz2009MolecularFraction=self%molecularFractionFunction            (s)
    else
       krumholz2009MolecularFraction=self%molecularFractionTable   %interpolate(s)
    end if
    return
  end function krumholz2009MolecularFraction
  
  subroutine krumholz2009SurfaceDensityFactors(self,node,radius,surfaceDensityGas,surfaceDensityGasDimensionless)
    !!{
    Compute surface density and related quantities needed for the \cite{krumholz_star_2009} star formation rate model.
    !!}
    use :: Coordinates               , only : coordinateCylindrical, assignment(=)
    use :: Galactic_Structure_Options, only : componentTypeDisk    , massTypeGaseous
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (starFormationRateSurfaceDensityDisksKrumholz2009), intent(inout) :: self
    type            (treeNode                                        ), intent(inout) :: node
    double precision                                                  , intent(in   ) :: radius
    double precision                                                  , intent(  out) :: surfaceDensityGas               , surfaceDensityGasDimensionless
    class           (massDistributionClass                           ), pointer       :: massDistribution_
    double precision                                                  , parameter     :: surfaceDensityTransition=85.0d12                                 !   M☉/Mpc²
    type            (coordinateCylindrical                           )                :: coordinates

    ! Get gas surface density.
    coordinates       =  [radius,0.0d0,0.0d0]
    massDistribution_ => node             %massDistribution(componentType=componentTypeDisk,massType=massTypeGaseous)
    surfaceDensityGas =  massDistribution_%surfaceDensity  (              coordinates                               )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Compute the cloud density factor.
    surfaceDensityGasDimensionless=+self%hydrogenMassFraction     &
         &                         *     surfaceDensityGas        &
         &                         /     surfaceDensityTransition
    return
  end subroutine krumholz2009SurfaceDensityFactors

  double precision function krumholz2009MolecularFractionFast(s)
    !!{
    Fast (but less accurate at low molecular fraction) fitting function from \cite{mckee_atomic--molecular_2010} for the
    molecular hydrogen fraction.
    !!}
    implicit none
    double precision, intent(in   ) :: s

    ! Check that s is below 2 - if it is, compute the molecular fraction, otherwise truncate to zero.
    if (s < 2.0d0) then
       krumholz2009MolecularFractionFast=1.0d0-0.75d0*s/(1.0d0+0.25d0*s)
    else
       krumholz2009MolecularFractionFast=0.0d0
    end if
    return
  end function krumholz2009MolecularFractionFast

  function krumholz2009Intervals(self,node,radiusInner,radiusOuter,intervalIsAnalytic,integralsAnalytic)
    !!{
    Returns intervals to use for integrating the \cite{krumholz_star_2009} star formation rate over a galactic disk.
    !!}
    implicit none
    class           (starFormationRateSurfaceDensityDisksKrumholz2009), intent(inout), target                      :: self
    double precision                                                                 , allocatable, dimension(:,:) :: krumholz2009Intervals
    type            (treeNode                                        ), intent(inout), target                      :: node
    double precision                                                  , intent(in   )                              :: radiusInner          , radiusOuter
    logical                                                           , intent(inout), allocatable, dimension(  :) :: intervalIsAnalytic
    double precision                                                  , intent(inout), allocatable, dimension(  :) :: integralsAnalytic
    double precision                                                                                               :: surfaceDensityGas    , surfaceDensityGasDimensionless, &
         &                                                                                                            radiusCritical       , radiusMaximum

    ! Check if we can assume a monotonic surface density.
    if (self%assumeMonotonicSurfaceDensity) then
       ! Ensure required factors are computed.
       call self%computeFactors(node)       
       ! First, if using the fast molecular fraction calculation, check if this truncates to zero somewhere in range.
       radiusMaximum=radiusOuter
       if (self%molecularFractionFast) then
          self_ => self
          node_ => node
          if (krumholz2009MolecularFractionRoot(radiusInner) <= 0.0d0) then
             ! The entire disk has no molecules, so the star formation rate is zero everywhere. Set zero intervals.
             allocate(krumholz2009Intervals(2,0))
             return
          else if (krumholz2009MolecularFractionRoot(radiusOuter) <= 0.0d0) then
             ! The outer regions of the disk have no molecules. Find the radius at which the molecular fraction drops to zero and
             ! set that as the maximum radius.
             if (self%radiusMaximumPrevious > 0.0d0) then
                radiusMaximum    =  self%finderMolecules%find(rootGuess=self%radiusMaximumPrevious)
             else
                radiusMaximum    =  self%finderMolecules%find(rootRange=[radiusInner,radiusOuter])
             end if
             self%radiusMaximumPrevious=radiusMaximum
          end if
       end if
       ! Test if the inner radius is below the surface density threshold.
       call self%surfaceDensityFactors(node,radiusInner,surfaceDensityGas,surfaceDensityGasDimensionless)
       if (surfaceDensityGasDimensionless <= 1.0d0) then
          ! The entire disk is below the critical surface density so use a single interval.
          allocate(krumholz2009Intervals(2,1))
          allocate(intervalIsAnalytic   (  1))
          intervalIsAnalytic   =.false.
          krumholz2009Intervals=reshape([radiusInner,radiusMaximum],[2,1])
       else
          ! Test the surface density at the outer radius.
          call self%surfaceDensityFactors(node,radiusMaximum,surfaceDensityGas,surfaceDensityGasDimensionless)
          if (surfaceDensityGasDimensionless >= 1.0d0) then
             ! Entire disk is above the critical surface density threshold so use a single interval.
             allocate(krumholz2009Intervals(2,1))
             allocate(intervalIsAnalytic   (  1))
             intervalIsAnalytic   =.false.
             krumholz2009Intervals=reshape([radiusInner,radiusMaximum],[2,1])
          else
             ! The disk transitions the critical surface density - attempt to locate the radius at which this happens and use two
             ! intervals split at this point.
             self_ => self
             node_ => node
             if (self%radiusCriticalPrevious > 0.0d0) then
                radiusCritical=self%finderCritical%find(rootGuess=self%radiusCriticalPrevious)
             else
                radiusCritical=self%finderCritical%find(rootRange=[radiusInner,radiusMaximum])
             end if
             self%radiusCriticalPrevious=radiusCritical
             allocate(krumholz2009Intervals(2,2))
             allocate(intervalIsAnalytic   (  2))
             intervalIsAnalytic   =.false.
             krumholz2009Intervals=reshape([radiusInner,radiusCritical,radiusCritical,radiusMaximum],[2,2])
          end if
       end if
    else
       ! Disk surface density can not be assumed to be monotonic - use a single interval.
       allocate(krumholz2009Intervals(2,1))
       allocate(intervalIsAnalytic   (  1))
       intervalIsAnalytic   =.false.
       krumholz2009Intervals=reshape([radiusInner,radiusOuter],[2,1])
    end if
    return
  end function krumholz2009Intervals
    
  double precision function krumholz2009MolecularFractionRoot(radius)
    !!{
    Function used in finding where the molecular fraction drops below a tiny threshold.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    double precision, parameter     :: molecularFractionTiny=1.0d-10
    double precision                :: surfaceDensityGas            , surfaceDensityGasDimensionless, &
         &                             sigmaMolecularComplex        , s                             , &
         &                             molecularFraction
    
    call self_%surfaceDensityFactors(node_,radius,surfaceDensityGas,surfaceDensityGasDimensionless)
    sigmaMolecularComplex            =+self_%sigmaMolecularComplexNormalization &
         &                            *      surfaceDensityGas
    s                                =+self_%sNormalization                     &
         &                            /      sigmaMolecularComplex
    molecularFraction                =+self_%molecularFraction    (s)
    krumholz2009MolecularFractionRoot=+molecularFraction                        &
         &                            -molecularFractionTiny
    return
  end function krumholz2009MolecularFractionRoot
    
  double precision function krumholz2009CriticalDensityRoot(radius)
    !!{
    Root function used in finding the radius in a disk where the surface density equals the critical surface density in the
    \cite{krumholz_star_2009} star formation rate model.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: surfaceDensityGas, surfaceDensityGasDimensionless

    call self_%surfaceDensityFactors(node_,radius,surfaceDensityGas,surfaceDensityGasDimensionless)
    krumholz2009CriticalDensityRoot=surfaceDensityGasDimensionless-1.0d0
    return
  end function krumholz2009CriticalDensityRoot

  logical function krumholz2009Unchanged(self,node)
    !!{
    Determine if the surface rate density of star formation is unchanged.
    !!}
    use :: Abundances_Structure, only : metallicityTypeLinearByMassSolar
    use :: Galacticus_Nodes    , only : nodeComponentDisk               , treeNode
    implicit none
    class           (starFormationRateSurfaceDensityDisksKrumholz2009), intent(inout) :: self
    type            (treeNode                                        ), intent(inout) :: node
    class           (nodeComponentDisk                               ), pointer       :: disk
    double precision                                                                  :: massGas       , radius
    type            (abundances                                      )                :: abundancesFuel

    disk    => node%disk   ()
    massGas =  disk%massGas()
    if (massGas > 0.0d0) then
       radius        =disk%radius       ()
       abundancesFuel=disk%abundancesGas()
       call abundancesFuel%massToMassFraction(massGas)
       if     (                                                                                                                                                             &
            &   massGas        == self%massGasPrevious                                                                                                                      &
            &  .and.                                                                                                                                                        &
            &   radius         == self%radiusPrevious                                                                                                                       &
            &  .and.                                                                                                                                                        &
            &   abundancesFuel%metallicity         (metallicityTypeLinearByMassSolar) == self%abundancesFuelPrevious%metallicity         (metallicityTypeLinearByMassSolar) &
            &  .and.                                                                                                                                                        &
            &   abundancesFuel%hydrogenMassFraction(                                ) == self%abundancesFuelPrevious%hydrogenMassFraction(                                ) &
            & ) then
          krumholz2009Unchanged      =.true.
       else
          krumholz2009Unchanged      =.false.
          self%massGasPrevious       =massGas
          self%radiusPrevious        =radius
          self%abundancesFuelPrevious=abundancesFuel
       end if
    else
       if (self%massGasPrevious == 0.0d0) then
          krumholz2009Unchanged=.true.
       else
          krumholz2009Unchanged=.false.
          self%massGasPrevious =0.0d0
       end if
    end if
    return
  end function krumholz2009Unchanged
