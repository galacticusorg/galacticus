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
  An implementation of accretion from the \gls{igm} onto halos accounting for the effects of isocurvature perturbations
  following the model of \cite{jessop_ripples_2026}.
  !!}

  use :: Cosmological_Density_Field     , only : criticalOverdensityClass
  use :: Cosmology_Parameters           , only : cosmologyParametersClass
  use :: Linear_Growth                  , only : linearGrowthClass
  use :: Power_Spectrum_Window_Functions, only : powerSpectrumWindowFunctionTopHat
  use :: Tables                         , only : table1DGeneric
  use :: Numerical_Interpolation        , only : interpolator

  !![
  <accretionHalo name="accretionHaloIsocurvature">
   <description>
    An accretion onto halos decorator class which reduces the mass accreted to account for the effects of isocurvature
    perturbations following the model of \cite{jessop_ripples_2026}. \gls{class} is used to compute perturbations in baryons and
    cold dark matter, $\delta_\mathrm{b}$ and $\delta_\mathrm{c}$ respectively, as a function of wavenumber. The correlation between
    these two is then computed using \citep[][eqn.~4]{jessop_ripples_2026}:
    \begin{equation}
     \alpha_0 = \frac{\int_0^\infty 4 \pi k^2 |\tilde{W}(k|M)|^2 \delta_\mathrm{bc}(k) \delta_\mathrm{c}(k) \mathrm{d}k}{\int_0^\infty 4 \pi k^2 |\tilde{W}(k|M)|^2 \delta_\mathrm{c}^2(k) \mathrm{d}k},
    \end{equation}
    where $\delta_\mathrm{bc} = \delta_\mathrm{b}-\delta_\mathrm{c}$, and $\tilde{W}(k|M)$ is the Fourier transform of a top-hat
    window function for mass scale $M$.

    The fraction of mass accreted into a halo is then reduced by a factor \citep[][eqn.~9]{jessop_ripples_2026}:
    \begin{equation}
     \frac{f_\mathrm{b}}{\bar{f}_\mathrm{b}} = 1 + \frac{(1-\bar{f}_\mathrm{b}) \delta_\mathrm{c}(t) \alpha_0}{D(t)},
    \end{equation}
    where $\bar{f}_\mathrm{b}$ is the universal baryon fraction, $\delta_\mathrm{c}(t)$ is the critical overdensity for halo
    collapse, and $D(t)$ is the linear growth factor.
   </description>
   <deepCopy>
     <functionClass variables="powerSpectrumWindowFunction_"/>
   </deepCopy>
   <stateStorable>
     <functionClass variables="powerSpectrumWindowFunction_"/>
   </stateStorable>
  </accretionHalo>
  !!]
  type, extends(accretionHaloClass) :: accretionHaloIsocurvature
     !!{
     A halo accretion class in which accretion is reduced to account for the effects of isocurvature perturbations
     following the model of \cite{jessop_ripples_2026}
     !!}
     private
     class           (accretionHaloClass               ), pointer :: accretionHalo_               => null()
     class           (cosmologyParametersClass         ), pointer :: cosmologyParameters_         => null()
     class           (criticalOverdensityClass         ), pointer :: criticalOverdensity_         => null()
     class           (linearGrowthClass                ), pointer :: linearGrowth_                => null()
     type            (powerSpectrumWindowFunctionTopHat), pointer :: powerSpectrumWindowFunction_ => null()
     double precision                                             :: wavenumberMaximum                     , fractionBaryonsUniversal, &
          &                                                          massMinimum                           , massMaximum
     logical                                                      :: wavenumberMaximumReached              , initialized
     integer                                                      :: countPerDecade
     type            (table1dGeneric                   )          :: perturbationsDarkMatter               , perturbationsBaryons
     type            (interpolator                     )          :: correlation
   contains
     !![
     <methods>
       <method method="fraction" description="Returns the fraction of baryons for halos of this mass, relative to the universal baryon fraction."/>
     </methods>
     !!]
     final     ::                              isocurvatureDestructor
     procedure :: fraction                  => isocurvatureFraction
     procedure :: branchHasBaryons          => isocurvatureBranchHasBaryons
     procedure :: accretionRate             => isocurvatureAccretionRate
     procedure :: accretedMass              => isocurvatureAccretedMass
     procedure :: failedAccretionRate       => isocurvatureFailedAccretionRate
     procedure :: failedAccretedMass        => isocurvatureFailedAccretedMass
     procedure :: accretionRateMetals       => isocurvatureAccretionRateMetals
     procedure :: accretedMassMetals        => isocurvatureAccretedMassMetals
     procedure :: failedAccretionRateMetals => isocurvatureFailedAccretionRateMetals
     procedure :: failedAccretedMassMetals  => isocurvatureFailedAccretedMassMetals
     procedure :: accretionRateChemicals    => isocurvatureAccretionRateChemicals
     procedure :: accretedMassChemicals     => isocurvatureAccretedMassChemicals
  end type accretionHaloIsocurvature

  interface accretionHaloIsocurvature
     !!{
     Constructors for the \refClass{accretionHaloIsocurvature} halo accretion class.
     !!}
     module procedure isocurvatureConstructorParameters
     module procedure isocurvatureConstructorInternal
  end interface accretionHaloIsocurvature

  ! Smallest maximum wavenumber to tabulate.
  double precision                           , parameter :: classWavenumberMaximumLimit=25000.0d0

  ! Submodule-scope variables used in integrands.
  class           (accretionHaloIsocurvature), pointer   :: self_
  double precision                                       :: massSmoothing_
  !$omp threadprivate(self_,massSmoothing_)
  
contains

  function isocurvatureConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the \refClass{accretionHaloIsocurvature} halo accretion class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (accretionHaloIsocurvature)                :: self
    type   (inputParameters          ), intent(inout) :: parameters
    class  (accretionHaloClass       ), pointer       :: accretionHalo_
    class  (cosmologyParametersClass ), pointer       :: cosmologyParameters_
    class  (criticalOverdensityClass ), pointer       :: criticalOverdensity_
    class  (linearGrowthClass        ), pointer       :: linearGrowth_
    integer                                           :: countPerDecade

    !![
    <inputParameter>
      <name>countPerDecade</name>
      <source>parameters</source>
      <defaultValue>100</defaultValue>
      <description>The number of points per decade of wavenumber to compute in the CLASS perturbations. A value of 0 allows CLASS to choose what it considers to be optimal spacing of wavenumbers.</description>
    </inputParameter>
    <objectBuilder class="accretionHalo"       name="accretionHalo_"       source="parameters"/>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="criticalOverdensity" name="criticalOverdensity_" source="parameters"/>
    <objectBuilder class="linearGrowth"        name="linearGrowth_"        source="parameters"/>
    !!]
    self=accretionHaloIsocurvature(countPerDecade,accretionHalo_,cosmologyParameters_,criticalOverdensity_,linearGrowth_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="accretionHalo_"      />
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="criticalOverdensity_"/>
    <objectDestructor name="linearGrowth_"       />
   !!]
    return
  end function isocurvatureConstructorParameters

  function isocurvatureConstructorInternal(countPerDecade,accretionHalo_,cosmologyParameters_,criticalOverdensity_,linearGrowth_) result(self)
    !!{
    Internal constructor for the \refClass{accretionHaloIsocurvature} halo accretion class.
    !!}
    use :: Cosmology_Parameters , only : hubbleUnitsLittleH
    implicit none
    type   (accretionHaloIsocurvature), target                :: self
    class  (accretionHaloClass       ), intent(in   ), target :: accretionHalo_
    class  (cosmologyParametersClass ), intent(in   ), target :: cosmologyParameters_
    class  (criticalOverdensityClass ), intent(in   ), target :: criticalOverdensity_
    class  (linearGrowthClass        ), intent(in   ), target :: linearGrowth_
    integer                           , intent(in   )         :: countPerDecade
    !![
    <constructorAssign variables="countPerDecade, *accretionHalo_, *cosmologyParameters_, *criticalOverdensity_, *linearGrowth_"/>
    !!]

    ! Set initialization state.
    self%initialized             =.false.
    self%massMinimum             =+huge(0.0d0)
    self%massMaximum             =-huge(0.0d0)
    ! Set maximum wavenumber.
    self%wavenumberMaximum       =+classWavenumberMaximumLimit                                        &
         &                        *self%cosmologyParameters_%hubbleConstant(units=hubbleUnitsLittleH)
    self%wavenumberMaximumReached=.false.
    ! Compute mean baryon fraction.
    self%fractionBaryonsUniversal=+self%cosmologyParameters_%OmegaBaryon() &
         &                        /self%cosmologyParameters_%OmegaMatter()
    ! Construct a window function.
    allocate(self%powerSpectrumWindowFunction_)
    !![
    <referenceConstruct owner="self" isResult="yes" object="powerSpectrumWindowFunction_">
      <constructor>
	powerSpectrumWindowFunctionTopHat(self%cosmologyParameters_)
      </constructor>
    </referenceConstruct>
    !!]
    return
  end function isocurvatureConstructorInternal

  subroutine isocurvatureDestructor(self)
    !!{
    Destructor for the \refClass{accretionHaloIsocurvature} halo accretion class.
    !!}
    implicit none
    type(accretionHaloIsocurvature), intent(inout) :: self

    !![
    <objectDestructor name="self%accretionHalo_"              />
    <objectDestructor name="self%cosmologyParameters_"        />
    <objectDestructor name="self%powerSpectrumWindowFunction_"/>
    <objectDestructor name="self%criticalOverdensity_"        />
    <objectDestructor name="self%linearGrowth_"               />
    !!]
    return
  end subroutine isocurvatureDestructor

  logical function isocurvatureBranchHasBaryons(self,node) result(branchHasBaryons)
    !!{
    Returns true if this branch can accrete any baryons.
    !!}
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodesBranch
    implicit none
    class(accretionHaloIsocurvature), intent(inout)          :: self
    type (treeNode                 ), intent(inout), target  :: node

    branchHasBaryons= self               %fraction        (node) > 0.0d0 &
         &           .and.                                               &
         &            self%accretionHalo_%branchHasBaryons(node)
    return
  end function isocurvatureBranchHasBaryons

  double precision function isocurvatureAccretionRate(self,node,accretionMode) result(rateAccretion)
    !!{
    Computes the baryonic accretion rate onto \mono{node}.
    !!}
    implicit none
    class(accretionHaloIsocurvature   ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode

    rateAccretion=+self               %fraction     (node              ) &
         &        *self%accretionHalo_%accretionRate(node,accretionMode)
    return
  end function isocurvatureAccretionRate

  double precision function isocurvatureAccretedMass(self,node,accretionMode) result(massAccreted)
    !!{
    Computes the mass of baryons accreted into \mono{node}.
    !!}
    implicit none
    class(accretionHaloIsocurvature   ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode
 
    massAccreted=+self               %fraction    (node              ) &
         &       *self%accretionHalo_%accretedMass(node,accretionMode)
    return
  end function isocurvatureAccretedMass

  double precision function isocurvatureFailedAccretionRate(self,node,accretionMode) result(rateAccretionFailed)
    !!{
    Computes the baryonic accretion rate onto \mono{node}.
    !!}
    implicit none
    class(accretionHaloIsocurvature   ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode

    rateAccretionFailed=+self               %fraction           (node              ) &
         &              *self%accretionHalo_%failedAccretionRate(node,accretionMode)
    return
  end function isocurvatureFailedAccretionRate

  double precision function isocurvatureFailedAccretedMass(self,node,accretionMode) result(massAccretedFailed)
    !!{
    Computes the mass of baryons accreted into \mono{node}.
    !!}
    implicit none
    class(accretionHaloIsocurvature   ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode

    massAccretedFailed=+self               %fraction          (node              ) &
         &             *self%accretionHalo_%failedAccretedMass(node,accretionMode)
    return
  end function isocurvatureFailedAccretedMass

  function isocurvatureAccretionRateMetals(self,node,accretionMode) result(rateAccretionMetals)
    !!{
    Computes the rate of mass of abundance accretion (in $M_\odot/$Gyr) onto \mono{node} from the intergalactic medium.
    !!}
    use :: Abundances_Structure, only : operator(*)
    implicit none
    type (abundances                  )                :: rateAccretionMetals
    class(accretionHaloIsocurvature   ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode

    rateAccretionMetals=+self               %fraction           (node              ) &
         &              *self%accretionHalo_%accretionRateMetals(node,accretionMode)
    return
  end function isocurvatureAccretionRateMetals

  function isocurvatureAccretedMassMetals(self,node,accretionMode) result(massAccretedMetals)
    !!{
    Computes the mass of abundances accreted (in $M_\odot$) onto \mono{node} from the intergalactic medium.
    !!}
    use :: Abundances_Structure, only : operator(*)
    implicit none
    type (abundances                  )                :: massAccretedMetals
    class(accretionHaloIsocurvature   ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode

    massAccretedMetals=+self               %fraction          (node              ) &
         &             *self%accretionHalo_%accretedMassMetals(node,accretionMode)
    return
  end function isocurvatureAccretedMassMetals

  function isocurvatureFailedAccretionRateMetals(self,node,accretionMode) result(rateAccretionMetalsFailed)
    !!{
    Computes the rate of failed mass of abundance accretion (in $M_\odot/$Gyr) onto \mono{node} from the intergalactic medium.
    !!}
    use :: Abundances_Structure, only : operator(*)
    implicit none
    type (abundances                  )                :: rateAccretionMetalsFailed
    class(accretionHaloIsocurvature   ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode
 
    rateAccretionMetalsFailed=+self               %fraction                 (node              ) &
         &                    *self%accretionHalo_%failedAccretionRateMetals(node,accretionMode)
    return
  end function isocurvatureFailedAccretionRateMetals

  function isocurvatureFailedAccretedMassMetals(self,node,accretionMode) result(massAccretedMetalsFailed)
    !!{
    Computes the mass of abundances that failed to accrete (in $M_\odot$) onto \mono{node} from the intergalactic medium.
    !!}
    use :: Abundances_Structure, only : operator(*)
    implicit none
    type (abundances                  )                :: massAccretedMetalsFailed
    class(accretionHaloIsocurvature   ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode
  
    massAccretedMetalsFailed=+self               %fraction                (node              ) &
         &                   *self%accretionHalo_%failedAccretedMassMetals(node,accretionMode)
    return
  end function isocurvatureFailedAccretedMassMetals
  
  function isocurvatureAccretionRateChemicals(self,node,accretionMode) result(rateAccretionChemicals)
    !!{
    Computes the rate of mass of chemicals accretion (in $M_\odot/$Gyr) onto \mono{node} from the intergalactic medium. Assumes a
    primordial mixture of hydrogen and helium and that accreted material is in collisional ionization equilibrium at the virial
    temperature.
    !!}
    use :: Chemical_Abundances_Structure, only : operator(*)
    implicit none
    type (chemicalAbundances          )                :: rateAccretionChemicals
    class(accretionHaloIsocurvature   ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode
  
    rateAccretionChemicals=+self               %fraction              (node              ) &
         &                 *self%accretionHalo_%accretionRateChemicals(node,accretionMode)
    return
  end function isocurvatureAccretionRateChemicals

  function isocurvatureAccretedMassChemicals(self,node,accretionMode) result(massAccretedChemicals)
    !!{
    Computes the mass of chemicals accreted (in $M_\odot$) onto \mono{node} from the intergalactic medium.
    !!}
    use :: Chemical_Abundances_Structure, only : operator(*)
    implicit none
    type (chemicalAbundances          )                :: massAccretedChemicals
    class(accretionHaloIsocurvature   ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode
 
    massAccretedChemicals=+self               %fraction             (node              ) &
         &                *self%accretionHalo_%accretedMassChemicals(node,accretionMode)
    return
  end function isocurvatureAccretedMassChemicals

  double precision function isocurvatureFraction(self,node) result(fraction)
    !!{    
    Computes the fraction of the universal baryon fraction present in the region from which the given \mono{node} is accreting,
    accounting for isocurvature perturbations following the model of \cite{jessop_ripples_2026}.
    !!}
    use :: Interfaces_CLASS     , only : Interface_CLASS_Perturbations
    use :: Galacticus_Nodes     , only : nodeComponentBasic
    use :: Numerical_Ranges     , only : Make_Range                   , rangeTypeLogarithmic
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (accretionHaloIsocurvature), intent(inout), target     :: self
    type            (treeNode                 ), intent(inout)             :: node
    integer                                    , parameter                 :: countPointsPerDecade =20
    double precision                           , parameter                 :: toleranceRelative    =1.0d-6, wavenumberLarge =1.0d3, &
         &                                                                    factorWavenumberLarge=1.0d3
    class           (nodeComponentBasic       )               , pointer    :: basic
    double precision                           , dimension(1)              :: redshifts
    double precision                           , dimension(:), allocatable :: mass                        , correlation
    double precision                                                       :: wavenumberMaximum           , massSmoothing         , &
         &                                                                    correlationDirect           , correlationCross
    logical                                                                :: makePerturbations           , remakeTable
    integer                                                                :: i                           , countPointsInTable
    type            (integrator               )                            :: integratorDirect            , integratorCross

    ! Determine if the tabulated correlation needs to be expanded.
    basic         => node %basic()
    massSmoothing =  basic%mass ()
    remakeTable   =  .false.
    if (massSmoothing < self%massMinimum) then
       self%massMinimum=massSmoothing/2.0d0
       remakeTable     =.true.
    end if
    if (massSmoothing > self%massMaximum) then
       self%massMaximum=massSmoothing*2.0d0
       remakeTable     =.true.
    end if
    if (remakeTable) then
       ! The tabulate correlation must be expanded. Next, determine if the tabulation of perturbations must be expanded.
       makePerturbations=.false.
       wavenumberMaximum=max(wavenumberLarge,factorWavenumberLarge/radiusLagrangian(self%massMinimum))
       if (self%initialized) then
          if     (                                                             &
               &   wavenumberMaximum > exp(self%perturbationsDarkMatter%x(-1)) &
               &  .and.                                                        &
               &   .not.self%wavenumberMaximumReached                          &
               & ) makePerturbations=.true.
       else
          makePerturbations=.true.
       end if
       if (makePerturbations) then
          ! The perturbations tables must be expanded. Call the CLASS interface to do this.
          redshifts=[0.0d0]
          call Interface_CLASS_Perturbations(                                                       &
               &                                                      self%cosmologyParameters_    , &
               &                                                           redshifts               , &
               &                                                           wavenumberMaximum       , &
               &                                                      self%wavenumberMaximum       , &
               &                                                      self%countPerDecade          , &
               &                             wavenumberMaximumReached=self%wavenumberMaximumReached, &
               &                             perturbationsDarkMatter =self%perturbationsDarkMatter , &
               &                             perturbationsBaryons    =self%perturbationsBaryons      &
               &                            )
       end if
       ! Construct the table of correlations.
       integratorDirect   =  integrator(integrandDirect,toleranceRelative=toleranceRelative)
       integratorCross    =  integrator(integrandCross ,toleranceRelative=toleranceRelative)
       self_              => self
       countPointsInTable =  int(log10(self%massMaximum/self%massMinimum)*dble(countPointsPerDecade))+1
       mass               =  Make_Range(self%massMinimum,self%massMaximum,countPointsInTable,rangeTypeLogarithmic)
       allocate(correlation(countPointsInTable))
       do i=1,countPointsInTable
          massSmoothing_      =mass(i)
          wavenumberMaximum   =max(wavenumberLarge,factorWavenumberLarge/radiusLagrangian(massSmoothing_))
          correlationDirect   =integratorDirect%integrate(0.0d0,wavenumberMaximum)
          correlationCross    =integratorCross %integrate(0.0d0,wavenumberMaximum)
          correlation      (i)=+correlationCross  &
               &               /correlationDirect
       end do
       ! Build the interpolator.
       self%correlation=interpolator(log(mass),correlation)
       self%initialized=.true.
    end if
    ! Evaluate the relative baryon fraction (eqn. 9 of Jessop et al.; 2026; arXiv:2512.02127)/
    fraction=+  1.0d0                                                           &
         &   +(                                                                 &
         &     +1.0d0                                                           &
         &     -self%fractionBaryonsUniversal                                   &
         &    )                                                                 &
         &   *self%criticalOverdensity_%value      (    basic%time         () ) &
         &   /self%linearGrowth_       %value      (    basic%time         () ) &
         &   *self%correlation         %interpolate(log(      massSmoothing  ))
    return

  contains

    double precision function radiusLagrangian(mass)
      !!{
      Compute the Lagrangian radius for a halo of the given \mono{mass}.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: mass

      radiusLagrangian=(                                             &
           &            +3.0d0                                       &
           &            /4.0d0                                       &
           &            /Pi                                          &
           &            *mass                                        &
           &            /self%cosmologyParameters_%densityCritical() &
           &            /self%cosmologyParameters_%OmegaMatter    () &
           &           )**(1.0d0/3.0d0)
      return
    end function radiusLagrangian
    
  end function isocurvatureFraction

  double precision function integrandDirect(wavenumber)
    !!{
    Integrand to compute the CDM--CDM autocorrelation:
    \begin{equation}
      \int \mathrm{d}k 4 \pi k^2 |W(k|M)|^2 \delta^2_\mathrm{c}(k).
    \end{equation}
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: wavenumber

    integrandDirect=+4.0d0                                                                             &
         &          *Pi                                                                                &
         &          *                                                   wavenumber                 **2 &
         &          *self_%powerSpectrumWindowFunction_%value      (    wavenumber ,massSmoothing_)**2 &
         &          *self_%perturbationsDarkMatter     %interpolate(log(wavenumber)               )**2
    return
  end function integrandDirect

  double precision function integrandCross(wavenumber)
    !!{
    Integrand to compute the baryon--CDM cross-correlation:
    \begin{equation}
      \int \mathrm{d}k 4 \pi k^2 |W(k|M)|^2 \delta_\mathrm{bc}(k) \delta_\mathrm{c}(k).
    \end{equation}
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: wavenumber

    integrandCross=+4.0d0                                                                               &
         &         *Pi                                                                                  &
         &         *                                                     wavenumber                 **2 &
         &         *  self_%powerSpectrumWindowFunction_%value      (    wavenumber ,massSmoothing_)**2 &
         &         *(                                                                                   &
         &           +self_%perturbationsBaryons        %interpolate(log(wavenumber)               )    &
         &           -self_%perturbationsDarkMatter     %interpolate(log(wavenumber)               )    &
         &         )                                                                                    &
         &         *  self_%perturbationsDarkMatter     %interpolate(log(wavenumber)               )
    return
  end function integrandCross
