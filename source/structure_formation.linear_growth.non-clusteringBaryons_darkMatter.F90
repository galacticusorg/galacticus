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
  An implementation of linear growth of cosmological structure in the limit where baryons do not cluster (i.e. small scales),
  and so has no wavenumber dependence. Also assumes no growth of radiation perturbations.
  !!}

  !![
  <linearGrowth name="linearGrowthNonClusteringBaryonsDarkMatter">
   <description>Linear growth of cosmological structure in the limit where baryons do not cluster (i.e. small scales), and so has no wavenumber dependence. Also assumes no growth of radiation perturbations.</description>
  </linearGrowth>
  !!]
  use :: Cosmology_Functions , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters, only : cosmologyParametersClass, hubbleUnitsTime
  use :: Tables              , only : table1D

  type, extends(linearGrowthClass) :: linearGrowthNonClusteringBaryonsDarkMatter
     !!{
     A linear growth of cosmological structure contrast class in the limit where baryons do not cluster (i.e. small scales),
     and so has no wavenumber dependence. Also assumes no growth of radiation perturbations.
     !!}
     private
     logical                                                 :: tableInitialized             =  .false.
     double precision                                        :: tableTimeMinimum                       , tableTimeMaximum, &
          &                                                     normalizationMatterDominated
     class           (table1D                 ), allocatable :: growthFactor
     class           (cosmologyParametersClass), pointer     :: cosmologyParameters_         => null()
     class           (cosmologyFunctionsClass ), pointer     :: cosmologyFunctions_          => null()
   contains
     !![
     <methods>
       <method description="Tabulate linear growth factor." method="retabulate" />
     </methods>
     !!]
     final     ::                                         nonClusteringBaryonsDarkMatterDestructor
     procedure :: value                                => nonClusteringBaryonsDarkMatterValue
     procedure :: logarithmicDerivativeExpansionFactor => nonClusteringBaryonsDarkMatterLogDerivativeExpansionFactor
     procedure :: logarithmicDerivativeWavenumber      => nonClusteringBaryonsDarkMatterLogDerivativeWavenumber
     procedure :: retabulate                           => nonClusteringBaryonsDarkMatterRetabulate
     procedure :: isWavenumberDependent                => nonClusteringBaryonsDarkMatterIsWavenumberDependent
  end type linearGrowthNonClusteringBaryonsDarkMatter

  interface linearGrowthNonClusteringBaryonsDarkMatter
     !!{
     Constructors for the \refClass{linearGrowthNonClusteringBaryonsDarkMatter} linear growth class.
     !!}
     module procedure nonClusteringBaryonsDarkMatterConstructorParameters
     module procedure nonClusteringBaryonsDarkMatterConstructorInternal
  end interface linearGrowthNonClusteringBaryonsDarkMatter

  ! Tolerance parameter used to ensure times do not exceed that at the Big Crunch.
  double precision, parameter :: timeToleranceRelative=1.0d-4

contains

  function nonClusteringBaryonsDarkMatterConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{linearGrowthNonClusteringBaryonsDarkMatter} linear growth class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (linearGrowthNonClusteringBaryonsDarkMatter)                :: self
    type (inputParameters                           ), intent(inout) :: parameters
    class(cosmologyParametersClass                  ), pointer       :: cosmologyParameters_
    class(cosmologyFunctionsClass                   ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    self=nonClusteringBaryonsDarkMatterConstructorInternal(cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function nonClusteringBaryonsDarkMatterConstructorParameters

  function nonClusteringBaryonsDarkMatterConstructorInternal(cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{linearGrowthNonClusteringBaryonsDarkMatter} linear growth class.
    !!}
    implicit none
    type            (linearGrowthNonClusteringBaryonsDarkMatter)                           :: self
    class           (cosmologyParametersClass                  ), target   , intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass                   ), target   , intent(in   ) :: cosmologyFunctions_
    double precision                                                                       :: timeBigCrunch
    !![
    <constructorAssign variables="*cosmologyParameters_, *cosmologyFunctions_"/>
    !!]

    self%tableInitialized=.false.
    self%tableTimeMinimum= 1.0d0
    self%tableTimeMaximum=20.0d0
    timeBigCrunch        =self%cosmologyFunctions_%timeBigCrunch()
    if (timeBigCrunch > 0.0d0) then
       ! A Big Crunch exists - avoid attempting to tabulate times beyond this epoch.
       if (self%tableTimeMinimum > timeBigCrunch) self%tableTimeMinimum= 0.5d0                       *timeBigCrunch
       if (self%tableTimeMaximum > timeBigCrunch) self%tableTimeMaximum=(1.0d0-timeToleranceRelative)*timeBigCrunch
    end if
    return
  end function nonClusteringBaryonsDarkMatterConstructorInternal

  subroutine nonClusteringBaryonsDarkMatterDestructor(self)
    !!{
    Destructor for the \refClass{linearGrowthNonClusteringBaryonsDarkMatter} linear growth class.
    !!}
    implicit none
    type (linearGrowthNonClusteringBaryonsDarkMatter), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    if (self%tableInitialized) then
       call self%growthFactor%destroy()
       deallocate(self%growthFactor)
    end if
    return
  end subroutine nonClusteringBaryonsDarkMatterDestructor

  subroutine nonClusteringBaryonsDarkMatterRetabulate(self,time)
    !!{
    Returns the linear growth factor $D(a)$ for expansion factor {\normalfont \ttfamily aExpansion}, normalized such that
    $D(1)=1$ for a nonClusteringBaryonsDarkMatter matter plus cosmological constant cosmology.
    !!}
    use :: Interface_GSL        , only : GSL_Success
    use :: Numerical_ODE_Solvers, only : odeSolver
    use :: Tables               , only : table1DLogarithmicLinear
    implicit none
    class           (linearGrowthNonClusteringBaryonsDarkMatter), intent(inout) :: self
    double precision                                            , intent(in   ) :: time
    double precision                                            , parameter     :: dominateFactor               =   1.0d+04
    double precision                                            , parameter     :: odeToleranceAbsolute         =   1.0d-10, odeToleranceRelative     =1.0d-10
    integer                                                     , parameter     :: growthTablePointsPerDecade   =1000
    double precision                                            , dimension(2)  :: growthFactorODEVariables
    logical                                                                     :: remakeTable
    integer                                                                     :: i
    double precision                                                            :: expansionFactorMatterDominant           , growthFactorDerivative           , &
         &                                                                         timeNow                                 , linearGrowthFactorPresent        , &
         &                                                                         timeMatterDominant                      , timePresent                      , &
         &                                                                         timeBigCrunch, exponent
    integer                                                                     :: growthTableNumberPoints
    type            (odeSolver                                 )                :: solver

    ! Check if we need to recompute our table.
    if (self%tableInitialized) then
       remakeTable=(                              &
            &        time < self%tableTimeMinimum &
            &       .or.                          &
            &        time > self%tableTimeMaximum &
            &      )
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       ! Find the present-day epoch.
       timePresent                 =     self%cosmologyFunctions_%cosmicTime           (                        1.0d0,collapsingPhase=self%cosmologyParameters_%HubbleConstant() < 0.0d0)
       ! Find epoch of matter-dark energy equality.
       expansionFactorMatterDominant=min(                                                                               &
            &                            self%cosmologyFunctions_%expansionFactor      (        self%tableTimeMinimum), &
            &                            self%cosmologyFunctions_%dominationEpochMatter(               dominateFactor)  &
            &                           )
       timeMatterDominant           =    self%cosmologyFunctions_%cosmicTime           (expansionFactorMatterDominant)
       ! Find minimum and maximum times to tabulate.
       self%tableTimeMinimum=min(self%tableTimeMinimum,min(timePresent,min(time/2.0,timeMatterDominant)      ))
       self%tableTimeMaximum=max(self%tableTimeMaximum,max(timePresent,max(time    ,timeMatterDominant)*2.0d0))
       timeBigCrunch        =self%cosmologyFunctions_%timeBigCrunch()
       if (timeBigCrunch > 0.0d0) then
          ! A Big Crunch exists - avoid attempting to tabulate times beyond this epoch.
          if (self%tableTimeMinimum > timeBigCrunch) self%tableTimeMinimum= 0.5d0                       *timeBigCrunch
          if (self%tableTimeMaximum > timeBigCrunch) self%tableTimeMaximum=(1.0d0-timeToleranceRelative)*timeBigCrunch
       end if
       ! Determine number of points to tabulate.
       growthTableNumberPoints=int(log10(self%tableTimeMaximum/self%tableTimeMinimum)*dble(growthTablePointsPerDecade))
       ! Destroy current table.
       if (allocated(self%growthFactor)) then
          call self%growthFactor%destroy()
          deallocate(self%growthFactor)
       end if
       ! Compute the exponent of the growth factor during the matter-dominated phase.
       exponent=(sqrt(1.0d0+24.0d0*(self%cosmologyParameters_%OmegaMatter()-self%cosmologyParameters_%OmegaBaryon())/self%cosmologyParameters_%OmegaMatter())-1.0d0)/6.0d0
       ! Create table.
       allocate(table1DLogarithmicLinear :: self%growthFactor)
       select type (growthFactor => self%growthFactor)
       type is (table1DLogarithmicLinear)
          call growthFactor%create(self%tableTimeMinimum,self%tableTimeMaximum,growthTableNumberPoints)
          ! Solve ODE to get corresponding expansion factors. Initialize with solution for matter dominated phase. We choose an
          ! initial growth factor of 1 (this is arbitrary as by definition the growth is linear so we can rescale to any
          ! value). Perturbations in the matter dominated phase grow as δ ∝ t^p, so the initial growth rate is dδ/dt = p δ/t.
          call growthFactor%populate(1.0d0,1)
          growthFactorDerivative=exponent/growthFactor%x(1)
          solver                =odeSolver(2_c_size_t,growthFactorODEs,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative)    
          do i=2,growthTableNumberPoints
             timeNow                    =growthFactor          %x(i-1)
             growthFactorODEVariables(1)=growthFactor          %y(i-1)
             growthFactorODEVariables(2)=growthFactorDerivative
             call solver      %solve   (timeNow,growthFactor%x(i),growthFactorODEVariables     )
             call growthFactor%populate(                          growthFactorODEVariables(1),i)
             growthFactorDerivative=growthFactorODEVariables(2)
          end do
          ! Normalize to growth factor of unity at present day.
          linearGrowthFactorPresent=growthFactor%interpolate(timePresent)
          call growthFactor%populate(reshape(growthFactor%ys(),[growthTableNumberPoints])/linearGrowthFactorPresent)
          ! Compute relative normalization factor such that growth factor behaves as expansion factor at early times.
          self%normalizationMatterDominated=+(                                                                &
               &                              +9.0d0                                                          &
               &                              *    self%cosmologyParameters_%OmegaMatter   (               )  &
               &                              /4.0d0                                                          &
               &                             )**(exponent/2.0d0)                                              &
               &                            *(                                                                &
               &                              +abs(self%cosmologyParameters_%HubbleConstant(hubbleUnitsTime)) &
               &                              *    growthFactor             %x             (              1)  &
               &                             )**exponent                                                      &
               &                            /      growthFactor             %y             (              1)
          self%tableInitialized=.true.
       end select
    end if
    return

  contains

    integer function growthFactorODEs(time,values,derivatives)
      !!{
      System of differential equations to solve for the growth factor.
      !!}
      double precision              , intent(in   ) :: time
      double precision, dimension(:), intent(in   ) :: values
      double precision, dimension(:), intent(  out) :: derivatives
      double precision                              :: expansionFactor

      expansionFactor   =+    self%cosmologyFunctions_%expansionFactor   (                           time)
      derivatives    (1)=+    values                                     (                              2)
      derivatives    (2)=+    1.5d0                                                                           &
           &             *    self%cosmologyFunctions_%expansionRate     (                expansionFactor)**2 &
           &             *    self%cosmologyFunctions_%omegaMatterEpochal(expansionFactor=expansionFactor)    &
           &             *    values                                     (                              1)    &
           &             *(                                                                                   &
           &              +(                                                                                  &
           &                + self%cosmologyParameters_%OmegaMatter      (                               )    &
           &                - self%cosmologyParameters_%OmegaBaryon      (                               )    &
           &             )                                                                                    &
           &             /    self%cosmologyParameters_%OmegaMatter      (                               )    &
           &             )                                                                                    &
           &             -    2.0d0                                                                           &
           &             *abs(self%cosmologyFunctions_%expansionRate     (                expansionFactor))   &
           &             *    values                                     (                              2)
      growthFactorODEs  = GSL_Success
    end function growthFactorODEs

  end subroutine nonClusteringBaryonsDarkMatterRetabulate

  double precision function nonClusteringBaryonsDarkMatterValue(self,time,expansionFactor,collapsing,normalize,component,wavenumber)
    !!{
    Return the linear growth factor at the given epoch.
    !!}
    implicit none
    class           (linearGrowthNonClusteringBaryonsDarkMatter), intent(inout)           :: self
    double precision                                            , intent(in   ), optional :: time      , expansionFactor
    logical                                                     , intent(in   ), optional :: collapsing
    type            (enumerationNormalizeType                  ), intent(in   ), optional :: normalize
    type            (enumerationComponentType                  ), intent(in   ), optional :: component
    double precision                                            , intent(in   ), optional :: wavenumber
    double precision                                                                      :: time_
    !![
    <optionalArgument name="normalize" defaultsTo="normalizePresentDay" />
    !!]
    !$GLC attributes unused :: component, wavenumber

    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Remake the table if necessary.
    call self%retabulate(time_)
    ! Interpolate to get the expansion factor.
    nonClusteringBaryonsDarkMatterValue=self%growthFactor%interpolate(time_)
    ! Normalize.
    select case (normalize_%ID)
    case (normalizeMatterDominated%ID)
       nonClusteringBaryonsDarkMatterValue=nonClusteringBaryonsDarkMatterValue*self%normalizationMatterDominated
    end select
    return
  end function nonClusteringBaryonsDarkMatterValue

  double precision function nonClusteringBaryonsDarkMatterLogDerivativeExpansionFactor(self,time,expansionFactor,collapsing,component,wavenumber)
    !!{
    Return the logarithmic gradient of linear growth factor with respect to expansion factor at the given epoch.
    !!}
    implicit none
    class           (linearGrowthNonClusteringBaryonsDarkMatter), intent(inout)           :: self
    double precision                                            , intent(in   ), optional :: time      , expansionFactor
    logical                                                     , intent(in   ), optional :: collapsing
    type            (enumerationComponentType                  ), intent(in   ), optional :: component
    double precision                                            , intent(in   ), optional :: wavenumber
    double precision                                                                      :: time_     , expansionFactor_
    !$GLC attributes unused :: component, wavenumber

    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_,expansionFactorOut=expansionFactor_)
    ! Remake the table if necessary.
    call self%retabulate(time_)
    ! Interpolate to get the expansion factor.
    nonClusteringBaryonsDarkMatterLogDerivativeExpansionFactor=+self%growthFactor       %interpolateGradient(time_           ) &
         &                                                     /self%growthFactor       %interpolate        (time_           ) &
         &                                                     /self%cosmologyFunctions_%expansionRate      (expansionFactor_)
    return
  end function nonClusteringBaryonsDarkMatterLogDerivativeExpansionFactor

  double precision function nonClusteringBaryonsDarkMatterLogDerivativeWavenumber(self,time,expansionFactor,collapsing,component,wavenumber)
    !!{
    Return the logarithmic gradient of linear growth factor with respect to wavenumber at the given epoch.
    !!}
    implicit none
    class           (linearGrowthNonClusteringBaryonsDarkMatter), intent(inout)           :: self
    double precision                                            , intent(in   ), optional :: time      , expansionFactor
    logical                                                     , intent(in   ), optional :: collapsing
    type            (enumerationComponentType                  ), intent(in   ), optional :: component
    double precision                                            , intent(in   ), optional :: wavenumber
    !$GLC attributes unused :: self, time, expansionFactor, collapsing, component, wavenumber

    ! No dependence on wavenumber.
    nonClusteringBaryonsDarkMatterLogDerivativeWavenumber=0.0d0
    return
  end function nonClusteringBaryonsDarkMatterLogDerivativeWavenumber

  logical function nonClusteringBaryonsDarkMatterIsWavenumberDependent(self,component)
    !!{
    Return false indicating that the growth function is not wavenumber-dependent.
    !!}
    implicit none
    class(linearGrowthNonClusteringBaryonsDarkMatter), intent(inout)           :: self
    type (enumerationComponentType                  ), intent(in   ), optional :: component
    !$GLC attributes unused :: self, component

    nonClusteringBaryonsDarkMatterIsWavenumberDependent=.false.
    return
  end function nonClusteringBaryonsDarkMatterIsWavenumberDependent
