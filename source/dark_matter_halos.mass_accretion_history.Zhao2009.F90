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
  An implementation of dark matter halo mass accretion histories using the \cite{zhao_accurate_2009} algorithm.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Cosmology_Functions       , only : cosmologyFunctions           , cosmologyFunctionsClass
  use :: Linear_Growth             , only : linearGrowth                 , linearGrowthClass

  !![
  <darkMatterHaloMassAccretionHistory name="darkMatterHaloMassAccretionHistoryZhao2009">
   <description>
    A dark matter halo mass accretion history class which uses the algorithm given by \cite{zhao_accurate_2009} to compute mass
    accretion histories. In particular, \cite{zhao_accurate_2009} give a fitting function for the quantity $\mathrm{d} \ln
    \sigma(M)/\mathrm{d} \ln \delta_\mathrm{c}(t)$ for the dimensionless growth rate in a mass accretion history at time $t$
    and halo mass $M$. This is converted to a dimensionful growth rate using
    \begin{equation}
     {\mathrm{d} M \over \mathrm{d} t} = \left({\mathrm{d} \ln \sigma(M) \over \mathrm{d} \ln M}\right)^{-1} \left({\mathrm{d}
     \delta_c(t) \over \mathrm{d} t}\right) \left( {M \over \delta_\mathrm{c}(t)} \right) \left({\mathrm{d} \ln \sigma(M) \over
     \mathrm{d} \ln \delta_\mathrm{c}(t)}\right).
    \end{equation}
    This differential equation is then solved numerically to find the mass accretion history.
   </description>
  </darkMatterHaloMassAccretionHistory>
  !!]
  type, extends(darkMatterHaloMassAccretionHistoryClass) :: darkMatterHaloMassAccretionHistoryZhao2009
     !!{
     A dark matter halo mass accretion historiy class using the \cite{zhao_accurate_2009} algorithm.
     !!}
     private
     class(criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class(cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class(linearGrowthClass            ), pointer :: linearGrowth_             => null()
     class(cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
   contains
     final     ::                      zhao2009Destructor
     procedure :: time              => zhao2009Time
     procedure :: massAccretionRate => zhao2009MassAccretionRate
  end type darkMatterHaloMassAccretionHistoryZhao2009

  interface darkMatterHaloMassAccretionHistoryZhao2009
     !!{
     Constructors for the \refClass{darkMatterHaloMassAccretionHistoryZhao2009} dark matter halo mass accretion history class.
     !!}
     module procedure zhao2009ConstructorParameters
     module procedure zhao2009ConstructorInternal
  end interface darkMatterHaloMassAccretionHistoryZhao2009

contains

  function zhao2009ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterHaloMassAccretionHistoryZhao2009} dark matter halo mass accretion history class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterHaloMassAccretionHistoryZhao2009)                :: self
    type (inputParameters                           ), intent(inout) :: parameters
    class(criticalOverdensityClass                  ), pointer       :: criticalOverdensity_
    class(cosmologicalMassVarianceClass             ), pointer       :: cosmologicalMassVariance_
    class(linearGrowthClass                         ), pointer       :: linearGrowth_
    class(cosmologyFunctionsClass                   ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    !!]
    self=darkMatterHaloMassAccretionHistoryZhao2009(criticalOverdensity_,cosmologicalMassVariance_,linearGrowth_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="linearGrowth_"            />
    <objectDestructor name="cosmologyFunctions_"      />
    !!]
    return
  end function zhao2009ConstructorParameters

  function zhao2009ConstructorInternal(criticalOverdensity_,cosmologicalMassVariance_,linearGrowth_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterHaloMassAccretionHistoryZhao2009} dark matter halo mass accretion history class.
    !!}
    implicit none
    type (darkMatterHaloMassAccretionHistoryZhao2009)                        :: self
    class(criticalOverdensityClass                  ), intent(in   ), target :: criticalOverdensity_
    class(cosmologicalMassVarianceClass             ), intent(in   ), target :: cosmologicalMassVariance_
    class(linearGrowthClass                         ), intent(in   ), target :: linearGrowth_
    class(cosmologyFunctionsClass                   ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="*criticalOverdensity_, *cosmologicalMassVariance_, *linearGrowth_, *cosmologyFunctions_"/>
    !!]

    return
  end function zhao2009ConstructorInternal

  subroutine zhao2009Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterHaloMassAccretionHistoryZhao2009} dark matter halo mass accretion history class.
    !!}
    implicit none
    type(darkMatterHaloMassAccretionHistoryZhao2009), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%linearGrowth_"            />
    <objectDestructor name="self%cosmologyFunctions_"      />
    !!]
    return
  end subroutine zhao2009Destructor

  double precision function zhao2009Time(self,node,mass)
    !!{
    Compute the time corresponding to {\normalfont \ttfamily mass} in the mass accretion history of {\normalfont \ttfamily
    node} using the algorithm of \cite{zhao_accurate_2009}.
    !!}
    use :: Error                , only : Error_Report
    use :: Galacticus_Nodes     , only : nodeComponentBasic, treeNode
    use :: Interface_GSL        , only : GSL_Success
    use :: Numerical_ODE_Solvers, only : odeSolver
    implicit none
    class           (darkMatterHaloMassAccretionHistoryZhao2009), intent(inout), target :: self
    type            (treeNode                                  ), intent(inout), target :: node
    double precision                                            , intent(in   )         :: mass
    class           (nodeComponentBasic                        ), pointer               :: basicBase
    double precision                                            , parameter             :: odeToleranceAbsolute          =1.0d-10, odeToleranceRelative =1.0d-10
    double precision                                            , dimension(1)          :: nowTime
    double precision                                                                    :: baseMass                              , baseTime                     , &
       &                                                                                   dSigmadMassLogarithmicObserved        , deltaCriticalObserved        , &
       &                                                                                   pObserved                             , sObserved                    , &
       &                                                                                   sigmaObserved                         , wObserved                    , &
       &                                                                                   currentMass
    type            (odeSolver                                 )                        :: solver

    ! Get properties of the base node.
    basicBase => node     %basic()
    baseMass  =  basicBase%mass ()
    baseTime  =  basicBase%time ()
    ! Trap cases where the mass occurs in the future.
    if (mass > baseMass) call Error_Report('specified mass is in the future'//{introspection:location})
    ! Calculate quantities which remain fixed through the ODE.
    ! Get σ(M) and its logarithmic derivative.
    call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(baseMass,baseTime,sigmaObserved,dSigmadMassLogarithmicObserved)
    ! Compute σ proxy.
    sObserved=sigmaObserved*10.0d0**dSigmadMassLogarithmicObserved ! Equation 8 from Zhao et al. (2009).
    ! Compute critical overdensities for collapse.
    deltaCriticalObserved=self%criticalOverdensity_%value(time=baseTime,mass=baseMass)
    ! Compute w factors.
    wObserved=deltaCriticalObserved/sObserved ! Equation 7 from Zhao et al. (2009).
    ! Compute p factors.
    pObserved=0.5d0*wObserved/(1.0d0+(0.25d0*wObserved)**6) ! Equation 11 from Zhao et al. (2009).
    ! Solve the ODE for the mass evolution.
    nowTime(1) =baseTime
    currentMass=baseMass
    solver     =odeSolver(1_c_size_t,growthRateODEs,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative)    
    call solver%solve(currentMass,mass,nowTime)
    ! Extract the time corresponding to the specified mass.
    zhao2009Time=nowTime(1)
    return

  contains

    integer function growthRateODEs(mass,nowTime,dNowTimedMass)
      !!{
      System of differential equations to solve for the growth rate.
      !!}
      implicit none
      double precision              , intent(in   ) :: mass
      double precision, dimension(:), intent(in   ) :: nowTime
      double precision, dimension(:), intent(  out) :: dNowTimedMass
      double precision                              :: dDeltaCriticaldtNow      , dSigmadDeltaCriticalLogarithmic, &
           &                                           dSigmadMassLogarithmicNow, deltaCriticalNow               , &
           &                                           pNow                     , sNow                           , &
           &                                           sigmaNow                 , wNow

      ! Trap unphysical cases.
      if (nowTime(1) <= 0.0d0 .or. nowTime(1) > baseTime .or. mass <= 0.0d0) then
         dNowTimedMass(1)=0.0d0
         growthRateODEs=GSL_Success
         return
      end if
      ! Get σ(M) and its logarithmic derivative.
      call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(mass,nowTime(1),sigmaNow,dSigmadMassLogarithmicNow)
      ! Compute σ proxy.
      sNow               =+sigmaNow                                                                                                                                          &
           &              *10.0d0**dSigmadMassLogarithmicNow                                                                                                                 &
           &              /                                                                  self%linearGrowth_      %value                               ( time=nowTime(1))   ! Equation 8 from Zhao et al. (2009).
      ! Compute critical overdensity for collapse.
      deltaCriticalNow   =+self%criticalOverdensity_%value       (time=nowTime(1),mass=mass)/self%linearGrowth_      %value                               ( time=nowTime(1))
      dDeltaCriticaldtNow=+self%criticalOverdensity_%gradientTime(time=nowTime(1),mass=mass)/self%linearGrowth_      %value                               ( time=nowTime(1)) &
           &              -self%criticalOverdensity_%value       (time=nowTime(1),mass=mass)/self%linearGrowth_      %value                               ( time=nowTime(1)) &
           &                                                                                *self%linearGrowth_      %logarithmicDerivativeExpansionFactor( time=nowTime(1)) &
           &                                                                                *self%cosmologyFunctions_%expansionRate                       (                  &
           &                                                                                 self%cosmologyFunctions_%expansionFactor                      (time=nowTime(1)) &
           &                                                                                                                                              )
      ! Compute w factor.
      wNow=deltaCriticalNow/sNow      ! Equation 7 from Zhao et al. (2009).
      ! Compute p factor.
      pNow=pObserved*max(0.0d0,1.0d0-log10(deltaCriticalNow/deltaCriticalObserved)*wObserved/0.272d0) ! Equation 10 from Zhao et al. (2009).
      ! Compute dimensionless growth rate.
      dSigmadDeltaCriticalLogarithmic=(wNow-pNow)/5.85d0 ! Equation 12 from Zhao et al. (2009).
      ! Convert to dimensionful units.
      dNowTimedMass(1)=(dSigmadMassLogarithmicNow/dDeltaCriticaldtNow)*(deltaCriticalNow/mass)/dSigmadDeltaCriticalLogarithmic
      ! Report success.
      growthRateODEs=GSL_Success
    end function growthRateODEs

  end function zhao2009Time

  double precision function zhao2009MassAccretionRate(self,node,time)
    !!{
    Compute the mass accretion rate at the given {\normalfont \ttfamily mass} in the mass accretion history of {\normalfont
    \ttfamily node} using the algorithm of \cite{zhao_accurate_2009}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (darkMatterHaloMassAccretionHistoryZhao2009), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    double precision                                            , intent(in   ) :: time
    !$GLC attributes unused :: self, node, time

    zhao2009MassAccretionRate=0.0d0
    call Error_Report('mass accretion rate is not implemented'//{introspection:location})
    return
  end function zhao2009MassAccretionRate
