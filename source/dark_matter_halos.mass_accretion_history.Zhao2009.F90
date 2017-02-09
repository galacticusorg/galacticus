!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

  !% An implementation of dark matter halo mass accretion histories using the \cite{zhao_accurate_2009} algorithm.

  !# <darkMatterHaloMassAccretionHistory name="darkMatterHaloMassAccretionHistoryZhao2009">
  !#  <description>Dark matter halo mass accretion histories using the \cite{zhao_accurate_2009} algorithm.</description>
  !# </darkMatterHaloMassAccretionHistory>

  type, extends(darkMatterHaloMassAccretionHistoryClass) :: darkMatterHaloMassAccretionHistoryZhao2009
     !% A dark matter halo mass accretion historiy class using the \cite{zhao_accurate_2009} algorithm.
     private
   contains
     procedure :: time => zhao2009Time
  end type darkMatterHaloMassAccretionHistoryZhao2009
  
contains
  
  double precision function zhao2009Time(self,node,mass)
    !% Compute the time corresponding to {\normalfont \ttfamily mass} in the mass accretion history of {\normalfont \ttfamily
    !% thisNode} using the algorithm of \cite{zhao_accurate_2009}.
    use FGSL
    use ODE_Solver
    use Galacticus_Error
    use Cosmological_Mass_Variance
    use Critical_Overdensities
    implicit none
    class           (darkMatterHaloMassAccretionHistoryZhao2009), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    double precision                                            , intent(in   ) :: mass
    class           (criticalOverdensityClass                  ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass             ), pointer       :: cosmologicalMassVariance_
    class           (nodeComponentBasic                        ), pointer       :: baseBasicComponent
    double precision                                            , parameter     :: odeToleranceAbsolute          =1.0d-10, odeToleranceRelative =1.0d-10
    double precision                                            , dimension(1)  :: nowTime
    double precision                                                            :: baseMass                              , baseTime                     , &
       &                                                                           dSigmadMassLogarithmicObserved        , deltaCriticalObserved        , &
       &                                                                           pObserved                             , sObserved                    , &
       &                                                                           sigmaObserved                         , wObserved                    , &
       &                                                                           currentMass
    type            (fgsl_odeiv_step                           )                :: odeStepper
    type            (fgsl_odeiv_control                        )                :: odeController
    type            (fgsl_odeiv_evolve                         )                :: odeEvolver
    type            (fgsl_odeiv_system                         )                :: odeSystem
    logical                                                                     :: odeReset                      =.true.
    !GCC$ attributes unused :: self
    
    ! Get required objects.
    criticalOverdensity_      => criticalOverdensity     ()
    cosmologicalMassVariance_ => cosmologicalMassVariance()
    ! Get properties of the base node.
    baseBasicComponent => node%basic()
    baseMass=baseBasicComponent%mass()
    baseTime=baseBasicComponent%time()

    ! Trap cases where the mass occurs in the future.
    if (mass > baseMass) call Galacticus_Error_Report('zhao2009Time','specified mass is in the future')

    ! Calculate quantities which remain fixed through the ODE.

    ! Get sigma(M) and its logarithmic derivative.
    call cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(baseMass,sigmaObserved,dSigmadMassLogarithmicObserved)

    ! Compute sigma proxy.
    sObserved=sigmaObserved*10.0d0**dSigmadMassLogarithmicObserved ! Equation 8 from Zhao et al. (2009).

    ! Compute critical overdensities for collapse.
    deltaCriticalObserved=criticalOverdensity_%value(time=baseTime,mass=baseMass)

    ! Compute w factors.
    wObserved=deltaCriticalObserved/sObserved ! Equation 7 from Zhao et al. (2009).

    ! Compute p factors.
    pObserved=0.5d0*wObserved/(1.0d0+(0.25d0*wObserved)**6) ! Equation 11 from Zhao et al. (2009).

    ! Solve the ODE for the mass evolution.
    nowTime(1) =baseTime
    currentMass=baseMass
    call ODE_Solve(odeStepper,odeController,odeEvolver,odeSystem,currentMass,mass,1,nowTime &
         &,growthRateODEs,odeToleranceAbsolute,odeToleranceRelative,reset=odeReset)
    odeReset=.true.

    ! Extract the time corresponding to the specified mass.
    zhao2009Time=nowTime(1)

    return

  contains
    
    integer function growthRateODEs(mass,nowTime,dNowTimedMass)
      !% System of differential equations to solve for the growth rate.
      use Critical_Overdensities
      implicit none
      double precision              , intent(in   ) :: mass
      double precision, dimension(:), intent(in   ) :: nowTime
      double precision, dimension(:), intent(  out) :: dNowTimedMass
      double precision                              :: dDeltaCriticaldtNow      , dSigmadDeltaCriticalLogarithmic, &
           &                                                 dSigmadMassLogarithmicNow, deltaCriticalNow               , &
           &                                                 pNow                     , sNow                           , &
           &                                                 sigmaNow                 , wNow
      
      ! Trap unphysical cases.
      if (nowTime(1) <= 0.0d0 .or. nowTime(1) > baseTime .or. mass <= 0.0d0) then
         dNowTimedMass(1)=0.0d0
         growthRateODEs=FGSL_Success
         return
      end if
      
      ! Get sigma(M) and its logarithmic derivative.
      call cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(mass,sigmaNow,dSigmadMassLogarithmicNow)
      
      ! Compute sigma proxy.
      sNow=sigmaNow*10.0d0**dSigmadMassLogarithmicNow ! Equation 8 from Zhao et al. (2009).
      
      ! Compute critical overdensity for collapse.
      deltaCriticalNow   =criticalOverdensity_%value       (time=nowTime(1),mass=mass)
      dDeltaCriticaldtNow=criticalOverdensity_%gradientTime(time=nowTime(1),mass=mass)
      
      ! Compute w factor.
      wNow=deltaCriticalNow/sNow      ! Equation 7 from Zhao et al. (2009).
      
      ! Compute p factor.
      pNow=pObserved*max(0.0d0,1.0d0-log10(deltaCriticalNow/deltaCriticalObserved)*wObserved/0.272d0) ! Equation 10 from Zhao et al. (2009).
      
      ! Compute dimensionless growth rate.
      dSigmadDeltaCriticalLogarithmic=(wNow-pNow)/5.85d0 ! Equation 12 from Zhao et al. (2009).
      
      ! Convert to dimensionful units.
      dNowTimedMass(1)=(dSigmadMassLogarithmicNow/dDeltaCriticaldtNow)*(deltaCriticalNow/mass)/dSigmadDeltaCriticalLogarithmic
      
      ! Report success.
      growthRateODEs=FGSL_Success
    end function growthRateODEs
    
  end function zhao2009Time
