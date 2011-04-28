!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements the \cite{zhao_accurate_2009} halo mass accretion algorithm.

module Dark_Matter_Halo_Mass_Accretion_Histories_Zhao2009
  !% Implements the \cite{zhao_accurate_2009} halo mass accretion algorithm.
  use Tree_Nodes
  use, intrinsic :: ISO_C_Binding
  use FGSL
  private
  public :: Dark_Matter_Mass_Accretion_Zhao2009_Initialize

  ! Module global variables used in ODE solving.
  double precision         :: baseMass,baseTime,sigmaObserved,dSigmadMassLogarithmicObserved,sObserved,wObserved,pObserved&
       &,deltaCriticalObserved
  !$omp threadprivate (baseMass,baseTime,sigmaObserved,dSigmadMassLogarithmicObserved,sObserved,wObserved,pObserved)
  !$omp threadprivate (deltaCriticalObserved)

  ! Variables used in the ODE solver.
  type(fgsl_odeiv_step)    :: odeStepper
  type(fgsl_odeiv_control) :: odeController
  type(fgsl_odeiv_evolve)  :: odeEvolver
  type(fgsl_odeiv_system)  :: odeSystem
  logical                  :: odeReset=.true. ! Ensure ODE variables will be reset on first call.
  !$omp threadprivate(odeStepper,odeController,odeEvolver,odeSystem,odeReset)

contains

  !# <darkMatterAccretionHistoryMethod>
  !#  <unitName>Dark_Matter_Mass_Accretion_Zhao2009_Initialize</unitName>
  !# </darkMatterAccretionHistoryMethod>
  subroutine Dark_Matter_Mass_Accretion_Zhao2009_Initialize(darkMatterAccretionHistoryMethod,Dark_Matter_Halo_Mass_Accretion_Time_Get)
    !% Initializes the ``Zhao 2009'' mass accretion history module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: darkMatterAccretionHistoryMethod
    procedure(double precision), pointer, intent(inout) :: Dark_Matter_Halo_Mass_Accretion_Time_Get
    
    if (darkMatterAccretionHistoryMethod == 'Zhao 2009') then
       ! Set procedure pointers.
       Dark_Matter_Halo_Mass_Accretion_Time_Get => Dark_Matter_Halo_Mass_Accretion_Time_Zhao2009       
    end if

    return
  end subroutine Dark_Matter_Mass_Accretion_Zhao2009_Initialize

  double precision function Dark_Matter_Halo_Mass_Accretion_Time_Zhao2009(baseNode,nodeMass)
    !% Compute the time corresponding to {\tt nodeMass} in the mass accretion history of {\tt thisNode} using the algorithm of
    !% \cite{zhao_accurate_2009}.
    use Tree_Nodes
    use ODE_Solver
    use Galacticus_Error
    use CDM_Power_Spectrum
    use Critical_Overdensity
    implicit none
    type(treeNode),   intent(inout), pointer :: baseNode
    double precision, intent(in)             :: nodeMass
    double precision, parameter              :: odeToleranceAbsolute=1.0d-10, odeToleranceRelative=1.0d-10
    double precision, dimension(1)           :: nowTime
    double precision                         :: currentMass
    type(c_ptr)                              :: parameterPointer

    ! Get properties of the base node.
    baseMass=Tree_Node_Mass(baseNode)
    baseTime=Tree_Node_Time(baseNode)

    ! Trap cases where the mass occurs in the future.
    if (nodeMass > baseMass) call Galacticus_Error_Report('Dark_Matter_Halo_Mass_Accretion_Time_Zhao2009','specified mass is in the future')

    ! Calculate quantities which remain fixed through the ODE.

    ! Get sigma(M) and its logarithmic derivative.
    call sigma_CDM_Plus_Logarithmic_Derivative(baseMass,sigmaObserved,dSigmadMassLogarithmicObserved)

    ! Compute sigma proxy.
    sObserved=sigmaObserved*10.0d0**dSigmadMassLogarithmicObserved ! Equation 8 from Zhao et al. (2009).

    ! Compute critical overdensities for collapse.
    deltaCriticalObserved=Critical_Overdensity_for_Collapse(baseTime)

    ! Compute w factors.
    wObserved=deltaCriticalObserved/sObserved ! Equation 7 from Zhao et al. (2009).

    ! Compute p factors.
    pObserved=0.5d0*wObserved/(1.0d0+(0.25d0*wObserved)**6) ! Equation 11 from Zhao et al. (2009).

    ! Solve the ODE for the mass evolution.
    nowTime(1) =baseTime
    currentMass=baseMass
    call ODE_Solve(odeStepper,odeController,odeEvolver,odeSystem,currentMass,nodeMass,1,nowTime &
         &,growthRateODEs,parameterPointer,odeToleranceAbsolute,odeToleranceRelative,reset=odeReset)
    odeReset=.true.

    ! Extract the time corresponding to the specified mass.    
    Dark_Matter_Halo_Mass_Accretion_Time_Zhao2009=nowTime(1)

    return
  end function Dark_Matter_Halo_Mass_Accretion_Time_Zhao2009

  function growthRateODEs(mass,nowTime,dNowTimedMass,parameterPointer) bind(c)
    !% System of differential equations to solve for the growth rate.
    use CDM_Power_Spectrum
    use Critical_Overdensity
    implicit none
    integer(c_int)                           :: growthRateODEs
    real(c_double), value                    :: mass
    real(c_double), dimension(1), intent(in) :: nowTime
    real(c_double), dimension(1)             :: dNowTimedMass
    type(c_ptr),    value                    :: parameterPointer
    real(c_double)                           :: sigmaNow,dSigmadMassLogarithmicNow,deltaCriticalNow&
         &,dDeltaCriticaldtNow,wNow,pNow ,dSigmadDeltaCriticalLogarithmic,sNow

    ! Trap unphysical cases.
    if (nowTime(1) <= 0.0d0 .or. nowTime(1) > baseTime .or. mass <= 0.0d0) then
       dNowTimedMass(1)=0.0d0
       growthRateODEs=FGSL_Success
      return
    end if

    ! Get sigma(M) and its logarithmic derivative.
    call sigma_CDM_Plus_Logarithmic_Derivative(mass,sigmaNow,dSigmadMassLogarithmicNow)

    ! Compute sigma proxy.
    sNow=sigmaNow*10.0d0**dSigmadMassLogarithmicNow ! Equation 8 from Zhao et al. (2009).

    ! Compute critical overdensity for collapse.
    deltaCriticalNow   =Critical_Overdensity_for_Collapse              (nowTime(1))
    dDeltaCriticaldtNow=Critical_Overdensity_for_Collapse_Time_Gradient(nowTime(1))

    ! Compute w factor.
    wNow=deltaCriticalNow/sNow      ! Equation 7 from Zhao et al. (2009).

    ! Compute p factor.
    pNow=pObserved*max(0.0d0,1.0d0-dlog10(deltaCriticalNow/deltaCriticalObserved)*wObserved/0.272d0) ! Equation 10 from Zhao et al. (2009).

    ! Compute dimensionless growth rate.
    dSigmadDeltaCriticalLogarithmic=(wNow-pNow)/5.85d0 ! Equation 12 from Zhao et al. (2009).

    ! Convert to dimensionful units.
    dNowTimedMass(1)=(dSigmadMassLogarithmicNow/dDeltaCriticaldtNow)*(deltaCriticalNow/mass)/dSigmadDeltaCriticalLogarithmic

    ! Report success.
    growthRateODEs=FGSL_Success
  end function growthRateODEs
  
end module Dark_Matter_Halo_Mass_Accretion_Histories_Zhao2009
