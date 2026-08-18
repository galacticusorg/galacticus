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

program Tests_Dark_Energy_CPL
  !!{RST
  Tests the CPL (Chevallier-Polarski-Linder) dark energy equation of state class. The dark energy exponent,
  :math:`\epsilon(a)`, is checked against a direct numerical integration of :math:`\d \ln \rho_\mathrm{DE} / \d \ln a =
  -3[1+w(a)]`, its derivative is checked against a finite difference of :math:`\epsilon(a)`, and the class is checked to reduce
  to the :galacticus-class:`cosmologyFunctionsMatterDarkEnergy` class when :math:`w_\mathrm{a}=0`.
  !!}
  use :: Cosmology_Functions  , only : cosmologyFunctionsMatterDarkEnergy, cosmologyFunctionsMatterDarkEnergyCPL
  use :: Cosmology_Parameters , only : cosmologyParametersSimple
  use :: Display              , only : displayVerbositySet               , verbosityLevelStandard
  use :: Numerical_Integration, only : integrator
  use :: Unit_Tests           , only : Assert                            , Unit_Tests_Begin_Group               , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  ! Expansion factors at which the equation of state is tested. These deliberately straddle a=1, where the expressions for the
  ! dark energy exponent and its derivative both have removable singularities, and include points just inside and just outside
  ! the range over which series expansions are used in their evaluation.
  double precision                                       , dimension(9), parameter :: expansionFactor          =[1.0d-2,1.0d-1,5.0d-1,9.0d-1,0.9999d0,1.0d0,1.0001d0,1.001d0,1.01d0]
  ! Values of (w₀,w_a) tested.
  double precision                                       , dimension(4), parameter :: equationOfState0         =[-1.0d0,-0.9d0,-1.1d0,-0.7d0]
  double precision                                       , dimension(4), parameter :: equationOfStateA         =[+0.0d0,+0.3d0,-0.4d0,+0.5d0]
  double precision                                                     , parameter :: toleranceRelative        =1.0d-9
  double precision                                                     , parameter :: expansionFactorOffset    =1.0d-6
  type            (cosmologyParametersSimple            )                          :: cosmologyParameters_
  type            (cosmologyFunctionsMatterDarkEnergyCPL )                         :: cosmologyFunctionsCPL_
  type            (cosmologyFunctionsMatterDarkEnergy   )                          :: cosmologyFunctionsJassal_
  character       (len=1024                             )                          :: message
  integer                                                                          :: iExpansion               , iEquationOfState
  double precision                                                                 :: exponentActual           , exponentExpected         , &
       &                                                                              derivativeActual         , derivativeExpected       , &
       &                                                                              exponentJassal           , equationOfState
  ! Equation of state parameters made available to the integrand function by host association.
  double precision                                                                 :: equationOfState0Integrand, equationOfStateAIntegrand

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("CPL dark energy equation of state")
  ! Define cosmological parameters.
  !![
  <referenceConstruct object="cosmologyParameters_">
   <constructor>
    cosmologyParametersSimple(                             &amp;
     &amp;                    OmegaMatter         =0.30d0, &amp;
     &amp;                    OmegaBaryon         =0.05d0, &amp;
     &amp;                    OmegaDarkEnergy     =0.70d0, &amp;
     &amp;                    temperatureCMB      =2.78d0, &amp;
     &amp;                    HubbleConstant      =7.00d1  &amp;
     &amp;                   )
   </constructor>
  </referenceConstruct>
  !!]
  ! Test the equation of state, the dark energy exponent, and its derivative for each set of (w₀,w_a).
  do iEquationOfState=1,size(equationOfState0)
     cosmologyFunctionsCPL_=cosmologyFunctionsMatterDarkEnergyCPL(                                                                 &
          &                                                       cosmologyParameters_       =cosmologyParameters_               , &
          &                                                       darkEnergyEquationOfStateW0=equationOfState0(iEquationOfState) , &
          &                                                       darkEnergyEquationOfStateWA=equationOfStateA(iEquationOfState)   &
          &                                                      )
     do iExpansion=1,size(expansionFactor)
        ! The equation of state must follow the CPL form exactly.
        equationOfState=cosmologyFunctionsCPL_%equationOfStateDarkEnergy(expansionFactor=expansionFactor(iExpansion))
        write (message,'(a,f5.2,a,f5.2,a,f8.4,a)') "w(a) [w₀=",equationOfState0(iEquationOfState),", w_a=",equationOfStateA(iEquationOfState),", a=",expansionFactor(iExpansion),"]"
        call Assert(                                      &
             &      trim(message)                       , &
             &      equationOfState                     , &
             &      +equationOfState0(iEquationOfState)   &
             &      +equationOfStateA(iEquationOfState)   &
             &      *(1.0d0-expansionFactor(iExpansion)), &
             &      relTol=toleranceRelative              &
             &     )
        ! The dark energy exponent, ε(a) ≡ ln[ρ_DE(a)/ρ_DE,0]/ln(a), must match a direct numerical integration of
        ! dlnρ_DE/dlna = -3[1+w(a)]. At a=1 the exponent is not defined by that expression, but tends to -3(1+w₀).
        exponentActual=cosmologyFunctionsCPL_%exponentDarkEnergy(expansionFactor=expansionFactor(iExpansion))
        if (expansionFactor(iExpansion) == 1.0d0) then
           exponentExpected=-3.0d0*(1.0d0+equationOfState0(iEquationOfState))
        else
           exponentExpected=+exponentIntegrated(equationOfState0(iEquationOfState),equationOfStateA(iEquationOfState),expansionFactor(iExpansion)) &
                &           /log(expansionFactor(iExpansion))
        end if
        write (message,'(a,f5.2,a,f5.2,a,f8.4,a)') "ε(a) [w₀=",equationOfState0(iEquationOfState),", w_a=",equationOfStateA(iEquationOfState),", a=",expansionFactor(iExpansion),"]"
        call Assert(trim(message),exponentActual,exponentExpected,relTol=toleranceRelative)
        ! The derivative of the exponent must match a finite difference of the exponent itself.
        derivativeActual  =cosmologyFunctionsCPL_%exponentDarkEnergyDerivative(expansionFactor=expansionFactor(iExpansion))
        derivativeExpected=+(                                                                                                              &
             &               +cosmologyFunctionsCPL_%exponentDarkEnergy(expansionFactor=expansionFactor(iExpansion)+expansionFactorOffset) &
             &               -cosmologyFunctionsCPL_%exponentDarkEnergy(expansionFactor=expansionFactor(iExpansion)-expansionFactorOffset) &
             &              )                                                                                                              &
             &             /2.0d0                                                                                                          &
             &             /expansionFactorOffset
        write (message,'(a,f5.2,a,f5.2,a,f8.4,a)') "dε/da [w₀=",equationOfState0(iEquationOfState),", w_a=",equationOfStateA(iEquationOfState),", a=",expansionFactor(iExpansion),"]"
        call Assert(trim(message),derivativeActual,derivativeExpected,relTol=1.0d-5)
     end do
  end do
  call Unit_Tests_End_Group()
  ! With w_a=0 the CPL form reduces to a constant equation of state, as does the Jassal form with w₁=0, so the two classes must
  ! agree exactly.
  call Unit_Tests_Begin_Group("Reduction to constant equation of state")
  do iEquationOfState=1,size(equationOfState0)
     cosmologyFunctionsCPL_   =cosmologyFunctionsMatterDarkEnergyCPL(                                                                &
          &                                                          cosmologyParameters_       =cosmologyParameters_              , &
          &                                                          darkEnergyEquationOfStateW0=equationOfState0(iEquationOfState), &
          &                                                          darkEnergyEquationOfStateWA=0.0d0                               &
          &                                                         )
     cosmologyFunctionsJassal_=cosmologyFunctionsMatterDarkEnergy   (                                                                &
          &                                                          cosmologyParameters_       =cosmologyParameters_              , &
          &                                                          darkEnergyEquationOfStateW0=equationOfState0(iEquationOfState), &
          &                                                          darkEnergyEquationOfStateW1=0.0d0                               &
          &                                                         )
     do iExpansion=1,size(expansionFactor)
        exponentActual=cosmologyFunctionsCPL_   %exponentDarkEnergy(expansionFactor=expansionFactor(iExpansion))
        exponentJassal=cosmologyFunctionsJassal_%exponentDarkEnergy(expansionFactor=expansionFactor(iExpansion))
        write (message,'(a,f5.2,a,f8.4,a)') "ε(a) matches Jassal form [w₀=",equationOfState0(iEquationOfState),", a=",expansionFactor(iExpansion),"]"
        call Assert(trim(message),exponentActual,exponentJassal,relTol=toleranceRelative)
     end do
     write (message,'(a,f5.2,a)') "cosmic age matches Jassal form [w₀=",equationOfState0(iEquationOfState),"]"
     call Assert(                                             &
          &      trim(message)                              , &
          &      cosmologyFunctionsCPL_   %cosmicTime(0.5d0), &
          &      cosmologyFunctionsJassal_%cosmicTime(0.5d0), &
          &      relTol=toleranceRelative                     &
          &     )
  end do
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

contains

  double precision function exponentIntegrated(equationOfState0_,equationOfStateA_,expansionFactor_) result(integral)
    !!{RST
    Return :math:`\ln [\rho_\mathrm{DE}(a)/\rho_\mathrm{DE,0}] = -3 \int_1^a [1+w(a^\prime)] \d \ln a^\prime` evaluated by direct
    numerical integration for the CPL equation of state.
    !!}
    implicit none
    double precision            , intent(in   ) :: equationOfState0_, equationOfStateA_, &
         &                                         expansionFactor_
    type            (integrator)                :: integrator_

    equationOfState0Integrand=equationOfState0_
    equationOfStateAIntegrand=equationOfStateA_
    integrator_=integrator(integrandExponent,toleranceRelative=1.0d-12)
    integral   =-3.0d0*integrator_%integrate(0.0d0,log(expansionFactor_))
    return
  end function exponentIntegrated

  double precision function integrandExponent(expansionFactorLogarithmic) result(integrand)
    !!{RST
    Integrand, :math:`1+w(a)`, used in evaluating the dark energy exponent by direct numerical integration.
    !!}
    implicit none
    double precision, intent(in   ) :: expansionFactorLogarithmic

    integrand=+1.0d0                                   &
         &    +equationOfState0Integrand               &
         &    +equationOfStateAIntegrand               &
         &    *(1.0d0-exp(expansionFactorLogarithmic))
    return
  end function integrandExponent

end program Tests_Dark_Energy_CPL
