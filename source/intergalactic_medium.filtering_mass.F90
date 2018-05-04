!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which provide utility functions related to calculations of the filtering mass.

module Intergalactic_Medium_Filtering_Masses
  !% Provides utility functions related to calculations of the filtering mass.
  private
  public :: Mass_Filtering_Early_Epoch, Mass_Filtering_ODE_Initial_Conditions, &
       &    Mass_Filtering_ODE_System
  
contains

  function Mass_Filtering_Early_Epoch(time,cosmologyParameters_,cosmologyFunctions_) result (massFiltering)
    !% Fitting function for the filtering mass at early epochs from \cite{naoz_formation_2007}. Checks for valid range of redshift
    !% and cosmology for the fit to be valid.
    use Cosmology_Parameters
    use Cosmology_Functions
    implicit none
    double precision                                          :: massFiltering
    double precision                          , intent(in   ) :: time
    class           (cosmologyParametersClass), intent(inout) :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), intent(inout) :: cosmologyFunctions_
    double precision                          , dimension(4)  :: coefficients
    double precision                                          :: expansionFactor
    
    ! Compute the expansion factor.
    expansionFactor=cosmologyFunctions_%expansionFactor    (time                                         )
    ! Get fitting function coefficients.
    coefficients   =Mass_Filtering_Early_Epoch_Coefficients(time,cosmologyParameters_,cosmologyFunctions_)
    ! Evaluate fitting function.
    massFiltering  =+exp(                                            &
         &               +coefficients(1)*(-log(expansionFactor))**3 &
         &               +coefficients(2)*(-log(expansionFactor))**2 &
         &               +coefficients(3)*(-log(expansionFactor))    &
         &               +coefficients(4)                            &
         &              )
    ! Compute the 
    return
  end function Mass_Filtering_Early_Epoch
  
  subroutine Mass_Filtering_ODE_Initial_Conditions(time,cosmologyParameters_,cosmologyFunctions_,linearGrowth_,massFilteringODEs,massFilteringScales)
    !% Compute initial conditions for a system of three variables used to solve for the evolution of the filtering mass. The ODE system to be solved is
    !% \begin{eqnarray}
    !%  \dot{y}_1 &=& y_2 \\
    !%  \dot{y}_2 &=& -2 (\dot{a}/a) D(t) (1+r_{\rm LSS}(t)) y_2 f_{\rm DM} {\rm k}_{\rm B} T(t)/\mu m_{\rm H} a^2 \\
    !%  \dot{y}_3 &=& 4 \pi^4 \bar{\rho}(t) \dot{k}_{\rm F}(t)/ k_{\rm F}^4(t)
    !% \end{eqnarray}
    !% with initial conditions
    !% \begin{eqnarray}
    !%  y_1 &=& D(t)/k_{\rm F}^2(t) \\
    !%  y_2 &=& \dot{y}_1 \\
    !%  y_3 &=& M_{\rm F}(t)
    !% \end{eqnarray}
    !% and where
    !% \begin{equation}
    !%  k_{\rm F}(t) = \pi / [M_{\rm F}(t) 3 / 4 \pi \bar{\rho}(t)]^{1/3}
    !% \end{equation},
    !% and $r_{\rm LSS}(t)$ is the function defined by \cite{naoz_formation_2007}.
    use Cosmology_Parameters
    use Cosmology_Functions
    use Linear_Growth
    use Numerical_Constants_Math
    implicit none
    double precision                          , intent(in   )                         :: time
    class           (cosmologyParametersClass), intent(inout)                         :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), intent(inout)                         :: cosmologyFunctions_
    class           (linearGrowthClass       ), intent(inout)                         :: linearGrowth_
    double precision                          , intent(  out), dimension(3)           :: massFilteringODEs
    double precision                          , intent(  out), dimension(3), optional :: massFilteringScales
    double precision                                         , dimension(4)           :: coefficients
    double precision                                                                  :: expansionFactor     , expansionRate      , &
         &                                                                               massFiltering       , wavenumberFiltering
    
    ! Get fitting function coefficients.
    coefficients          =Mass_Filtering_Early_Epoch_Coefficients(time,cosmologyParameters_,cosmologyFunctions_)
    ! Compute expansion factor and rate.
    expansionFactor       =cosmologyFunctions_%expansionFactor    (time                                         )
    expansionRate         =cosmologyFunctions_%expansionRate      (expansionFactor                              )
    ! Get the filtering mass at the initial time.
    massFiltering         =Mass_Filtering_Early_Epoch             (time,cosmologyParameters_,cosmologyFunctions_)
    ! Find the corresponding filtering wavenumber.
    wavenumberFiltering   =+Pi                                       &
         &                 /(                                        &
         &                   +massFiltering                          &
         &                   *3.0d0                                  &
         &                   /4.0d0                                  &
         &                   /Pi                                     &
         &                   /cosmologyParameters_%OmegaMatter    () &
         &                   /cosmologyParameters_%densityCritical() &
         &                  )**(1.0d0/3.0d0)
    ! Evaluate the three ODE variables at the initial time.
    massFilteringODEs  (1)=+linearGrowth_%value                               (time) &
         &                 /wavenumberFiltering**2
    massFilteringODEs  (2)=+linearGrowth_%value                               (time) &
         &                 /time                                                     &
         &                 *linearGrowth_%logarithmicDerivativeExpansionFactor(time) &
         &                 /wavenumberFiltering**2                                   &
         &                 +2.0d0                                                    &
         &                 /3.0d0                                                    &
         &                 *linearGrowth_%value                               (time) &
         &                 /wavenumberFiltering**2                                   &
         &                 *(                                                        &
         &                   -3.0d0                                                  &
         &                   *coefficients(1)                                        &
         &                   *log(expansionFactor)**2                                &
         &                   +2.0d0                                                  &
         &                   *coefficients(2)                                        &
         &                   *log(expansionFactor)                                   &
         &                   -coefficients(3)                                        &
         &                  )                                                        &
         &                 *expansionRate
    massFilteringODEs  (3)=+massFiltering
    ! Evaluate suitable absolute tolerance scales for the ODE variables.
    if (present(massFilteringScales)) then
       massFilteringScales(1)=+linearGrowth_%value(time)                              &
            &                 *(                                                      &
            &                   +massFiltering                                        &
            &                   *3.0d0                                                &
            &                   /(                                                    &
            &                     +4.0d0                                              &
            &                     *Pi                                                 &
            &                     *cosmologyParameters_%OmegaMatter    ()             &
            &                     *cosmologyParameters_%densityCritical()             &
            &                     )                                                   &
            &                  )**(2.0d0/3.0d0)                                       &
            &                 *Pi**2
       massFilteringScales(2)=+massFilteringScales(1)                               &
            &                 *cosmologyParameters_%HubbleConstant(hubbleUnitsTime)
       massFilteringScales(3)=+massFiltering
    end if
    return
  end subroutine Mass_Filtering_ODE_Initial_Conditions

  function Mass_Filtering_ODE_System(cosmologyParameters_,cosmologyFunctions_,linearGrowth_,time,massParticleMean,temperature,massFilteringODEs) result (massFilteringODEsRateOfChange)
    !% Compute the rates of change of the filtering mass ODE system.
    use Cosmology_Parameters
    use Cosmology_Functions
    use Linear_Growth
    use Numerical_Constants_Math
    use Numerical_Constants_Astronomical
    implicit none
    double precision                                         , dimension(3) :: massFilteringODEsRateOfChange
    class           (cosmologyParametersClass), intent(inout)               :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), intent(inout)               :: cosmologyFunctions_
    class           (linearGrowthClass       ), intent(inout)               :: linearGrowth_
    double precision                          , intent(in   )               :: time                           , massParticleMean   , &
         &                                                                     temperature
    double precision                          , intent(in   ), dimension(3) :: massFilteringODEs
    double precision                                                        :: darkMatterFraction             , wavenumberFiltering, &
         &                                                                     wavenumberFilteringRateOfChange

    ! Compute dark matter mass fraction.
    darkMatterFraction              =1.0d0-cosmologyParameters_%OmegaBaryon()/cosmologyParameters_%OmegaMatter()
    ! Evaluate filtering mass composite terms.
    massFilteringODEsRateOfChange(1)=massFilteringODEs(2)
    if (massParticleMean > 0.0d0) then
       massFilteringODEsRateOfChange(2)=-2.0d0                                                                            &
            &                           *cosmologyFunctions_%expansionRate(cosmologyFunctions_%expansionFactor (time))    &
            &                           *massFilteringODEs(2)                                                             &
            &                           +darkMatterFraction                                                               &
            &                           /                                  cosmologyFunctions_%expansionFactor (time)**2  &
            &                           *boltzmannsConstant                                                               &
            &                           *temperature                                                                      & 
            &                           /massParticleMean                                                                 &
            &                           *linearGrowth_%value                                                   (time)     &
            &                           *(                                                                                &
            &                             +1.0d0                                                                          &
            &                             +rLSS(                                                                          &
            &                                                              cosmologyParameters_%OmegaMatter    (    )   , &
            &                                                              cosmologyFunctions_ %expansionFactor(time)     &
            &                                  )                                                                          &
            &                            )                                                                                &
            &                           *gigayear  **2                                                                    &
            &                           /megaparsec**2
    else
       massFilteringODEsRateOfChange(2)=0.0d0
    end if
    if (massFilteringODEs(1) > 0.0d0) then
       wavenumberFiltering             =+sqrt(                           &
            &                                 +linearGrowth_%value(time) &
            &                                 /massFilteringODEs  (1   ) &
            &                                )
       wavenumberFilteringRateOfChange =+0.5d0                                                                                                 &
            &                           /sqrt(                                                                                                 &
            &                                 +linearGrowth_%value                                                                     (time)  &
            &                                 /massFilteringODEs                                                                       (1   )  &
            &                                )                                                                                                 &
            &                           *(                                                                                                     &
            &                             +linearGrowth_      %logarithmicDerivativeExpansionFactor                                    (time)  &
            &                             *cosmologyFunctions_%expansionRate                       (cosmologyFunctions_%expansionFactor(time)) &
            &                             *linearGrowth_%value                                                                         (time)  &
            &                             *massFilteringODEs                                                                           (1   )  &
            &                             -linearGrowth_%value                                                                         (time)  &
            &                             *massFilteringODEs                                                                           (2   )  &
            &                            )                                                                                                     &
            &                           /massFilteringODEs(1)**2
       massFilteringODEsRateOfChange(3)=-4.0d0                                     &
            &                           *Pi                                    **4 &
            &                           *cosmologyParameters_%densityCritical()    &
            &                           *cosmologyParameters_%OmegaMatter    ()    & 
            &                           /wavenumberFiltering                   **4 &
            &                           *wavenumberFilteringRateOfChange
    else
       massFilteringODEsRateOfChange(3)=+0.0d0
    end if
    return
  end function Mass_Filtering_ODE_System

  function Mass_Filtering_Early_Epoch_Coefficients(time,cosmologyParameters_,cosmologyFunctions_) result (coefficients)
    !% Return the coefficients of the fitting function for the filtering mass at early epochs from
    !% \cite{naoz_formation_2007}. Checks for valid range of redshift and cosmology for the fit to be valid.
    use Cosmology_Parameters
    use Cosmology_Functions
    use Galacticus_Error
    implicit none
    double precision                          , dimension(4)  :: coefficients
    double precision                          , intent(in   ) :: time
    class           (cosmologyParametersClass), intent(inout) :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), intent(inout) :: cosmologyFunctions_
    double precision                                          :: omegaMatter         , expansionFactor, &
         &                                                       redshift

    ! Extract matter density and redshift.
    omegaMatter    =cosmologyParameters_%OmegaMatter                (               )
    expansionFactor=cosmologyFunctions_ %expansionFactor            (time           )
    redshift       =cosmologyFunctions_ %redshiftFromExpansionFactor(expansionFactor)
    ! Validate input.
    if (omegaMatter < 0.25d0 .or. omegaMatter >   0.40d0) call Galacticus_Warn('Mass_Filtering_Early_Epoch: matter density outside validated range of fitting function; 0.25 ≤ Ωₘ ≤ 0.40')
    if (redshift    < 7.00d0 .or. redshift    > 150.00d0) call Galacticus_Warn('Mass_Filtering_Early_Epoch: redshift outside validated range of fitting function; 7 ≤ z ≤ 150'           )
    ! Evaluate fitting function.
    coefficients(1)=-0.38d0*(omegaMatter**2)+ 0.41d0*omegaMatter- 0.16d0
    coefficients(2)=+3.30d0*(omegaMatter**2)- 3.38d0*omegaMatter+ 1.15d0
    coefficients(3)=-9.64d0*(omegaMatter**2)+ 9.75d0*omegaMatter- 2.37d0
    coefficients(4)=+9.80d0*(omegaMatter**2)-10.68d0*omegaMatter+11.60d0
    return
  end function Mass_Filtering_Early_Epoch_Coefficients

  double precision function rLSS(omegaMatter,expansionFactor)
    !% Evaluate the $r_\mathrm{LSS}$ parameter of \cite{naoz_formation_2007} using their fitting formula.
    implicit none
    double precision, intent(in   ) :: omegaMatter     , expansionFactor
    double precision                :: rLSSCoefficient1, rLSSCoefficient2, rLSSCoefficient3

    rLSSCoefficient1=1.0d-4*(-1.99d0*(omegaMatter**2)+2.41d0*omegaMatter+0.21d0)
    rLSSCoefficient2=1.0d-3*(+6.37d0*(omegaMatter**2)-6.99d0*omegaMatter-1.76d0)
    rLSSCoefficient3=1.0d-2*(-1.83d0*(omegaMatter**2)+2.40d0*omegaMatter-0.54d0)
    ! Note that the coefficients for the different exponents of expansion factor are reversed from that given in Naoz &
    ! Barkana. Without this change the fit for rLSS does not work.
    rLSS            =+rLSSCoefficient1/expansionFactor**1.5d0 &
         &           +rLSSCoefficient2/expansionFactor        &
         &           +rLSSCoefficient3 
    return
  end function rLSS

end module Intergalactic_Medium_Filtering_Masses
