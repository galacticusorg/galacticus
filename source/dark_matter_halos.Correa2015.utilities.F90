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
Contains a module which implements utility functions for the \cite{correa_accretion_2015} dark matter halo models.
!!}

module Dark_Matter_Halos_Correa2015
  !!{
  Implements utility functions for the \cite{correa_accretion_2015} dark matter halo models.
  !!}
  implicit none
  private
  public :: Dark_Matter_Halo_Correa2015_Fit_Parameters

contains

  subroutine Dark_Matter_Halo_Correa2015_Fit_Parameters(mass,expansionFactor,cosmologyFunctions_,linearGrowth_,cosmologicalMassVariance_,aTilde,bTilde)
    !!{
    Computes fitting function parameters for the \cite{correa_accretion_2015} dark matter halo models.
    !!}
    use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
    use :: Cosmology_Functions       , only : cosmologyFunctionsClass
    use :: Linear_Growth             , only : linearGrowthClass
    use :: Numerical_Constants_Math  , only : Pi
    implicit none
    double precision                               , intent(in   ) :: mass                             , expansionFactor
    double precision                               , intent(  out) :: aTilde                           , bTilde
    class           (cosmologyFunctionsClass      ), intent(inout) :: cosmologyFunctions_
    class           (linearGrowthClass            ), intent(inout) :: linearGrowth_
    class           (cosmologicalMassVarianceClass), intent(inout) :: cosmologicalMassVariance_
    double precision                               , parameter     :: deltaCritical            =1.686d0
    double precision                                               :: redshiftFormation                , q              , &
         &                                                            sigma                            , sigmaQ         , &
         &                                                            f

    redshiftFormation=+1.8837d0                & ! Correa et al. eqn. 6
         &            +0.0237d0*log10(mass)    &
         &            -0.0064d0*log10(mass)**2
    q     =4.137d0/redshiftFormation**0.9476d0   ! Correa et al. eqn. 5.
    sigma =cosmologicalMassVariance_%rootVariance(mass  ,cosmologyFunctions_%cosmicTime(1.0d0))
    sigmaQ=cosmologicalMassVariance_%rootVariance(mass/q,cosmologyFunctions_%cosmicTime(1.0d0))
    f     =1.0d0/sqrt(sigmaQ**2-sigma**2)
    aTilde=+f                                                                                     & ! Correa et al. eqn. 2.
         & *(                                                                                     &
         &   +1.0d0                                                                               &
         &   -sqrt(2.0d0/Pi)                                                                      &
         &   *deltaCritical                                                                       &
         &   *                                                                   expansionFactor  &
         &   *linearGrowth_%logarithmicDerivativeExpansionFactor(expansionFactor=expansionFactor) &
         &   /linearGrowth_%value                               (expansionFactor=expansionFactor) &
         &  )
    bTilde=-f ! Correa et al. eqn. 3.
    return
  end subroutine Dark_Matter_Halo_Correa2015_Fit_Parameters

end module Dark_Matter_Halos_Correa2015
