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
  Implements calculations of satellite merging times using the \cite{lacey_merger_1993} method with a parameterization of
  orbital parameters designed to fit the results of \cite{tormen_rise_1997} as described by \cite{cole_hierarchical_2000}.
  !!}

  !![
  <satelliteMergingTimescales name="satelliteMergingTimescalesLaceyCole1993Tormen">
   <description>
    A satellite merging timescale class which computes merging timescales using the dynamical friction calculation of
    \cite{lacey_merger_1993} with a parameterization of orbital parameters designed to fit the results of
    \cite{tormen_rise_1997} as described by \cite{cole_hierarchical_2000}. Timescales are multiplied by the value of the
    {\normalfont \ttfamily mergingTimescaleMultiplier} input parameter. Specifically, the merging time is taken to be:
    \begin{equation}
     \tau_\mathrm{merge} = {f_\tau \Phi \tau_\mathrm{dynamical} \over 2 B(1)} { M_\mathrm{host}/M_\mathrm{satellite} \over \ln
     (M_\mathrm{host}/M_\mathrm{satellite})}
    \end{equation}
    where $f_\tau=${\normalfont \ttfamily mergingTimescaleMultiplier}, $\tau_\mathrm{dynamical}$ is the dynamical time of the
    host halo and $B(x)=\hbox{erf}(x)-2 x \exp(x)/\sqrt{\Pi}$. The orbital factor $\Phi \equiv \epsilon^{0.78}
    (R_\mathrm{c}/R_\mathrm{virial})^2$ is drawn at random from a log-normal distribution with median $-0.14$ and dispersion
    $0.26$ as found by \cite{cole_hierarchical_2000}.
   </description>
  </satelliteMergingTimescales>
  !!]
  type, extends(satelliteMergingTimescalesLaceyCole1993) :: satelliteMergingTimescalesLaceyCole1993Tormen
     !!{
     A class implementing the \cite{cole_hierarchical_2000} method for satellite merging timescales.
     !!}
     private
   contains
     procedure :: timeUntilMerging => laceyCole1993TormenTimeUntilMerging
  end type satelliteMergingTimescalesLaceyCole1993Tormen

  interface satelliteMergingTimescalesLaceyCole1993Tormen
     !!{
     Constructors for the \cite{cole_hierarchical_2000} merging timescale class.
     !!}
     module procedure laceyCole1993TormenConstructorParameters
     module procedure laceyCole1993TormenConstructorInternal
  end interface satelliteMergingTimescalesLaceyCole1993Tormen

contains

  function laceyCole1993TormenConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{cole_hierarchical_2000} merging timescale class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (satelliteMergingTimescalesLaceyCole1993Tormen)                :: self
    type (inputParameters                              ), intent(inout) :: parameters

    self%satelliteMergingTimescalesLaceyCole1993=satelliteMergingTimescalesLaceyCole1993(parameters)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function laceyCole1993TormenConstructorParameters

  function laceyCole1993TormenConstructorInternal(timescaleMultiplier,darkMatterHaloScale_) result(self)
    !!{
    Constructor for the \cite{cole_hierarchical_2000} merging timescale class.
    !!}
    implicit none
    type            (satelliteMergingTimescalesLaceyCole1993Tormen)                        :: self
    double precision                                               , intent(in   )         :: timescaleMultiplier
    class           (darkMatterHaloScaleClass                     ), intent(in   ), target :: darkMatterHaloScale_

    self%satelliteMergingTimescalesLaceyCole1993=satelliteMergingTimescalesLaceyCole1993(timescaleMultiplier,darkMatterHaloScale_)
    return
  end function laceyCole1993TormenConstructorInternal

  double precision function laceyCole1993TormenTimeUntilMerging(self,node,orbit)
    !!{
    Return the timescale for merging satellites using the \cite{lacey_merger_1993} method with a parameterization of orbital
    parameters designed to fit the results of \cite{tormen_rise_1997} as described by \cite{cole_hierarchical_2000}.
    !!}
    implicit none
    class           (satelliteMergingTimescalesLaceyCole1993Tormen), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    type            (keplerOrbit                                  ), intent(inout) :: orbit
    double precision                                               , parameter     :: orbitalFactorDistributionSigma=+0.26d0                   !   Cole et al. (2000).
    double precision                                               , parameter     :: orbitalFactorDistributionMean =-0.14d0                   !   Cole et al. (2000).
    double precision                                                               :: log10OrbitalFactor                    , randomDeviate, &
         &                                                                            orbitalFactor
    !$GLC attributes unused :: orbit

    ! Compute the orbital factor - selected at random from a lognormal distribution.
    randomDeviate     =+orbitalFactorDistributionSigma                              &
         &             *node%hostTree%randomNumberGenerator_%standardNormalSample()
    log10OrbitalFactor=+orbitalFactorDistributionMean                               &
         &             +randomDeviate
    orbitalFactor     =+10.0d0**log10OrbitalFactor
    ! Compute the timescale.
    laceyCole1993TormenTimeUntilMerging=orbitalFactor*self%timeUntilMergingMassDependence(node)
    return
  end function laceyCole1993TormenTimeUntilMerging
