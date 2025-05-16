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

!+    Contributions to this file made by: Sachi Weerasooriya

  !!{
  Implementation of a class for the distribution of hydrogen density in a HII region which assumes a log-normal distribution.
  !!}

  use :: Statistics_Distributions, only : distributionFunction1DLogNormal

  !![
  <hiiRegionDensityDistribution name="hiiRegionDensityDistributionLogNormal">
   <description>
    A class for the  distribution of hydrogen density in a HII region in which the distribution is a lognormal, specifically:
    \begin{equation}
     p(n_\mathrm{H}) = \left\{ \begin{array}{ll} \frac{1}{\sqrt{2\pi} n_\mathrm{H} \sigma} \exp\left(-\frac{1}{2}\left[\frac{\log(n_\mathrm{H})-\log(n_\mathrm{H,0})}{\sigma}\right]^2\right) &amp; \hbox{ if } n_\mathrm{H,min} &lt; n_\mathrm{H} &lt; n_\mathrm{H,max}, \\ 0 &amp; \hbox{ otherwise.} \end{array} \right.
    \end{equation} 
   </description>
  </hiiRegionDensityDistribution>
  !!]
  type, extends(hiiRegionDensityDistributionClass) :: hiiRegionDensityDistributionLogNormal
     !!{
     A class for the  distribution of hydrogen density in a HII region in which the distribution is a lognormal.
     !!}
     private
     type            (distributionFunction1DLogNormal) :: distribution
     double precision                                  :: densityHydrogenMinimum  , densityHydrogenMaximum, &
          &                                               densityHydrogenReference, sigma
     
   contains
     procedure :: cumulativeDensityDistribution => lognormalCumulativeDensityDistribution
  end type hiiRegionDensityDistributionLogNormal

  interface hiiRegionDensityDistributionLogNormal
     !!{
     Constructors for the \refClass{hiiRegionDensityDistributionLogNormal} HII region density distribution class.
     !!}
     module procedure logNormalConstructorParameters
     module procedure logNormalConstructorInternal
  end interface hiiRegionDensityDistributionLogNormal

contains

  function logNormalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the hiiRegionDensityDistributionLogNormal class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hiiRegionDensityDistributionLogNormal)                :: self
    type            (inputParameters                    ), intent(inout)   :: parameters
    double precision                                                       :: densityHydrogenMinimum  , densityHydrogenMaximum, &
         &                                                                    densityHydrogenReference, sigma    
    !![
    <inputParameter>
      <name>densityHydrogenReference</name>
      <defaultValue>250.0d0</defaultValue>
      <defaultSource>Fit to FIRE II data from \cite{yang_efficient_2023}, provided by S. Yang (private communication).</defaultSource>
      <description>The parameter $n_\mathrm{H,0}$ in the log normal distribution of HII region densities.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>sigma</name>
      <defaultValue>1.8d0</defaultValue>
      <defaultSource>Fit to FIRE II data from \cite{yang_efficient_2023}, provided by S. Yang (private communication).</defaultSource>
      <description>The parameter $\sigma$ in the log normal distribution of HII region densities.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>densityHydrogenMinimum</name>
      <defaultValue>1.0d1</defaultValue>
      <defaultSource>Fit to FIRE II data from \cite{yang_efficient_2023}, provided by S. Yang (private communication).</defaultSource>
      <description>Minimum value of hydrogen density in HII regions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>densityHydrogenMaximum</name>
      <defaultValue>1.0d5</defaultValue>
      <defaultSource>Fit to FIRE II data from \cite{yang_efficient_2023}, provided by S. Yang (private communication).</defaultSource>
      <description>Maximum value of hydrogen density in HII regions.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=hiiRegionDensityDistributionLogNormal(densityHydrogenReference,sigma,densityHydrogenMinimum,densityHydrogenMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function logNormalConstructorParameters

  function logNormalConstructorInternal(densityHydrogenReference,sigma,densityHydrogenMinimum,densityHydrogenMaximum) result(self)
    !!{
    Internal constructor for the \refClass{hiiRegionDensityDistributionLogNormal} class.
    !!}
    implicit none
    type            (hiiRegionDensityDistributionLogNormal)                :: self
    double precision                                       , intent(in   ) :: densityHydrogenMinimum  , densityHydrogenMaximum, &
         &                                                                    densityHydrogenReference, sigma
    !![
    <constructorAssign variables="densityHydrogenReference, sigma, densityHydrogenMinimum, densityHydrogenMaximum"/>
    !!]
    
    self%distribution=distributionFunction1DLogNormal(x0=densityHydrogenReference,sigma=sigma,limitLower=densityHydrogenMinimum,limitUpper=densityHydrogenMaximum)
    return
  end function logNormalConstructorInternal

  double precision function lognormalCumulativeDensityDistribution(self,densityHydrogenMinimum,densityHydrogenMaximum) result(weightDensityHydrogen)
    !!{
    Compute the cumulative distribution function of the hydrogen density in HII regions for a log-normal distribution.
    !!}
    implicit none
    class           (hiiRegionDensityDistributionLogNormal), intent(inout) :: self    
    double precision                                       , intent(in   ) :: densityHydrogenMinimum, densityHydrogenMaximum
    
    ! Compute the distribution of HII regions for a log normal function.
    weightDensityHydrogen=+self%distribution%cumulative(densityHydrogenMaximum) &
         &                -self%distribution%cumulative(densityHydrogenMinimum)
    
    return
  end function lognormalCumulativeDensityDistribution

