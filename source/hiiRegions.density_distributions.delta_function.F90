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

!+    Contributions to this file made by: Sachi Weerasooriya

  !!{
  Implements a class for the distribution of hydrogen density in a HII region in which the distribution is a delta function.
  !!}
  
  !![
  <hiiRegionDensityDistribution name="hiiRegionDensityDistributionDeltaFunction">
   <description>
    A class for the distribution of hydrogen density in a HII region in which the distribution is a delta function. Specifically:
    \begin{equation}
    p(n_\mathrm{H}) = \delta(n_\mathrm{H} - n_\mathrm{H,0}),
    \end{equation}
    where $n_\mathrm{H,0}=${\normalfont \ttfamily [densityHydrogen]}.
   </description>
  </hiiRegionDensityDistribution>
  !!]
  type, extends(hiiRegionDensityDistributionClass) :: hiiRegionDensityDistributionDeltaFunction
     !!{
     A class for the distribution of hydrogen density in a HII region in which the distribution is a delta function.
     !!}
     private
     double precision :: densityHydrogen
   contains
     procedure :: cumulativeDensityDistribution => deltaFunctionCumulativeDensityDistribution
  end type hiiRegionDensityDistributionDeltaFunction

  interface hiiRegionDensityDistributionDeltaFunction
     !!{
     Constructors for the \refClass{hiiRegionDensityDistributionDeltaFunction} HII region density distribution class.
     !!}
     module procedure deltaFunctionConstructorParameters
     module procedure deltaFunctionConstructorInternal
  end interface hiiRegionDensityDistributionDeltaFunction

contains

  function deltaFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{hiiRegionDensityDistributionDeltaFunction} class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hiiRegionDensityDistributionDeltaFunction)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    double precision                                                           :: densityHydrogen

    !![
    <inputParameter>
      <name>densityHydrogen</name>
      <defaultValue>100.0d0</defaultValue>
      <description>The density of hydrogen, $n_\mathrm{H}$, in HII regions (in units of cm$^{-3}$).</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=hiiRegionDensityDistributionDeltaFunction(densityHydrogen)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function deltaFunctionConstructorParameters

  function deltaFunctionConstructorInternal(densityHydrogen) result(self)
    !!{
    Internal constructor for the \refClass{hiiRegionDensityDistributionDeltaFunction} class.
    !!}
    
    implicit none
    type            (hiiRegionDensityDistributionDeltaFunction)                :: self
    double precision                                           , intent(in   ) :: densityHydrogen
    !![
    <constructorAssign variables="densityHydrogen"/>
    !!]

    return
  end function deltaFunctionConstructorInternal

  double precision function deltaFunctionCumulativeDensityDistribution(self,densityHydrogenMinimum,densityHydrogenMaximum) result(distributionFunction)
    !!{
    Compute the cumulative distribution function of the hydrogen density in HII regions. A delta-function distribution is assumed.
    !!}
    implicit none
    class           (hiiRegionDensityDistributionDeltaFunction), intent(inout) :: self
    double precision                                           , intent(in   ) :: densityHydrogenMinimum, densityHydrogenMaximum
    
    if     (                                                &
         &   self%densityHydrogen >= densityHydrogenMinimum &
         &  .and.                                           &
         &   self%densityHydrogen <= densityHydrogenMaximum &
         & ) then
       distributionFunction=+1.0d0
    else
       distributionFunction=+0.0d0
    end if    
    return
  end function deltaFunctionCumulativeDensityDistribution

