!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
!+    Contributions to this file made by: Sachi Weerasooriya, Andrew Benson

  !!{
  Implementation of a class for the distribution of hydrogen density in a HII region.
  !!}
  use :: Statistics_Distributions, only : distributionFunction1DLogNormal
  !![
  <hiiRegionDensityDistribution name="hiiRegionDensityDistributionLogNormal">
   <description>
    A function class that calculates the distribution of hydrogen density based on a lognormal distribution.
    \begin{equation}
     if(nH_min \less n_H \less nH_max)
        weight_nH=exp(-½[{ln(n_H)-ln(n₀)}/σ]²)
     else
        weight_nH=0

    \end{equation}
 
   </description>
  </hiiRegionDensityDistribution>
  !!]
  type, extends(hiiRegionDensityDistributionClass) :: hiiRegionDensityDistributionLogNormal
     !!{
     Implementation of a log normal density distribution
     !!}
     private
     type            (distributionFunction1DLogNormal)                :: distribution
     double precision     :: densityHydrogenMin, densityHydrogenMax, x0, sigma
     
   contains
     procedure :: cumulativeDensityDistribution => hiiRegionDensityDistributionLogNormal
  end type hiiRegionDensityDistributionLogNormal

  interface hiiRegionDensityDistributionLogNormal
     !!{
     Constructors for the hiiRegionDensityDistributionLogNormal for star HII region density distribution class.
     !!}
     module procedure hiiRegionDensityDistributionLogNormalConstructorParameters
     module procedure hiiRegionDensityDistributionLogNormalConstructorInternal
  end interface hiiRegionDensityDistributionLogNormal

contains

  function hiiRegionDensityDistributionLogNormalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the hiiRegionDensityDistributionLogNormal class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hiiRegionDensityDistributionLogNormal)                :: self
    type            (inputParameters                    ), intent(inout)   :: parameters
    double precision                                                       :: densityHydrogenMin, densityHydrogenMax, x0, sigma    
    !![
    <inputParameter>
      <name>x0</name>
      <defaultValue>250.0d0</defaultValue>
      <defaultSource>Fit to FIRE II data from Yang+</defaultSource>
      <description>x0 of log normal distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>sigma</name>
      <defaultValue>1.8d0</defaultValue>
      <defaultSource>Fit to FIRE II data from Yang+</defaultSource>
      <description>Standard deviation</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>densityHydrogenMin</name>
      <defaultValue>10d0</defaultValue>
      <defaultSource></defaultSource>
      <description>Minimum value of hydrogen density in HII region</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>densityHydrogenMax</name>
      <defaultValue>1d5</defaultValue>
      <defaultSource></defaultSource>
      <description>Maximum value of hydrogen density in HII region</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=hiiRegionDensityDistributionLogNormal(x0,sigma,densityHydrogenMin,densityHydrogenMax)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function hiiRegionDensityDistributionLogNormalConstructorParameters

  function hiiRegionDensityDistributionLogNormalConstructorInternal(x0,sigma,densityHydrogenMin, densityHydrogenMax) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily hiiRegionDensityDistributionLogNormal} class.
    !!}
    !use :: Statistics_Distributions, only : distributionFunction1DLogNormal
    implicit none
    type            (hiiRegionDensityDistributionLogNormal)          :: self
    !type            (distributionFunction1DLogNormal)                :: distribution 
    double precision                                 , intent(in   ) :: densityHydrogenMin, densityHydrogenMax, x0, sigma, 
    !![
    <constructorAssign variables="densityHydrogenMin, densityHydrogenMax, x0, sigma"/>
    !!]
    self%distribution=distributionFunction1DLogNormal(x0=x0,sigma=sigma,limitLower=densityHydrogenMin,limitUpper=densityHydrogenMax)
    return
  end function hiiRegionDensityDistributionLogNormalConstructorInternal

  double precision function hiiRegionDensityDistributionLogNormal(self,densityHydrogenMin, densityHydrogenMax) result(weightDensityHydrogen)
    !!{
    \begin{equation}
    if(n_H,min \less n_H \less n_H,max)
       weight_nH= exp(-½[{ln(n_H)-ln(n₀)}/σ]²) 
    Else
       weight_nH=0
    \end{equation}
    !!}
    implicit none
    class           (hiiRegionDensityDistributionLogNormal), intent(inout) :: self    
    double precision, intent(in)                                                       :: densityHydrogenMin, densityHydrogenMax         
    
    ! Compute the distribution of HII regions for a log normal function.
    weightDensityHydrogen=self%distribution%cumulative(densityHydrogenMax) - self%distribution%cumulative(densityHydrogenMin) 
    
    return
  end function hiiRegionDensityDistributionLogNormal

