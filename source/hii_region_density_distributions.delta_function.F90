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

  !!{
  Implementation of a class for the distribution of hydrogen density in a HII region.
  !!}
  !![
  <hiiRegionDensityDistribution name="hiiRegionDensityDistributionDeltaFunction">
   <description>
    A function class that calculates the distribution of hydrogen density.
    \begin{equation}
     if(nH_min \less n_H \less nH_max)
        weight_nH=1
     else
        weight_nH=0

    \end{equation}
 
   </description>
  </hiiRegionDensityDistribution>
  !!]
  type, extends(hiiRegionDensityDistributionClass) :: hiiRegionDensityDistributionDeltaFunction
     !!{
     Implementation of a density distribution to a delta function
     !!}
     private
     double precision, intent(in   )     :: densityHydrogen, densityHydrogenMin, densityHydrogenMax
     
   contains
     procedure :: hiiRegionDensityDistribution => hiiRegionDensityDistributionDeltaFunction
  end type hiiRegionDensityDistributionDeltaFunction

  interface hiiRegionDensityDistributionDeltaFunction
     !!{
     Constructors for the hiiRegionDensityDistributionDeltaFunction for star HII region density distribution class.
     !!}
     module procedure hiiRegionDensityDistributionDeltaFunctionConstructorParameters
     module procedure hiiRegionDensityDistributionDeltaFunctionConstructorInternal
  end interface hiiRegionDensityDistributionDeltaFunction

contains

  function hiiRegionDensityDistributionDeltaFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the hiiRegionDensityDistributionDeltaFunction class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hiiRegionDensityDistributionDeltaFunction)                :: self
    type            (inputParameters                    ), intent(inout).      :: parameters
    self=hiiRegionDensityDistributionDeltaFunction(densityHydrogen,densityHydrogenMin, densityHydrogenMax)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function hiiRegionDensityDistributionDeltaFunctionConstructorParameters

  function hiiRegionDensityDistributionDeltaFunctionConstructorInternal(densityHydrogen) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily hiiRegionDensityDistributionDeltaFunction} class.
    !!}
    
    implicit none
    type            (hiiRegionDensityDistributionDeltaFunction)          :: self
    double precision                                     , intent(in   ) :: densityHydrogen,densityHydrogenMin, densityHydrogenMax
   

    !![
    <constructorAssign variables="densityHydrogen,densityHydrogenMin, densityHydrogenMax"/>
    !!]
    
    
    
    return
  end function hiiRegionDensityDistributionDeltaFunctionConstructorInternal

  double precision function hiiRegionDensityDistributionDeltaFunction(self,densityHydrogen,densityHydrogenMin, densityHydrogenMax) result(weightDensityHydrogen)
    !!{
   
    \begin{equation}
    if(n_H,min \less n_H \less n_H,max)
       n_H=1
    Else
       n_H=0
    \end{equation}
    !!}

  
    implicit none
    class  (hiiRegionDensityDistributionDeltaFunction), intent(inout) :: self
   
    double precision                parameter     :: densityHydrogen,densityHydrogenMin, densityHydrogenMax        
    
    ! Compute the distribution of HII regions for a delta function.
    if (densityHydrogenMin .lt. self%densityHydrogen .lt. densityHydrogenMax) then
	weightDensityHydrogen=1
    else
        weightDensityHydrogen=0
    end if
    
    return
  end function hiiRegionDensityDistributionDeltaFunction

