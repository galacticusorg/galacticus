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

  !!{
  An implementation of the dark matter halo spin distribution which assumes a
  $\delta$-function distribution.
  !!}

  !![
  <haloSpinDistribution name="haloSpinDistributionDeltaFunction">
   <description>
    A halo spin distribution class in which the spin is drawn from a $\delta$-function distribution, $P(\lambda) =
    \delta(\lambda-\lambda_0)$, where $\lambda_0=${\normalfont \ttfamily [spin]}, i.e. a fixed value of spin equal to
    $\lambda_0$ is returned.
   </description>
  </haloSpinDistribution>
  !!]
  type, extends(haloSpinDistributionClass) :: haloSpinDistributionDeltaFunction
     !!{
     A dark matter halo spin distribution concentration class which assumes a
     $\delta$-function distribution.
     !!}
     private
     double precision :: spin
   contains
     procedure :: sample       => deltaFunctionSample
     procedure :: distribution => deltaFunctionDistribution
  end type haloSpinDistributionDeltaFunction

  interface haloSpinDistributionDeltaFunction
     !!{
     Constructors for the \refClass{haloSpinDistributionDeltaFunction} dark matter halo spin
     distribution class.
     !!}
     module procedure deltaFunctionConstructorParameters
     module procedure deltaFunctionConstructorInternal
  end interface haloSpinDistributionDeltaFunction

contains

  function deltaFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloSpinDistributionDeltaFunction} dark matter halo spin
    distribution class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloSpinDistributionDeltaFunction)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    double precision                                                   :: spin

    !![
    <inputParameter>
      <name>spin</name>
      <source>parameters</source>
      <defaultValue>0.03687d0</defaultValue>
      <defaultSource>\citep{bett_spin_2007}</defaultSource>
      <description>The fixed value of spin in a $\delta$-function spin distribution.</description>
    </inputParameter>
    !!]
    self=haloSpinDistributionDeltaFunction(spin)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function deltaFunctionConstructorParameters

  function deltaFunctionConstructorInternal(spin) result(self)
    !!{
    Internal constructor for the \refClass{haloSpinDistributionDeltaFunction} dark matter halo spin
    distribution class.
    !!}
    implicit none
    type            (haloSpinDistributionDeltaFunction)                :: self
    double precision                                   , intent(in   ) :: spin
    !![
    <constructorAssign variables="spin"/>
    !!]
    
    return
  end function deltaFunctionConstructorInternal

  double precision function deltaFunctionSample(self,node)
    !!{
    Sample from a $\delta$-function spin parameter distribution for the given {\normalfont
    \ttfamily node}.
    !!}
    implicit none
    class(haloSpinDistributionDeltaFunction), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node
    !$GLC attributes unused :: node

    deltaFunctionSample=self%spin
    return
  end function deltaFunctionSample

  double precision function deltaFunctionDistribution(self,node)
    !!{
    Return the spin parameter distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(haloSpinDistributionDeltaFunction), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    deltaFunctionDistribution=0.0d0
    call Error_Report('distribution function can not be evaluated'//{introspection:location})
    return
  end function deltaFunctionDistribution
