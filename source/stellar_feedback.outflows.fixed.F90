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
  Implementation of a fixed fraction stellar feedback model.
  !!}

  !![
  <stellarFeedbackOutflows name="stellarFeedbackOutflowsFixed">
   <description>
  !!]

  !![
    A stellar feedback outflow class in which the outflow rate is fixed. Specifically,
    \begin{equation}
     \dot{M}_\mathrm{outflow} = f_\mathrm{outflow} {\dot{E} \over E_\mathrm{canonical}},
    \end{equation}
    where $f_\mathrm{outflow}=${\normalfont \ttfamily [fraction]} is the fraction of the star formation rate that goes into
    outflow, $\dot{E}$ is the rate of energy input from stellar populations and $E_\mathrm{canonical}$ is the total energy
    input by a canonical stellar population normalized to $1 M_\odot$ after infinite time.
   </description>
  </stellarFeedbackOutflows>
  !!]
  type, extends(stellarFeedbackOutflowsClass) :: stellarFeedbackOutflowsFixed
     !!{
     Implementation of a fixed fraction stellar feedback model.
     !!}
     private
     double precision :: fraction
   contains
     procedure :: outflowRate => fixedOutflowRate
  end type stellarFeedbackOutflowsFixed

  interface stellarFeedbackOutflowsFixed
     !!{
     Constructors for the fixed fraction stellar feedback class.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface stellarFeedbackOutflowsFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the fixed fraction stellar feedback class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (stellarFeedbackOutflowsFixed)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    double precision                                              :: fraction

    !![
    <inputParameter>
      <name>fraction</name>
      <source>parameters</source>
      <defaultValue>0.01d0</defaultValue>
      <description>The ratio of outflow rate to star formation rate in disks.</description>
    </inputParameter>
    !!]
    self=stellarFeedbackOutflowsFixed(fraction)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(fraction) result(self)
    !!{
    Internal constructor for the fixed stellar feedback class.
    !!}
    implicit none
    type            (stellarFeedbackOutflowsFixed)                :: self
    double precision                              , intent(in   ) :: fraction

    !![
    <constructorAssign variables="fraction"/>
    !!]
    return
  end function fixedConstructorInternal

  subroutine fixedOutflowRate(self,component,rateStarFormation,rateEnergyInput,rateOutflowEjective,rateOutflowExpulsive)
    !!{
    Returns the outflow rate (in $M_\odot$ Gyr$^{-1}$) due to stellar feedback in the given {\normalfont \ttfamily
   component}. Assumes a fixed ratio of outflow rate to star formation rate.
   !!}
    use :: Stellar_Feedback, only : feedbackEnergyInputAtInfinityCanonical
    implicit none
    class           (stellarFeedbackOutflowsFixed), intent(inout) :: self
    class           (nodeComponent               ), intent(inout) :: component
    double precision                              , intent(in   ) :: rateEnergyInput    , rateStarFormation
    double precision                              , intent(  out) :: rateOutflowEjective, rateOutflowExpulsive
    !$GLC attributes unused :: component, rateStarFormation

    rateOutflowEjective =+self%fraction                          &
         &               *rateEnergyInput                        &
         &               /feedbackEnergyInputAtInfinityCanonical
    rateOutflowExpulsive=+0.0d0
    return
  end subroutine fixedOutflowRate
