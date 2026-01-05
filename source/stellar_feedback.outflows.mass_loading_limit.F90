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
  Implementation of a stellar feedback model which limits to a maximum mass-loading factor.
  !!}
  
  !![
  <stellarFeedbackOutflows name="stellarFeedbackOutflowsMassLoadingLimit">
    <description>
    A stellar feedback model which limits to a maximum mass-loading factor. The outflow rate will be
    \begin{equation}
    \dot{M}_\mathrm{out} = \mathrm{tanh} (\beta / \beta_\mathrm{max}) \beta_\mathrm{max} \dot{M}_\star,
    \end{equation}
    where $\dot{M}_\star$ is the star formation rate, $\beta$ is the mass-loading factor of the decorated
    model, and $\beta_\mathrm{max}=${\normalfont \ttfamily [factorMassLoadingMaximum]} is the maximum mass
    loading factor.
    </description>
  </stellarFeedbackOutflows>
  !!]
  type, extends(stellarFeedbackOutflowsClass) :: stellarFeedbackOutflowsMassLoadingLimit
     !!{
     Implementation of a stellar feedback model which limits to a maximum mass-loading factor.
     !!}
     private
     class           (stellarFeedbackOutflowsClass), pointer :: stellarFeedbackOutflows_ => null()
     double precision                                        :: factorMassLoadingMaximum
   contains
     final     ::                massLoadingLimitDestructor
     procedure :: outflowRate => massLoadingLimitOutflowRate
  end type stellarFeedbackOutflowsMassLoadingLimit

  interface stellarFeedbackOutflowsMassLoadingLimit
     !!{
     Constructors for the rate-limiting stellar feedback class.
     !!}
     module procedure massLoadingLimitConstructorParameters
     module procedure massLoadingLimitConstructorInternal
  end interface stellarFeedbackOutflowsMassLoadingLimit

contains

  function massLoadingLimitConstructorParameters(parameters) result(self)
    !!{
    Constructor for the rate-limiting stellar feedback class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (stellarFeedbackOutflowsMassLoadingLimit)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (stellarFeedbackOutflowsClass           ), pointer       :: stellarFeedbackOutflows_
    double precision                                                         :: factorMassLoadingMaximum

    !![
    <inputParameter>
      <name>factorMassLoadingMaximum</name>
      <source>parameters</source>
      <description>The maximum mass loading factor for outflows due to stellar feedback.</description>
    </inputParameter>
    <objectBuilder class="stellarFeedbackOutflows" name="stellarFeedbackOutflows_" source="parameters"/>
    !!]
    self=stellarFeedbackOutflowsMassLoadingLimit(factorMassLoadingMaximum,stellarFeedbackOutflows_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="stellarFeedbackOutflows_"/>
    !!]
    return
  end function massLoadingLimitConstructorParameters

  function massLoadingLimitConstructorInternal(factorMassLoadingMaximum,stellarFeedbackOutflows_) result(self)
    !!{
    Internal constructor for the massLoadingLimit stellar feedback class.
    !!}
    implicit none
    type            (stellarFeedbackOutflowsMassLoadingLimit)                        :: self
    class           (stellarFeedbackOutflowsClass           ), intent(in   ), target :: stellarFeedbackOutflows_
    double precision                                         , intent(in   )         :: factorMassLoadingMaximum
    !![
    <constructorAssign variables="factorMassLoadingMaximum, *stellarFeedbackOutflows_"/>
    !!]
    
    return
  end function massLoadingLimitConstructorInternal

  subroutine massLoadingLimitDestructor(self)
    !!{
    Internal constructor for the massLoadingLimit stellar feedback class.
    !!}
    implicit none
    type(stellarFeedbackOutflowsMassLoadingLimit), intent(inout) :: self
    
    !![
    <objectDestructor name="self%stellarFeedbackOutflows_"/>
    !!]
    return
  end subroutine massLoadingLimitDestructor

  subroutine massLoadingLimitOutflowRate(self,component,rateStarFormation,rateEnergyInput,rateOutflowEjective,rateOutflowExpulsive)
    !!{
    Limits the outflow rate from another stellar feedback class such that the mass loading factor never exceeds a given
    value.
    !!}
    implicit none
    class           (stellarFeedbackOutflowsMassLoadingLimit), intent(inout) :: self
    class           (nodeComponent                          ), intent(inout) :: component
    double precision                                         , intent(in   ) :: rateEnergyInput    , rateStarFormation
    double precision                                         , intent(  out) :: rateOutflowEjective, rateOutflowExpulsive
    double precision                                                         :: factorMassLoading  , factorReduction

    if (rateStarFormation > 0.0d0) then
       call self%stellarFeedbackOutflows_%outflowRate(component,rateStarFormation,rateEnergyInput,rateOutflowEjective,rateOutflowExpulsive)
       factorMassLoading   =+(                      &
            &                 +rateOutflowEjective  &
            &                 +rateOutflowExpulsive &
            &                 )                     &
            &               /  rateStarFormation
       factorReduction     =+tanh(     factorMassLoading       /self%factorMassLoadingMaximum) &
            &               *     self%factorMassLoadingMaximum/      factorMassLoading
       
       rateOutflowEjective =+rateOutflowEjective *factorReduction
       rateOutflowExpulsive=+rateOutflowExpulsive*factorReduction       
    else
       rateOutflowEjective =+0.0d0
       rateOutflowExpulsive=+0.0d0
    end if
    return
  end subroutine massLoadingLimitOutflowRate
