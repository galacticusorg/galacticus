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
  Implementation of a posterior sampling differential evolution proposal size temperature exponent class in which the exponent
  is fixed.
  !!}

  !![
  <posteriorSampleDffrntlEvltnPrpslSzTmpExp name="posteriorSampleDffrntlEvltnPrpslSzTmpExpFixed">
   <description>
    This class uses a fixed $\alpha=${\normalfont \ttfamily [alpha]}.
   </description>
  </posteriorSampleDffrntlEvltnPrpslSzTmpExp>
  !!]
  type, extends(posteriorSampleDffrntlEvltnPrpslSzTmpExpClass) :: posteriorSampleDffrntlEvltnPrpslSzTmpExpFixed
     !!{
     Implementation of a posterior sampling differential evolution proposal size class in which the exponent is fixed.
     !!}
     private
     double precision :: exponentValue
   contains
     procedure :: exponent => fixedExponent
  end type posteriorSampleDffrntlEvltnPrpslSzTmpExpFixed

  interface posteriorSampleDffrntlEvltnPrpslSzTmpExpFixed
     !!{
     Constructors for the \refClass{posteriorSampleDffrntlEvltnPrpslSzTmpExpFixed} posterior sampling differential evolution random jump class.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface posteriorSampleDffrntlEvltnPrpslSzTmpExpFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleDffrntlEvltnPrpslSzTmpExpFixed} posterior sampling differential evolution random jump class which builds
    the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleDffrntlEvltnPrpslSzTmpExpFixed)                 :: self
    type            (inputParameters                              ), intent(inout)  :: parameters
    double precision                                                                :: exponentValue

    !![
    <inputParameter>
      <name>exponentValue</name>
      <description>The exponent of temperature.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=posteriorSampleDffrntlEvltnPrpslSzTmpExpFixed(exponentValue)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(exponentValue) result(self)
    !!{
    Internal constructor for the \refClass{posteriorSampleDffrntlEvltnPrpslSzTmpExpFixed} posterior sampling differential evolution proposal size
    temperature exponent class.
    !!}
    implicit none
    type            (posteriorSampleDffrntlEvltnPrpslSzTmpExpFixed)                :: self
    double precision                                               , intent(in   ) :: exponentValue
    !![
    <constructorAssign variables="exponentValue"/>
    !!]

    return
  end function fixedConstructorInternal

  double precision function fixedExponent(self,temperedStates,temperatures,simulationState,simulationConvergence)
    !!{
    Return the fixed differential evolution proposal size temperature exponent.
    !!}
    implicit none
    class           (posteriorSampleDffrntlEvltnPrpslSzTmpExpFixed), intent(inout)               :: self
    class           (posteriorSampleStateClass                    ), intent(inout)               :: simulationState
    class           (posteriorSampleStateClass                    ), intent(inout), dimension(:) :: temperedStates
    double precision                                               , intent(in   ), dimension(:) :: temperatures
    class           (posteriorSampleConvergenceClass              ), intent(inout)               :: simulationConvergence
    !$GLC attributes unused :: temperedStates, temperatures, simulationState, simulationConvergence

    fixedExponent=self%exponentValue
    return
  end function fixedExponent

