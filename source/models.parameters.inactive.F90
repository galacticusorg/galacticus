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
  Implementation of an inactive model parameter class.
  !!}

  !![
  <modelParameter name="modelParameterInactive">
   <description>An inactive model parameter class.</description>
  </modelParameter>
  !!]
  type, extends(modelParameterClass) :: modelParameterInactive
     !!{
     Implementation of an inactive model parameter class.
     !!}
     private
     type(varying_string) :: name_
   contains
     procedure :: logPrior           => inactiveLogPrior
     procedure :: name               => inactiveName
     procedure :: priorSample        => inactivePriorSample
     procedure :: priorInvert        => inactivePriorInvert
     procedure :: priorMinimum       => inactivePriorMinimum
     procedure :: priorMaximum       => inactivePriorMaximum
     procedure :: randomPerturbation => inactiveRandomPerturbation
     procedure :: map                => inactiveMap
     procedure :: unmap              => inactiveUnmap
  end type modelParameterInactive

  interface modelParameterInactive
     !!{
     Constructors for the \refClass{modelParameterInactive} 1D distribution function class.
     !!}
     module procedure inactiveConstructorParameters
     module procedure inactiveConstructorInternal
  end interface modelParameterInactive

contains

  function inactiveConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{modelParameterInactive} model parameter class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(modelParameterInactive)                :: self
    type(inputParameters       ), intent(inout) :: parameters
    type (varying_string       )                :: name

    !![
    <inputParameter>
      <name>name</name>
      <description>The name of the parameter.</description>
      <defaultValue>var_str('')</defaultValue>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=modelParameterInactive(name)
     !![
     <inputParametersValidate source="parameters"/>
     !!]
   return
  end function inactiveConstructorParameters

  function inactiveConstructorInternal(name_) result(self)
    !!{
    Internal constructor for the \refClass{modelParameterInactive} model parameter class.
    !!}
    implicit none
    type (modelParameterInactive)                :: self
    type (varying_string        ), intent(in   ) :: name_
    !![
    <constructorAssign variables="name_"/>
    !!]

    return
  end function inactiveConstructorInternal

  function inactiveName(self)
    !!{
    Return the name of this parameter.
    !!}
    implicit none
    type (varying_string        )                :: inactiveName
    class(modelParameterInactive), intent(inout) :: self

    inactiveName=self%name_
    return
  end function inactiveName

  double precision function inactiveLogPrior(self,x)
    !!{
    Return the log-prior on this parameter.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (modelParameterInactive), intent(inout) :: self
    double precision                        , intent(in   ) :: x
    !$GLC attributes unused :: self, x

    inactiveLogPrior=0.0d0
    call Error_Report('parameter is inactive'//{introspection:location})
    return
  end function inactiveLogPrior

  double precision function inactivePriorSample(self)
    !!{
    Sample from the of this parameter.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(modelParameterInactive), intent(inout) :: self
    !$GLC attributes unused :: self

    inactivePriorSample=0.0d0
    call Error_Report('parameter is inactive'//{introspection:location})
    return
  end function inactivePriorSample

  double precision function inactivePriorInvert(self,f)
    !!{
    Invert the prior, returning the parameter value given the cumulative probability.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (modelParameterInactive), intent(inout) :: self
    double precision                        , intent(in   ) :: f
    !$GLC attributes unused :: self, f

    inactivePriorInvert=0.0d0
    call Error_Report('parameter is inactive'//{introspection:location})
    return
  end function inactivePriorInvert

  double precision function inactivePriorMinimum(self)
    !!{
    Return the minimum value for which the prior is non-zero.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(modelParameterInactive), intent(inout) :: self
    !$GLC attributes unused :: self

    inactivePriorMinimum=0.0d0
    call Error_Report('parameter is inactive'//{introspection:location})
    return
  end function inactivePriorMinimum

  double precision function inactivePriorMaximum(self)
    !!{
    Return the maximum value for which the prior is non-zero.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(modelParameterInactive), intent(inout) :: self
    !$GLC attributes unused :: self

    inactivePriorMaximum=0.0d0
    call Error_Report('parameter is inactive'//{introspection:location})
    return
  end function inactivePriorMaximum

  double precision function inactiveRandomPerturbation(self)
    !!{
    Return a random perturbation to this parameter.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(modelParameterInactive), intent(inout) :: self
    !$GLC attributes unused :: self

    inactiveRandomPerturbation=0.0d0
    call Error_Report('parameter is inactive'//{introspection:location})
    return
  end function inactiveRandomPerturbation

  double precision function inactiveMap(self,x)
    !!{
    Map this parameter.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (modelParameterInactive), intent(inout) :: self
    double precision                        , intent(in   ) :: x
    !$GLC attributes unused :: self, x

    inactiveMap=0.0d0
    call Error_Report('parameter is inactive'//{introspection:location})
    return
  end function inactiveMap

  double precision function inactiveUnmap(self,x)
    !!{
    Unmap this parameter.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (modelParameterInactive), intent(inout) :: self
    double precision                        , intent(in   ) :: x
    !$GLC attributes unused :: self, x

    inactiveUnmap=0.0d0
    call Error_Report('parameter is inactive'//{introspection:location})
    return
  end function inactiveUnmap

