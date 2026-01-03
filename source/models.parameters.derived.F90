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
  Implementation of a model parameter class in which the parameter value is derived from other parameters.
  !!}

  !![
  <modelParameter name="modelParameterDerived">
   <description>A model parameter class in which the parameter value is derived from other parameters.</description>
  </modelParameter>
  !!]
  type, extends(modelParameterInactive) :: modelParameterDerived
     !!{
     Implementation of a model parameter class in which the parameter value is derived from other parameters.
     !!}
     private
     type(varying_string) :: definition_
   contains
     !![
     <methods>
       <method description="Return the definition for this parameter." method="definition" />
     </methods>
     !!]
     procedure :: definition => derivedDefinition
  end type modelParameterDerived

  interface modelParameterDerived
     !!{
     Constructors for the \refClass{modelParameterDerived} 1D distribution function class.
     !!}
     module procedure derivedConstructorParameters
     module procedure derivedConstructorInternal
  end interface modelParameterDerived

contains

  function derivedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{modelParameterDerived} model parameter class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(modelParameterDerived)                :: self
    type(inputParameters      ), intent(inout) :: parameters
    type(varying_string       )                :: name      , definition

    !![
    <inputParameter>
      <name>name</name>
      <description>The name of the parameter.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>definition</name>
      <description>The definition of the parameter.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=modelParameterDerived(name,definition)
     !![
     <inputParametersValidate source="parameters"/>
     !!]
   return
  end function derivedConstructorParameters

  function derivedConstructorInternal(name_,definition_) result(self)
    !!{
    Internal constructor for the \refClass{modelParameterDerived} model parameter class.
    !!}
    implicit none
    type(modelParameterDerived)                :: self
    type(varying_string       ), intent(in   ) :: name_, definition_
    !![
    <constructorAssign variables="name_, definition_"/>
    !!]

    return
  end function derivedConstructorInternal

  function derivedDefinition(self)
    !!{
    Return the definition of this parameter.
    !!}
    implicit none
    type (varying_string       )                :: derivedDefinition
    class(modelParameterDerived), intent(inout) :: self

    derivedDefinition=self%definition_
    return
  end function derivedDefinition

