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

  !!{RST
  Implementation of a model parameter class in which the parameter value is derived from other parameters.
  !!}

  !![
  <modelParameter name="modelParameterDerived" docformat="rst">
   <description>
   A model parameter class in which the parameter value is not sampled directly but instead derived algorithmically from other model parameters via a mathematical expression. The parameter name is set by ``[name]`` and the derivation formula by ``[definition]``, enabling computed quantities to be tracked as named parameters during posterior sampling.
   </description>
  </modelParameter>
  !!]
  type, extends(modelParameterInactive) :: modelParameterDerived
     !!{RST
     Implementation of a model parameter class in which the parameter value is derived from other parameters.
     !!}
     private
     type   (varying_string) :: definition_
     logical                 :: isInteger_
   contains
     !![
     <methods>
       <method description="Return the definition for this parameter."   method="definition"/>
       <method description="Return true if the parameter is an integer." method="isInteger" />
     </methods>
     !!]
     procedure :: definition => derivedDefinition
     procedure :: isInteger  => derivedIsInteger
  end type modelParameterDerived

  interface modelParameterDerived
     !!{RST
     Constructors for the :galacticus-class:`modelParameterDerived` model parameter class.
     !!}
     module procedure derivedConstructorParameters
     module procedure derivedConstructorInternal
  end interface modelParameterDerived

contains

  function derivedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`modelParameterDerived` model parameter class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (modelParameterDerived)                :: self
    type   (inputParameters      ), intent(inout) :: parameters
    type   (varying_string       )                :: name      , definition
    logical                                       :: isInteger

    !![
    <inputParameter docformat="rst">
      <name>name</name>
      <description>
      The name of the derived model parameter as it appears in the Galacticus parameter file and posterior sampling output; used to identify this computed quantity in the parameter chain.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>definition</name>
      <description>
      The mathematical expression or formula used to derive this parameter's value from other model parameters during posterior sampling.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>isInteger</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true the derived value is treated as an integer, otherwise it is considered to be a float.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=modelParameterDerived(name,definition,isInteger)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function derivedConstructorParameters

  function derivedConstructorInternal(name_,definition_,isInteger_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`modelParameterDerived` model parameter class.
    !!}
    implicit none
    type   (modelParameterDerived)                :: self
    type   (varying_string       ), intent(in   ) :: name_     , definition_
    logical                       , intent(in   ) :: isInteger_
    !![
    <constructorAssign variables="name_, definition_, isInteger_"/>
    !!]

    return
  end function derivedConstructorInternal

  function derivedDefinition(self)
    !!{RST
    Return the definition of this parameter.
    !!}
    implicit none
    type (varying_string       )                :: derivedDefinition
    class(modelParameterDerived), intent(inout) :: self

    derivedDefinition=self%definition_
    return
  end function derivedDefinition

  logical function derivedIsInteger(self)
    !!{RST
    Return true if the derived parameter is to be treated as an integer.
    !!}
    implicit none
    class(modelParameterDerived), intent(inout) :: self

    derivedIsInteger=self%isInteger_
    return
  end function derivedIsInteger

