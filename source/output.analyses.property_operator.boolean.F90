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
Implements a boolean analysis property operator class.
!!}

  !![
  <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorBoolean">
   <description>A boolean analysis property operator class, specifically $x \rightarrow 0$ if $x=0$, and $x \rightarrow 1$ otherwise..</description>
  </outputAnalysisPropertyOperator>
  !!]
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorBoolean
     !!{
     A boolean property operator class.
     !!}
     private
     logical :: preciseZero
   contains
     procedure :: operate => booleanOperate
  end type outputAnalysisPropertyOperatorBoolean

  interface outputAnalysisPropertyOperatorBoolean
     !!{
     Constructors for the \refClass{outputAnalysisPropertyOperatorBoolean} output analysis class.
     !!}
     module procedure booleanConstructorParameters
     module procedure booleanConstructorInternal
  end interface outputAnalysisPropertyOperatorBoolean

contains

  function booleanConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisPropertyOperatorBoolean} output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(outputAnalysisPropertyOperatorBoolean)                :: self
    type(inputParameters                      ), intent(inout) :: parameters
    logical :: preciseZero
    
    !![
    <inputParameter>
      <name>preciseZero</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, then input value of 0 will be mapped to precisely 0. Otherwise, they are mapped to the smallest representable non-zero value, {\normalfont \ttfamily epsilon(0.0d0)}. This is useful since a precise 0 is treated differently by the \refClass{outputAnalysisMeanFunction1D} class.</description>
    </inputParameter>
    !!]
    self=outputAnalysisPropertyOperatorBoolean(preciseZero)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function booleanConstructorParameters

  function booleanConstructorInternal(preciseZero) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisPropertyOperatorBoolean} output analysis property operator class.
    !!}
    implicit none
    type   (outputAnalysisPropertyOperatorBoolean)                :: self
    logical                                       , intent(in   ) :: preciseZero
    !![
    <constructorAssign variables="preciseZero"/>   
    !!]
   
    return
  end function booleanConstructorInternal

  double precision function booleanOperate(self,propertyValue,node,propertyType,outputIndex)
    !!{
    Implement an boolean output analysis property operator.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (outputAnalysisPropertyOperatorBoolean    ), intent(inout)           :: self
    double precision                                           , intent(in   )           :: propertyValue
    type            (treeNode                                 ), intent(inout), optional :: node
    type            (enumerationOutputAnalysisPropertyTypeType), intent(inout), optional :: propertyType
    integer         (c_size_t                                 ), intent(in   ), optional :: outputIndex
    !$GLC attributes unused :: propertyType, outputIndex, node

    if (propertyValue == 0.0d0) then
       ! Optionally set the result to either precisely zero, or the smallest possible non-zero value.
       if (self%preciseZero) then
          booleanOperate=        0.0d0
       else
          booleanOperate=epsilon(0.0d0)
       end if
    else
       booleanOperate=1.0d0
    end if
    return
  end function booleanOperate
