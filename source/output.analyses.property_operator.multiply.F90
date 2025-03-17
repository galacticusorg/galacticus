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
Implements a multiplication analysis property operator class.
!!}

  !![
  <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorMultiply">
   <description>A high-pass filter analysis property operator class.</description>
  </outputAnalysisPropertyOperator>
  !!]
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorMultiply
     !!{
     A multiplication property operator class.
     !!}
     private
     double precision :: multiplier
   contains
     procedure :: operate => multiplyOperate
  end type outputAnalysisPropertyOperatorMultiply

  interface outputAnalysisPropertyOperatorMultiply
     !!{
     Constructors for the ``multiply'' output analysis class.
     !!}
     module procedure multiplyConstructorParameters
     module procedure multiplyConstructorInternal
  end interface outputAnalysisPropertyOperatorMultiply

contains

  function multiplyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``multiply'' output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(outputAnalysisPropertyOperatorMultiply)                :: self
    type(inputParameters                       ), intent(inout) :: parameters
    double precision                                            :: multiplier

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>multiplier</name>
      <source>parameters</source>
      <variable>multiplier</variable>
      <description>Multiplying factor.</description>
    </inputParameter>
    !!]
    self=outputAnalysisPropertyOperatorMultiply(multiplier)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function multiplyConstructorParameters

  function multiplyConstructorInternal(multiplier) result (self)
    !!{
    Internal constructor for the ``multiply'' output analysis distribution operator class.
    !!}
    implicit none
    type            (outputAnalysisPropertyOperatorMultiply)                :: self
    double precision                                        , intent(in   ) :: multiplier
    !![
    <constructorAssign variables="multiplier"/>
    !!]

    return
  end function multiplyConstructorInternal

  double precision function multiplyOperate(self,propertyValue,node,propertyType,outputIndex)
    !!{
    Implement an multiply output analysis property operator.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (outputAnalysisPropertyOperatorMultiply   ), intent(inout)           :: self
    double precision                                           , intent(in   )           :: propertyValue
    type            (treeNode                                 ), intent(inout), optional :: node
    type            (enumerationOutputAnalysisPropertyTypeType), intent(inout), optional :: propertyType
    integer         (c_size_t                                 ), intent(in   ), optional :: outputIndex
    !$GLC attributes unused :: propertyType, outputIndex, node

    multiplyOperate=+propertyValue   &
         &          *self%multiplier
    return
  end function multiplyOperate
