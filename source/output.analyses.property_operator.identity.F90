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
Implements an identity output analysis property operator class.
!!}

  !![
  <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorIdentity">
   <description>An identity output analysis property operator class.</description>
  </outputAnalysisPropertyOperator>
  !!]
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorIdentity
     !!{
     An identity output property operator class.
     !!}
     private
   contains
     procedure :: operate  => identityOperate
  end type outputAnalysisPropertyOperatorIdentity

  interface outputAnalysisPropertyOperatorIdentity
     !!{
     Constructors for the {\normalfont \ttfamily identity} output analysis class.
     !!}
     module procedure identityConstructorParameters
  end interface outputAnalysisPropertyOperatorIdentity

contains

  function identityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily identity} output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(outputAnalysisPropertyOperatorIdentity)                :: self
    type(inputParameters                       ), intent(inout) :: parameters

    self=outputAnalysisPropertyOperatorIdentity()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function identityConstructorParameters

  double precision function identityOperate(self,propertyValue,node,propertyType,outputIndex)
    !!{
    Implement an identity output analysis property operator.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (outputAnalysisPropertyOperatorIdentity   ), intent(inout)           :: self
    double precision                                           , intent(in   )           :: propertyValue
    type            (treeNode                                 ), intent(inout), optional :: node
    type            (enumerationOutputAnalysisPropertyTypeType), intent(inout), optional :: propertyType
    integer         (c_size_t                                 ), intent(in   ), optional :: outputIndex
    !$GLC attributes unused :: self, outputIndex, propertyType, node

    identityOperate=propertyValue
    return
  end function identityOperate
