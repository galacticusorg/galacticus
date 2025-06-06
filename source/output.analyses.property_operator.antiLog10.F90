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
Implements an anti-$\log_{10}()$ output analysis property operator class.
!!}

  !![
  <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorAntiLog10">
   <description>An anti-$\log_{10}()$ output analysis property operator class.</description>
  </outputAnalysisPropertyOperator>
  !!]
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorAntiLog10
     !!{
     An anti-$\log_{10}()$ output property operator class.
     !!}
     private
   contains
     procedure :: operate => antiLog10Operate
  end type outputAnalysisPropertyOperatorAntiLog10

  interface outputAnalysisPropertyOperatorAntiLog10
     !!{
     Constructors for the \refClass{outputAnalysisPropertyOperatorAntiLog10} output analysis class.
     !!}
     module procedure antiLog10ConstructorParameters
  end interface outputAnalysisPropertyOperatorAntiLog10

contains

  function antiLog10ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisPropertyOperatorAntiLog10} output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(outputAnalysisPropertyOperatorAntiLog10)                :: self
    type(inputParameters                        ), intent(inout) :: parameters

    self=outputAnalysisPropertyOperatorAntiLog10()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function antiLog10ConstructorParameters

  double precision function antiLog10Operate(self,propertyValue,node,propertyType,outputIndex)
    !!{
    Implement an antiLog10 output analysis property operator.
    !!}
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear, outputAnalysisPropertyTypeLog10, outputAnalysisPropertyTypeUnknown
    implicit none
    class           (outputAnalysisPropertyOperatorAntiLog10  ), intent(inout)           :: self
    double precision                                           , intent(in   )           :: propertyValue
    type            (treeNode                                 ), intent(inout), optional :: node
    type            (enumerationOutputAnalysisPropertyTypeType), intent(inout), optional :: propertyType
    integer         (c_size_t                                 ), intent(in   ), optional :: outputIndex
    !$GLC attributes unused :: self, outputIndex, node

    antiLog10Operate=10.0d0**propertyValue
    ! Change the property type.
    if (present(propertyType)) then
       if (propertyType == outputAnalysisPropertyTypeLog10) then
          propertyType=outputAnalysisPropertyTypeLinear
       else
          propertyType=outputAnalysisPropertyTypeUnknown
       end if
    end if
    return
  end function antiLog10Operate
