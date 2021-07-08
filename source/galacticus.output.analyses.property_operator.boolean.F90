!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements a boolean analysis property operator class.
!!}

  !![
  <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorBoolean">
   <description>A boolean analysis property operator class, specifically $x \rightarrow x/|x|$, that is, the operator maintains the sign of the input while normalizing the magnitude to unity (or zero for zero input).</description>
  </outputAnalysisPropertyOperator>
  !!]
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorBoolean
     !!{
     A boolean property operator class.
     !!}
     private
   contains
     procedure :: operate => booleanOperate
  end type outputAnalysisPropertyOperatorBoolean

  interface outputAnalysisPropertyOperatorBoolean
     !!{
     Constructors for the ``boolean'' output analysis class.
     !!}
     module procedure booleanConstructorParameters
  end interface outputAnalysisPropertyOperatorBoolean

contains

  function booleanConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``boolean'' output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(outputAnalysisPropertyOperatorBoolean)                :: self
    type(inputParameters                      ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=outputAnalysisPropertyOperatorBoolean()
    return
  end function booleanConstructorParameters

  double precision function booleanOperate(self,propertyValue,node,propertyType,outputIndex)
    !!{
    Implement an boolean output analysis property operator.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (outputAnalysisPropertyOperatorBoolean), intent(inout)           :: self
    double precision                                       , intent(in   )           :: propertyValue
    type            (treeNode                             ), intent(inout), optional :: node
    integer                                                , intent(inout), optional :: propertyType
    integer         (c_size_t                             ), intent(in   ), optional :: outputIndex
    !$GLC attributes unused :: self, propertyType, outputIndex, node

    if (propertyValue == 0.0d0) then
       booleanOperate=     0.0d0
    else
       booleanOperate=sign(1.0d0,propertyValue)
    end if
    return
  end function booleanOperate
