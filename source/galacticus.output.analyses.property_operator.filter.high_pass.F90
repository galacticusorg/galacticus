!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements a high-pass filter analysis property operator class.

  !# <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorFilterHighPass">
  !#  <description>A high-pass filter analysis property operator class.</description>
  !# </outputAnalysisPropertyOperator>
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorFilterHighPass
     !% A high-pass filter property operator class.
     private
     double precision :: filterThreshold
   contains
     procedure :: operate => filterHighPassOperate
  end type outputAnalysisPropertyOperatorFilterHighPass

  interface outputAnalysisPropertyOperatorFilterHighPass
     !% Constructors for the ``filterHighPass'' output analysis class.
     module procedure filterHighPassConstructorParameters
     module procedure filterHighPassConstructorInternal
  end interface outputAnalysisPropertyOperatorFilterHighPass

contains

  function filterHighPassConstructorParameters(parameters) result(self)
    !% Constructor for the ``filterHighPass'' output analysis property operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(outputAnalysisPropertyOperatorFilterHighPass)                :: self
    type(inputParameters                             ), intent(inout) :: parameters
    double precision                                                  :: filterThreshold

    ! Check and read parameters.
    !# <inputParameter>
    !#   <name>filterThreshold</name>
    !#   <source>parameters</source>
    !#   <variable>filterThreshold</variable>
    !#   <description>Threshold for the high-pass filter distribution operator.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    self=outputAnalysisPropertyOperatorFilterHighPass(filterThreshold)
    !# <inputParametersValidate source="parameters"/>
    return
  end function filterHighPassConstructorParameters

  function filterHighPassConstructorInternal(filterThreshold) result (self)
    !% Internal constructor for the ``filterHighPass'' output analysis distribution operator class.
    implicit none
    type            (outputAnalysisPropertyOperatorFilterHighPass)                :: self
    double precision                                              , intent(in   ) :: filterThreshold
    !# <constructorAssign variables="filterThreshold"/>

    return
  end function filterHighPassConstructorInternal

  double precision function filterHighPassOperate(self,propertyValue,node,propertyType,outputIndex)
    !% Implement an filterHighPass output analysis property operator.
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (outputAnalysisPropertyOperatorFilterHighPass), intent(inout)           :: self
    double precision                                              , intent(in   )           :: propertyValue
    type            (treeNode                                    ), intent(inout), optional :: node
    integer                                                       , intent(inout), optional :: propertyType
    integer         (c_size_t                                    ), intent(in   ), optional :: outputIndex
    !$GLC attributes unused :: propertyType, outputIndex, node

    if (propertyValue > self%filterThreshold) then
       filterHighPassOperate=propertyValue
    else
       filterHighPassOperate=0.0d0
    end if
    return
  end function filterHighPassOperate
