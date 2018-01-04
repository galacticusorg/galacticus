!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% Contains a module which implements a property operator class which transforms into various definitions of metallicity.
  
  !# <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorMetallicity">
  !#  <description>A property operator class in which the property value is replaced with an integral over a metallicity distribution between given limits, using the property value at the mean of the distribution.</description>
  !# </outputAnalysisPropertyOperator>
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorMetallicity
     !% A metallicity property operator class.
     private
     integer          :: type
     double precision :: valueSolar
   contains
     procedure :: operate => metallicityOperate
  end type outputAnalysisPropertyOperatorMetallicity

  interface outputAnalysisPropertyOperatorMetallicity
     !% Constructors for the ``metallicity'' output analysis class.
     module procedure metallicityConstructorParameters
     module procedure metallicityConstructorInternal
  end interface outputAnalysisPropertyOperatorMetallicity

contains

  function metallicityConstructorParameters(parameters) result(self)
    !% Constructor for the ``metallicity'' output analysis property operator class which takes a parameter set as input.
    use Abundances_Structure
    use Input_Parameters
    implicit none
    type            (outputAnalysisPropertyOperatorMetallicity)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    type            (varying_string                           )                :: type
    double precision                                                           :: valueSolar

    ! Check and read parameters.
    !# <inputParameter>
    !#   <name>type</name>
    !#   <source>parameters</source>
    !#   <description>The type of metallicity to transform into.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>valueSolar</name>
    !#   <source>parameters</source>
    !#   <description>The value at Solar metallicity</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    self=outputAnalysisPropertyOperatorMetallicity(enumerationMetallicityTypeEncode(char(type),includesPrefix=.false.),valueSolar)
    !# <inputParametersValidate source="parameters"/>
    return
  end function metallicityConstructorParameters

  function metallicityConstructorInternal(type,valueSolar) result (self)
    !% Internal constructor for the ``metallicity'' output analysis distribution operator class.
    implicit none
    type            (outputAnalysisPropertyOperatorMetallicity)                :: self
    integer                                                    , intent(in   ) :: type
    double precision                                           , intent(in   ) :: valueSolar
    !# <constructorAssign variables="type, valueSolar"/>

    return
  end function metallicityConstructorInternal

  double precision function metallicityOperate(self,propertyValue,node,propertyType,outputIndex)
    !% Implement an metallicity output analysis property operator.
    use, intrinsic :: ISO_C_Binding
    use               Abundances_Structure
    use               Galacticus_Error
    use               Numerical_Constants_Astronomical
    implicit none
    class           (outputAnalysisPropertyOperatorMetallicity), intent(inout)           :: self
    double precision                                           , intent(in   )           :: propertyValue
    type            (treeNode                                 ), intent(inout), optional :: node
    integer                                                    , intent(inout), optional :: propertyType
    integer         (c_size_t                                 ), intent(in   ), optional :: outputIndex
    !GCC$ attributes unused :: propertyType, outputIndex, node

    select case (self%type)
    case (metallicityTypeLogarithmicByNumberSolarPlus12)
       if (propertyValue > 0.0d0) then
          metallicityOperate=log10(propertyValue/metallicitySolar)+self%valueSolar
       else
          metallicityOperate=-huge(0.0d0)
       end if
    case default
       metallicityOperate=0.0d0
       call Galacticus_Error_Report('unsupported metallicity type'//{introspection:location})
    end select
    return
  end function metallicityOperate
