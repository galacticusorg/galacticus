!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements a conversion of ISM mass to HI mass analysis property operator class.

  use Output_Analysis_Molecular_Ratios
  
  !# <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorHIMass">
  !#  <description>A conversion of ISM mass to HI mass analysis property operator class.</description>
  !# </outputAnalysisPropertyOperator>
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorHIMass
     !% A conversion of ISM mass to HI mass property operator class.
     private
     class(outputAnalysisMolecularRatioClass), pointer :: outputAnalysisMolecularRatio_
   contains
     procedure :: operate => hiMassOperate
  end type outputAnalysisPropertyOperatorHIMass

  interface outputAnalysisPropertyOperatorHIMass
     !% Constructors for the ``hiMass'' output analysis class.
     module procedure hiMassConstructorParameters
     module procedure hiMassConstructorInternal
  end interface outputAnalysisPropertyOperatorHIMass

contains

  function hiMassConstructorParameters(parameters) result(self)
    !% Constructor for the ``hiMass'' output analysis property operator class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type (outputAnalysisPropertyOperatorHIMass)                :: self
    type (inputParameters                     ), intent(inout) :: parameters
    class(outputAnalysisMolecularRatioClass   ), pointer       :: outputAnalysisMolecularRatio_

    ! Check and read parameters.
    !# <objectBuilder class="outputAnalysisMolecularRatio" name="outputAnalysisMolecularRatio_" source="parameters" />
    self=outputAnalysisPropertyOperatorHIMass(outputAnalysisMolecularRatio_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function hiMassConstructorParameters

  function hiMassConstructorInternal(outputAnalysisMolecularRatio_) result (self)
    !% Internal constructor for the ``hiMass'' output analysis distribution operator class.
    use Input_Parameters
    implicit none
    type (outputAnalysisPropertyOperatorHIMass)                        :: self
    class(outputAnalysisMolecularRatioClass   ), intent(in   ), target :: outputAnalysisMolecularRatio_
    !# <constructorAssign variables="*outputAnalysisMolecularRatio_"/>

    return
  end function hiMassConstructorInternal

  double precision function hiMassOperate(self,propertyValue,node,propertyType,outputIndex)
    !% Implement an hiMass output analysis property operator.
    use, intrinsic :: ISO_C_Binding
    use            :: Galacticus_Error
    use            :: Numerical_Constants_Astronomical
    implicit none
    class           (outputAnalysisPropertyOperatorHIMass), intent(inout)           :: self
    double precision                                      , intent(in   )           :: propertyValue
    type            (treeNode                            ), intent(inout), optional :: node
    integer                                               , intent(inout), optional :: propertyType
    integer         (c_size_t                            ), intent(in   ), optional :: outputIndex
    !GCC$ attributes unused :: propertyType, outputIndex

    if (.not.present(node)) call Galacticus_Error_Report('node must be provided'//{introspection:location})
    hiMassOperate=+propertyValue                                                &
         &        *self%outputAnalysisMolecularRatio_%ratio(propertyValue,node)
    return
  end function hiMassOperate
