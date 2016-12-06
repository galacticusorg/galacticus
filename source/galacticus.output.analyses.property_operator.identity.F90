!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements an identity output analysis property operator class.

  !# <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorIdentity">
  !#  <description>An identity output analysis property operator class.</description>
  !# </outputAnalysisPropertyOperator>
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorIdentity
     !% An identity output property operator class.
     private
   contains
     procedure :: operate  => identityOperate
  end type outputAnalysisPropertyOperatorIdentity

  interface outputAnalysisPropertyOperatorIdentity
     !% Constructors for the ``identity'' output analysis class.
     module procedure identityConstructorParameters
  end interface outputAnalysisPropertyOperatorIdentity

contains

  function identityConstructorParameters(parameters)
    !% Constructor for the ``identity'' output analysis property operateor class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(outputAnalysisPropertyOperatorIdentity)                :: identityConstructorParameters
    type(inputParameters                       ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    identityConstructorParameters=outputAnalysisPropertyOperatorIdentity()
    return
  end function identityConstructorParameters

  double precision function identityOperate(self,propertyValue)
    !% Implement an identity output analysis property operator.
    implicit none
    class           (outputAnalysisPropertyOperatorIdentity), intent(inout) :: self
    double precision                                        , intent(in   ) :: propertyValue
    !GCC$ attributes unused :: self

    identityOperate=propertyValue
    return
  end function identityOperate
