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

!% Contains a module which implements an log10 output analysis property operator class.

  !# <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorLog10">
  !#  <description>An log10 output analysis property operator class.</description>
  !# </outputAnalysisPropertyOperator>
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorLog10
     !% An log10 output property operator class.
     private
   contains
     procedure :: operate => log10Operate
  end type outputAnalysisPropertyOperatorLog10

  interface outputAnalysisPropertyOperatorLog10
     !% Constructors for the ``log10'' output analysis class.
     module procedure log10ConstructorParameters
  end interface outputAnalysisPropertyOperatorLog10

contains

  function log10ConstructorParameters(parameters)
    !% Constructor for the ``log10'' output analysis property operator class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(outputAnalysisPropertyOperatorLog10)                :: log10ConstructorParameters
    type(inputParameters                    ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    log10ConstructorParameters=outputAnalysisPropertyOperatorLog10()
    return
  end function log10ConstructorParameters

  double precision function log10Operate(self,propertyValue)
    !% Implement an log10 output analysis property operator.
    implicit none
    class           (outputAnalysisPropertyOperatorLog10), intent(inout) :: self
    double precision                                     , intent(in   ) :: propertyValue
    !GCC$ attributes unused :: self

    if (propertyValue > 0.0d0) then
       log10Operate=log10(propertyValue)
    else
       log10Operate=-huge(1.0d0)
    end if
    return
  end function log10Operate
