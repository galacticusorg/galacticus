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

!% Contains a module which implements a null output analysis class.

  !# <outputAnalysis name="outputAnalysisNull">
  !#  <description>A null output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisClass) :: outputAnalysisNull
     !% A null output analysis class.
     private
   contains
     procedure :: analyze  => nullAnalyze
     procedure :: finalize => nullFinalize
  end type outputAnalysisNull

  interface outputAnalysisNull
     !% Constructors for the ``null'' output analysis class.
     module procedure nullConstructorParameters
  end interface outputAnalysisNull

contains

  function nullConstructorParameters(parameters)
    !% Constructor for the ``null'' output analysis class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(outputAnalysisNull)                :: nullConstructorParameters
    type(inputParameters   ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    nullConstructorParameters=outputAnalysisNull()
    return
  end function nullConstructorParameters

  subroutine nullAnalyze(self,node,iOutput)
    !% Implement a null output analysis.
    implicit none
    class  (outputAnalysisNull), intent(inout) :: self
    type   (treeNode          ), intent(inout) :: node
    integer(c_size_t          ), intent(in   ) :: iOutput
    !GCC$ attributes unused :: self, node, iOutput
    
    return
  end subroutine nullAnalyze

  subroutine nullFinalize(self)
    !% Implement a null output analysis finalization.
    implicit none
    class(outputAnalysisNull), intent(inout) :: self
    !GCC$ attributes unused :: self
    
    return
  end subroutine nullFinalize
