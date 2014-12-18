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

  !% Implements a null modifier of halo mass samples.
 
  !# <massFunctionSamplingModifier name="massFunctionSamplingModifierNull">
  !#  <description>Makes no modification of halo mass samples.</description>
  !# </massFunctionSamplingModifier>

  type, extends(massFunctionSamplingModifierClass) :: massFunctionSamplingModifierNull
     !% A class implementing halo mass sample modification which performs no modification.
     private
   contains
     procedure :: modify => nullModify
  end type massFunctionSamplingModifierNull

  interface massFunctionSamplingModifierNull
     !% Constructors for the null halo mass sample modifier class.
     module procedure nullDefaultConstructor
  end interface massFunctionSamplingModifierNull

contains

  function nullDefaultConstructor()
    !% Default constructor for the null halo mass sample modifier class.
    use Galacticus_Display
    use Input_Parameters
    implicit none
    type(massFunctionSamplingModifierNull) :: nullDefaultConstructor

    return
  end function nullDefaultConstructor

  subroutine nullModify(self,treeHaloMass,treeBaseTime)
    !% Perform no modification of a halo mass sample.
    implicit none
    class           (massFunctionSamplingModifierNull)                           , intent(inout) :: self
    double precision                                  , allocatable, dimension(:), intent(inout) :: treeHaloMass
    double precision                                                             , intent(in   ) :: treeBaseTime
    
    return
  end subroutine nullModify
  
