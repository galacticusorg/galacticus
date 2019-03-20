!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implements a galactic structure initial radius class in which the initial and final radii are equal.
  
  !# <galacticStructureRadiiInitial name="galacticStructureRadiiInitialStatic">
  !#  <description>A galactic structure initial radius class in which the initial and final radii are equal.</description>
  !# </galacticStructureRadiiInitial>
  type, extends(galacticStructureRadiiInitialClass) :: galacticStructureRadiiInitialStatic
     !% A galactic structure initial radius class in which the initial and final radii are equal.
     private
   contains
     procedure :: radius           => staticRadius
     procedure :: radiusDerivative => staticRadiusDerivative
  end type galacticStructureRadiiInitialStatic

  interface galacticStructureRadiiInitialStatic
     !% Constructors for the {\normalfont \ttfamily static} galactic structure initial radius class.
     module procedure staticConstructorParameters
  end interface galacticStructureRadiiInitialStatic

contains

  function staticConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily static} galactic structure initial radius class which takes a parameter list as
    !% input.
    use Input_Parameters
    implicit none
    type(galacticStructureRadiiInitialStatic)                :: self
    type(inputParameters                    ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters
    
    self=galacticStructureRadiiInitialStatic()
    return
  end function staticConstructorParameters

  double precision function staticRadius(self,node,radius)
    !% Compute the initial radius in the dark matter halo assuming the halo is static.
    implicit none
    class           (galacticStructureRadiiInitialStatic), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    double precision                                     , intent(in   ) :: radius
    !GCC$ attributes unused :: self, node, radius

    staticRadius=radius
    return
  end function staticRadius

  double precision function staticRadiusDerivative(self,node,radius)
    !% Compute the derivative of the initial radius in the dark matter halo assuming the halo is static.
    implicit none
    class           (galacticStructureRadiiInitialStatic), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    double precision                                     , intent(in   ) :: radius
    !GCC$ attributes unused :: self, node, radius

    staticRadiusDerivative=1.0d0
    return
  end function staticRadiusDerivative
