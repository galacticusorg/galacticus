!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a fixed covering fraction for the dusty torus around an AGN.
  !!}

  !![
  <blackHoleTorusCoveringFraction name="blackHoleTorusCoveringFractionFixed" docformat="rst">
   <description>
   A covering fraction of the dusty torus around an AGN which is the same, ``[coveringFraction]``, for every black
   hole, whatever its luminosity.
   </description>
  </blackHoleTorusCoveringFraction>
  !!]
  type, extends(blackHoleTorusCoveringFractionClass) :: blackHoleTorusCoveringFractionFixed
     !!{RST
     A fixed covering fraction for the dusty torus around an AGN.
     !!}
     private
     double precision :: coveringFraction_
   contains
     procedure :: coveringFraction => fixedCoveringFraction
  end type blackHoleTorusCoveringFractionFixed

  interface blackHoleTorusCoveringFractionFixed
     !!{RST
     Constructors for the :galacticus-class:`blackHoleTorusCoveringFractionFixed` torus covering fraction class.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface blackHoleTorusCoveringFractionFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`blackHoleTorusCoveringFractionFixed` torus covering fraction class which
    takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (blackHoleTorusCoveringFractionFixed)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    double precision                                                     :: coveringFraction

    !![
    <inputParameter docformat="rst">
      <name>coveringFraction</name>
      <defaultValue>0.5d0</defaultValue>
      <description>
      The fraction of the sky, as seen from the accretion disk, covered by the dusty torus.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=blackHoleTorusCoveringFractionFixed(coveringFraction)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(coveringFraction_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`blackHoleTorusCoveringFractionFixed` torus covering fraction class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (blackHoleTorusCoveringFractionFixed)                :: self
    double precision                                     , intent(in   ) :: coveringFraction_
    !![
    <constructorAssign variables="coveringFraction_"/>
    !!]

    if (coveringFraction_ < 0.0d0 .or. coveringFraction_ > 1.0d0) call Error_Report('`coveringFraction` must lie between zero and one'//{introspection:location})
    return
  end function fixedConstructorInternal

  double precision function fixedCoveringFraction(self,node) result(coveringFraction)
    !!{RST
    Return the fixed covering fraction of the dusty torus.
    !!}
    implicit none
    class(blackHoleTorusCoveringFractionFixed), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node
    !$GLC attributes unused :: node

    coveringFraction=self%coveringFraction_
    return
  end function fixedCoveringFraction
