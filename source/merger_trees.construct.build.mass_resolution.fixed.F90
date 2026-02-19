!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  An implementation of a merger tree builder mass resolution which assumes a fixed resolution.
  !!}

  !![
  <mergerTreeMassResolution name="mergerTreeMassResolutionFixed">
   <description>
    A merger tree mass resolution class which assumes a fixed mass resolution of {\normalfont \ttfamily [massResolution]} for
    all merger trees.
   </description>
  </mergerTreeMassResolution>
  !!]
  type, extends(mergerTreeMassResolutionClass) :: mergerTreeMassResolutionFixed
     !!{
     A merger tree mass resolution class which assumes a fixed mass resolution.
     !!}
     private
     double precision :: massResolution
   contains
     procedure :: resolution => fixedResolution
  end type mergerTreeMassResolutionFixed

  interface mergerTreeMassResolutionFixed
     !!{
     Constructors for the \refClass{mergerTreeMassResolutionFixed} merger tree resolution class.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface mergerTreeMassResolutionFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeMassResolutionFixed} merger tree building mass resolution class which reads parameters from a
    provided parameter list.
    !!}
    implicit none
    type            (mergerTreeMassResolutionFixed)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    double precision                                               :: massResolution

    !![
    <inputParameter>
      <name>massResolution</name>
      <source>parameters</source>
      <defaultValue>5.0d9</defaultValue>
      <description>The mass resolution to use when building merger trees.</description>
    </inputParameter>
    !!]
    self=mergerTreeMassResolutionFixed(massResolution)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(massResolution) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeMassResolutionFixed} merger tree building mass resolution class.
    !!}
    implicit none
    type            (mergerTreeMassResolutionFixed)                :: self
    double precision                               , intent(in   ) :: massResolution
    !![
    <constructorAssign variables="massResolution"/>
    !!]
    
    return
  end function fixedConstructorInternal

  double precision function fixedResolution(self,tree)
    !!{
    Returns a fixed mass resolution to use when building merger trees.
    !!}
    implicit none
    class(mergerTreeMassResolutionFixed), intent(inout) :: self
    type (mergerTree                   ), intent(in   ) :: tree
    !$GLC attributes unused :: tree

    fixedResolution=self%massResolution
    return
  end function fixedResolution
