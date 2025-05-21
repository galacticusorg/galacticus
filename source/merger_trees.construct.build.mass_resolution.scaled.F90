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
  An implementation of a merger tree builder mass resolution which assumes a resolution which scales with tree mass.
  !!}

  !![
  <mergerTreeMassResolution name="mergerTreeMassResolutionScaled">
   <description>
    A merger tree mass resolution class which computes the mass resolution to be min({\normalfont \ttfamily
    [massResolutionMaximum]},max({\normalfont \ttfamily [massResolutionMinimum]},{\normalfont \ttfamily
    [massResolutionFractional]}$\times M_\mathrm{base}$)), where $M_\mathrm{base}$ is the base mass of the merger tree.
   </description>
  </mergerTreeMassResolution>
  !!]
  type, extends(mergerTreeMassResolutionClass) :: mergerTreeMassResolutionScaled
     !!{
     A merger tree mass resolution class which assumes a mass resolution that scales with tree mass.
     !!}
     private
     double precision :: massResolutionMinimum   , massResolutionMaximum, &
          &              massResolutionFractional
   contains
     procedure :: resolution => scaledResolution
  end type mergerTreeMassResolutionScaled

  interface mergerTreeMassResolutionScaled
     !!{
     Constructors for the \refClass{mergerTreeMassResolutionScaled} merger tree resolution class.
     !!}
     module procedure scaledConstructorParameters
     module procedure scaledConstructorInternal
  end interface mergerTreeMassResolutionScaled

contains

  function scaledConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeMassResolutionScaled} merger tree building mass resolution class which reads parameters from a
    provided parameter list.
    !!}
    implicit none
    type            (mergerTreeMassResolutionScaled)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    double precision                                                :: massResolutionMinimum    , massResolutionMaximum, &
         &                                                             massResolutionFractional

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>massResolutionMinimum</name>
      <source>parameters</source>
      <defaultValue>5.0d9</defaultValue>
      <description>The minimum mass resolution to use when building merger trees.</description>
    </inputParameter>
    <inputParameter>
      <name>massResolutionMaximum</name>
      <source>parameters</source>
      <defaultValue>huge(0.0d0)</defaultValue>
      <description>The maximum mass resolution to use when building merger trees.</description>
    </inputParameter>
    <inputParameter>
      <name>massResolutionFractional</name>
      <source>parameters</source>
      <defaultValue>1.0d-3</defaultValue>
      <description>The fraction of the tree's root node mass to be used for the mass resolution when building merger trees.</description>
    </inputParameter>
    !!]
    self=mergerTreeMassResolutionScaled(massResolutionMinimum,massResolutionMaximum,massResolutionFractional)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function scaledConstructorParameters

  function scaledConstructorInternal(massResolutionMinimum,massResolutionMaximum,massResolutionFractional) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeMassResolutionScaled} merger tree building mass resolution class.
    !!}
    implicit none
    type            (mergerTreeMassResolutionScaled)                :: self
    double precision                                , intent(in   ) :: massResolutionMinimum    , massResolutionMaximum, &
         &                                                             massResolutionFractional
    !![
    <constructorAssign variables="massResolutionMinimum, massResolutionMaximum, massResolutionFractional"/>
    !!]

    return
  end function scaledConstructorInternal

  double precision function scaledResolution(self,tree)
    !!{
    Returns a scaled mass resolution to use when building merger trees.
    !!}
    use :: Galacticus_Nodes, only : mergerTree, nodeComponentBasic
    implicit none
    class(mergerTreeMassResolutionScaled), intent(inout) :: self
    type (mergerTree                    ), intent(in   ) :: tree
    class(nodeComponentBasic            ), pointer       :: basicBase

    basicBase        => tree%nodeBase%basic()
    scaledResolution =  max(                                                    &
         &                      self%massResolutionMinimum                    , &
         &                  min(                                                &
         &                      self%massResolutionMaximum                    , &
         &                      self%massResolutionFractional*basicBase%mass()  &
         &                     )                                                &
         &                 )
    return
  end function scaledResolution
