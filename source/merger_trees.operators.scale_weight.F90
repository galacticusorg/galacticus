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

!!{
Implements a merger tree operator that scales the weight of each tree by a fixed factor.
!!}

  !![
  <mergerTreeOperator name="mergerTreeOperatorScaleWeight">
   <description>A merger tree operator that scales the weight of each tree by a fixed factor.</description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorScaleWeight
     !!{
     A merger tree operator that scales the weight of each tree by a fixed factor.
     !!}
     private
     double precision :: scaleFactor
   contains
     procedure :: operatePreInitialization => scaleFactorOperatePreInitialization
  end type mergerTreeOperatorScaleWeight

  interface mergerTreeOperatorScaleWeight
     !!{
     Constructors for the scaleWeight merger tree operator class.
     !!}
     module procedure scaleWeightConstructorParameters
     module procedure scaleWeightConstructorInternal
  end interface mergerTreeOperatorScaleWeight

contains

  function scaleWeightConstructorParameters(parameters) result(self)
    !!{
    Constructor for the scaleWeight merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (mergerTreeOperatorScaleWeight)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    double precision                                               :: scaleFactor
    
    !![
    <inputParameter>
      <name>scaleFactor</name>
      <description>The factor by which to scale merger tree weights.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=mergerTreeOperatorScaleWeight(scaleFactor)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function scaleWeightConstructorParameters

  function scaleWeightConstructorInternal(scaleFactor) result(self)
    !!{
    Internal constructor for the scaleWeight merger tree operator class.
    !!}
    implicit none
    type            (mergerTreeOperatorScaleWeight)                :: self
    double precision                               , intent(inout) :: scaleFactor
    !![
    <constructorAssign variables="scaleFactor"/>
    !!]
    
    return
  end function scaleWeightConstructorInternal

  subroutine scaleFactorOperatePreInitialization(self,tree)
    !!{
    Scale the weight of a merger tree by a fixed factor.
    !!}
    implicit none
    class(mergerTreeOperatorScaleWeight), intent(inout), target :: self
    type (mergerTree                   ), intent(inout), target :: tree

    tree%volumeWeight=+tree%volumeWeight &
         &            *self%scaleFactor
    return
  end subroutine scaleFactorOperatePreInitialization
