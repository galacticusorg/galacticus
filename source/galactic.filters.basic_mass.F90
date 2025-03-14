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
Implements a galactic high-pass filter for the default ``basic'' halo mass.
!!}

  !![
  <galacticFilter name="galacticFilterBasicMass">
   <description>
   A high-pass filter for basic mass. Halos with a basic mass mass greater than or equal to a fixed threshold,
   $M_0=${\normalfont \ttfamily [massThreshold]}.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterBasicMass
     !!{
     A galactic high-pass filter class for basic mass.
     !!}
     private
     double precision :: massThreshold
   contains
     procedure :: passes => basicMassPasses
  end type galacticFilterBasicMass

  interface galacticFilterBasicMass
     !!{
     Constructors for the ``basicMass'' galactic filter class.
     !!}
     module procedure basicMassConstructorParameters
     module procedure basicMassConstructorInternal
  end interface galacticFilterBasicMass

contains

  function basicMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``basicMass'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterBasicMass)                :: self
    type            (inputParameters        ), intent(inout) :: parameters
    double precision                                         :: massThreshold

    !![
    <inputParameter>
      <name>massThreshold</name>
      <source>parameters</source>
      <description>The parameter $M_0$ (in units of $M_\odot$) appearing in the basic mass threshold for the basic mass galactic filter class.</description>
    </inputParameter>
    !!]
    self=galacticFilterBasicMass(massThreshold)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function basicMassConstructorParameters

  function basicMassConstructorInternal(massThreshold) result(self)
    !!{
    Internal constructor for the ``basicMass'' galactic filter class.
    !!}
    implicit none
    type            (galacticFilterBasicMass)                :: self
    double precision                         , intent(in   ) :: massThreshold
    !![
    <constructorAssign variables="massThreshold"/>
    !!]
    
    return
  end function basicMassConstructorInternal

  logical function basicMassPasses(self,node)
    !!{
    Implement a  basic mass high-pass galactic filter.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(galacticFilterBasicMass), intent(inout)         :: self
    type (treeNode               ), intent(inout), target :: node
    class(nodeComponentBasic     ), pointer               :: basic

    basic           => node %basic        ()
    basicMassPasses =  basic%mass         () &
         &             >=                    &
         &             self %massThreshold
    return
  end function basicMassPasses
