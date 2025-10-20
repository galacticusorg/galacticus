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
Implements a galactic high-pass filter for total black hole mass.
!!}

  !![
  <galacticFilter name="galacticFilterBlackHoleMass">
   <description>
   A galactic high-pass filter for black hole mass. Galaxies with a central black hole mass greater than or equal
   to a fixed threshold, $M_{\bullet,0}=${\normalfont \ttfamily [massThreshold]}, are passed.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterBlackHoleMass
     !!{
     A galactic high-pass filter class for black hole mass.
     !!}
     private
     double precision :: massThreshold
   contains
     procedure :: passes => blackHoleMassPasses
  end type galacticFilterBlackHoleMass

  interface galacticFilterBlackHoleMass
     !!{
     Constructors for the \refClass{galacticFilterBlackHoleMass} galactic filter class.
     !!}
     module procedure blackHoleMassConstructorParameters
     module procedure blackHoleMassConstructorInternal
  end interface galacticFilterBlackHoleMass

contains

  function blackHoleMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterBlackHoleMass} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterBlackHoleMass)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    double precision                                             :: massThreshold

    !![
    <inputParameter>
      <name>massThreshold</name>
      <source>parameters</source>
      <description>The parameter $M_0$ (in units of $M_\odot$) appearing in the black hole mass threshold for the black hole mass galactic filter class.</description>
    </inputParameter>
    !!]
    self=galacticFilterBlackHoleMass(massThreshold)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function blackHoleMassConstructorParameters

  function blackHoleMassConstructorInternal(massThreshold) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterBlackHoleMass} galactic filter class.
    !!}
    implicit none
    type            (galacticFilterBlackHoleMass)                :: self
    double precision                             , intent(in   ) :: massThreshold
    !![
    <constructorAssign variables="massThreshold"/>
    !!]
    
    return
  end function blackHoleMassConstructorInternal

  logical function blackHoleMassPasses(self,node) result(passes)
    !!{
    Implement a black hole mass high-pass galactic filter.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    class(galacticFilterBlackHoleMass), intent(inout)          :: self
    type (treeNode                   ), intent(inout), target  :: node
    class(nodeComponentBlackHole     )               , pointer :: blackHole

    blackHole =>  node     %blackHole    ()
    passes    =   blackHole%mass         () &
         &       >=                         &
         &        self     %massThreshold
    return
  end function blackHoleMassPasses
