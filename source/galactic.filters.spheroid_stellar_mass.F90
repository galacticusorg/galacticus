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
Implements a galactic high-pass filter for spheroid stellar mass.
!!}

  !![
  <galacticFilter name="galacticFilterSpheroidStellarMass">
   <description>
   A galactic high-pass filter for stellar mass. Galaxies with a spheroid stellar mass greater than or equal
   to a fixed threshold, $M_{\star,0}=${\normalfont \ttfamily [massThreshold]}.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterSpheroidStellarMass
     !!{
     A galactic high-pass filter class for spheroid stellar mass.
     !!}
     private
     double precision :: massThreshold
   contains
     procedure :: passes => spheroidStellarMassPasses
  end type galacticFilterSpheroidStellarMass

  interface galacticFilterSpheroidStellarMass
     !!{
     Constructors for the ``spheroidStellarMass'' galactic filter class.
     !!}
     module procedure spheroidStellarMassConstructorParameters
     module procedure spheroidStellarMassConstructorInternal
  end interface galacticFilterSpheroidStellarMass

contains

  function spheroidStellarMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``spheroidStellarMass'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterSpheroidStellarMass)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    double precision                                                   :: massThreshold

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>massThreshold</name>
      <source>parameters</source>
      <description>The parameter $M_0$ (in units of $M_\odot$) appearing in the stellar mass threshold for the spheroid stellar mass galactic filter class.</description>
    </inputParameter>
    !!]
    self=galacticFilterSpheroidStellarMass(massThreshold)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function spheroidStellarMassConstructorParameters

  function spheroidStellarMassConstructorInternal(massThreshold) result(self)
    !!{
    Internal constructor for the ``spheroidStellarMass'' galactic filter class.
    !!}
    implicit none
    type            (galacticFilterSpheroidStellarMass)                :: self
    double precision                                   , intent(in   ) :: massThreshold
    !![
    <constructorAssign variables="massThreshold"/>
    !!]

    return
  end function spheroidStellarMassConstructorInternal

  logical function spheroidStellarMassPasses(self,node)
    !!{
    Implement a  stellar mass high-pass galactic filter.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    class           (galacticFilterSpheroidStellarMass), intent(inout)         :: self
    type            (treeNode                         ), intent(inout), target :: node
    class           (nodeComponentSpheroid            ), pointer               :: spheroid
    double precision                                                           :: spheroidStellarMass

    spheroid                  =>  node    %spheroid   ()
    spheroidStellarMass       =  +spheroid%massStellar()
    spheroidStellarMassPasses =   spheroidStellarMass    &
         &                       >=                      &
         &                        self%massThreshold
    return
  end function spheroidStellarMassPasses
