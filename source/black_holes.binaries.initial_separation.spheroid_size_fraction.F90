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
Implements a class for black hole binary initial separation in which the radius is a fixed fraction of the scale radius of the
larger of the host and satellite spheroids.
!!}

  !![
  <blackHoleBinaryInitialSeparation name="blackHoleBinaryInitialSeparationSpheroidRadiusFraction">
   <description>
    A black hole binary initial separation class that assumes that the initial separation of the binary is equal to a fixed
    fraction {\normalfont \ttfamily [spheroidRadiusFraction]} of the larger of the spheroid scale radii of the two merging
    galaxies.
   </description>
  </blackHoleBinaryInitialSeparation>
  !!]
  type, extends(blackHoleBinaryInitialSeparationClass) :: blackHoleBinaryInitialSeparationSpheroidRadiusFraction
     !!{
     A black hole binary initial separation class in which the radius is a fixed fraction of the scale radius of the larger of the host and satellite spheroids.
     !!}
     private
     double precision :: spheroidRadiusFraction
   contains
     procedure :: separationInitial => spheroidRadiusFractionSeparationInitial
  end type blackHoleBinaryInitialSeparationSpheroidRadiusFraction

  interface blackHoleBinaryInitialSeparationSpheroidRadiusFraction
     !!{
     Constructors for the \refClass{blackHoleBinaryInitialSeparationSpheroidRadiusFraction} black hole binary initial radius class.
     !!}
     module procedure spheroidRadiusFractionConstructorParameters
     module procedure spheroidRadiusFractionConstructorInternal
  end interface blackHoleBinaryInitialSeparationSpheroidRadiusFraction

contains

  function spheroidRadiusFractionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{blackHoleBinaryInitialSeparationSpheroidRadiusFraction} black hole binary recoiled class which takes a parameter list as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (blackHoleBinaryInitialSeparationSpheroidRadiusFraction)                :: self
    type            (inputParameters                                       ), intent(inout) :: parameters
    double precision                                                                        :: spheroidRadiusFraction

    !![
    <inputParameter>
      <name>spheroidRadiusFraction</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The fraction of the spheroid radius at which merging black holes will be initially placed.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=blackHoleBinaryInitialSeparationSpheroidRadiusFraction(spheroidRadiusFraction)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function spheroidRadiusFractionConstructorParameters

  function spheroidRadiusFractionConstructorInternal(spheroidRadiusFraction) result(self)
    !!{
    Constructor for the \refClass{blackHoleBinaryInitialSeparationSpheroidRadiusFraction} black hole binary recoiled class which takes a parameter list as
    input.
    !!}
    implicit none
    type            (blackHoleBinaryInitialSeparationSpheroidRadiusFraction)                :: self
    double precision                                                        , intent(in   ) :: spheroidRadiusFraction
    !![
    <constructorAssign variables="spheroidRadiusFraction"/>
    !!]

    return
  end function spheroidRadiusFractionConstructorInternal

  double precision function spheroidRadiusFractionSeparationInitial(self,node,nodeHost)
    !!{
    Returns an initial separation for a binary black holes that is a fixed fraction of the scale radius of the larger of the
    host and satellite spheroids.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    class(blackHoleBinaryInitialSeparationSpheroidRadiusFraction), intent(inout), target :: self
    type (treeNode                                              ), intent(inout), target :: nodeHost    , node
    class(nodeComponentSpheroid                                 ), pointer               :: spheroidHost, spheroid

    spheroid                                =>  node    %spheroid()
    spheroidHost                            =>  nodeHost%spheroid()
    spheroidRadiusFractionSeparationInitial =  +self%spheroidRadiusFraction &
         &                                     *max(                        &
         &                                          spheroid    %radius(),  &
         &                                          spheroidHost%radius()   &
         &                                         )
    return
  end function spheroidRadiusFractionSeparationInitial
