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
Implements a property extractor class for bar instability timescales.
!!}

  use :: Galactic_Dynamics_Bar_Instabilities, only : galacticDynamicsBarInstabilityClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorTimescaleBarInstability">
   <description>A property extractor class for bar instability timescales.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorTimescaleBarInstability
     !!{
     A property extractor class for bar instability timescales.
     !!}
     private
     class(galacticDynamicsBarInstabilityClass), pointer :: galacticDynamicsBarInstability_ => null()
   contains
     final     ::                timescaleBarInstabilityDestructor
     procedure :: extract     => timescaleBarInstabilityExtract
     procedure :: name        => timescaleBarInstabilityName
     procedure :: description => timescaleBarInstabilityDescription
     procedure :: unitsInSI   => timescaleBarInstabilityUnitsInSI
  end type nodePropertyExtractorTimescaleBarInstability

  interface nodePropertyExtractorTimescaleBarInstability
     !!{
     Constructors for the \refClass{nodePropertyExtractorTimescaleBarInstability} output analysis class.
     !!}
     module procedure timescaleBarInstabilityConstructorParameters
     module procedure timescaleBarInstabilityConstructorInternal
  end interface nodePropertyExtractorTimescaleBarInstability

contains

  function timescaleBarInstabilityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorTimescaleBarInstability} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorTimescaleBarInstability)                :: self
    type (inputParameters                             ), intent(inout) :: parameters
    class(galacticDynamicsBarInstabilityClass         ), pointer       :: galacticDynamicsBarInstability_

    !![
    <objectBuilder class="galacticDynamicsBarInstability" name="galacticDynamicsBarInstability_" source="parameters"/>
    !!]
    self=nodePropertyExtractorTimescaleBarInstability(galacticDynamicsBarInstability_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticDynamicsBarInstability_"/>
    !!]
    return
  end function timescaleBarInstabilityConstructorParameters

  function timescaleBarInstabilityConstructorInternal(galacticDynamicsBarInstability_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorTimescaleBarInstability} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorTimescaleBarInstability)                        :: self
    class(galacticDynamicsBarInstabilityClass         ), intent(in   ), target :: galacticDynamicsBarInstability_
    !![
    <constructorAssign variables="*galacticDynamicsBarInstability_"/>
    !!]

    return
  end function timescaleBarInstabilityConstructorInternal

  subroutine timescaleBarInstabilityDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorTimescaleBarInstability} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorTimescaleBarInstability), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticDynamicsBarInstability_"/>
    !!]
    return
  end subroutine timescaleBarInstabilityDestructor

  double precision function timescaleBarInstabilityExtract(self,node,instance)
    !!{
    Implement extraction of the bar instability timescale.
    !!}
    implicit none
    class           (nodePropertyExtractorTimescaleBarInstability), intent(inout), target   :: self
    type            (treeNode                                    ), intent(inout), target   :: node
    type            (multiCounter                                ), intent(inout), optional :: instance
    double precision                                                                        :: fractionAngularMomentumRetainedDisk, fractionAngularMomentumRetainedSpheroid, &
         &                                                                                     torqueSpecificExternal
    !$GLC attributes unused :: instance

    call self%galacticDynamicsBarInstability_%timescale(node,timescaleBarInstabilityExtract,torqueSpecificExternal,fractionAngularMomentumRetainedDisk,fractionAngularMomentumRetainedSpheroid)
    return
  end function timescaleBarInstabilityExtract

  function timescaleBarInstabilityName(self)
    !!{
    Return the name of the bar instability timescale property.
    !!}
    implicit none
    type (varying_string                              )                :: timescaleBarInstabilityName
    class(nodePropertyExtractorTimescaleBarInstability), intent(inout) :: self
    !$GLC attributes unused :: self

    timescaleBarInstabilityName=var_str('timescaleBarInstabilityDisk')
    return
  end function timescaleBarInstabilityName

  function timescaleBarInstabilityDescription(self)
    !!{
    Return a description of the bar instability timescale property.
    !!}
    implicit none
    type (varying_string                              )                :: timescaleBarInstabilityDescription
    class(nodePropertyExtractorTimescaleBarInstability), intent(inout) :: self
    !$GLC attributes unused :: self

    timescaleBarInstabilityDescription=var_str('The timescale for bar instability in the disk [Gyr].')
    return
  end function timescaleBarInstabilityDescription

  double precision function timescaleBarInstabilityUnitsInSI(self)
    !!{
    Return the units of the bar instability timescale property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    class(nodePropertyExtractorTimescaleBarInstability), intent(inout) :: self
    !$GLC attributes unused :: self

    timescaleBarInstabilityUnitsInSI=gigaYear
    return
  end function timescaleBarInstabilityUnitsInSI


