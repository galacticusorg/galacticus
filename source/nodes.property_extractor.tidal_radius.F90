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
Implements a tidal radius property extractor class.
!!}

  use :: Satellite_Tidal_Stripping_Radii, only : satelliteTidalStrippingRadiusClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRadiusTidal">
   <description>
    A tidal radius property extractor class. Extracts the tidal radius in the halo in Mpc.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRadiusTidal
     !!{
     A tidal radius property extractor class.
     !!}
     private
     class(satelliteTidalStrippingRadiusClass), pointer :: satelliteTidalStrippingRadius_ => null()
   contains
     final     ::                radiusTidalDestructor
     procedure :: extract     => radiusTidalExtract
     procedure :: name        => radiusTidalName
     procedure :: description => radiusTidalDescription
     procedure :: unitsInSI   => radiusTidalUnitsInSI
  end type nodePropertyExtractorRadiusTidal

  interface nodePropertyExtractorRadiusTidal
     !!{
     Constructors for the \refClass{nodePropertyExtractorRadiusTidal} output analysis class.
     !!}
     module procedure radiusTidalConstructorParameters
     module procedure radiusTidalConstructorInternal
  end interface nodePropertyExtractorRadiusTidal

contains

  function radiusTidalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorRadiusTidal} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorRadiusTidal  )                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(satelliteTidalStrippingRadiusClass), pointer       :: satelliteTidalStrippingRadius_

    !![
    <objectBuilder class="satelliteTidalStrippingRadius" name="satelliteTidalStrippingRadius_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRadiusTidal(satelliteTidalStrippingRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteTidalStrippingRadius_"/>
    !!]
    return
  end function radiusTidalConstructorParameters

  function radiusTidalConstructorInternal(satelliteTidalStrippingRadius_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorRadiusTidal} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorRadiusTidal  )                        :: self
    class(satelliteTidalStrippingRadiusClass), intent(in   ), target :: satelliteTidalStrippingRadius_
    !![
    <constructorAssign variables="*satelliteTidalStrippingRadius_"/>
    !!]

    return
  end function radiusTidalConstructorInternal

  subroutine radiusTidalDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorRadiusTidal} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRadiusTidal), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteTidalStrippingRadius_"/>
    !!]
    return
  end subroutine radiusTidalDestructor

  double precision function radiusTidalExtract(self,node,instance)
    !!{
    Implement a tidal radius property extractor.
    !!}
    implicit none
    class(nodePropertyExtractorRadiusTidal), intent(inout), target   :: self
    type (treeNode                        ), intent(inout), target   :: node
    type (multiCounter                    ), intent(inout), optional :: instance
    !$GLC attributes unused :: instance

    radiusTidalExtract=self%satelliteTidalStrippingRadius_%radius(node)
    return
  end function radiusTidalExtract

  function radiusTidalName(self)
    !!{
    Return the name of the tidal radius property.
    !!}
    implicit none
    type (varying_string                  )                :: radiusTidalName
    class(nodePropertyExtractorRadiusTidal), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusTidalName=var_str('satelliteRadiusTidal')
    return
  end function radiusTidalName

  function radiusTidalDescription(self)
    !!{
    Return a description of the tidal radius property.
    !!}
    implicit none
    type (varying_string                  )                :: radiusTidalDescription
    class(nodePropertyExtractorRadiusTidal), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusTidalDescription=var_str('Tidal radius in the halo [Mpc].')
    return
  end function radiusTidalDescription

  double precision function radiusTidalUnitsInSI(self)
    !!{
    Return the units of the tidal radius property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorRadiusTidal), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusTidalUnitsInSI=megaParsec
    return
  end function radiusTidalUnitsInSI


