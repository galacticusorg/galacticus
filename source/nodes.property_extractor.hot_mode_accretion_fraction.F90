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
Implements a hot mode accretion fraction rate property extractor class.
!!}

  use :: Accretion_Halos, only : accretionHalo, accretionHaloClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorFractionAccretionHotMode">
   <description>A hot mode accretion fraction property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorFractionAccretionHotMode
     !!{
     A hot mode accretion fraction property extractor class.
     !!}
     private
     class(accretionHaloClass), pointer :: accretionHalo_ => null()
   contains
     final     ::                fractionAccretionHotModeDestructor
     procedure :: extract     => fractionAccretionHotModeExtract
     procedure :: name        => fractionAccretionHotModeName
     procedure :: description => fractionAccretionHotModeDescription
     procedure :: unitsInSI   => fractionAccretionHotModeUnitsInSI
  end type nodePropertyExtractorFractionAccretionHotMode

  interface nodePropertyExtractorFractionAccretionHotMode
     !!{
     Constructors for the {\normalfont \ttfamily fractionAccretionHotMode} output analysis class.
     !!}
     module procedure fractionAccretionHotModeConstructorParameters
     module procedure fractionAccretionHotModeConstructorInternal
  end interface nodePropertyExtractorFractionAccretionHotMode

contains

  function fractionAccretionHotModeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily fractionAccretionHotMode} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorFractionAccretionHotMode)                :: self
    type (inputParameters                              ), intent(inout) :: parameters
    class(accretionHaloClass                           ), pointer       :: accretionHalo_

    !![
    <objectBuilder class="accretionHalo" name="accretionHalo_" source="parameters"/>
    !!]
    self=nodePropertyExtractorFractionAccretionHotMode(accretionHalo_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="accretionHalo_"/>
    !!]
    return
  end function fractionAccretionHotModeConstructorParameters

  function fractionAccretionHotModeConstructorInternal(accretionHalo_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily fractionAccretionHotMode} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorFractionAccretionHotMode)                        :: self
    class(accretionHaloClass                           ), intent(in   ), target :: accretionHalo_
    !![
    <constructorAssign variables="*accretionHalo_"/>
    !!]

    return
  end function fractionAccretionHotModeConstructorInternal

  subroutine fractionAccretionHotModeDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily fractionAccretionHotMode} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorFractionAccretionHotMode), intent(inout) :: self

    !![
    <objectDestructor name="self%accretionHalo_"/>
    !!]
    return
  end subroutine fractionAccretionHotModeDestructor

  double precision function fractionAccretionHotModeExtract(self,node,instance)
    !!{
    Implement a {\normalfont \ttfamily fractionAccretionHotMode} property extractor.
    !!}
    use :: Accretion_Halos, only : accretionModeHot, accretionModeTotal
    implicit none
    class           (nodePropertyExtractorFractionAccretionHotMode), intent(inout), target   :: self
    type            (treeNode                                     ), intent(inout), target   :: node
    type            (multiCounter                                 ), intent(inout), optional :: instance
    double precision                                                                         :: accretionRateHot, accretionRateTotal
    !$GLC attributes unused :: instance

    accretionRateHot  =self%accretionHalo_%accretionRate(node,accretionModeHot  )
    accretionRateTotal=self%accretionHalo_%accretionRate(node,accretionModeTotal)
    if (accretionRateTotal /= 0.0d0) then
       fractionAccretionHotModeExtract=accretionRateHot/accretionRateTotal
    else
       fractionAccretionHotModeExtract=0.0d0
    end if
    return
  end function fractionAccretionHotModeExtract

  function fractionAccretionHotModeName(self)
    !!{
    Return the name of the {\normalfont \ttfamily fractionAccretionHotMode} property.
    !!}
    implicit none
    type (varying_string                               )                :: fractionAccretionHotModeName
    class(nodePropertyExtractorFractionAccretionHotMode), intent(inout) :: self
    !$GLC attributes unused :: self

    fractionAccretionHotModeName=var_str('haloAccretionHotModeFraction')
    return
  end function fractionAccretionHotModeName

  function fractionAccretionHotModeDescription(self)
    !!{
    Return a description of the {\normalfont \ttfamily fractionAccretionHotMode} property.
    !!}
    implicit none
    type (varying_string                               )                :: fractionAccretionHotModeDescription
    class(nodePropertyExtractorFractionAccretionHotMode), intent(inout) :: self
    !$GLC attributes unused :: self

    fractionAccretionHotModeDescription=var_str('Fraction of halo accretion rate occuring via the hot mode.')
    return
  end function fractionAccretionHotModeDescription

  double precision function fractionAccretionHotModeUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily fractionAccretionHotMode} property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorFractionAccretionHotMode), intent(inout) :: self
    !$GLC attributes unused :: self

    fractionAccretionHotModeUnitsInSI=0.0d0
    return
  end function fractionAccretionHotModeUnitsInSI


