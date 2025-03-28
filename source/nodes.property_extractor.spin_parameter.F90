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
Implements a spin parameter output analysis property extractor class.
!!}
  
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorSpin">
   <description>A spin parameter output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorSpin
     !!{
     A spin parameter property extractor output analysis class.
     !!}
     private
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::                spinDestructor
     procedure :: extract     => spinExtract
     procedure :: name        => spinName
     procedure :: description => spinDescription
     procedure :: unitsInSI   => spinUnitsInSI
  end type nodePropertyExtractorSpin

  interface nodePropertyExtractorSpin
     !!{
     Constructors for the {\normalfont \ttfamily spin} output property extractor class.
     !!}
     module procedure spinConstructorParameters
     module procedure spinConstructorInternal
  end interface nodePropertyExtractorSpin

contains

  function spinConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily spin} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorSpin)                :: self
    type (inputParameters          ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass ), pointer       :: darkMatterHaloScale_
    !$GLC attributes unused :: parameters

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodePropertyExtractorSpin(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function spinConstructorParameters

  function spinConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily spin} output analysis property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorSpin)                        :: self
    class(darkMatterHaloScaleClass ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    return
  end function spinConstructorInternal

  subroutine spinDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily spin} output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorSpin), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine spinDestructor

  double precision function spinExtract(self,node,instance)
    !!{
    Implement a spin output property extractor.
    !!}
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum_Scale
    use :: Galacticus_Nodes      , only : nodeComponentSpin
    implicit none
    class(nodePropertyExtractorSpin), intent(inout), target   :: self
    type (treeNode                 ), intent(inout), target   :: node
    type (multiCounter             ), intent(inout), optional :: instance
    class(nodeComponentSpin        ), pointer                 :: spin
    !$GLC attributes unused :: self, instance

    spin        =>  node               %spin()
    spinExtract =  +spin%angularMomentum    ()                                              &
         &         /Dark_Matter_Halo_Angular_Momentum_Scale(node,self%darkmatterHaloScale_)
    return
  end function spinExtract


  function spinName(self)
    !!{
    Return the name of the spin property.
    !!}
    implicit none
    type (varying_string           )                :: spinName
    class(nodePropertyExtractorSpin), intent(inout) :: self
    !$GLC attributes unused :: self

    spinName=var_str('spinParameter')
    return
  end function spinName

  function spinDescription(self)
    !!{
    Return a description of the spin property.
    !!}
    implicit none
    type (varying_string           )                :: spinDescription
    class(nodePropertyExtractorSpin), intent(inout) :: self
    !$GLC attributes unused :: self

    spinDescription=var_str('The spin parameter of the dark matter halos.')
    return
  end function spinDescription

  double precision function spinUnitsInSI(self)
    !!{
    Return the units of the spin property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorSpin), intent(inout) :: self
    !$GLC attributes unused :: self

    spinUnitsInSI=0.0d0
    return
  end function spinUnitsInSI
