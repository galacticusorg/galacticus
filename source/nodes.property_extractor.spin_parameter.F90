!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements a spin parameter output analysis property extractor class.
!!}

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
   contains
     procedure :: extract     => spinExtract
     procedure :: type        => spinType
     procedure :: name        => spinName
     procedure :: description => spinDescription
     procedure :: unitsInSI   => spinUnitsInSI
  end type nodePropertyExtractorSpin

  interface nodePropertyExtractorSpin
     !!{
     Constructors for the ``spin'' output property extractor class.
     !!}
     module procedure spinConstructorParameters
  end interface nodePropertyExtractorSpin

contains

  function spinConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``spin'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorSpin)                :: self
    type (inputParameters          ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    ! Build the object.
    self=nodePropertyExtractorSpin()
    return
  end function spinConstructorParameters

  double precision function spinExtract(self,node,instance)
    !!{
    Implement a spin output property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpin, treeNode
    implicit none
    class(nodePropertyExtractorSpin), intent(inout)           :: self
    type (treeNode                 ), intent(inout), target   :: node
    type (multiCounter             ), intent(inout), optional :: instance
    class(nodeComponentSpin        ), pointer                 :: spin
    !$GLC attributes unused :: self, instance

    spin        => node%spin()
    spinExtract =  spin%spin()
    return
  end function spinExtract

  integer function spinType(self)
    !!{
    Return the type of the spin parameter property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorSpin), intent(inout) :: self
    !$GLC attributes unused :: self

    spinType=outputAnalysisPropertyTypeLinear
    return
  end function spinType

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
