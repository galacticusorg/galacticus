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
Implements a dark matter profile scale radius output analysis property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorDarkMatterProfileShapeParameter">
   <description>A  dark matter profile scale radius output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorDarkMatterProfileShapeParameter
     !!{
     A dark matter profile scale radius output property extractor class.
     !!}
     private
   contains
     procedure :: extract     => darkMatterProfileShapeParameterExtract
     procedure :: name        => darkMatterProfileShapeParameterName
     procedure :: description => darkMatterProfileShapeParameterDescription
     procedure :: unitsInSI   => darkMatterProfileShapeParameterUnitsInSI
  end type nodePropertyExtractorDarkMatterProfileShapeParameter

  interface nodePropertyExtractorDarkMatterProfileShapeParameter
     !!{
     Constructors for the {\normalfont \ttfamily darkMatterProfileShapeParameter} output analysis class.
     !!}
     module procedure darkMatterProfileShapeParameterConstructorParameters
  end interface nodePropertyExtractorDarkMatterProfileShapeParameter

contains

  function darkMatterProfileShapeParameterConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily darkMatterProfileShapeParameter} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorDarkMatterProfileShapeParameter)                :: self
    type (inputParameters                                  ), intent(inout) :: parameters
    
    self=nodePropertyExtractorDarkMatterProfileShapeParameter()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function darkMatterProfileShapeParameterConstructorParameters

  double precision function darkMatterProfileShapeParameterExtract(self,node,instance)
    !!{
    Implement a {\normalfont \ttfamily darkMatterProfileShapeParameter} output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, treeNode
    implicit none
    class(nodePropertyExtractorDarkMatterProfileShapeParameter), intent(inout), target   :: self
    type (treeNode                                            ), intent(inout), target   :: node
    type (multiCounter                                        ), intent(inout), optional :: instance
    class(nodeComponentDarkMatterProfile                      ), pointer                 :: darkMatterProfile
    !$GLC attributes unused :: self, instance

    darkMatterProfile                      => node             %darkMatterProfile()
    darkMatterProfileShapeParameterExtract =  darkMatterProfile%shape            ()
    return
  end function darkMatterProfileShapeParameterExtract


  function darkMatterProfileShapeParameterName(self)
    !!{
    Return the name of the {\normalfont \ttfamily darkMatterProfileShapeParameter} property.
    !!}
    implicit none
    type (varying_string                                      )                :: darkMatterProfileShapeParameterName
    class(nodePropertyExtractorDarkMatterProfileShapeParameter), intent(inout) :: self
    !$GLC attributes unused :: self

    darkMatterProfileShapeParameterName=var_str('darkMatterProfileShapeParameter')
    return
  end function darkMatterProfileShapeParameterName

  function darkMatterProfileShapeParameterDescription(self)
    !!{
    Return a description of the {\normalfont \ttfamily darkMatterProfileShapeParameter} property.
    !!}
    implicit none
    type (varying_string                                      )                :: darkMatterProfileShapeParameterDescription
    class(nodePropertyExtractorDarkMatterProfileShapeParameter), intent(inout) :: self
    !$GLC attributes unused :: self

    darkMatterProfileShapeParameterDescription=var_str('The shape parameter of the dark matter profile.')
    return
  end function darkMatterProfileShapeParameterDescription

  double precision function darkMatterProfileShapeParameterUnitsInSI(self)
    !!{
    Return the units of the darkMatterProfileShapeParameter property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorDarkMatterProfileShapeParameter), intent(inout) :: self
    !$GLC attributes unused :: self

    darkMatterProfileShapeParameterUnitsInSI=0.0d0
    return
  end function darkMatterProfileShapeParameterUnitsInSI
