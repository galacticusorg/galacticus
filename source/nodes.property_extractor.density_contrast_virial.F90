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
Implements a virial density contrast output analysis property extractor class.
!!}

  use :: Virial_Density_Contrast, only : virialDensityContrastClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorDensityContrastVirial">
   <description>A virial density contrast output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorDensityContrastVirial
     !!{
     A virial density contrast extractor output analysis class. Note that the density contrast is defined here at the time at
     which is was last isolated.
     !!}
     private
     class(virialDensityContrastClass), pointer :: virialDensityContrast_ => null()
   contains
     final     ::                densityContrastVirialDestructor
     procedure :: extract     => densityContrastVirialExtract
     procedure :: name        => densityContrastVirialName
     procedure :: description => densityContrastVirialDescription
     procedure :: unitsInSI   => densityContrastVirialUnitsInSI
  end type nodePropertyExtractorDensityContrastVirial

  interface nodePropertyExtractorDensityContrastVirial
     !!{
     Constructors for the \refClass{nodePropertyExtractorDensityContrastVirial} output analysis class.
     !!}
     module procedure densityContrastVirialConstructorParameters
     module procedure densityContrastVirialConstructorInternal
  end interface nodePropertyExtractorDensityContrastVirial

contains

  function densityContrastVirialConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorDensityContrastVirial} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorDensityContrastVirial)                :: self
    type (inputParameters                           ), intent(inout) :: parameters
    class(virialDensityContrastClass                ), pointer       :: virialDensityContrast_

    !![
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    !!]
    self=nodePropertyExtractorDensityContrastVirial(virialDensityContrast_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="virialDensityContrast_"/>
    !!]
    return
  end function densityContrastVirialConstructorParameters

  function densityContrastVirialConstructorInternal(virialDensityContrast_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorDensityContrastVirial} output analysis property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorDensityContrastVirial)                        :: self
    class(virialDensityContrastClass                ), intent(in   ), target :: virialDensityContrast_
    !![
    <constructorAssign variables="*virialDensityContrast_"/>
    !!]

    return
  end function densityContrastVirialConstructorInternal

  subroutine densityContrastVirialDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorDensityContrastVirial} output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorDensityContrastVirial), intent(inout) :: self

    !![
    <objectDestructor name="self%virialDensityContrast_"/>
    !!]
    return
  end subroutine densityContrastVirialDestructor

  double precision function densityContrastVirialExtract(self,node,instance) result(densityContrast)
    !!{
    Implement a densityContrastVirial output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodePropertyExtractorDensityContrastVirial), intent(inout), target   :: self
    type (treeNode                                  ), intent(inout), target   :: node
    type (multiCounter                              ), intent(inout), optional :: instance
    class(nodeComponentBasic                        ), pointer                 :: basic
    !$GLC attributes unused :: instance
    
    basic           => node%basic                                 (                               &
         &                                                        )
    densityContrast =  self%virialDensityContrast_%densityContrast(                               &
         &                                                         mass=basic%mass            (), &
         &                                                         time=basic%timeLastIsolated()  &
         &                                                        )
    return
  end function densityContrastVirialExtract

  function densityContrastVirialName(self) result(name)
    !!{
    Return the name of the densityContrastVirial property.
    !!}
    implicit none
    type (varying_string                            )                :: name
    class(nodePropertyExtractorDensityContrastVirial), intent(inout) :: self
    !$GLC attributes unused :: self

    name=var_str('densityContrastVirial')
    return
  end function densityContrastVirialName

  function densityContrastVirialDescription(self) result(description)
    !!{
    Return a description of the densityContrastVirial property.
    !!}
    implicit none
    type (varying_string                            )                :: description
    class(nodePropertyExtractorDensityContrastVirial), intent(inout) :: self
    !$GLC attributes unused :: self

    description=var_str('The virial density contrast of the halo.')
    return
  end function densityContrastVirialDescription

  double precision function densityContrastVirialUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the densityContrastVirial property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorDensityContrastVirial), intent(inout) :: self
    !$GLC attributes unused :: self

    unitsInSI=1.0d0
    return
  end function densityContrastVirialUnitsInSI
