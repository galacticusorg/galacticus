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
Implements an output analysis property extractor class that extracts the bound mass.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassBound">
   <description>An output analysis property extractor class that extracts the bound mass.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMassBound
     !!{
     A property extractor output analysis class that extracts the bound mass.
     !!}
     private
   contains
     procedure :: extract     => massBoundExtract
     procedure :: name        => massBoundName
     procedure :: description => massBoundDescription
     procedure :: unitsInSI   => massBoundUnitsInSI
  end type nodePropertyExtractorMassBound

  interface nodePropertyExtractorMassBound
     !!{
     Constructors for the ``massBound'' output analysis class.
     !!}
     module procedure massBoundConstructorParameters
  end interface nodePropertyExtractorMassBound

contains

  function massBoundConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``massBound'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorMassBound)                :: self
    type (inputParameters               ), intent(inout) :: parameters
    
    self=nodePropertyExtractorMassBound()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massBoundConstructorParameters

  double precision function massBoundExtract(self,node,instance)
    !!{
    Implement a massBound output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class(nodePropertyExtractorMassBound), intent(inout), target   :: self
    type (treeNode                      ), intent(inout), target   :: node
    type (multiCounter                  ), intent(inout), optional :: instance
    class(nodeComponentSatellite        ), pointer                 :: satellite
    !$GLC attributes unused :: self, instance

    satellite        => node     %satellite()
    massBoundExtract =  satellite%boundMass()
    return
  end function massBoundExtract


  function massBoundName(self)
    !!{
    Return the name of the massBound property.
    !!}
    implicit none
    type (varying_string                )                :: massBoundName
    class(nodePropertyExtractorMassBound), intent(inout) :: self
    !$GLC attributes unused :: self

    massBoundName=var_str('massBound')
    return
  end function massBoundName

  function massBoundDescription(self)
    !!{
    Return a description of the massBound property.
    !!}
    implicit none
    type (varying_string                )                :: massBoundDescription
    class(nodePropertyExtractorMassBound), intent(inout) :: self
    !$GLC attributes unused :: self

    massBoundDescription=var_str('The bound mass of the node.')
    return
  end function massBoundDescription

  double precision function massBoundUnitsInSI(self)
    !!{
    Return the units of the massBound property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(nodePropertyExtractorMassBound), intent(inout) :: self
    !$GLC attributes unused :: self

    massBoundUnitsInSI=massSolar
    return
  end function massBoundUnitsInSI
