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
Implements an output analysis property extractor class that extracts the basic mass.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassBasic">
   <description>An output analysis property extractor class that extracts the basic mass.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMassBasic
     !!{
     A property extractor output analysis class that extracts the basic mass.
     !!}
     private
   contains
     procedure :: extract     => massBasicExtract
     procedure :: name        => massBasicName
     procedure :: description => massBasicDescription
     procedure :: unitsInSI   => massBasicUnitsInSI
  end type nodePropertyExtractorMassBasic

  interface nodePropertyExtractorMassBasic
     !!{
     Constructors for the ``massBasic'' output analysis class.
     !!}
     module procedure massBasicConstructorParameters
  end interface nodePropertyExtractorMassBasic

contains

  function massBasicConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``massBasic'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorMassBasic)                :: self
    type (inputParameters               ), intent(inout) :: parameters
    
    self=nodePropertyExtractorMassBasic()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massBasicConstructorParameters

  double precision function massBasicExtract(self,node,instance)
    !!{
    Implement a massBasic output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(nodePropertyExtractorMassBasic), intent(inout), target   :: self
    type (treeNode                      ), intent(inout), target   :: node
    type (multiCounter                  ), intent(inout), optional :: instance
    class(nodeComponentBasic            ), pointer                 :: basic
    !$GLC attributes unused :: self, instance

    basic            => node %basic()
    massBasicExtract =  basic%mass ()
    return
  end function massBasicExtract


  function massBasicName(self)
    !!{
    Return the name of the massBasic property.
    !!}
    implicit none
    type (varying_string                )                :: massBasicName
    class(nodePropertyExtractorMassBasic), intent(inout) :: self
    !$GLC attributes unused :: self

    massBasicName=var_str('massBasic')
    return
  end function massBasicName

  function massBasicDescription(self)
    !!{
    Return a description of the massBasic property.
    !!}
    implicit none
    type (varying_string                )                :: massBasicDescription
    class(nodePropertyExtractorMassBasic), intent(inout) :: self
    !$GLC attributes unused :: self

    massBasicDescription=var_str('The basic mass of the node.')
    return
  end function massBasicDescription

  double precision function massBasicUnitsInSI(self)
    !!{
    Return the units of the massBasic property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(nodePropertyExtractorMassBasic), intent(inout) :: self
    !$GLC attributes unused :: self

    massBasicUnitsInSI=massSolar
    return
  end function massBasicUnitsInSI
