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
Implements a cosmic time output analysis property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorTime">
   <description>A cosmic time output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorTime
     !!{
     A cosmic time property extractor output analysis class.
     !!}
     private
   contains
     procedure :: extract     => timeExtract
     procedure :: name        => timeName
     procedure :: description => timeDescription
     procedure :: unitsInSI   => timeUnitsInSI
  end type nodePropertyExtractorTime

  interface nodePropertyExtractorTime
     !!{
     Constructors for the ``time'' output analysis class.
     !!}
     module procedure timeConstructorParameters
  end interface nodePropertyExtractorTime

contains

  function timeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``time'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorTime)                :: self
    type (inputParameters          ), intent(inout) :: parameters
    
    self=nodePropertyExtractorTime()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function timeConstructorParameters

  double precision function timeExtract(self,node,instance)
    !!{
    Implement a time output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(nodePropertyExtractorTime), intent(inout), target   :: self
    type (treeNode                 ), intent(inout), target   :: node
    type (multiCounter             ), intent(inout), optional :: instance
    class(nodeComponentBasic       ), pointer                 :: basic
    !$GLC attributes unused :: self, instance

    basic       => node %basic()
    timeExtract =  basic%time ()
    return
  end function timeExtract


  function timeName(self)
    !!{
    Return the name of the time property.
    !!}
    implicit none
    type (varying_string           )                :: timeName
    class(nodePropertyExtractorTime), intent(inout) :: self
    !$GLC attributes unused :: self

    timeName=var_str('time')
    return
  end function timeName

  function timeDescription(self)
    !!{
    Return a description of the time property.
    !!}
    implicit none
    type (varying_string           )                :: timeDescription
    class(nodePropertyExtractorTime), intent(inout) :: self
    !$GLC attributes unused :: self

    timeDescription=var_str('The cosmic time at which the node exists.')
    return
  end function timeDescription

  double precision function timeUnitsInSI(self)
    !!{
    Return the units of the time property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    class(nodePropertyExtractorTime), intent(inout) :: self
    !$GLC attributes unused :: self

    timeUnitsInSI=gigaYear
    return
  end function timeUnitsInSI
