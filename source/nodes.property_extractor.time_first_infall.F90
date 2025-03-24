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
  <nodePropertyExtractor name="nodePropertyExtractorTimeFirstInfall">
   <description>A time of first infall property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorTimeFirstInfall
     !!{
     A time of first infall property extractor output analysis class.
     !!}
     private
     integer :: timeFirstInfallID
   contains
     procedure :: extract     => timeFirstInfallExtract
     procedure :: name        => timeFirstInfallName
     procedure :: description => timeFirstInfallDescription
     procedure :: unitsInSI   => timeFirstInfallUnitsInSI
  end type nodePropertyExtractorTimeFirstInfall

  interface nodePropertyExtractorTimeFirstInfall
     !!{
     Constructors for the {\normalfont \ttfamily timeFirstInfall} output analysis class.
     !!}
     module procedure timeFirstInfallConstructorParameters
     module procedure timeFirstInfallConstructorInternal
  end interface nodePropertyExtractorTimeFirstInfall

contains

  function timeFirstInfallConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily timeFirstInfall} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorTimeFirstInfall)                :: self
    type(inputParameters                     ), intent(inout) :: parameters
    
    self=nodePropertyExtractorTimeFirstInfall()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function timeFirstInfallConstructorParameters

  function timeFirstInfallConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily timeFirstInfall} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorTimeFirstInfall) :: self

    !![
    <addMetaProperty component="basic" name="timeFirstInfall" id="self%timeFirstInfallID" isEvolvable="no" isCreator="no"/>
    !!]
    return
  end function timeFirstInfallConstructorInternal

  double precision function timeFirstInfallExtract(self,node,instance)
    !!{
    Implement a time of first infall property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(nodePropertyExtractorTimeFirstInfall), intent(inout), target   :: self
    type (treeNode                            ), intent(inout), target   :: node
    type (multiCounter                        ), intent(inout), optional :: instance
    class(nodeComponentBasic                  ), pointer                 :: basic
    !$GLC attributes unused :: self, instance

    basic       => node %basic                    (                      )
    timeFirstInfallExtract =  basic%floatRank0MetaPropertyGet(self%timeFirstInfallID)
    return
  end function timeFirstInfallExtract

  function timeFirstInfallName(self)
    !!{
    Return the name of the time of first infall property.
    !!}
    implicit none
    type (varying_string                      )                :: timeFirstInfallName
    class(nodePropertyExtractorTimeFirstInfall), intent(inout) :: self
    !$GLC attributes unused :: self

    timeFirstInfallName=var_str('timeFirstInfall')
    return
  end function timeFirstInfallName

  function timeFirstInfallDescription(self)
    !!{
    Return a description of the time of first infall property.
    !!}
    implicit none
    type (varying_string                      )                :: timeFirstInfallDescription
    class(nodePropertyExtractorTimeFirstInfall), intent(inout) :: self
    !$GLC attributes unused :: self

    timeFirstInfallDescription=var_str('The cosmic time at which the node first infell into another node.')
    return
  end function timeFirstInfallDescription

  double precision function timeFirstInfallUnitsInSI(self)
    !!{
    Return the units of the time of first infall property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    class(nodePropertyExtractorTimeFirstInfall), intent(inout) :: self
    !$GLC attributes unused :: self

    timeFirstInfallUnitsInSI=gigaYear
    return
  end function timeFirstInfallUnitsInSI
