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

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassBertschinger">
   <description>
     A node property extractor which extracts the Bertschinger mass of the halo.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMassBertschinger
     !!{
     A property extractor which extracts the Bertschinger mass of the halo.
     !!}
     private
     integer :: massBertschingerID
   contains
     procedure :: extract     => massBertschingerExtract
     procedure :: name        => massBertschingerName
     procedure :: description => massBertschingerDescription
     procedure :: unitsInSI   => massBertschingerUnitsInSI
  end type nodePropertyExtractorMassBertschinger

  interface nodePropertyExtractorMassBertschinger
     !!{
     Constructors for the \refClass{nodePropertyExtractorMassBertschinger} output extractor class.
     !!}
     module procedure massBertschingerConstructorParameters
     module procedure massBertschingerConstructorInternal
  end interface nodePropertyExtractorMassBertschinger

contains

  function massBertschingerConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorMassBertschinger} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorMassBertschinger)                :: self
    type(inputParameters                      ), intent(inout) :: parameters
    
    self=nodePropertyExtractorMassBertschinger()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massBertschingerConstructorParameters

  function massBertschingerConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorMassBertschinger} output extractor property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorMassBertschinger) :: self
    
    !![
    <addMetaProperty component="basic" name="massBertschinger" id="self%massBertschingerID" isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function massBertschingerConstructorInternal

  double precision function massBertschingerExtract(self,node,instance)
    !!{
    Implement an output extractor for the mass of gas cooled out of the \gls{cgm}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodePropertyExtractorMassBertschinger), intent(inout), target   :: self
    type (treeNode                             ), intent(inout), target   :: node
    type (multiCounter                         ), intent(inout), optional :: instance
    class(nodeComponentBasic                   )               , pointer  :: basic
    !$GLC attributes unused :: instance

    basic                   => node %basic                    (                       )
    massBertschingerExtract =  basic%floatRank0MetaPropertyGet(self%massBertschingerID)
    return
  end function massBertschingerExtract

  function massBertschingerName(self)
    !!{
    Return the names of the {\normalfont \ttfamily massBertschinger} property.
    !!}
    implicit none
    type (varying_string                        )               :: massBertschingerName
    class(nodePropertyExtractorMassBertschinger), intent(inout) :: self
    !$GLC attributes unused :: self

    massBertschingerName=var_str('massBertschinger')
    return
  end function massBertschingerName

  function massBertschingerDescription(self)
    !!{
    Return the description of the {\normalfont \ttfamily massBertschinger} property.
    !!}
    implicit none
    type (varying_string                       )                :: massBertschingerDescription
    class(nodePropertyExtractorMassBertschinger), intent(inout) :: self
    !$GLC attributes unused :: self

    massBertschingerDescription=var_str('Bertschinger mass of the halo [M☉].')
    return
  end function massBertschingerDescription

  double precision function massBertschingerUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily massBertschinger} property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(nodePropertyExtractorMassBertschinger), intent(inout) :: self
    !$GLC attributes unused :: self

    massBertschingerUnitsInSI=massSolar
    return
  end function massBertschingerUnitsInSI
