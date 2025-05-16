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
  <nodePropertyExtractor name="nodePropertyExtractorMassCooled">
   <description>
     A node property extractor which extracts the mass of gas cooled out of the \gls{cgm}. If the parameter {\normalfont \ttfamily
     [resetAfterExtract]}$=${\normalfont \ttfamily true} then the cooled mass is reset to zero after extraction.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMassCooled
     !!{
     A property extractor which extracts the mass of gas cooled out of the \gls{cgm}.
     !!}
     private
     integer :: massCooledID
     logical :: resetAfterExtract
   contains
     procedure :: extract     => massCooledExtract
     procedure :: name        => massCooledName
     procedure :: description => massCooledDescription
     procedure :: unitsInSI   => massCooledUnitsInSI
  end type nodePropertyExtractorMassCooled

  interface nodePropertyExtractorMassCooled
     !!{
     Constructors for the \refClass{nodePropertyExtractorMassCooled} output extractor class.
     !!}
     module procedure massCooledConstructorParameters
     module procedure massCooledConstructorInternal
  end interface nodePropertyExtractorMassCooled

contains

  function massCooledConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorMassCooled} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorMassCooled)                :: self
    type   (inputParameters                ), intent(inout) :: parameters
    logical                                                 :: resetAfterExtract
    
    !![
    <inputParameter>
      <name>resetAfterExtract</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, the mass of gas cooled is reset to zero after being extracted.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodePropertyExtractorMassCooled(resetAfterExtract)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massCooledConstructorParameters

  function massCooledConstructorInternal(resetAfterExtract) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorMassCooled} output extractor property extractor class.
    !!}
    implicit none
    type   (nodePropertyExtractorMassCooled)                :: self
    logical                                 , intent(in   ) :: resetAfterExtract
    !![
    <constructorAssign variables="resetAfterExtract"/>
    !!]
    
    !![
    <addMetaProperty component="hotHalo" name="massCooled" id="self%massCooledID" isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function massCooledConstructorInternal

  double precision function massCooledExtract(self,node,instance)
    !!{
    Implement an output extractor for the mass of gas cooled out of the \gls{cgm}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo
    implicit none
    class(nodePropertyExtractorMassCooled), intent(inout), target   :: self
    type (treeNode                       ), intent(inout), target   :: node
    type (multiCounter                   ), intent(inout), optional :: instance
    class(nodeComponentHotHalo           )               , pointer  :: hotHalo
    !$GLC attributes unused :: instance

    hotHalo => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! Spheroid does not yet exist.
       massCooledExtract=0.0d0
    class default
       massCooledExtract=hotHalo%floatRank0MetaPropertyGet(self%massCooledID)
       if (self%resetAfterExtract) call hotHalo%floatRank0MetaPropertySet(self%massCooledID,0.0d0)
    end select
    return
  end function massCooledExtract

  function massCooledName(self)
    !!{
    Return the names of the {\normalfont \ttfamily massCooled} property.
    !!}
    implicit none
    type (varying_string                 )                :: massCooledName
    class(nodePropertyExtractorMassCooled), intent(inout) :: self
    !$GLC attributes unused :: self

    massCooledName=var_str('massCooledCGM')
    return
  end function massCooledName

  function massCooledDescription(self)
    !!{
    Return the description of the {\normalfont \ttfamily massCooled} property.
    !!}
    implicit none
    type (varying_string                 )                :: massCooledDescription
    class(nodePropertyExtractorMassCooled), intent(inout) :: self
    !$GLC attributes unused :: self

    massCooledDescription=var_str('Mass of gas cooled out of the CGM [M☉].')
    return
  end function massCooledDescription

  double precision function massCooledUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily massCooled} property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(nodePropertyExtractorMassCooled), intent(inout) :: self
    !$GLC attributes unused :: self

    massCooledUnitsInSI=massSolar
    return
  end function massCooledUnitsInSI
