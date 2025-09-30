!e Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
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
  Implements a property extractor class for the core mass of the \gls{fdm} soliton.
  !!}
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorSolitonMassCore">
   <description>
    A property extractor class for the core mass of the \gls{fdm} soliton.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorSolitonMassCore
     !!{
     A property extractor class for the core mass of the \gls{fdm} soliton.
     !!}
     private
     integer :: massCoreID
   contains
     procedure :: extract     => solitonMassCoreExtract
     procedure :: name        => solitonMassCoreName
     procedure :: description => solitonMassCoreDescription
     procedure :: unitsInSI   => solitonMassCoreUnitsInSI
  end type nodePropertyExtractorSolitonMassCore

  interface nodePropertyExtractorSolitonMassCore
     !!{
     Constructors for the \refClass{nodePropertyExtractorSolitonMassCore} property extractor class.
     !!}
     module procedure solitonMassCoreConstructorParameters
     module procedure solitonMassCoreConstructorInternal
  end interface nodePropertyExtractorSolitonMassCore

contains

  function solitonMassCoreConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorSolitonMassCore} property extractor class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorSolitonMassCore)                :: self
    type(inputParameters                     ), intent(inout) :: parameters

    self=nodePropertyExtractorSolitonMassCore()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function solitonMassCoreConstructorParameters

  function solitonMassCoreConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorSolitonMassCore} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorSolitonMassCore) :: self
    
    !![
    <addMetaProperty component="darkMatterProfile" name="massCore" id="self%massCoreID" isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function solitonMassCoreConstructorInternal

  double precision function solitonMassCoreExtract(self,node,instance) result(massCore)
    !!{
    Implement a {\normalfont \ttfamily solitonMassCore} property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(nodePropertyExtractorSolitonMassCore), intent(inout), target   :: self
    type (treeNode                            ), intent(inout), target   :: node
    type (multiCounter                        ), intent(inout), optional :: instance
    class(nodeComponentDarkMatterProfile      )               , pointer  :: darkMatterProfile
    !$GLC attributes unused :: instance

    darkMatterProfile => node%darkMatterProfile()
    select type (darkMatterProfile)
    type is (nodeComponentDarkMatterProfile)
       ! Dark matter profile does not exist.
       massCore=0.0d0
    class default
       massCore=darkMatterProfile%floatRank0MetaPropertyGet(self%massCoreID)
    end select
    return
  end function solitonMassCoreExtract

  function solitonMassCoreName(self) result(name)
    !!{
    Return the name of the {\normalfont \ttfamily solitonMassCore} property.
    !!}
    implicit none
    type (varying_string                      )                :: name
    class(nodePropertyExtractorSolitonMassCore), intent(inout) :: self
    !$GLC attributes unused :: self
    
    name=var_str('solitonMassCore')
    return
  end function solitonMassCoreName

  function solitonMassCoreDescription(self) result(description)
    !!{
    Return a description of the {\normalfont \ttfamily solitonMassCore} property.
    !!}
    implicit none
    type (varying_string                      )                :: description
    class(nodePropertyExtractorSolitonMassCore), intent(inout) :: self

    description=var_str('The solitonic core mass of the FDM halo, in units of M☉.')
    return
  end function solitonMassCoreDescription

 double precision function solitonMassCoreUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily SolitonMassCore} property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(nodePropertyExtractorSolitonMassCore), intent(inout) :: self

    unitsInSI=massSolar
    return
  end function solitonMassCoreUnitsInSI
  
