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
  Implements a property extractor class for the properties of the \gls{fdm} soliton in the \cite{chan_diversity_2022} model.
  !!}
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorSoliton">
   <description>
    A property extractor class for the properties of the \gls{fdm} soliton in the \cite{chan_diversity_2022} model.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorSoliton
     !!{
     A property extractor class for the velocity dispersion at a set of radii.
     !!}
     private
     integer   :: massCoreID, massHaloID, massParticleID, zeta0ID, zetazID, expansionFactorID
   contains
     procedure :: elementCount       => solitonElementCount
     procedure :: extract            => solitonExtract
     procedure :: names              => solitonNames
     procedure :: descriptions       => solitonDescriptions
     procedure :: unitsInSI          => solitonUnitsInSI
  end type nodePropertyExtractorSoliton

  interface nodePropertyExtractorSoliton
     !!{
     Constructors for the \refClass{nodePropertyExtractorSoliton} output analysis class.
     !!}
     module procedure solitonConstructorParameters
     module procedure solitonConstructorInternal
  end interface nodePropertyExtractorSoliton

contains

  function solitonConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorSoliton} property extractor class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorSoliton)                :: self
    type(inputParameters             ), intent(inout) :: parameters

    self=nodePropertyExtractorSoliton()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function solitonConstructorParameters

  function solitonConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorSoliton} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorSoliton) :: self
    
    !![
    <addMetaProperty component="darkMatterProfile" name="massCore"        id="self%massCoreID"        isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="zeta0"           id="self%zeta0ID"           isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="zetaz"           id="self%zetazID"           isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="expansionFactor" id="self%expansionFactorID" isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="massParticle"    id="self%massParticleID"    isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="massHalo"        id="self%massHaloID"        isEvolvable="no"  isCreator="no"/>
    !!]
    return
  end function solitonConstructorInternal

  integer function solitonElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily soliton} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorSoliton), intent(inout) :: self
    double precision                              , intent(in   ) :: time
    !$GLC attributes unused :: time

    solitonElementCount=6
    return
  end function solitonElementCount

  function solitonExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily soliton} property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    double precision                                , dimension(:) , allocatable :: solitonExtract
    class           (nodePropertyExtractorSoliton  ), intent(inout), target      :: self
    type            (treeNode                      ), intent(inout), target      :: node
    double precision                                , intent(in   )              :: time
    type            (multiCounter                  ), intent(inout), optional    :: instance
    class           (nodeComponentDarkMatterProfile)               , pointer     :: darkMatterProfile

    !$GLC attributes unused :: time, instance

    allocate(solitonExtract(6))
    darkMatterProfile => node%darkMatterProfile()
    select type (darkMatterProfile)
    type is (nodeComponentDarkMatterProfile)
       ! Dark matter profile does not exist.
      solitonExtract=[        & 
        &             0.0d0,  &
        &             0.0d0,  &
        &             0.0d0,  &
        &             0.0d0,  &
        &             0.0d0,  &
        &             0.0d0   &
        &            ]
    class default
      solitonExtract=[                                                                     &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%massCoreID)       , &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%zeta0ID)          , &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%zetazID)          , &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%expansionFactorID), &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%massParticleID)   , &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%massHaloID)         &
       &             ]    
    end select
    return
  end function solitonExtract

  subroutine solitonNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily soliton} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorSoliton), intent(inout)                             :: self
    double precision                              , intent(in   )                             :: time
    type            (varying_string              ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time
    
    allocate(names(6))
    names(1)=var_str('solitonMassCore')
    names(2)=var_str('solitonZeta0')
    names(3)=var_str('solitonZetaZ')
    names(4)=var_str('solitonExpansionFactor')
    names(5)=var_str('solitonMassParticle')
    names(6)=var_str('solitonMassHalo')
    return
  end subroutine solitonNames

  subroutine solitonDescriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily soliton} property.
    !!}
    implicit none
    class           (nodePropertyExtractorSoliton), intent(inout)                             :: self
    double precision                              , intent(in   )                             :: time
    type            (varying_string              ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(6))
    descriptions(1) = var_str('The solitonic core mass of the FDM halo, in units of M☉.')
    descriptions(2) = var_str('The normalization constant ζ₀ of the soliton profile, dimensionless.')
    descriptions(3) = var_str('The effective ζ(z) factor controlling redshift evolution of the soliton–halo relation, dimensionless.')
    descriptions(4) = var_str('The cosmological expansion factor a = 1/(1+z), dimensionless.')
    descriptions(5) = var_str('The mass of an individual FDM particle, in units of eV.')
    descriptions(6) = var_str('The virial mass of the host dark matter halo, in units of M☉.')
    return
  end subroutine solitonDescriptions

  function solitonUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily Soliton} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    use :: Numerical_Constants_Units       , only : electronVolt
    implicit none
    double precision                              , allocatable  , dimension(:) :: solitonUnitsInSI
    class           (nodePropertyExtractorSoliton), intent(inout)               :: self
    double precision                              , intent(in   )               :: time
    !$GLC attributes unused :: time

    allocate(solitonUnitsInSI(6))
    solitonUnitsInSI=[              &
         &            massSolar   , &
         &            1.0d0       , &
         &            1.0d0       , &
         &            1.0d0       , &
         &            electronVolt, &
         &            massSolar     &
         &           ]
    return
  end function solitonUnitsInSI
  
