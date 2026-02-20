!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  <nodePropertyExtractor name="nodePropertyExtractorSoliton">
   <description>
    A property extractor class for the properties of the \gls{fdm} soliton.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorSoliton
     !!{
     A property extractor class for the properties of the \gls{fdm} soliton.
     !!}
     private
     integer :: massCoreNormalID, massCoreID, densityCoreID, radiusCoreID, radiusSolitonID, zetaID
   contains
     procedure :: elementCount => solitonElementCount
     procedure :: extract      => solitonExtract
     procedure :: names        => solitonNames
     procedure :: descriptions => solitonDescriptions
     procedure :: unitsInSI    => solitonUnitsInSI
  end type nodePropertyExtractorSoliton

  interface nodePropertyExtractorSoliton
     !!{
     Constructors for the \refClass{nodePropertyExtractorSoliton} property extractor class.
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
    <addMetaProperty component="darkMatterProfile" name="solitonMassCoreNormal" id="self%massCoreNormalID" isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonMassCore"       id="self%massCoreID"       isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonDensityCore"    id="self%densityCoreID"    isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonRadiusCore"     id="self%radiusCoreID"     isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonRadiusSoliton"  id="self%radiusSolitonID"  isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonZeta"           id="self%zetaID"           isEvolvable="no"  isCreator="no"/>
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
    !$GLC attributes unused :: self, time

    solitonElementCount=6
    return
  end function solitonElementCount

  function solitonExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily soliton} property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    double precision                     , dimension(: ), allocatable :: solitonExtract
    class(nodePropertyExtractorSoliton  ), intent(inout), target      :: self
    type (treeNode                      ), intent(inout), target      :: node
    double precision                     , intent(in   )              :: time
    type (multiCounter                  ), intent(inout), optional    :: instance
    class(nodeComponentDarkMatterProfile)               , pointer     :: darkMatterProfile
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
      solitonExtract=[                                                                    &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%massCoreNormalID), &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%massCoreID      ), &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%densityCoreID   ), &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%radiusCoreID    ), &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%radiusSolitonID ), &
       &              darkMatterProfile%floatRank0MetaPropertyGet(self%zetaID          )  &
       &             ]    
    end select
    return
  end function solitonExtract

  subroutine solitonNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily soliton} property.
    !!}
    implicit none
    class(nodePropertyExtractorSoliton), intent(inout)                             :: self
    double precision                   , intent(in   )                             :: time
    type            (varying_string   ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(6))
    names(1)=var_str('solitonMassCoreNormal')
    names(2)=var_str('solitonMassCore'      )
    names(3)=var_str('solitonDensityCore'   )
    names(4)=var_str('solitonRadiusCore'    )
    names(5)=var_str('solitonRadiusSoliton' )
    names(6)=var_str('solitonZetazOverZeta0')
    return
  end subroutine solitonNames

  subroutine solitonDescriptions(self,time,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily soliton} property.
    !!}
    implicit none
    class(nodePropertyExtractorSoliton), intent(inout)                             :: self
    double precision                   , intent(in   )                             :: time
    type            (varying_string   ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(6))
    descriptions(1)=var_str('The solitonic core mass of the FDM halo (without scatter), in units of M☉.')
    descriptions(2)=var_str('The solitonic core mass of the FDM halo, in units of M☉.'                  )
    descriptions(3)=var_str('The soliton central density of the FDM halo, in units of [M☉/Mpc³].'       )
    descriptions(4)=var_str('The soliton core radius of the FDM halo, in units of Mpc.'                 )
    descriptions(5)=var_str('The soliton transition radius of the FDM halo, in units of Mpc.'           )
    descriptions(6)=var_str('The ratio of density contrast at redshift z and 0, dimensionless.'         )
    return
  end subroutine solitonDescriptions

  function solitonUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily Soliton} property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    implicit none
    double precision                   , allocatable  , dimension(:) :: solitonUnitsInSI
    class(nodePropertyExtractorSoliton), intent(inout)               :: self
    double precision                                 , intent(in   )              :: time
    !$GLC attributes unused :: self, time

    allocate(solitonUnitsInSI(6))
    solitonUnitsInSI=[                         &
         &            massSolar              , &
         &            massSolar              , &
         &            massSolar/megaParsec**3, &
         &            megaParsec             , &
         &            megaParsec             , &
         &            1.0d0                    &
         &           ]
    return
  end function solitonUnitsInSI
  
