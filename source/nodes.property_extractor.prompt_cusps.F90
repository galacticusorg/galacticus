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
  Implements a property extractor class for the parameters of prompt cusps in the \cite{delos_cusp-halo_2025} model.
  !!}
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorPromptCusps">
   <description>
    A property extractor class for the properties of the nuclear star cluster at the moment of the black hole formation.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorPromptCusps
     !!{
     A property extractor class for the velocity dispersion at a set of radii.
     !!}
     private
     integer :: promptCuspMassID, promptCuspAmplitudeID, &
          &     promptCuspNFWYID
   contains
     procedure :: elementCount       => promptCuspsElementCount
     procedure :: extract            => promptCuspsExtract
     procedure :: names              => promptCuspsNames
     procedure :: descriptions       => promptCuspsDescriptions
     procedure :: unitsInSI          => promptCuspsUnitsInSI
  end type nodePropertyExtractorPromptCusps

  interface nodePropertyExtractorPromptCusps
     !!{
     Constructors for the \refClass{nodePropertyExtractorPromptCusps} output analysis class.
     !!}
     module procedure promptCuspsConstructorParameters
     module procedure promptCuspsConstructorInternal
  end interface nodePropertyExtractorPromptCusps

contains

  function promptCuspsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorPromptCusps} property extractor class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorPromptCusps)                :: self
    type(inputParameters                 ), intent(inout) :: parameters

    self=nodePropertyExtractorPromptCusps()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function promptCuspsConstructorParameters

  function promptCuspsConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorPromptCusps} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorPromptCusps) :: self
    
    !![
    <addMetaProperty component="darkMatterProfile" name="promptCuspAmplitude" id="self%promptCuspAmplitudeID" isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspMass"      id="self%promptCuspMassID"      isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspNFWY"      id="self%promptCuspNFWYID"      isEvolvable="no" isCreator="no"/>
    !!]
    return
  end function promptCuspsConstructorInternal

  integer function promptCuspsElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily promptCusps} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorPromptCusps), intent(inout) :: self
    double precision                                  , intent(in   ) :: time
    !$GLC attributes unused :: time

    promptCuspsElementCount=3
    return
  end function promptCuspsElementCount

  function promptCuspsExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily promptCusps} property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    double precision                                  , dimension(:) , allocatable :: promptCuspsExtract
    class           (nodePropertyExtractorPromptCusps), intent(inout), target      :: self
    type            (treeNode                        ), intent(inout), target      :: node
    double precision                                  , intent(in   )              :: time
    type            (multiCounter                    ), intent(inout), optional    :: instance
    class           (nodeComponentDarkMatterProfile  )               , pointer     :: darkMatterProfile

    !$GLC attributes unused :: time, instance

    allocate(promptCuspsExtract(3))
    darkMatterProfile => node%darkMatterProfile()
    select type (darkMatterProfile)
    type is (nodeComponentDarkMatterProfile)
       ! Dark matter profile does not exist.
      promptCuspsExtract=[       & 
        &                 0.0d0, &
        &                 0.0d0, &
        &                 0.0d0  &
        &                ]
    class default
      promptCuspsExtract=[                                                                         &
       &                  darkMatterProfile%floatRank0MetaPropertyGet(self%promptCuspAmplitudeID), &
       &                  darkMatterProfile%floatRank0MetaPropertyGet(self%promptCuspMassID     ), &
       &                  darkMatterProfile%floatRank0MetaPropertyGet(self%promptCuspNFWYID     )  &
       &                 ]
    end select
    return
  end function promptCuspsExtract

  subroutine promptCuspsNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily promptCusps} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorPromptCusps), intent(inout)                             :: self
    double precision                                  , intent(in   )                             :: time
    type            (varying_string                  ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time
    
    allocate(names(3))
    names(1)=var_str('darkMatterProfilePromptCuspAmplitude')
    names(2)=var_str('darkMatterProfilePromptCuspMass'     )
    names(3)=var_str('darkMatterProfilePromptCuspNFWY'     )
    return
  end subroutine promptCuspsNames

  subroutine promptCuspsDescriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily promptCusps} property.
    !!}
    implicit none
    class           (nodePropertyExtractorPromptCusps), intent(inout)                             :: self
    double precision                                  , intent(in   )                             :: time
    type            (varying_string                  ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(3))
    descriptions(1)=var_str('The amplitude of the prompt cusp, A, in units of M☉/Mpc^1.5.')
    descriptions(2)=var_str('The mass of the prompt cusp, in units of M☉.'                )
    descriptions(3)=var_str('The y-parameter of the cusp-NFW profile.'                    )
    return
  end subroutine promptCuspsDescriptions

  function promptCuspsUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily PromptCusps} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    implicit none
    double precision                                  , allocatable  , dimension(:) :: promptCuspsUnitsInSI
    class           (nodePropertyExtractorPromptCusps), intent(inout)               :: self
    double precision                                  , intent(in   )               :: time
    !$GLC attributes unused :: time

    allocate(promptCuspsUnitsInSI(3))
    promptCuspsUnitsInSI=[                             &
         &                massSolar/megaParsec**1.5d0, &
         &                massSolar                  , &
         &                1.0d0                        &
         &               ]
    return
  end function promptCuspsUnitsInSI
  
