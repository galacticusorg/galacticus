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
Implements a stellar mass-weighted morphology output analysis property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassStellarMorphology">
   <description>A stellar mass-weighted morphology output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMassStellarMorphology
     !!{
     A stellar mass output analysis class.
     !!}
     private
   contains
     procedure :: extract     => massStellarMorphologyExtract
     procedure :: name        => massStellarMorphologyName
     procedure :: description => massStellarMorphologyDescription
     procedure :: unitsInSI   => massStellarMorphologyUnitsInSI
  end type nodePropertyExtractorMassStellarMorphology

  interface nodePropertyExtractorMassStellarMorphology
     !!{
     Constructors for the {\normalfont \ttfamily massStellarMorphology} output analysis class.
     !!}
     module procedure massStellarMorphologyConstructorParameters
  end interface nodePropertyExtractorMassStellarMorphology

contains

  function massStellarMorphologyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily massStellarMorphology} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorMassStellarMorphology)                :: self
    type (inputParameters                           ), intent(inout) :: parameters

   
    self=nodePropertyExtractorMassStellarMorphology()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massStellarMorphologyConstructorParameters

  double precision function massStellarMorphologyExtract(self,node,instance)
    !!{
    Implement a stellar mass-weighted morphology output analysis.
    !!}
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Galactic_Structure_Options, only : componentTypeDisk    , componentTypeSpheroid, massTypeStellar
    implicit none
    class           (nodePropertyExtractorMassStellarMorphology), intent(inout), target   :: self
    type            (treeNode                                  ), intent(inout), target   :: node
    type            (multiCounter                              ), intent(inout), optional :: instance
    class           (massDistributionClass                     )               , pointer  :: massDistributionDisk, massDistributionSpheroid
    double precision                                                                      :: massStellarDisk     , massStellarSpheroid
    !$GLC attributes unused :: self, instance

    massDistributionDisk     => node                    %massDistribution(massType=massTypeStellar,componentType=componentTypeDisk    )
    massDistributionSpheroid => node                    %massDistribution(massType=massTypeStellar,componentType=componentTypeSpheroid)
    massStellarDisk          =  massDistributionDisk    %massTotal       (                                                            )
    massStellarSpheroid      =  massDistributionSpheroid%massTotal       (                                                            )
    if (massStellarDisk+massStellarSpheroid > 0.0d0) then
       massStellarMorphologyExtract=+  massStellarSpheroid &
            &                       /(                     &
            &                         +massStellarDisk     &
            &                         +massStellarSpheroid &
            &                       )
    else
       massStellarMorphologyExtract=+0.0d0
    end if
    !![
    <objectDestructor name="massDistributionDisk"    />
    <objectDestructor name="massDistributionSpheroid"/>
    !!]
    return
  end function massStellarMorphologyExtract


  function massStellarMorphologyName(self)
    !!{
    Return the name of the massStellarMorphology property.
    !!}
    implicit none
    type (varying_string                            )                :: massStellarMorphologyName
    class(nodePropertyExtractorMassStellarMorphology), intent(inout) :: self
    !$GLC attributes unused :: self

    massStellarMorphologyName=var_str('morphologyMassStellarWeighted')
    return
  end function massStellarMorphologyName

  function massStellarMorphologyDescription(self)
    !!{
    Return a description of the massStellarMorphology property.
    !!}
    implicit none
    type (varying_string                            )                :: massStellarMorphologyDescription
    class(nodePropertyExtractorMassStellarMorphology), intent(inout) :: self
    !$GLC attributes unused :: self

    massStellarMorphologyDescription=var_str('The stellar mass-weighted mophology (bulge-to-total ratio) of the galaxy.')
    return
  end function massStellarMorphologyDescription

  double precision function massStellarMorphologyUnitsInSI(self)
    !!{
    Return the units of the massStellarMorphology property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorMassStellarMorphology), intent(inout) :: self
    !$GLC attributes unused :: self

    massStellarMorphologyUnitsInSI=0.0d0
    return
  end function massStellarMorphologyUnitsInSI
