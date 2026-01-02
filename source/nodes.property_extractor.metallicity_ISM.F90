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
Implements an ISM metallicity output analysis property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMetallicityISM">
   <description>An ISM metallicity output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMetallicityISM
     !!{
     An ISM metallicity output analysis property extractor class.
     !!}
     private
     integer          :: indexElement
     character(len=3) :: element
   contains
     procedure :: extract     => metallicityISMExtract
     procedure :: name        => metallicityISMName
     procedure :: description => metallicityISMDescription
     procedure :: unitsInSI   => metallicityISMUnitsInSI
  end type nodePropertyExtractorMetallicityISM

  interface nodePropertyExtractorMetallicityISM
     !!{
     Constructors for the \refClass{nodePropertyExtractorMetallicityISM} output analysis class.
     !!}
     module procedure metallicityISMConstructorParameters
     module procedure metallicityISMConstructorInternal
  end interface nodePropertyExtractorMetallicityISM

contains

  function metallicityISMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorMetallicityISM} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Abundances_Structure, only : Abundances_Index_From_Name
    use :: Input_Parameters    , only : inputParameter            , inputParameters
    use :: Error               , only : Error_Report
    implicit none
    type     (nodePropertyExtractorMetallicityISM)                :: self
    type     (inputParameters                    ), intent(inout) :: parameters
    integer                                                       :: indexElement
    character(len=3                              )                :: element

    !![
    <inputParameter>
      <name>element</name>
      <source>parameters</source>
      <description>The atomic symbol for the element to use to define metallicity.</description>
    </inputParameter>
    !!]
    indexElement=Abundances_Index_From_Name(element)
    if (indexElement < 0) call Error_Report('element "'//trim(element)//'" is not being tracked'//{introspection:location})
    self=nodePropertyExtractorMetallicityISM(indexElement)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function metallicityISMConstructorParameters

  function metallicityISMConstructorInternal(indexElement) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorMetallicityISM} output analysis property extractor class.
    !!}
    use :: Abundances_Structure, only : Abundances_Names
    implicit none
    type(nodePropertyExtractorMetallicityISM)                :: self
    integer                                  , intent(in   ) :: indexElement
    !![
    <constructorAssign variables="indexElement"/>
    !!]

    self%element=Abundances_Names(indexElement)
    return
  end function metallicityISMConstructorInternal

  double precision function metallicityISMExtract(self,node,instance)
    !!{
    Extracts the metallicity (defined as the mass ratio of a specified element to hydrogen) in the ISM.
    !!}
    use :: Abundances_Structure, only : abundances
    use :: Galacticus_Nodes    , only : nodeComponentDisk, nodeComponentSpheroid, treeNode
    implicit none
    class           (nodePropertyExtractorMetallicityISM), intent(inout), target       :: self
    type            (treeNode                           ), intent(inout), target       :: node
    type            (multiCounter                       ), intent(inout), optional     :: instance
    class           (nodeComponentDisk                  ), pointer                     :: disk
    class           (nodeComponentSpheroid              ), pointer                     :: spheroid
    double precision                                     , allocatable  , dimension(:) :: massElementsDisk, massElementsSpheroid
    type            (abundances                         )                              :: abundancesDisk  , abundancesSpheroid
    integer                                                                            :: countElements
    !$GLC attributes unused :: self, instance

    disk               => node    %disk         ()
    spheroid           => node    %spheroid     ()
    abundancesDisk     =  disk    %abundancesGas()
    abundancesSpheroid =  spheroid%abundancesGas()
    if (disk%massGas() > 0.0d0 .or. spheroid%massGas() > 0.0d0) then
       countElements=abundancesDisk%serializeCount()
       allocate(massElementsDisk    (countElements))
       allocate(massElementsSpheroid(countElements))
       call abundancesDisk    %serialize(massElementsDisk    )
       call abundancesSpheroid%serialize(massElementsSpheroid)
       metallicityISMExtract=+(                                                              &
            &                  +massElementsDisk    (self%indexElement)                      &
            &                  +massElementsSpheroid(self%indexElement)                      &
            &                 )                                                              &
            &                /(                                                              &
            &                  +disk    %massGas()*abundancesDisk    %hydrogenMassFraction() &
            &                  +spheroid%massGas()*abundancesSpheroid%hydrogenMassFraction() &
            &                 )
    else
       metallicityISMExtract=0.0d0
    end if
    return
  end function metallicityISMExtract


  function metallicityISMName(self)
    !!{
    Return the name of the metallicityISM property.
    !!}
    use :: Abundances_Structure, only : Abundances_Names
    implicit none
    type (varying_string                     )                :: metallicityISMName
    class(nodePropertyExtractorMetallicityISM), intent(inout) :: self
    !$GLC attributes unused :: self

    metallicityISMName=var_str('metallicityISM')//Abundances_Names(self%indexElement)
    return
  end function metallicityISMName

  function metallicityISMDescription(self)
    !!{
    Return a description of the metallicityISM property.
    !!}
    use :: Abundances_Structure, only : Abundances_Names
    implicit none
    type (varying_string                     )                :: metallicityISMDescription
    class(nodePropertyExtractorMetallicityISM), intent(inout) :: self
    !$GLC attributes unused :: self

    metallicityISMDescription=var_str('The metallicity ')//Abundances_Names(self%indexElement)//'/H (by mass) of the interstellar medium.'
    return
  end function metallicityISMDescription

  double precision function metallicityISMUnitsInSI(self)
    !!{
    Return the units of the metallicityISM property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorMetallicityISM), intent(inout) :: self
    !$GLC attributes unused :: self

    metallicityISMUnitsInSI=0.0d0
    return
  end function metallicityISMUnitsInSI
