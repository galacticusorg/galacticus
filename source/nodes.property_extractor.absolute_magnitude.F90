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
Implements a node property extractor class for absolute magnitudes.
!!}

  use :: Galactic_Structure_Options, only : enumerationComponentTypeType

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMagnitudesAbsolute">
   <description>
    A node property extractor which extracts stellar absolute magnitudes in all available bands.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorMagnitudesAbsolute
     !!{
     A node property extractor which extracts stellar absolute magnitudes in all available bands.
     !!}
     private
     type(enumerationComponentTypeType) :: component
   contains
     procedure :: elementCount => magnitudesAbsoluteElementCount
     procedure :: extract      => magnitudesAbsoluteExtract
     procedure :: names        => magnitudesAbsoluteNames
     procedure :: descriptions => magnitudesAbsoluteDescriptions
     procedure :: unitsInSI    => magnitudesAbsoluteUnitsInSI
  end type nodePropertyExtractorMagnitudesAbsolute

  interface nodePropertyExtractorMagnitudesAbsolute
     !!{
     Constructors for the \refClass{nodePropertyExtractorMagnitudesAbsolute} output analysis class.
     !!}
     module procedure magnitudesAbsoluteConstructorParameters
     module procedure magnitudesAbsoluteConstructorInternal
  end interface nodePropertyExtractorMagnitudesAbsolute

contains

  function magnitudesAbsoluteConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorMagnitudesAbsolute} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode
    implicit none
    type(nodePropertyExtractorMagnitudesAbsolute)                :: self
    type(inputParameters                        ), intent(inout) :: parameters
    type(varying_string                         )                :: component
    
    !![
    <inputParameter>
      <name>component</name>
      <source>parameters</source>
      <description>The component from which to extract star formation rate.</description>
    </inputParameter>
    !!]
    self=nodePropertyExtractorMagnitudesAbsolute(enumerationComponentTypeEncode(char(component),includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
   return
  end function magnitudesAbsoluteConstructorParameters

  function magnitudesAbsoluteConstructorInternal(component) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorMagnitudesAbsolute} property extractor class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorMagnitudesAbsolute)                :: self
    type(enumerationComponentTypeType           ), intent(in   ) :: component
    !![
    <constructorAssign variables="component"/>
    !!]

    return
  end function magnitudesAbsoluteConstructorInternal
  
  integer function magnitudesAbsoluteElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily magnitudesAbsolute} property extractor class.
    !!}
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    class           (nodePropertyExtractorMagnitudesAbsolute), intent(inout) :: self
    double precision                                         , intent(in   ) :: time
    !$GLC attributes unused :: self

    magnitudesAbsoluteElementCount=unitStellarLuminosities%luminosityOutputCount(time)
    return
  end function magnitudesAbsoluteElementCount

  function magnitudesAbsoluteExtract(self,node,time,instance) result(magnitudes)
    !!{
    Implement a {\normalfont \ttfamily magnitudesAbsolute} property extractor.
    !!}
    use :: Galacticus_Nodes              , only : nodeComponentDisk  , nodeComponentSpheroid
    use :: Galactic_Structure_Options    , only : componentTypeDisk  , componentTypeSpheroid
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities, unitStellarLuminosities
    use :: Error                         , only : Error_Report
    implicit none
    double precision                                         , dimension(:) , allocatable :: magnitudes
    class           (nodePropertyExtractorMagnitudesAbsolute), intent(inout), target      :: self
    type            (treeNode                               ), intent(inout), target      :: node
    double precision                                         , intent(in   )              :: time
    type            (multiCounter                           ), intent(inout), optional    :: instance
    class           (nodeComponentDisk                      )               , pointer     :: disk
    class           (nodeComponentSpheroid                  )               , pointer     :: spheroid
    type            (stellarLuminosities                    )                             :: luminosities
    integer                                                                               :: i           , j
    !$GLC attributes unused :: instance
    
    allocate(magnitudes(unitStellarLuminosities%luminosityOutputCount(time)))
    select case (self%component%ID)
    case (componentTypeDisk    %ID)
       disk         => node    %disk                   ()
       luminosities =  disk    %luminositiesStellar    ()
    case (componentTypeSpheroid%ID)
       spheroid     => node    %spheroid               ()
       luminosities =  spheroid%luminositiesStellar    ()
    case default
       luminosities =           unitStellarLuminosities
       call Error_Report("only 'disk' and 'spheroid' components are supported"//{introspection:location})
    end select
    ! Convert luminosities to magnitudes.
    j=0
    do i=1,unitStellarLuminosities%luminosityCount()
       if (unitStellarLuminosities%isOutput(i,time)) then
          j=j+1
          magnitudes(j)=luminosities%luminosity(i)
       end if
    end do
    where (magnitudes > 0.0d0)
       ! Where the luminosity is positive, compute the corresponding absolute magnitude.
       magnitudes=-2.5d0             &
            &     *log10(magnitudes)
    elsewhere
       ! For non-positive luminosities, return the faintest absolute magnitude possible.
       magnitudes=+huge (0.0d0     )
    end where
    return
  end function magnitudesAbsoluteExtract

  subroutine magnitudesAbsoluteNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily magnitudesAbsolute} properties.
    !!}
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    use :: Galactic_Structure_Options    , only : enumerationComponentTypeDecode
    implicit none
    class           (nodePropertyExtractorMagnitudesAbsolute), intent(inout)                             :: self
    double precision                                         , intent(in   )                             :: time
    type            (varying_string                         ), intent(inout), dimension(:) , allocatable :: names
    integer                                                                                              :: i    , j
    !$GLC attributes unused :: self

    allocate(names(unitStellarLuminosities%luminosityOutputCount(time)))
    j=0
    do i=1,unitStellarLuminosities%luminosityCount()
       if (unitStellarLuminosities%isOutput(i,time)) then
          j=j+1
          names(j)=enumerationComponentTypeDecode(self%component)//'MagnitudeAbsoluteStellar:'//unitStellarLuminosities%name(i)
       end if
    end do
    return
  end subroutine magnitudesAbsoluteNames

  subroutine magnitudesAbsoluteDescriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily magnitudesAbsolute} property extractor class.
    !!}
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    class           (nodePropertyExtractorMagnitudesAbsolute), intent(inout)                             :: self
    double precision                                         , intent(in   )                             :: time
    type            (varying_string                         ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(unitStellarLuminosities%luminosityOutputCount(time)))
    descriptions=var_str('Stellar absolute magnitude [AB system].')
    return
  end subroutine magnitudesAbsoluteDescriptions

  function magnitudesAbsoluteUnitsInSI(self,time) result(unitsInSI)
    !!{
    Return the units of absolute magnitude in the SI system.
    !!}
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    double precision                                         , allocatable  , dimension(:) :: unitsInSI
    class           (nodePropertyExtractorMagnitudesAbsolute), intent(inout)               :: self
    double precision                                         , intent(in   )               :: time
    !$GLC attributes unused :: self

    allocate(unitsInSI(unitStellarLuminosities%luminosityOutputCount(time)))
    unitsInSI=0.0d0
    return
  end function magnitudesAbsoluteUnitsInSI

