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
Implements a node property extractor class for apparent magnitudes.
!!}

  use :: Galactic_Structure_Options, only : enumerationComponentTypeType
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorMagnitudesApparent">
   <description>
    A node property extractor which extracts stellar apparent magnitudes in all available bands.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorMagnitudesApparent
     !!{
     A node property extractor which extracts stellar apparent magnitudes in all available bands.
     !!}
     private
     class(cosmologyFunctionsClass     ), pointer :: cosmologyFunctions_ => null()
     type (enumerationComponentTypeType)          :: component
   contains
     final     ::                 magnitudesApparentDestructor
     procedure :: elementCount => magnitudesApparentElementCount
     procedure :: extract      => magnitudesApparentExtract
     procedure :: names        => magnitudesApparentNames
     procedure :: descriptions => magnitudesApparentDescriptions
     procedure :: unitsInSI    => magnitudesApparentUnitsInSI
  end type nodePropertyExtractorMagnitudesApparent

  interface nodePropertyExtractorMagnitudesApparent
     !!{
     Constructors for the {\normalfont \ttfamily magnitudesApparent} output analysis class.
     !!}
     module procedure magnitudesApparentConstructorParameters
     module procedure magnitudesApparentConstructorInternal
  end interface nodePropertyExtractorMagnitudesApparent

contains

  function magnitudesApparentConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily magnitudesApparent} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode
    implicit none
    type (nodePropertyExtractorMagnitudesApparent)                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    type (varying_string                         )                :: component
    
    !![
    <inputParameter>
      <name>component</name>
      <source>parameters</source>
      <description>The component from which to extract star formation rate.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=nodePropertyExtractorMagnitudesApparent(enumerationComponentTypeEncode(char(component),includesPrefix=.false.),cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
   return
  end function magnitudesApparentConstructorParameters

  function magnitudesApparentConstructorInternal(component,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily magnitudesApparent} property extractor class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorMagnitudesApparent)                        :: self
    type (enumerationComponentTypeType           ), intent(in   )         :: component
    class(cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="component, *cosmologyFunctions_"/>
    !!]

    return
  end function magnitudesApparentConstructorInternal
  
  subroutine magnitudesApparentDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily magnitudesApparent} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorMagnitudesApparent), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine magnitudesApparentDestructor

  integer function magnitudesApparentElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily magnitudesApparent} property extractor class.
    !!}
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    class           (nodePropertyExtractorMagnitudesApparent), intent(inout) :: self
    double precision                                         , intent(in   ) :: time
    !$GLC attributes unused :: self

    magnitudesApparentElementCount=unitStellarLuminosities%luminosityOutputCount(time)
    return
  end function magnitudesApparentElementCount

  function magnitudesApparentExtract(self,node,time,instance) result(magnitudes)
    !!{
    Implement a {\normalfont \ttfamily magnitudesApparent} property extractor.
    !!}
    use :: Galacticus_Nodes              , only : nodeComponentDisk  , nodeComponentSpheroid
    use :: Galactic_Structure_Options    , only : componentTypeDisk  , componentTypeSpheroid
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities, unitStellarLuminosities
    use :: Error                         , only : Error_Report
    implicit none
    double precision                                         , dimension(:) , allocatable :: magnitudes
    class           (nodePropertyExtractorMagnitudesApparent), intent(inout), target      :: self
    type            (treeNode                               ), intent(inout), target      :: node
    double precision                                         , intent(in   )              :: time
    type            (multiCounter                           ), intent(inout), optional    :: instance
    class           (nodeComponentDisk                      )               , pointer     :: disk
    class           (nodeComponentSpheroid                  )               , pointer     :: spheroid
    type            (stellarLuminosities                    )                             :: luminosities
    integer                                                                               :: i              , j
    double precision                                                                      :: distanceModulus, distanceLuminosity
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
    ! Find the distance modulus, including the extra -2.5log₁₀(1+z) term that accounts for compression of photon frequencies due to redshifting.
    distanceLuminosity=self%cosmologyFunctions_%distanceLuminosity(time)
    if (distanceLuminosity > 0.0d0) then
       distanceModulus=+25.0d0                                                          &
            &          + 5.0d0*log10(                         distanceLuminosity      ) &
            &          + 2.5d0*log10(self%cosmologyFunctions_%expansionFactor   (time))
       ! Convert luminosities to magnitudes.
       j=0
       do i=1,unitStellarLuminosities%luminosityCount()
          if (unitStellarLuminosities%isOutput(i,time)) then
             j=j+1
             magnitudes(j)=luminosities%luminosity(i)
          end if
       end do
       where (magnitudes > 0.0d0)
          ! Where the luminosity is positive, compute the corresponding apparent magnitude.
          magnitudes=-2.5d0             &
               &     *log10(magnitudes) &
               &     +distanceModulus
       elsewhere
          ! For non-positive luminosities, return the faintest apparent magnitude possible.
          magnitudes=+huge (0.0d0     )
       end where
    else
       ! Luminosity distance is zero, so set apparent magnitudes to the brightest magnitude possible.
       magnitudes   =-huge (0.0d0     )
    end if
    return
  end function magnitudesApparentExtract

  subroutine magnitudesApparentNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily magnitudesApparent} properties.
    !!}
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    use :: Galactic_Structure_Options    , only : enumerationComponentTypeDecode
    implicit none
    class           (nodePropertyExtractorMagnitudesApparent), intent(inout)                             :: self
    double precision                                         , intent(in   )                             :: time
    type            (varying_string                         ), intent(inout), dimension(:) , allocatable :: names
    integer                                                                                              :: i    , j
    !$GLC attributes unused :: self

    allocate(names(unitStellarLuminosities%luminosityOutputCount(time)))
    j=0
    do i=1,unitStellarLuminosities%luminosityCount()
       if (unitStellarLuminosities%isOutput(i,time)) then
          j=j+1
          names(j)=enumerationComponentTypeDecode(self%component)//'MagnitudeApparentStellar:'//unitStellarLuminosities%name(i)
       end if
    end do
    return
  end subroutine magnitudesApparentNames

  subroutine magnitudesApparentDescriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily magnitudesApparent} property extractor class.
    !!}
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    class           (nodePropertyExtractorMagnitudesApparent), intent(inout)                             :: self
    double precision                                         , intent(in   )                             :: time
    type            (varying_string                         ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(unitStellarLuminosities%luminosityOutputCount(time)))
    descriptions=var_str('Stellar apparent magnitude [AB system].')
    return
  end subroutine magnitudesApparentDescriptions

  function magnitudesApparentUnitsInSI(self,time) result(unitsInSI)
    !!{
    Return the units of apparent magnitude in the SI system.
    !!}
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    double precision                                         , allocatable  , dimension(:) :: unitsInSI
    class           (nodePropertyExtractorMagnitudesApparent), intent(inout)               :: self
    double precision                                         , intent(in   )               :: time
    !$GLC attributes unused :: self

    allocate(unitsInSI(unitStellarLuminosities%luminosityOutputCount(time)))
    unitsInSI=0.0d0
    return
  end function magnitudesApparentUnitsInSI

