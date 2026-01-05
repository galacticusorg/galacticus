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

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScale, darkMatterHaloScaleClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorVirialProperties">
   <description>
    A node property extractor which extracts the following quantities related to the virialized region of each node:
    \begin{description}
     \item [{\normalfont \ttfamily nodeVirialRadius}] The virial radius (following whatever definition of virial overdensity is
     specified by the virial density contrast (see \refPhysics{virialDensityContrast}) in units of Mpc;
     \item [{\normalfont \ttfamily nodeVirialVelocity}] The circular velocity at the virial radius (in km/s).
    \end{description}
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorVirialProperties
     !!{
     A property extractor which extracts virialProperties properties.
     !!}
     private
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::                 virialPropertiesDestructor
     procedure :: elementCount => virialPropertiesElementCount
     procedure :: extract      => virialPropertiesExtract
     procedure :: names        => virialPropertiesNames
     procedure :: descriptions => virialPropertiesDescriptions
     procedure :: unitsInSI    => virialPropertiesUnitsInSI
  end type nodePropertyExtractorVirialProperties

  interface nodePropertyExtractorVirialProperties
     !!{
     Constructors for the \refClass{nodePropertyExtractorVirialProperties} output extractor class.
     !!}
     module procedure virialPropertiesConstructorParameters
     module procedure virialPropertiesConstructorInternal
  end interface nodePropertyExtractorVirialProperties

contains

  function virialPropertiesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorVirialProperties} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorVirialProperties)                :: self
    type (inputParameters                      ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass             ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=nodePropertyExtractorVirialProperties(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    !!]
    return
  end function virialPropertiesConstructorParameters

  function virialPropertiesConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorVirialProperties} output extractor property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorVirialProperties)                        :: self
    class(darkMatterHaloScaleClass             ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    return
  end function virialPropertiesConstructorInternal

  subroutine virialPropertiesDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorVirialProperties} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorVirialProperties), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine virialPropertiesDestructor

  integer function virialPropertiesElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily virialProperties} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorVirialProperties), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    virialPropertiesElementCount=3
    return
  end function virialPropertiesElementCount

  function virialPropertiesExtract(self,node,time,instance)
    !!{
    Implement a virialProperties output extractor.
    !!}
    implicit none
    double precision                                       , dimension(:) , allocatable :: virialPropertiesExtract
    class           (nodePropertyExtractorVirialProperties), intent(inout), target      :: self
    type            (treeNode                             ), intent(inout), target      :: node
    double precision                                       , intent(in   )              :: time
    type            (multiCounter                         ), intent(inout), optional    :: instance
    !$GLC attributes unused :: time, instance

    allocate(virialPropertiesExtract(3))
    virialPropertiesExtract=[                                                   &
         &                   self%darkMatterHaloScale_%radiusVirial     (node), &
         &                   self%darkMatterHaloScale_%velocityVirial   (node), &
         &                   self%darkMatterHaloScale_%temperatureVirial(node)  &
         &                  ]
    return
  end function virialPropertiesExtract

  subroutine virialPropertiesNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily virialProperties} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorVirialProperties), intent(inout)                             :: self
    double precision                                       , intent(in   )                             :: time
    type            (varying_string                       ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(3))
    names(1)=var_str('darkMatterOnlyRadiusVirial'     )
    names(2)=var_str('darkMatterOnlyVelocityVirial'   )
    names(3)=var_str('darkMatterOnlyTemperatureVirial')
    return
  end subroutine virialPropertiesNames

  subroutine virialPropertiesDescriptions(self,time,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily virialProperties} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorVirialProperties), intent(inout)                             :: self
    double precision                                       , intent(in   )                             :: time
    type            (varying_string                       ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(3))
    descriptions(1)=var_str('Virial radius of the dark matter only node [Mpc].'   )
    descriptions(2)=var_str('Virial velocity of the dark matter only node [km/s].')
    descriptions(3)=var_str('Virial temperature of the dark matter only node [K].')
    return
  end subroutine virialPropertiesDescriptions

  function virialPropertiesUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily virialProperties} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                       , dimension(:) , allocatable :: virialPropertiesUnitsInSI
    class           (nodePropertyExtractorVirialProperties), intent(inout)              :: self
    double precision                                       , intent(in   )              :: time
    !$GLC attributes unused :: self, time

    allocate(virialPropertiesUnitsInSI(3))
    virialPropertiesUnitsInSI=[            &
         &                     megaParsec, &
         &                     kilo      , &
         &                     1.0d0       &
         &                    ]
    return
  end function virialPropertiesUnitsInSI

