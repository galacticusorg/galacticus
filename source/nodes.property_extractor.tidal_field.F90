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

!+    Contributions to this file made by: Andrew Benson, Charles Gannon.

!!{
Implements a tidal field property extractor class.
!!}

  use :: Satellites_Tidal_Fields, only : satelliteTidalFieldClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorTidalField">
   <description> 
    A property extractor class which extracts the radial component of the tidal tensor in units of $\mathrm{Gyr^{-2}}$.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorTidalField
     !!{
      A property extractor class which extracts the radial component of the tidal tensor in units of $\mathrm{Gyr^{-2}}$.
     !!}
     class(satelliteTidalFieldClass), pointer :: satelliteTidalField_ => null()
   contains
     final     ::                tidalFieldDestructor
     procedure :: extract     => tidalFieldExtract
     procedure :: name        => tidalFieldName
     procedure :: description => tidalFieldDescription
     procedure :: unitsInSI   => tidalFieldUnitsInSI
  end type nodePropertyExtractorTidalField

  interface nodePropertyExtractorTidalField
     !!{
     Constructors for the \refClass{nodePropertyExtractorTidalField} output analysis class.
     !!}
     module procedure tidalFieldConstructorParameters
     module procedure tidalFieldConstructorInternal
  end interface nodePropertyExtractorTidalField

contains

  function tidalFieldConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorTidalField} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorTidalField)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(satelliteTidalFieldClass       ), pointer       :: satelliteTidalField_
    
    !![
    <objectBuilder class="satelliteTidalField" name="satelliteTidalField_" source="parameters"/>
    !!]
    self=nodePropertyExtractorTidalField(satelliteTidalField_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteTidalField_"/>
    !!]
    return
  end function tidalFieldConstructorParameters

  function tidalFieldConstructorInternal(satelliteTidalField_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorTidalField} class.
    !!}
    implicit none
    type (nodePropertyExtractorTidalField)                        :: self
    class(satelliteTidalFieldClass       ), intent(in   ), target :: satelliteTidalField_
    !![
    <constructorAssign variables="*satelliteTidalField_"/>
    !!]

    return
  end function tidalFieldConstructorInternal

  subroutine tidalFieldDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorTidalField} class.
    !!}
    implicit none
    type(nodePropertyExtractorTidalField), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteTidalField_"/>
    !!]
    return
  end subroutine tidalFieldDestructor

  double precision function tidalFieldExtract(self,node,instance)
    !!{
    Return the radial part of the tidal tensor for satellite halos.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear, megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class(nodePropertyExtractorTidalField), intent(inout), target   :: self
    type (treeNode                       ), intent(inout), target   :: node
    type (multiCounter                   ), intent(inout), optional :: instance 
    !$GLC attributes unused :: instance

    if (node%isSatellite()) then
       ! Get the tidal field and convert to Gyr⁻² units.
       tidalFieldExtract=+self%satelliteTidalField_%tidalTensorRadial(node,includeCentrifugalAcceleration=.true.) &
            &            *gigaYear  **2                                                                           &
            &            *kilo      **2                                                                           &
            &            /megaParsec**2
    else
       ! For isolated halos, always return zero tidal field.
       tidalFieldExtract=+0.0d0
    end if
    return
  end function tidalFieldExtract

  function tidalFieldName(self)
    !!{
    Return the name of the tidal radius property.
    !!}
    implicit none
    type (varying_string                 )                :: tidalFieldName
    class(nodePropertyExtractorTidalField), intent(inout) :: self
    !$GLC attributes unused :: self

    tidalFieldName=var_str('satelliteTidalField')
    return
  end function tidalFieldName

  function tidalFieldDescription(self)
    !!{
    Return a description of the tidal radius property.
    !!}
    implicit none
    type (varying_string                 )                :: tidalFieldDescription
    class(nodePropertyExtractorTidalField), intent(inout) :: self
    !$GLC attributes unused :: self

    tidalFieldDescription=var_str('Tidal field in the halo Gyr⁻².')
    return
  end function tidalFieldDescription

  double precision function tidalFieldUnitsInSI(self)
    !!{
    Return the units of the tidal radius property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    class(nodePropertyExtractorTidalField), intent(inout) :: self
    !$GLC attributes unused :: self

    tidalFieldUnitsInSI=1.0d0/gigaYear**2
    return
  end function tidalFieldUnitsInSI


