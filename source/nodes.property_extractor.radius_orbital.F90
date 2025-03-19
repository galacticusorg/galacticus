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
Implements an orbital radius output analysis property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRadiusOrbital">
   <description>An orbital radius output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRadiusOrbital
     !!{
     An orbital radius property extractor output analysis class.
     !!}
     private
   contains
     procedure :: extract     => radiusOrbitalExtract
     procedure :: name        => radiusOrbitalName
     procedure :: description => radiusOrbitalDescription
     procedure :: unitsInSI   => radiusOrbitalUnitsInSI
  end type nodePropertyExtractorRadiusOrbital

  interface nodePropertyExtractorRadiusOrbital
     !!{
     Constructors for the ``radiusOrbital'' output analysis class.
     !!}
     module procedure radiusOrbitalConstructorParameters
  end interface nodePropertyExtractorRadiusOrbital

contains

  function radiusOrbitalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``radiusOrbital'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorRadiusOrbital)                :: self
    type(inputParameters                   ), intent(inout) :: parameters

    self=nodePropertyExtractorRadiusOrbital()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function radiusOrbitalConstructorParameters

  double precision function radiusOrbitalExtract(self,node,instance)
    !!{
    Implement a radiusOrbital output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    use :: Vectors         , only : Vector_Magnitude
    implicit none
    class           (nodePropertyExtractorRadiusOrbital), intent(inout), target   :: self
    type            (treeNode                          ), intent(inout), target   :: node
    type            (multiCounter                      ), intent(inout), optional :: instance
    type            (treeNode                          ), pointer                 :: nodeWork
    class           (nodeComponentSatellite            ), pointer                 :: satellite
    double precision                                    , dimension(3)            :: position
    !$GLC attributes unused :: self, instance

    position =  0.0d0
    nodeWork => node
    ! Walk up through all host halos of this node, accumulating position offsets from the host node center.
    do while (associated(nodeWork))
       satellite =>  nodeWork %satellite()
       position  =  +          position    &
            &       +satellite%position ()
       if (nodeWork%isSatellite()) then
          nodeWork => nodeWork%parent
       else
          nodeWork => null()
       end if
    end do
    radiusOrbitalExtract=Vector_Magnitude(position)
    return
  end function radiusOrbitalExtract


  function radiusOrbitalName(self)
    !!{
    Return the name of the radiusOrbital property.
    !!}
    implicit none
    type (varying_string                    )                :: radiusOrbitalName
    class(nodePropertyExtractorRadiusOrbital), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusOrbitalName=var_str('radiusOrbitalCurrent')
    return
  end function radiusOrbitalName

  function radiusOrbitalDescription(self)
    !!{
    Return a description of the radiusOrbital property.
    !!}
    implicit none
    type (varying_string                    )                :: radiusOrbitalDescription
    class(nodePropertyExtractorRadiusOrbital), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusOrbitalDescription=var_str('The orbital radius of the halo.')
    return
  end function radiusOrbitalDescription

  double precision function radiusOrbitalUnitsInSI(self)
    !!{
    Return the units of the radiusOrbital property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorRadiusOrbital), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusOrbitalUnitsInSI=megaParsec
    return
  end function radiusOrbitalUnitsInSI
