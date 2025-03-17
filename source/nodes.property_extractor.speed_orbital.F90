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
Implements an orbital speed output analysis property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorSpeedOrbital">
   <description>An orbital speed output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorSpeedOrbital
     !!{
     An orbital speed property extractor output analysis class.
     !!}
     private
   contains
     procedure :: extract     => speedOrbitalExtract
     procedure :: name        => speedOrbitalName
     procedure :: description => speedOrbitalDescription
     procedure :: unitsInSI   => speedOrbitalUnitsInSI
  end type nodePropertyExtractorSpeedOrbital

  interface nodePropertyExtractorSpeedOrbital
     !!{
     Constructors for the ``speedOrbital'' output analysis class.
     !!}
     module procedure speedOrbitalConstructorParameters
  end interface nodePropertyExtractorSpeedOrbital

contains

  function speedOrbitalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``speedOrbital'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorSpeedOrbital)                :: self
    type(inputParameters                  ), intent(inout) :: parameters

    self=nodePropertyExtractorSpeedOrbital()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function speedOrbitalConstructorParameters

  double precision function speedOrbitalExtract(self,node,instance)
    !!{
    Implement a speedOrbital output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    use :: Vectors         , only : Vector_Magnitude
    implicit none
    class           (nodePropertyExtractorSpeedOrbital), intent(inout), target   :: self
    type            (treeNode                         ), intent(inout), target   :: node
    type            (multiCounter                     ), intent(inout), optional :: instance
    type            (treeNode                         ), pointer                 :: nodeWork
    class           (nodeComponentSatellite           ), pointer                 :: satellite
    double precision                                   , dimension(3)            :: velocity
    !$GLC attributes unused :: self, instance

    velocity =  0.0d0
    nodeWork => node
    ! Walk up through all host halos of this node, accumulating velocity offsets from the host node center.
    do while (associated(nodeWork))
       satellite =>  nodeWork %satellite()
       velocity  =  +          velocity    &
            &       +satellite%velocity ()
       if (nodeWork%isSatellite()) then
          nodeWork => nodeWork%parent
       else
          nodeWork => null()
       end if
    end do
    speedOrbitalExtract=Vector_Magnitude(velocity)
    return
  end function speedOrbitalExtract


  function speedOrbitalName(self)
    !!{
    Return the name of the speedOrbital property.
    !!}
    implicit none
    type (varying_string                   )                :: speedOrbitalName
    class(nodePropertyExtractorSpeedOrbital), intent(inout) :: self
    !$GLC attributes unused :: self

    speedOrbitalName=var_str('speedOrbitalCurrent')
    return
  end function speedOrbitalName

  function speedOrbitalDescription(self)
    !!{
    Return a description of the speedOrbital property.
    !!}
    implicit none
    type (varying_string                   )                :: speedOrbitalDescription
    class(nodePropertyExtractorSpeedOrbital), intent(inout) :: self
    !$GLC attributes unused :: self

    speedOrbitalDescription=var_str('The orbital speed of the halo.')
    return
  end function speedOrbitalDescription

  double precision function speedOrbitalUnitsInSI(self)
    !!{
    Return the units of the speedOrbital property in the SI system.
    !!}
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    class(nodePropertyExtractorSpeedOrbital), intent(inout) :: self
    !$GLC attributes unused :: self

    speedOrbitalUnitsInSI=kilo
    return
  end function speedOrbitalUnitsInSI
