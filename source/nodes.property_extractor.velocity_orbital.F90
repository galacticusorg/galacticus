!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements an orbital velocity output analysis property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorVelocityOrbital">
   <description>An orbital velocity output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorVelocityOrbital
     !!{
     An orbital velocity property extractor output analysis class.
     !!}
     private
   contains
     procedure :: extract     => velocityOrbitalExtract
     procedure :: type        => velocityOrbitalType
     procedure :: name        => velocityOrbitalName
     procedure :: description => velocityOrbitalDescription
     procedure :: unitsInSI   => velocityOrbitalUnitsInSI
  end type nodePropertyExtractorVelocityOrbital

  interface nodePropertyExtractorVelocityOrbital
     !!{
     Constructors for the ``velocityOrbital'' output analysis class.
     !!}
     module procedure velocityOrbitalConstructorParameters
  end interface nodePropertyExtractorVelocityOrbital

contains

  function velocityOrbitalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``velocityOrbital'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorVelocityOrbital)                :: self
    type(inputParameters                     ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=nodePropertyExtractorVelocityOrbital()
    return
  end function velocityOrbitalConstructorParameters

  double precision function velocityOrbitalExtract(self,node,instance)
    !!{
    Implement a velocityOrbital output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    use :: Vectors         , only : Vector_Magnitude
    implicit none
    class           (nodePropertyExtractorVelocityOrbital), intent(inout)           :: self
    type            (treeNode                            ), intent(inout), target   :: node
    type            (multiCounter                        ), intent(inout), optional :: instance
    type            (treeNode                            ), pointer                 :: nodeWork
    class           (nodeComponentSatellite              ), pointer                 :: satellite
    double precision                                      , dimension(3)            :: velocity
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
    velocityOrbitalExtract=Vector_Magnitude(velocity)
    return
  end function velocityOrbitalExtract

  integer function velocityOrbitalType(self)
    !!{
    Return the type of the halo mass property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorVelocityOrbital), intent(inout) :: self
    !$GLC attributes unused :: self

    velocityOrbitalType=outputAnalysisPropertyTypeLinear
    return
  end function velocityOrbitalType

  function velocityOrbitalName(self)
    !!{
    Return the name of the velocityOrbital property.
    !!}
    implicit none
    type (varying_string                      )                :: velocityOrbitalName
    class(nodePropertyExtractorVelocityOrbital), intent(inout) :: self
    !$GLC attributes unused :: self

    velocityOrbitalName=var_str('velocityOrbitalCurrent')
    return
  end function velocityOrbitalName

  function velocityOrbitalDescription(self)
    !!{
    Return a description of the velocityOrbital property.
    !!}
    implicit none
    type (varying_string                      )                :: velocityOrbitalDescription
    class(nodePropertyExtractorVelocityOrbital), intent(inout) :: self
    !$GLC attributes unused :: self

    velocityOrbitalDescription=var_str('The orbital velocity of the halo.')
    return
  end function velocityOrbitalDescription

  double precision function velocityOrbitalUnitsInSI(self)
    !!{
    Return the units of the velocityOrbital property in the SI system.
    !!}
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    class(nodePropertyExtractorVelocityOrbital), intent(inout) :: self
    !$GLC attributes unused :: self

    velocityOrbitalUnitsInSI=kilo
    return
  end function velocityOrbitalUnitsInSI
