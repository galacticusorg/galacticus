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
Implements a massHost property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassHost">
   <description>A host halo mass property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMassHost
     !!{
     A host halo mass property extractor class.
     !!}
     private
   contains
     procedure :: extract     => massHostExtract
     procedure :: name        => massHostName
     procedure :: description => massHostDescription
     procedure :: unitsInSI   => massHostUnitsInSI
  end type nodePropertyExtractorMassHost

  interface nodePropertyExtractorMassHost
     !!{
     Constructors for the \refClass{nodePropertyExtractorMassHost} output analysis class.
     !!}
     module procedure massHostConstructorParameters
  end interface nodePropertyExtractorMassHost

contains

  function massHostConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorMassHost} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorMassHost)                :: self
    type(inputParameters              ), intent(inout) :: parameters

    self=nodePropertyExtractorMassHost()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massHostConstructorParameters

  double precision function massHostExtract(self,node,instance)
    !!{
    Implement a last isolated redshift output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (nodePropertyExtractorMassHost), intent(inout), target   :: self
    type            (treeNode                     ), intent(inout), target   :: node
    type            (multiCounter                 ), intent(inout), optional :: instance
    type            (treeNode                     ), pointer                 :: nodeHost
    class           (nodeComponentBasic           ), pointer                 :: basicHost
    !$GLC attributes unused :: self, instance

    if (node%isSatellite()) then
       nodeHost        => node     %parent
       basicHost       => nodeHost %basic ()
       massHostExtract =  basicHost%mass  ()
    else
       massHostExtract =  -1.0d0
    end if
    return
  end function massHostExtract

  function massHostName(self)
    !!{
    Return the name of the host halo mass property.
    !!}
    implicit none
    type (varying_string               )                :: massHostName
    class(nodePropertyExtractorMassHost), intent(inout) :: self
    !$GLC attributes unused :: self

    massHostName=var_str('satelliteHostHaloMass')
    return
  end function massHostName

  function massHostDescription(self)
    !!{
    Return a description of the host halo mass property.
    !!}
    implicit none
    type (varying_string               )                :: massHostDescription
    class(nodePropertyExtractorMassHost), intent(inout) :: self
    !$GLC attributes unused :: self

    massHostDescription=var_str("Mass of the satellite's host halo [M☉].")
    return
  end function massHostDescription

  double precision function massHostUnitsInSI(self)
    !!{
    Return the units of the host halo mass in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(nodePropertyExtractorMassHost), intent(inout) :: self
    !$GLC attributes unused :: self

    massHostUnitsInSI=massSolar
    return
  end function massHostUnitsInSI

