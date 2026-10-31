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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Implements a disk heating property extractor class.
!!}

  use :: Galactic_Disk_Heatings, only : galacticDiskHeatingClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorDiskHeating" docformat="rst">
   <description>
   Extracts the dimensionless heating parameter, :math:`s`, of the disk (see :galacticus-class:`galacticDiskHeatingClass`), such
   that a ring of the disk at radius :math:`R` is spread in radius with a dispersion :math:`s R^2/R_\mathrm{d}`. Heating is
   applied only to satellites, so zero is returned for all other nodes, and for nodes with no disk. The heating parameter is
   computed using the :galacticus-class:`galacticDiskHeatingClass` given to this extractor (or found by searching upward from
   it in the parameter file), which should therefore be that used by the disk component. The value returned is that given by
   the :galacticus-class:`galacticDiskHeatingClass`, before any limit imposed when it is applied to the mass distribution of the
   disk (see :galacticus-class:`massDistributionExponentialDiskTidal`).
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorDiskHeating
     !!{RST
     A disk heating property extractor class.
     !!}
     private
     class(galacticDiskHeatingClass), pointer :: galacticDiskHeating_ => null()
   contains
     final     ::                diskHeatingDestructor
     procedure :: extract     => diskHeatingExtract
     procedure :: name        => diskHeatingName
     procedure :: description => diskHeatingDescription
     procedure :: unitsInSI   => diskHeatingUnitsInSI
     procedure :: units       => diskHeatingUnits
  end type nodePropertyExtractorDiskHeating

  interface nodePropertyExtractorDiskHeating
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorDiskHeating` property extractor class.
     !!}
     module procedure diskHeatingConstructorParameters
     module procedure diskHeatingConstructorInternal
  end interface nodePropertyExtractorDiskHeating

contains

  function diskHeatingConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorDiskHeating` property extractor class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorDiskHeating)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(galacticDiskHeatingClass        ), pointer       :: galacticDiskHeating_

    !![
    <objectBuilder class="galacticDiskHeating" name="galacticDiskHeating_" source="parameters"/>
    !!]
    self=nodePropertyExtractorDiskHeating(galacticDiskHeating_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticDiskHeating_"/>
    !!]
    return
  end function diskHeatingConstructorParameters

  function diskHeatingConstructorInternal(galacticDiskHeating_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorDiskHeating` property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorDiskHeating)                        :: self
    class(galacticDiskHeatingClass        ), intent(in   ), target :: galacticDiskHeating_
    !![
    <constructorAssign variables="*galacticDiskHeating_"/>
    !!]

    return
  end function diskHeatingConstructorInternal

  subroutine diskHeatingDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodePropertyExtractorDiskHeating` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorDiskHeating), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticDiskHeating_"/>
    !!]
    return
  end subroutine diskHeatingDestructor

  double precision function diskHeatingExtract(self,node,instance)
    !!{RST
    Implement a disk heating property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSatellite
    implicit none
    class(nodePropertyExtractorDiskHeating), intent(inout), target   :: self
    type (treeNode                        ), intent(inout), target   :: node
    type (multiCounter                    ), intent(inout), optional :: instance
    class(nodeComponentDisk               ), pointer                 :: disk
    class(nodeComponentSatellite          ), pointer                 :: satellite
    !$GLC attributes unused :: instance

    diskHeatingExtract=0.0d0
    if (.not.self%galacticDiskHeating_%isActive() .or. .not.node%isSatellite()) return
    disk => node%disk()
    select type (disk)
    type is (nodeComponentDisk)
       ! Disk does not yet exist.
    class default
       satellite          => node                     %satellite(                                                                )
       diskHeatingExtract =  self%galacticDiskHeating_%heating  (satellite%tidalHeatingNormalized(),disk%radius(),disk%velocity())
    end select
    return
  end function diskHeatingExtract

  function diskHeatingName(self)
    !!{RST
    Return the name of the disk heating property.
    !!}
    implicit none
    type (varying_string                  )                :: diskHeatingName
    class(nodePropertyExtractorDiskHeating), intent(inout) :: self
    !$GLC attributes unused :: self

    diskHeatingName=var_str('diskTidalHeating')
    return
  end function diskHeatingName

  function diskHeatingDescription(self)
    !!{RST
    Return a description of the disk heating property.
    !!}
    implicit none
    type (varying_string                  )                :: diskHeatingDescription
    class(nodePropertyExtractorDiskHeating), intent(inout) :: self
    !$GLC attributes unused :: self

    diskHeatingDescription=var_str('Dimensionless tidal heating parameter, s, of the disk, such that a ring at radius R is spread in radius with dispersion s R²/R_d.')
    return
  end function diskHeatingDescription

  double precision function diskHeatingUnitsInSI(self)
    !!{RST
    Return the units of the disk heating property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorDiskHeating), intent(inout) :: self
    !$GLC attributes unused :: self

    diskHeatingUnitsInSI=1.0d0
    return
  end function diskHeatingUnitsInSI

  function diskHeatingUnits(self) result(units)
    !!{RST
    Return the units of the disk heating property.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type (unitType                        )                :: units
    class(nodePropertyExtractorDiskHeating), intent(inout) :: self
    !$GLC attributes unused :: self

    units=unitType(1.0d0)
    return
  end function diskHeatingUnits
