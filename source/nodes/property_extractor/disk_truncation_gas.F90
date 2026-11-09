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

!+    Contributions to this file made by: Claude.

!!{RST
Implements a property extractor class that extracts the degree of truncation of the gas of the disk.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorDiskTruncationGas" docformat="rst">
   <description>
   A node property extractor which extracts the degree of truncation, :math:`x_\mathrm{gas}=-\log f_\mathrm{gas}`, of the gas
   of the disk, where :math:`f_\mathrm{gas}` is the fraction of the mass of its untruncated mass distribution which the gas
   retains (see :galacticus-class:`galacticComponentTruncationClass`). This is increased by tidal stripping (see
   :galacticus-class:`nodeOperatorTidalMassLossDisks`) and ram pressure stripping (see
   :galacticus-class:`nodeOperatorRamPressureMassLossDisks`) when the disk is truncated. Zero is returned if the disk does not
   exist, or if its gas is not truncated. (Small negative values, which the ODE solver can produce as gas arrives in the disk,
   are returned as zero, as they are treated as zero wherever the truncation is used.)
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorDiskTruncationGas
     !!{RST
     A property extractor which extracts the degree of truncation of the gas of the disk.
     !!}
     private
     integer :: truncationGasID
   contains
     procedure :: extract     => diskTruncationGasExtract
     procedure :: name        => diskTruncationGasName
     procedure :: description => diskTruncationGasDescription
     procedure :: unitsInSI   => diskTruncationGasUnitsInSI
     procedure :: units       => diskTruncationGasUnits
  end type nodePropertyExtractorDiskTruncationGas

  interface nodePropertyExtractorDiskTruncationGas
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorDiskTruncationGas` property extractor class.
     !!}
     module procedure diskTruncationGasConstructorParameters
     module procedure diskTruncationGasConstructorInternal
  end interface nodePropertyExtractorDiskTruncationGas

contains

  function diskTruncationGasConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorDiskTruncationGas` property extractor class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorDiskTruncationGas)                :: self
    type(inputParameters                       ), intent(inout) :: parameters

    self=nodePropertyExtractorDiskTruncationGas()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function diskTruncationGasConstructorParameters

  function diskTruncationGasConstructorInternal() result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorDiskTruncationGas` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorDiskTruncationGas) :: self

    !![
    <addMetaProperty component="disk" name="truncationGas" id="self%truncationGasID" isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function diskTruncationGasConstructorInternal

  double precision function diskTruncationGasExtract(self,node,instance)
    !!{RST
    Implement a ``diskTruncationGas`` property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, treeNode
    implicit none
    class(nodePropertyExtractorDiskTruncationGas), intent(inout), target   :: self
    type (treeNode                              ), intent(inout), target   :: node
    type (multiCounter                          ), intent(inout), optional :: instance
    class(nodeComponentDisk                     ), pointer                 :: disk
    !$GLC attributes unused :: instance

    diskTruncationGasExtract=0.0d0
    disk => node%disk()
    select type (disk)
    type is (nodeComponentDisk)
       ! Disk does not yet exist.
    class default
       if (disk%floatRank0MetaPropertyIsCreated(self%truncationGasID)) &
            & diskTruncationGasExtract=max(0.0d0,disk%floatRank0MetaPropertyGet(self%truncationGasID))
    end select
    return
  end function diskTruncationGasExtract

  function diskTruncationGasName(self)
    !!{RST
    Return the name of the ``diskTruncationGas`` property.
    !!}
    implicit none
    type (varying_string                        )                :: diskTruncationGasName
    class(nodePropertyExtractorDiskTruncationGas), intent(inout) :: self
    !$GLC attributes unused :: self

    diskTruncationGasName=var_str('diskTruncationGas')
    return
  end function diskTruncationGasName

  function diskTruncationGasDescription(self)
    !!{RST
    Return a description of the ``diskTruncationGas`` property.
    !!}
    implicit none
    type (varying_string                        )                :: diskTruncationGasDescription
    class(nodePropertyExtractorDiskTruncationGas), intent(inout) :: self
    !$GLC attributes unused :: self

    diskTruncationGasDescription=var_str('Degree of truncation, -log(f), of the gas of the disk, where f is the fraction of its untruncated mass retained.')
    return
  end function diskTruncationGasDescription

  double precision function diskTruncationGasUnitsInSI(self)
    !!{RST
    Return the units of the ``diskTruncationGas`` property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorDiskTruncationGas), intent(inout) :: self
    !$GLC attributes unused :: self

    diskTruncationGasUnitsInSI=1.0d0
    return
  end function diskTruncationGasUnitsInSI

  function diskTruncationGasUnits(self) result(units)
    !!{RST
    Return the units of the ``diskTruncationGas`` property.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type (unitType                              )                :: units
    class(nodePropertyExtractorDiskTruncationGas), intent(inout) :: self
    !$GLC attributes unused :: self

    units=unitType(1.0d0)
    return
  end function diskTruncationGasUnits
