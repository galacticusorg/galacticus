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
  Implements a property extractor class for the ram pressure history of satellites.
  !!}

  use :: Hot_Halo_Ram_Pressure_Forces, only : hotHaloRamPressureForceClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRamPressureHistory" docformat="rst">
   <description>
   A property extractor class which extracts the ram pressure history of satellites, as accumulated by the
   :galacticus-class:`nodeOperatorRamPressureHistory` node operator: the (decaying) running maximum of the ram pressure,
   ``satelliteRamPressureMaximum``, and its (decaying) time integral, ``satelliteRamPressureTimeIntegrated``. The instantaneous
   ram pressure, ``satelliteRamPressure``, computed using the :galacticus-class:`hotHaloRamPressureForceClass`, is also
   extracted. All three are zero for non-satellite nodes. The :galacticus-class:`nodeOperatorRamPressureHistory` node operator
   must be included in the model.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorRamPressureHistory
     !!{RST
     A property extractor class for the ram pressure history of satellites.
     !!}
     private
     class  (hotHaloRamPressureForceClass), pointer :: hotHaloRamPressureForce_    => null()
     integer                                        :: ramPressureTimeIntegratedID          , ramPressureMaximumID
   contains
     final     ::                 ramPressureHistoryDestructor
     procedure :: elementCount => ramPressureHistoryElementCount
     procedure :: extract      => ramPressureHistoryExtract
     procedure :: names        => ramPressureHistoryNames
     procedure :: descriptions => ramPressureHistoryDescriptions
     procedure :: unitsInSI    => ramPressureHistoryUnitsInSI
     procedure :: units        => ramPressureHistoryUnits
  end type nodePropertyExtractorRamPressureHistory

  interface nodePropertyExtractorRamPressureHistory
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorRamPressureHistory` property extractor class.
     !!}
     module procedure ramPressureHistoryConstructorParameters
     module procedure ramPressureHistoryConstructorInternal
  end interface nodePropertyExtractorRamPressureHistory

contains

  function ramPressureHistoryConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorRamPressureHistory` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorRamPressureHistory)                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(hotHaloRamPressureForceClass           ), pointer       :: hotHaloRamPressureForce_

    !![
    <objectBuilder class="hotHaloRamPressureForce" name="hotHaloRamPressureForce_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRamPressureHistory(hotHaloRamPressureForce_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="hotHaloRamPressureForce_"/>
    !!]
    return
  end function ramPressureHistoryConstructorParameters

  function ramPressureHistoryConstructorInternal(hotHaloRamPressureForce_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorRamPressureHistory` property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorRamPressureHistory)                        :: self
    class(hotHaloRamPressureForceClass           ), intent(in   ), target :: hotHaloRamPressureForce_
    !![
    <constructorAssign variables="*hotHaloRamPressureForce_"/>
    <addMetaProperty component="satellite" name="ramPressureTimeIntegrated" id="self%ramPressureTimeIntegratedID" isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="satellite" name="ramPressureMaximum"        id="self%ramPressureMaximumID"        isEvolvable="yes" isCreator="no"/>
    !!]

    return
  end function ramPressureHistoryConstructorInternal

  subroutine ramPressureHistoryDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodePropertyExtractorRamPressureHistory` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRamPressureHistory), intent(inout) :: self

    !![
    <objectDestructor name="self%hotHaloRamPressureForce_"/>
    !!]
    return
  end subroutine ramPressureHistoryDestructor

  integer function ramPressureHistoryElementCount(self,time)
    !!{RST
    Return the number of elements in the ram pressure history property extractor.
    !!}
    implicit none
    class           (nodePropertyExtractorRamPressureHistory), intent(inout) :: self
    double precision                                         , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    ramPressureHistoryElementCount=3
    return
  end function ramPressureHistoryElementCount

  function ramPressureHistoryExtract(self,node,time,instance)
    !!{RST
    Implement a ram pressure history property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite, treeNode
    implicit none
    double precision                                         , dimension(:) , allocatable :: ramPressureHistoryExtract
    class           (nodePropertyExtractorRamPressureHistory), intent(inout), target      :: self
    type            (treeNode                               ), intent(inout), target      :: node
    double precision                                         , intent(in   )              :: time
    type            (multiCounter                           ), intent(inout), optional    :: instance
    class           (nodeComponentSatellite                 ), pointer                    :: satellite
    !$GLC attributes unused :: time, instance

    allocate(ramPressureHistoryExtract(3))
    if (node%isSatellite()) then
       satellite                 => node%satellite()
       ramPressureHistoryExtract =  [                                                                                           &
            &                        satellite                    %floatRank0MetaPropertyGet(self%ramPressureMaximumID       ), &
            &                        satellite                    %floatRank0MetaPropertyGet(self%ramPressureTimeIntegratedID), &
            &                        self%hotHaloRamPressureForce_%force                    (node                            )  &
            &                       ]
    else
       ramPressureHistoryExtract =  0.0d0
    end if
    return
  end function ramPressureHistoryExtract

  subroutine ramPressureHistoryNames(self,time,names)
    !!{RST
    Return the names of the ram pressure history properties.
    !!}
    implicit none
    class           (nodePropertyExtractorRamPressureHistory), intent(inout)                             :: self
    double precision                                         , intent(in   )                             :: time
    type            (varying_string                         ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(3))
    names(1)=var_str('satelliteRamPressureMaximum'       )
    names(2)=var_str('satelliteRamPressureTimeIntegrated')
    names(3)=var_str('satelliteRamPressure'              )
    return
  end subroutine ramPressureHistoryNames

  subroutine ramPressureHistoryDescriptions(self,time,descriptions)
    !!{RST
    Return the descriptions of the ram pressure history properties.
    !!}
    implicit none
    class           (nodePropertyExtractorRamPressureHistory), intent(inout)                             :: self
    double precision                                         , intent(in   )                             :: time
    type            (varying_string                         ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(3))
    descriptions(1)=var_str('The (decaying) running maximum of the ram pressure experienced by the satellite [M☉ (km/s)² Mpc⁻³].'  )
    descriptions(2)=var_str('The (decaying) time integral of the ram pressure experienced by the satellite [M☉ (km/s)² Mpc⁻³ Gyr].')
    descriptions(3)=var_str('The instantaneous ram pressure experienced by the satellite [M☉ (km/s)² Mpc⁻³].'                      )
    return
  end subroutine ramPressureHistoryDescriptions

  function ramPressureHistoryUnitsInSI(self,time)
    !!{RST
    Return the units of the ram pressure history properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear, massSolar, megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                         , dimension(:) , allocatable :: ramPressureHistoryUnitsInSI
    class           (nodePropertyExtractorRamPressureHistory), intent(inout)              :: self
    double precision                                         , intent(in   )              :: time
    double precision                                                                      :: unitsPressure
    !$GLC attributes unused :: self, time

    unitsPressure=massSolar*kilo**2/megaParsec**3
    allocate(ramPressureHistoryUnitsInSI(3))
    ramPressureHistoryUnitsInSI=[unitsPressure,unitsPressure*gigaYear,unitsPressure]
    return
  end function ramPressureHistoryUnitsInSI

  function ramPressureHistoryUnits(self,time) result(units)
    !!{RST
    Return the units of the ram pressure history properties.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type            (unitType                               ), dimension(:), allocatable :: units
    class           (nodePropertyExtractorRamPressureHistory), intent(inout)             :: self
    double precision                                         , intent(in   )             :: time
    double precision                                         , dimension(:), allocatable :: unitsInSI

    unitsInSI=self%unitsInSI(time)
    allocate(units(3))
    units(1)=unitType(unitsInSI(1),description='M☉ (km/s)² Mpc⁻³'    ,quantity='solMass*km^2/(s^2*Mpc^3)'    )
    units(2)=unitType(unitsInSI(2),description='M☉ (km/s)² Mpc⁻³ Gyr',quantity='solMass*km^2*Gyr/(s^2*Mpc^3)')
    units(3)=unitType(unitsInSI(3),description='M☉ (km/s)² Mpc⁻³'    ,quantity='solMass*km^2/(s^2*Mpc^3)'    )
    return
  end function ramPressureHistoryUnits
