!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements satellite orbital extrema property extractor class.

  !# <nodePropertyExtractor name="nodePropertyExtractorSatelliteOrbitalExtrema">
  !#  <description>A satellite orbital extrema luminosity property extractor class.</description>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorSatelliteOrbitalExtrema
     !% A satellite orbital extrema property extractor class.
     private
     integer :: offsetPericenter , offsetApocenter , &
          &     elementCount_
     logical :: extractPericenter, extractApocenter
   contains
     procedure :: elementCount => satelliteOrbitalExtremaElementCount
     procedure :: extract      => satelliteOrbitalExtremaExtract
     procedure :: names        => satelliteOrbitalExtremaNames
     procedure :: descriptions => satelliteOrbitalExtremaDescriptions
     procedure :: unitsInSI    => satelliteOrbitalExtremaUnitsInSI
     procedure :: type         => satelliteOrbitalExtremaType
  end type nodePropertyExtractorSatelliteOrbitalExtrema

  interface nodePropertyExtractorSatelliteOrbitalExtrema
     !% Constructors for the ``satelliteOrbitalExtrema'' output analysis class.
     module procedure satelliteOrbitalExtremaConstructorParameters
     module procedure satelliteOrbitalExtremaConstructorInternal
  end interface nodePropertyExtractorSatelliteOrbitalExtrema

contains

  function satelliteOrbitalExtremaConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily satelliteOrbitalExtrema} property extractor class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorSatelliteOrbitalExtrema)                :: self
    type   (inputParameters                             ), intent(inout) :: parameters
    logical                                                              :: extractPericenter, extractApocenter

    !# <inputParameter>
    !#   <name>extractPericenter</name>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Specifies whether or not satellite orbital pericenter data (radius, velocity) should be extracted.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>extractApocenter</name>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Specifies whether or not satellite orbital apocenter data (radius, velocity) should be extracted.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    self=nodePropertyExtractorSatelliteOrbitalExtrema(extractPericenter,extractApocenter)
    !# <inputParametersValidate source="parameters"/>
    return
  end function satelliteOrbitalExtremaConstructorParameters

  function satelliteOrbitalExtremaConstructorInternal(extractPericenter, extractApocenter) result(self)
    !% Internal constructor for the {\normalfont \ttfamily satelliteOrbitalExtrema} property extractor class.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type (nodePropertyExtractorSatelliteOrbitalExtrema)                :: self
    logical                                            , intent(in   ) :: extractPericenter, extractApocenter
    !# <constructorAssign variables="extractPericenter, extractApocenter"/>

    self%elementCount_=0
    if (extractPericenter) then
       self%offsetPericenter=self%elementCount_+1
       self%elementCount_   =self%elementCount_+2
    end if
    if (extractApocenter ) then
       self%offsetApocenter =self%elementCount_+1
       self%elementCount_   =self%elementCount_+2
    end if
    if (self%elementCount_ == 0) call Galacticus_Error_Report('no properties selected for extraction'//{introspection:location})
    return
  end function satelliteOrbitalExtremaConstructorInternal

  integer function satelliteOrbitalExtremaElementCount(self,time)
    !% Return the number of elements in the {\normalfont \ttfamily satelliteOrbitalExtrema} property extractor class.
    implicit none
    class           (nodePropertyExtractorSatelliteOrbitalExtrema), intent(inout) :: self
    double precision                                              , intent(in   ) :: time
    !$GLC attributes unused :: time

    satelliteOrbitalExtremaElementCount=self%elementCount_
    return
  end function satelliteOrbitalExtremaElementCount

  function satelliteOrbitalExtremaExtract(self,node,time,instance)
    !% Implement a {\normalfont \ttfamily satelliteOrbitalExtrema} property extractor
    use :: Galacticus_Nodes, only : nodeComponentSatellite                          , treeNode
    use :: Kepler_Orbits   , only : keplerOrbit
    use :: Satellite_Orbits, only : Satellite_Orbit_Extremum_Phase_Space_Coordinates, extremumApocenter, extremumPericenter
    implicit none
    double precision                                              , dimension(:) , allocatable :: satelliteOrbitalExtremaExtract
    class           (nodePropertyExtractorSatelliteOrbitalExtrema), intent(inout), target      :: self
    type            (treeNode                                    ), intent(inout), target      :: node
    double precision                                              , intent(in   )              :: time
    type            (multiCounter                                ), intent(inout), optional    :: instance
    type            (treeNode                                    ), pointer                    :: nodeHost
    class           (nodeComponentSatellite                      ), pointer                    :: satellite
    type            (keplerOrbit                                 )                             :: orbit
    double precision                                                                           :: radiusOrbital                 , velocityOrbital
    !$GLC attributes unused :: time, instance

    allocate(satelliteOrbitalExtremaExtract(self%elementCount_))
    if (self%extractPericenter) then
       if (node%isSatellite()) then
          nodeHost  => node     %parent
          satellite => node     %satellite  ()
          orbit     =  satellite%virialOrbit()
          call Satellite_Orbit_Extremum_Phase_Space_Coordinates(nodeHost,orbit,extremumPericenter,radiusOrbital,velocityOrbital)
       else
          radiusOrbital  =0.0d0
          velocityOrbital=0.0d0
       end if
       satelliteOrbitalExtremaExtract(self%offsetPericenter:self%offsetPericenter+1)=[radiusOrbital,velocityOrbital]
    end if
    if (self%extractApocenter) then
       if (node%isSatellite()) then
          nodeHost  => node     %parent
          satellite => node     %satellite  ()
          orbit     =  satellite%virialOrbit()
          call Satellite_Orbit_Extremum_Phase_Space_Coordinates(nodeHost,orbit,extremumApocenter ,radiusOrbital,velocityOrbital)
       else
          radiusOrbital  =0.0d0
          velocityOrbital=0.0d0
       end if
       satelliteOrbitalExtremaExtract(self%offsetApocenter :self%offsetApocenter+1 )=[radiusOrbital,velocityOrbital]
    end if
    return
  end function satelliteOrbitalExtremaExtract

  function satelliteOrbitalExtremaNames(self,time)
    !% Return the name of the {\normalfont \ttfamily satelliteOrbitalExtrema} property.
    implicit none
    type            (varying_string                              ), dimension(:) , allocatable :: satelliteOrbitalExtremaNames
    class           (nodePropertyExtractorSatelliteOrbitalExtrema), intent(inout)              :: self
    double precision                                              , intent(in   )              :: time
    !$GLC attributes unused :: time

    allocate(satelliteOrbitalExtremaNames(self%elementCount_))
    if (self%extractPericenter) satelliteOrbitalExtremaNames(self%offsetPericenter:self%offsetPericenter+1)=[var_str('satellitePericenterRadius'),var_str('satellitePericenterVelocity')]
    if (self%extractApocenter ) satelliteOrbitalExtremaNames(self%offsetApocenter :self%offsetApocenter +1)=[var_str('satelliteApocenterRadius'),var_str('satelliteApocenterVelocity')]
    return
  end function satelliteOrbitalExtremaNames

  function satelliteOrbitalExtremaDescriptions(self,time)
    !% Return a description of the satelliteOrbitalExtrema property.
    implicit none
    type            (varying_string                              ), dimension(:) , allocatable :: satelliteOrbitalExtremaDescriptions
    class           (nodePropertyExtractorSatelliteOrbitalExtrema), intent(inout)              :: self
    double precision                                              , intent(in   )              :: time
    !$GLC attributes unused :: time

    allocate(satelliteOrbitalExtremaDescriptions(self%elementCount_))
    if (self%extractPericenter) satelliteOrbitalExtremaDescriptions(self%offsetPericenter:self%offsetPericenter+1)=[var_str('Pericenteric radius of satellite orbit [Mpc].'),var_str('Pericenteric velocity of satellite orbit [km/s].')]
    if (self%extractApocenter ) satelliteOrbitalExtremaDescriptions(self%offsetApocenter :self%offsetApocenter +1)=[var_str('Apocenteric radius of satellite orbit [Mpc].' ),var_str('Apocenteric velocity of satellite orbit [km/s].' )]
    return
  end function satelliteOrbitalExtremaDescriptions

  function satelliteOrbitalExtremaUnitsInSI(self,time)
    !% Return the units of the {\normalfont \ttfamily satelliteOrbitalExtrema} property in the SI system.
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                              , allocatable  , dimension(:) :: satelliteOrbitalExtremaUnitsInSI
    class           (nodePropertyExtractorSatelliteOrbitalExtrema), intent(inout)               :: self
    double precision                                              , intent(in   )               :: time
    !$GLC attributes unused :: time

    allocate(satelliteOrbitalExtremaUnitsInSI(self%elementCount_))
    if (self%extractPericenter) satelliteOrbitalExtremaUnitsInSI(self%offsetPericenter:self%offsetPericenter+1)=[megaParsec,kilo]
    if (self%extractApocenter ) satelliteOrbitalExtremaUnitsInSI(self%offsetApocenter :self%offsetApocenter +1)=[megaParsec,kilo]
    return
  end function satelliteOrbitalExtremaUnitsInSI

  integer function satelliteOrbitalExtremaType(self)
    !% Return the type of the {\normalfont \ttfamily satelliteOrbitalExtrema} property.
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorSatelliteOrbitalExtrema), intent(inout) :: self
    !$GLC attributes unused :: self

    satelliteOrbitalExtremaType=outputAnalysisPropertyTypeLinear
    return
  end function satelliteOrbitalExtremaType
