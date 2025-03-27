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
Implements an ISM mass output analysis property extractor class.
!!}

  !![
  <enumeration>
   <name>satelliteStatusDiscriminator</name>
   <description>Enumeration of possible discriminators for satellite orphan status.</description>
   <visibility>private</visibility>
   <encodeFunction>yes</encodeFunction>
   <entry label="boundMass"/>
   <entry label="position" />
  </enumeration>
  !!]

  !![
  <nodePropertyExtractor name="nodePropertyExtractorSatelliteStatus">
   <description>An ISM mass output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorSatelliteStatus
     !!{
     A stellar mass output analysis class.
     !!}
     private
     type(enumerationSatelliteStatusDiscriminatorType) :: discriminator
   contains
     procedure :: extract     => satelliteStatusExtract
     procedure :: name        => satelliteStatusName
     procedure :: description => satelliteStatusDescription
  end type nodePropertyExtractorSatelliteStatus

  interface nodePropertyExtractorSatelliteStatus
     !!{
     Constructors for the {\normalfont \ttfamily satelliteStatus} output analysis class.
     !!}
     module procedure satelliteStatusConstructorParameters
     module procedure satelliteStatusConstructorInternal
  end interface nodePropertyExtractorSatelliteStatus

contains

  function satelliteStatusConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily satelliteStatus} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorSatelliteStatus)                :: self
    type(inputParameters                     ), intent(inout) :: parameters
    type(varying_string                      )                :: discriminator

    !![
    <inputParameter>
      <name>discriminator</name>
      <defaultValue>var_str('boundMass')</defaultValue>
      <description>Specifies whether bound mass or position history will be used to determine satellite orphan status.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodePropertyExtractorSatelliteStatus(enumerationSatelliteStatusDiscriminatorEncode(char(discriminator),includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function satelliteStatusConstructorParameters

  function satelliteStatusConstructorInternal(discriminator) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily satelliteStatus} node property extractor class.
    !!}
    use :: Error           , only : Component_List          , Error_Report
    use :: Galacticus_Nodes, only : defaultPositionComponent, defaultSatelliteComponent
    implicit none
    type(nodePropertyExtractorSatelliteStatus       )                :: self
    type(enumerationSatelliteStatusDiscriminatorType), intent(in   ) :: discriminator
    !![
    <constructorAssign variables="discriminator"/>
    !!]

    select case (discriminator%ID)
    case (satelliteStatusDiscriminatorBoundMass%ID)
       if     (                                                                                                             &
            &  .not.                                                                                                        &
            &       (                                                                                                       &
            &        defaultSatelliteComponent% boundMassHistoryIsGettable()                                                &
            &       )                                                                                                       &
            & ) call Error_Report                                                                                           &
            &        (                                                                                                      &
            &         'this method requires that boundMassHistory property must be gettable for the satellite component.'// &
            &         Component_List(                                                                                       &
            &                        'satellite'                                                                         ,  &
            &                        defaultSatelliteComponent%boundMassHistoryAttributeMatch(requireGettable=.true.)       &
            &                       )                                                                                    // &
            &         {introspection:location}                                                                              &
            &        )
    case (satelliteStatusDiscriminatorPosition %ID)
       if     (                                                                                                             &
            &  .not.                                                                                                        &
            &       (                                                                                                       &
            &        defaultPositionComponent %  positionHistoryIsGettable()                                                &
            &       )                                                                                                       &
            & ) call Error_Report                                                                                           &
            &        (                                                                                                      &
            &         'this method requires that positionHistory property must be gettable for the position component.'  // &
            &         Component_List(                                                                                       &
            &                        'position'                                                                          ,  &
            &                        defaultPositionComponent%positionHistoryAttributeMatch  (requireGettable=.true.)       &
            &                       )                                                                                    // &
            &         {introspection:location}                                                                              &
            &        )
    end select
    return
  end function satelliteStatusConstructorInternal

  function satelliteStatusExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily satelliteStatus} node property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentPosition, nodeComponentSatellite, treeNode
    use :: Histories       , only : history
    implicit none
    integer         (kind_int8                           )                          :: satelliteStatusExtract
    class           (nodePropertyExtractorSatelliteStatus), intent(inout)           :: self
    type            (treeNode                            ), intent(inout), target   :: node
    double precision                                      , intent(in   )           :: time
    type            (multiCounter                        ), intent(inout), optional :: instance
    class           (nodeComponentBasic                  ), pointer                 :: basic
    class           (nodeComponentSatellite              ), pointer                 :: satellite
    class           (nodeComponentPosition               ), pointer                 :: position
    type            (history                             )                          :: discriminatorHistory
    !$GLC attributes unused :: instance, time

    if (node%isSatellite()) then
       basic => node%basic()
       select case (self%discriminator%ID)
       case (satelliteStatusDiscriminatorBoundMass%ID)
          satellite            => node     %satellite       ()
          discriminatorHistory =  satellite%boundMassHistory()
       case (satelliteStatusDiscriminatorPosition %ID)
          position             => node     %position        ()
          discriminatorHistory =  position %positionHistory ()
       end select
       if (discriminatorHistory%exists()) then
          if (basic%time() > discriminatorHistory%time(size(discriminatorHistory%time))) then
             ! Orphaned satellite.
             satelliteStatusExtract=2
          else
             ! Non-orphaned satellite.
             satelliteStatusExtract=1
          end if
       else
          ! Orphaned satellite.
          satelliteStatusExtract=2
       end if
    else
       ! Not a satellite.
       satelliteStatusExtract=0
    end if
    return
  end function satelliteStatusExtract


  function satelliteStatusName(self)
    !!{
    Return the name of the satelliteStatus property.
    !!}
    implicit none
    type (varying_string                      )                :: satelliteStatusName
    class(nodePropertyExtractorSatelliteStatus), intent(inout) :: self
    !$GLC attributes unused :: self

    satelliteStatusName=var_str('satelliteStatus')
    return
  end function satelliteStatusName

  function satelliteStatusDescription(self)
    !!{
    Return a description of the satelliteStatus property.
    !!}
    implicit none
    type (varying_string                      )                :: satelliteStatusDescription
    class(nodePropertyExtractorSatelliteStatus), intent(inout) :: self
    !$GLC attributes unused :: self

    satelliteStatusDescription=var_str('Satellite status flag (0=not a satellite; 1=satellite with halo; 2=orphaned satellite).')
    return
  end function satelliteStatusDescription


