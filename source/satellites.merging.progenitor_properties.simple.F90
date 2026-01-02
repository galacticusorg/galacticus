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

  !!{
  Implements a merger progenitor properties class which uses a simple calculation.
  !!}

  use :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass

  !![
  <mergerProgenitorProperties name="mergerProgenitorPropertiesSimple">
   <description>A merger progenitor properties class which uses a simple calculation.</description>
  </mergerProgenitorProperties>
  !!]
  type, extends(mergerProgenitorPropertiesClass) :: mergerProgenitorPropertiesSimple
     !!{
     A merger progenitor properties class which uses a simple calculation.
     !!}
     private
     class(mergerMassMovementsClass), pointer :: mergerMassMovements_ => null()
   contains
     final     ::        simpleDestructor
     procedure :: get => simpleGet
  end type mergerProgenitorPropertiesSimple

  interface mergerProgenitorPropertiesSimple
     !!{
     Constructors for the \refClass{mergerProgenitorPropertiesSimple} merger progenitor properties class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface mergerProgenitorPropertiesSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerProgenitorPropertiesSimple} merger progenitor properties class which takes a parameter list as input.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Error_Report        , Component_List
    use :: Galacticus_Nodes, only : defaultDiskComponent, defaultSpheroidComponent
    use :: Input_Parameters, only : inputParameter      , inputParameters
    implicit none
    type (mergerProgenitorPropertiesSimple)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(mergerMassMovementsClass        ), pointer       :: mergerMassMovements_

    ! Ensure that required methods are supported.
    if     (                                                                                                                                          &
         &  .not.                                                                                                                                     &
         &       (                                                                                                                                    &
         &        defaultDiskComponent    %    massStellarIsGettable().and.                                                                           &
         &        defaultDiskComponent    %        massGasIsGettable().and.                                                                           &
         &        defaultDiskComponent    % halfMassRadiusIsGettable()                                                                                &
         &  )                                                                                                                                         &
         & ) call Error_Report                                                                                                                        &
         &        (                                                                                                                                   &
         &         'this method requires that massStellar, massGas, and halfMassRadius properties must all be gettable for the disk component.'    // &
         &         Component_List(                                                                                                                    &
         &                        'disk'                                                                                                           ,  &
         &                        defaultDiskComponent    %    massStellarAttributeMatch(requireGettable=.true.).intersection.                        &
         &                        defaultDiskComponent    %        massGasAttributeMatch(requireGettable=.true.).intersection.                        &
         &                        defaultDiskComponent    % halfMassRadiusAttributeMatch(requireGettable=.true.)                                      &
         &                       )                                                                                                                 // &
         &         {introspection:location}                                                                                                           &
         &        )
    if     (                                                                                                                                          &
         &  .not.                                                                                                                                     &
         &       (                                                                                                                                    &
         &        defaultSpheroidComponent%    massStellarIsGettable().and.                                                                           &
         &        defaultSpheroidComponent%        massGasIsGettable().and.                                                                           &
         &        defaultSpheroidComponent% halfMassRadiusIsGettable()                                                                                &
         &  )                                                                                                                                         &
         & ) call Error_Report                                                                                                                        &
         &        (                                                                                                                                   &
         &         'this method requires that massStellar, massGas, and halfMassRadius properties must all be gettable for the spheroid component.'// &
         &         Component_List(                                                                                                                    &
         &                        'spheroid'                                                                                                       ,  &
         &                        defaultSpheroidComponent%    massStellarAttributeMatch(requireGettable=.true.).intersection.                        &
         &                        defaultSpheroidComponent%        massGasAttributeMatch(requireGettable=.true.).intersection.                        &
         &                        defaultSpheroidComponent% halfMassRadiusAttributeMatch(requireGettable=.true.)                                      &
         &                       )                                                                                                                 // &
         &         {introspection:location}                                                                                                           &
         &        )
    !![
    <objectBuilder class="mergerMassMovements" name="mergerMassMovements_" source="parameters"/>
    !!]
    self=mergerProgenitorPropertiesSimple(mergerMassMovements_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerMassMovements_"/>
    !!]
    return
  end function simpleConstructorParameters

 function simpleConstructorInternal(mergerMassMovements_) result(self)
    !!{
    Internal constructor for the \refClass{mergerProgenitorPropertiesSimple} merger progenitor properties class.
    !!}
    implicit none
    type (mergerProgenitorPropertiesSimple)                        :: self
    class(mergerMassMovementsClass        ), intent(in   ), target :: mergerMassMovements_
    !![
    <constructorAssign variables="*mergerMassMovements_"/>
    !!]

    return
  end function simpleConstructorInternal

  subroutine simpleDestructor(self)
    !!{
    Destructor for the \refClass{mergerProgenitorPropertiesSimple} merger progenitor properties class.
    !!}
    implicit none
    type(mergerProgenitorPropertiesSimple), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerMassMovements_"/>
    !!]
    return
  end subroutine simpleDestructor

  subroutine simpleGet(self,nodeSatellite,nodeHost,massSatellite,massHost,massSpheroidSatellite,massSpheroidHost,massSpheroidHostPreMerger,radiusSatellite,radiusHost,factorAngularMomentum,massSpheroidRemnant,massGasSpheroidRemnant)
    !!{
    Computes various properties of the progenitor galaxies useful for calculations of merger remnant sizes.
    !!}
    use :: Galactic_Structure_Options      , only : massTypeGalactic
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentDisk    , nodeComponentSpheroid    , treeNode
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerDisk, destinationMergerSpheroid, destinationMergerUnmoved, enumerationDestinationMergerType
    use :: Mass_Distributions              , only : massDistributionClass
    implicit none
    class           (mergerProgenitorPropertiesSimple), intent(inout), target :: self
    type            (treeNode                        ), intent(inout), target :: nodeSatellite            , nodeHost
    double precision                                  , intent(  out)         :: factorAngularMomentum    , massHost                 , &
         &                                                                       radiusHost               , massSpheroidHost         , &
         &                                                                       massSpheroidHostPreMerger, massGasSpheroidRemnant   , &
         &                                                                       massSpheroidRemnant      , massSatellite            , &
         &                                                                       radiusSatellite          , massSpheroidSatellite
    class           (massDistributionClass           ), pointer               :: massDistributionHost     , massDistributionSatellite
    class           (nodeComponentDisk               ), pointer               :: diskHost                 , diskSatelite
    class           (nodeComponentSpheroid           ), pointer               :: spheroidHost             , spheroidSatellite
    type            (enumerationDestinationMergerType)                        :: destinationGasSatellite  , destinationGasHost       , &
         &                                                                       destinationStarsHost     , destinationStarsSatellite
    logical                                                                   :: mergerIsMajor

    ! Find how mass is moved by the merger.
    call self%mergerMassMovements_%get(nodeSatellite,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    ! Get the disk and spheroid components of host and satellite.
    diskHost          => nodeHost     %disk    ()
    spheroidHost      => nodeHost     %spheroid()
    diskSatelite      => nodeSatellite%disk    ()
    spheroidSatellite => nodeSatellite%spheroid()
    ! Find the baryonic masses of the two galaxies.
    massDistributionSatellite => nodeSatellite            %massDistribution(massType=massTypeGalactic)
    massDistributionHost      => nodeHost                 %massDistribution(massType=massTypeGalactic)
    massSatellite             =  massDistributionSatellite%massTotal       (                         )
    massHost                  =  massDistributionHost     %massTotal       (                         )
    ! Find the masses of material that will end up in the spheroid component of the remnant.
    select case (destinationGasHost%ID)
    case (destinationMergerSpheroid%ID)
       massSpheroidHost      =+spheroidHost%massGas()+diskHost    %massGas       ()
       massGasSpheroidRemnant=+spheroidHost%massGas()+diskHost    %massGas       ()
       massSpheroidRemnant   =+spheroidHost%massGas()+diskHost    %massGas       ()
       radiusHost            =+spheroidHost%massGas()*spheroidHost%halfMassRadius() &
            &                 +diskHost    %massGas()*diskHost    %halfMassRadius()
    case (destinationMergerDisk   %ID)
       massSpheroidHost      =+0.0d0
       massGasSpheroidRemnant=+0.0d0
       massSpheroidRemnant   =+0.0d0
    case (destinationMergerUnmoved%ID)
       massSpheroidHost      =+spheroidHost%massGas()
       massGasSpheroidRemnant=+spheroidHost%massGas()
       massSpheroidRemnant   =+spheroidHost%massGas()
       radiusHost            =+spheroidHost%massGas()*spheroidHost%halfMassRadius()
    case default
       call Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    select case (destinationStarsHost%ID)
    case (destinationMergerSpheroid%ID)
       massSpheroidHost     =+massSpheroidHost       +spheroidHost%massStellar()+    diskHost%massStellar   ()
       massSpheroidRemnant  =+massSpheroidRemnant    +spheroidHost%massStellar()+    diskHost%massStellar   ()
       radiusHost           =+radiusHost             +spheroidHost%massStellar()*spheroidHost%halfMassRadius() &
            &                                        +    diskHost%massStellar()*    diskHost%halfMassRadius()
    case (destinationMergerDisk   %ID)
       massSpheroidHost     =+massSpheroidHost
    case (destinationMergerUnmoved%ID)
       massSpheroidHost     =+massSpheroidHost       +spheroidHost%massStellar()
       massSpheroidRemnant  =+massSpheroidRemnant    +spheroidHost%massStellar()
       radiusHost           =+radiusHost             +spheroidHost%massStellar()*spheroidHost%halfMassRadius()
    case default
       call Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    select case (destinationGasSatellite%ID)
    case (destinationMergerSpheroid%ID)
       massSpheroidSatellite =                       +spheroidSatellite%massGas()+     diskSatelite%massGas       ()
       massGasSpheroidRemnant=+massGasSpheroidRemnant+spheroidSatellite%massGas()+     diskSatelite%massGas       ()
       massSpheroidRemnant   =+massSpheroidRemnant   +spheroidSatellite%massGas()+     diskSatelite%massGas       ()
       radiusSatellite       =                       +spheroidSatellite%massGas()*spheroidSatellite%halfMassRadius() &
            &                                        +     diskSatelite%massGas()*     diskSatelite%halfMassRadius()
    case (destinationMergerDisk   %ID)
       massSpheroidSatellite =+0.0d0
    case (destinationMergerUnmoved%ID)
       massSpheroidSatellite =                       +spheroidSatellite%massGas()
       massGasSpheroidRemnant=+massGasSpheroidRemnant+spheroidSatellite%massGas()
       massSpheroidRemnant   =+massSpheroidRemnant   +spheroidSatellite%massGas()
       radiusSatellite       =                       +spheroidSatellite%massGas()*spheroidSatellite%halfMassRadius()
    case default
       call Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    select case (destinationStarsSatellite%ID)
    case (destinationMergerSpheroid%ID)
       massSpheroidSatellite=+massSpheroidSatellite+spheroidSatellite%massStellar()+     diskSatelite%massStellar   ()
       massSpheroidRemnant  =+massSpheroidRemnant  +spheroidSatellite%massStellar()+     diskSatelite%massStellar   ()
       radiusSatellite      =+radiusSatellite      +spheroidSatellite%massStellar()*spheroidSatellite%halfMassRadius() &
            &                                      +     diskSatelite%massStellar()*     diskSatelite%halfMassRadius()
    case (destinationMergerDisk   %ID)
       massSpheroidSatellite=+massSpheroidSatellite
    case (destinationMergerUnmoved%ID)
       massSpheroidSatellite=+massSpheroidSatellite+spheroidSatellite%massStellar()
       massSpheroidRemnant  =+massSpheroidRemnant  +spheroidSatellite%massStellar()
       radiusSatellite      =+radiusSatellite      +spheroidSatellite%massStellar()*spheroidSatellite%halfMassRadius()
    case default
       call Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    ! Compute the half-mass radii of the material that will end up in the remnant spheroid.
    if (massSpheroidHost      > 0.0d0) radiusHost     =radiusHost     /massSpheroidHost
    if (massSpheroidSatellite > 0.0d0) radiusSatellite=radiusSatellite/massSpheroidSatellite
    ! Compute the angular momentum factor.
    factorAngularMomentum=1.0d0
    ! Compute the mass of the host spheroid before the merger.
    massSpheroidHostPreMerger=spheroidHost%massStellar()+spheroidHost%massGas()
    ! Clean up.
    !![
    <objectDestructor name="massDistributionSatellite"/>
    <objectDestructor name="massDistributionHost     "/>
    !!]
    return
  end subroutine simpleGet
