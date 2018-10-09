!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% Implements a merger progenitor properties class which uses a simple calculation.

  !# <mergerProgenitorProperties name="mergerProgenitorPropertiesSimple">
  !#  <description>A merger progenitor properties class which uses a simple calculation.</description>
  !# </mergerProgenitorProperties>
  type, extends(mergerProgenitorPropertiesClass) :: mergerProgenitorPropertiesSimple
     !% A merger progenitor properties class which uses a simple calculation.
     private
   contains
     procedure :: get => simpleGet
  end type mergerProgenitorPropertiesSimple

  interface mergerProgenitorPropertiesSimple
     !% Constructors for the {\normalfont \ttfamily simple} merger progenitor properties class.
     module procedure simpleConstructorParameters
  end interface mergerProgenitorPropertiesSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily simple} merger progenitor properties class which takes a parameter list as input.
    use Galacticus_Error
    use Array_Utilities
    use Input_Parameters
    implicit none
    type(mergerProgenitorPropertiesSimple)                :: self
    type(inputParameters                 ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters
    
    ! Ensure that required methods are supported.
    if     (                                                                                                                                          &
         &  .not.                                                                                                                                     &
         &       (                                                                                                                                    &
         &        defaultDiskComponent    %    massStellarIsGettable().and.                                                                           &
         &        defaultDiskComponent    %        massGasIsGettable().and.                                                                           &
         &        defaultDiskComponent    % halfMassRadiusIsGettable()                                                                                &
         &  )                                                                                                                                         &
         & ) call Galacticus_Error_Report                                                                                                             &
         &        (                                                                                                                                   &
         &         'this method requires that massStellar, massGas, and halfMassRadius properties must all be gettable for the disk component.'    // &
         &         Galacticus_Component_List(                                                                                                         &
         &                                   'disk'                                                                                                ,  &
         &                                   defaultDiskComponent    %    massStellarAttributeMatch(requireGettable=.true.).intersection.             &
         &                                   defaultDiskComponent    %        massGasAttributeMatch(requireGettable=.true.).intersection.             &
         &                                   defaultDiskComponent    % halfMassRadiusAttributeMatch(requireGettable=.true.)                           &
         &                                  )                                                                                                      // &
         &         {introspection:location}                                                                                                           &
         &        )
    if     (                                                                                                                                          &
         &  .not.                                                                                                                                     &
         &       (                                                                                                                                    &
         &        defaultSpheroidComponent%    massStellarIsGettable().and.                                                                           &
         &        defaultSpheroidComponent%        massGasIsGettable().and.                                                                           &
         &        defaultSpheroidComponent% halfMassRadiusIsGettable()                                                                                &
         &  )                                                                                                                                         &
         & ) call Galacticus_Error_Report                                                                                                             &
         &        (                                                                                                                                   &
         &         'this method requires that massStellar, massGas, and halfMassRadius properties must all be gettable for the spheroid component.'// &
         &         Galacticus_Component_List(                                                                                                         &
         &                                   'spheroid'                                                                                            ,  &
         &                                   defaultSpheroidComponent%    massStellarAttributeMatch(requireGettable=.true.).intersection.             &
         &                                   defaultSpheroidComponent%        massGasAttributeMatch(requireGettable=.true.).intersection.             &
         &                                   defaultSpheroidComponent% halfMassRadiusAttributeMatch(requireGettable=.true.)                           &
         &                                  )                                                                                                      // &
         &         {introspection:location}                                                                                                           &
         &        )
    self=mergerProgenitorPropertiesSimple()
    return
  end function simpleConstructorParameters

  subroutine simpleGet(self,nodeSatellite,nodeHost,massSatellite,massHost,massSpheroidSatellite,massSpheroidHost,massSpheroidHostPreMerger,radiusSatellite,radiusHost,factorAngularMomentum,massSpheroidRemnant,massGasSpheroidRemnant)
    !% Computes various properties of the progenitor galaxies useful for calculations of merger remnant sizes.
    use Galactic_Structure_Radii
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Satellite_Merging_Mass_Movements_Descriptors
    use Numerical_Constants_Physical
    use Galacticus_Error
    implicit none
    class           (mergerProgenitorPropertiesSimple), intent(inout)         :: self
    type            (treeNode                        ), intent(inout), target :: nodeSatellite            , nodeHost
    double precision                                  , intent(  out)         :: factorAngularMomentum    , massHost              , &
         &                                                                       radiusHost               , massSpheroidHost      , &
         &                                                                       massSpheroidHostPreMerger, massGasSpheroidRemnant, &
         &                                                                       massSpheroidRemnant      , massSatellite         , &
         &                                                                       radiusSatellite          , massSpheroidSatellite
    class           (nodeComponentDisk               ), pointer               :: diskHost                 , diskSatelite
    class           (nodeComponentSpheroid           ), pointer               :: spheroidHost             , spheroidSatellite
    !GCC$ attributes unused :: self
    
    ! Get the disk and spheroid components of host and satellite.
    diskHost          => nodeHost     %disk    ()
    spheroidHost      => nodeHost     %spheroid()
    diskSatelite      => nodeSatellite%disk    ()
    spheroidSatellite => nodeSatellite%spheroid()
    ! Solve for the radii of the host and satellite nodes, to ensure they are computed and up to date.
    call Galactic_Structure_Radii_Solve(nodeHost     )
    call Galactic_Structure_Radii_Solve(nodeSatellite)
    ! Find the baryonic masses of the two galaxies.
    massSatellite=Galactic_Structure_Enclosed_Mass(nodeSatellite,massType=massTypeGalactic)
    massHost     =Galactic_Structure_Enclosed_Mass(nodeHost     ,massType=massTypeGalactic)
    ! Find the masses of material that will end up in the spheroid component of the remnant.
    select case (thisHostGasMovesTo)
    case (movesToSpheroid)
       massSpheroidHost      =+spheroidHost%massGas()+diskHost    %massGas       ()
       massGasSpheroidRemnant=+spheroidHost%massGas()+diskHost    %massGas       ()
       massSpheroidRemnant   =+spheroidHost%massGas()+diskHost    %massGas       ()
       radiusHost            =+spheroidHost%massGas()*spheroidHost%halfMassRadius() &
            &                 +diskHost    %massGas()*diskHost    %halfMassRadius()
    case (movesToDisk)
       massSpheroidHost      =+0.0d0
       massGasSpheroidRemnant=+0.0d0
       massSpheroidRemnant   =+0.0d0
    case (doesNotMove)
       massSpheroidHost      =+spheroidHost%massGas()
       massGasSpheroidRemnant=+spheroidHost%massGas()
       massSpheroidRemnant   =+spheroidHost%massGas()
       radiusHost            =+spheroidHost%massGas()*spheroidHost%halfMassRadius()
    case default
       call Galacticus_Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    select case (thisHostStarsMoveTo)
    case (movesToSpheroid)
       massSpheroidHost     =+massSpheroidHost       +spheroidHost%massStellar()+    diskHost%massStellar   ()
       massSpheroidRemnant  =+massSpheroidRemnant    +spheroidHost%massStellar()+    diskHost%massStellar   ()
       radiusHost           =+radiusHost             +spheroidHost%massStellar()*spheroidHost%halfMassRadius() &
            &                                        +    diskHost%massStellar()*    diskHost%halfMassRadius()
    case (movesToDisk)
       massSpheroidHost     =+massSpheroidHost
    case (doesNotMove)
       massSpheroidHost     =+massSpheroidHost       +spheroidHost%massStellar()
       massSpheroidRemnant  =+massSpheroidRemnant    +spheroidHost%massStellar()
       radiusHost           =+radiusHost             +spheroidHost%massStellar()*spheroidHost%halfMassRadius()
    case default
       call Galacticus_Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    select case (thisMergerGasMovesTo)
    case (movesToSpheroid)
       massSpheroidSatellite =                       +spheroidSatellite%massGas()+     diskSatelite%massGas       ()
       massGasSpheroidRemnant=+massGasSpheroidRemnant+spheroidSatellite%massGas()+     diskSatelite%massGas       ()
       massSpheroidRemnant   =+massSpheroidRemnant   +spheroidSatellite%massGas()+     diskSatelite%massGas       ()
       radiusSatellite       =                       +spheroidSatellite%massGas()*spheroidSatellite%halfMassRadius() &
            &                                        +     diskSatelite%massGas()*     diskSatelite%halfMassRadius()
    case (movesToDisk)
       massSpheroidSatellite =+0.0d0
    case (doesNotMove)
       massSpheroidSatellite =                       +spheroidSatellite%massGas()
       massGasSpheroidRemnant=+massGasSpheroidRemnant+spheroidSatellite%massGas()
       massSpheroidRemnant   =+massSpheroidRemnant   +spheroidSatellite%massGas()
       radiusSatellite       =                       +spheroidSatellite%massGas()*spheroidSatellite%halfMassRadius()
    case default
       call Galacticus_Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    select case (thisMergerStarsMoveTo)
    case (movesToSpheroid)
       massSpheroidSatellite=+massSpheroidSatellite+spheroidSatellite%massStellar()+     diskSatelite%massStellar   ()
       massSpheroidRemnant  =+massSpheroidRemnant  +spheroidSatellite%massStellar()+     diskSatelite%massStellar   ()
       radiusSatellite      =+radiusSatellite      +spheroidSatellite%massStellar()*spheroidSatellite%halfMassRadius() &
            &                                      +     diskSatelite%massStellar()*     diskSatelite%halfMassRadius()
    case (movesToDisk)
       massSpheroidSatellite=+massSpheroidSatellite
    case (doesNotMove)
       massSpheroidSatellite=+massSpheroidSatellite+spheroidSatellite%massStellar()
       massSpheroidRemnant  =+massSpheroidRemnant  +spheroidSatellite%massStellar()
       radiusSatellite      =+radiusSatellite      +spheroidSatellite%massStellar()*spheroidSatellite%halfMassRadius()
    case default
       call Galacticus_Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    ! Compute the half-mass radii of the material that will end up in the remnant spheroid.
    if (massSpheroidHost      > 0.0d0) radiusHost     =radiusHost     /massSpheroidHost
    if (massSpheroidSatellite > 0.0d0) radiusSatellite=radiusSatellite/massSpheroidSatellite
    ! Compute the angular momentum factor.
    factorAngularMomentum=1.0d0
    ! Compute the mass of the host spheroid before the merger.
    massSpheroidHostPreMerger=spheroidHost%massStellar()+spheroidHost%massGas()
    return
  end subroutine simpleGet
