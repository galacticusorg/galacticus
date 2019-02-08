!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implements a merger progenitor properties class which uses the algorithm of \cite{cole_hierarchical_2000}.

  use Satellite_Merging_Mass_Movements

  !# <mergerProgenitorProperties name="mergerProgenitorPropertiesCole2000">
  !#  <description>A merger progenitor properties class which uses the algorithm of \cite{cole_hierarchical_2000}.</description>
  !# </mergerProgenitorProperties>
  type, extends(mergerProgenitorPropertiesClass) :: mergerProgenitorPropertiesCole2000
     !% A merger progenitor properties class which uses the algorithm of \cite{cole_hierarchical_2000}.
     private
     class(mergerMassMovementsClass), pointer :: mergerMassMovements_ => null()
   contains
     final     ::        cole2000Destructor
     procedure :: get => cole2000Get
  end type mergerProgenitorPropertiesCole2000

  interface mergerProgenitorPropertiesCole2000
     !% Constructors for the {\normalfont \ttfamily cole2000} merger progenitor properties class.
     module procedure cole2000ConstructorParameters
     module procedure cole2000ConstructorInternal
  end interface mergerProgenitorPropertiesCole2000

  ! Module global variables used in root finding.
  type            (treeNode  ), pointer :: cole2000Node
  integer                               :: cole2000DestinationGas, cole2000DestinationStars
  double precision                      :: cole2000HalfMass
  !$omp threadprivate(cole2000Node,cole2000DestinationGas,cole2000DestinationStars,cole2000HalfMass)

contains

  function cole2000ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily cole2000} merger progenitor properties class which takes a parameter list as input.
    use Galacticus_Nodes , only : defaultDiskComponent, defaultSpheroidComponent
    use Galacticus_Error
    use Array_Utilities
    use Input_Parameters
    implicit none
    type (mergerProgenitorPropertiesCole2000)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(mergerMassMovementsClass          ), pointer       :: mergerMassMovements_
    
    ! Ensure that required methods are supported.
    if     (                                                                                                                                                           &
         &  .not.                                                                                                                                                      &
         &       (                                                                                                                                                     &
         &        defaultDiskComponent    %    massStellarIsGettable().and.                                                                                            &
         &        defaultDiskComponent    %        massGasIsGettable().and.                                                                                            &
         &        defaultDiskComponent    % halfMassRadiusIsGettable().and.                                                                                            &
         &        defaultDiskComponent    %angularMomentumIsGettable()                                                                                                 &
         &  )                                                                                                                                                          &
         & ) call Galacticus_Error_Report                                                                                                                              &
         &        (                                                                                                                                                    &
         &         'this method requires that massStellar, massGas, halfMassRadius, and angularMomentum properties must all be gettable for the disk component.'    // &
         &         Galacticus_Component_List(                                                                                                                          &
         &                                   'disk'                                                                                                                 ,  &
         &                                   defaultDiskComponent    %    massStellarAttributeMatch(requireGettable=.true.).intersection.                              &
         &                                   defaultDiskComponent    %        massGasAttributeMatch(requireGettable=.true.).intersection.                              &
         &                                   defaultDiskComponent    % halfMassRadiusAttributeMatch(requireGettable=.true.).intersection.                              &
         &                                   defaultDiskComponent    %angularMomentumAttributeMatch(requireGettable=.true.)                                            &
         &                                  )                                                                                                                       // &
         &         {introspection:location}                                                                                                                            &
         &        )
    if     (                                                                                                                                                           &
         &  .not.                                                                                                                                                      &
         &       (                                                                                                                                                     &
         &        defaultSpheroidComponent%    massStellarIsGettable().and.                                                                                            &
         &        defaultSpheroidComponent%        massGasIsGettable().and.                                                                                            &
         &        defaultSpheroidComponent% halfMassRadiusIsGettable().and.                                                                                            &
         &        defaultSpheroidComponent%angularMomentumIsGettable()                                                                                                 &
         &  )                                                                                                                                                          &
         & ) call Galacticus_Error_Report                                                                                                                              &
         &        (                                                                                                                                                    &
         &         'this method requires that massStellar, massGas, halfMassRadius, and angularMomentum properties must all be gettable for the spheroid component.'// &
         &         Galacticus_Component_List(                                                                                                                          &
         &                                   'spheroid'                                                                                                             ,  &
         &                                   defaultSpheroidComponent%    massStellarAttributeMatch(requireGettable=.true.).intersection.                              &
         &                                   defaultSpheroidComponent%        massGasAttributeMatch(requireGettable=.true.).intersection.                              &
         &                                   defaultSpheroidComponent% halfMassRadiusAttributeMatch(requireGettable=.true.).intersection.                              &
         &                                   defaultSpheroidComponent%angularMomentumAttributeMatch(requireGettable=.true.)                                            &
         &                                  )                                                                                                                       // &
         &         {introspection:location}                                                                                                                            &
         &        )
    !# <objectBuilder class="mergerMassMovements" name="mergerMassMovements_" source="parameters"/>
    self=mergerProgenitorPropertiesCole2000(mergerMassMovements_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="mergerMassMovements_"/>
    return
  end function cole2000ConstructorParameters
  
 function cole2000ConstructorInternal(mergerMassMovements_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily cole2000} merger progenitor properties class.
    implicit none
    type (mergerProgenitorPropertiesCole2000)                        :: self
    class(mergerMassMovementsClass          ), intent(in   ), target :: mergerMassMovements_
    !# <constructorAssign variables="*mergerMassMovements_"/>
    
    return
  end function cole2000ConstructorInternal

  subroutine cole2000Destructor(self)
    !% Destructor for the {\normalfont \ttfamily cole2000} merger progenitor properties class.
    implicit none
    type (mergerProgenitorPropertiesCole2000) :: self

    !# <objectDestructor name="self%mergerMassMovements_"/>
    return
  end subroutine cole2000Destructor
  
  subroutine cole2000Get(self,nodeSatellite,nodeHost,massSatellite,massHost,massSpheroidSatellite,massSpheroidHost,massSpheroidHostPreMerger,radiusSatellite,radiusHost,factorAngularMomentum,massSpheroidRemnant,massGasSpheroidRemnant)
    !% Computes various properties of the progenitor galaxies useful for calculations of merger remnant sizes.
    use Galacticus_Nodes                  , only : nodeComponentDisk, nodeComponentSpheroid
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Galacticus_Error
    use Root_Finder
    implicit none
    class           (mergerProgenitorPropertiesCole2000), intent(inout)         :: self
    type            (treeNode                          ), intent(inout), target :: nodeSatellite                  , nodeHost
    double precision                                    , intent(  out)         :: factorAngularMomentum          , massHost                         , &
         &                                                                         radiusHost                     , massSpheroidHost                 , &
         &                                                                         massSpheroidHostPreMerger      , massGasSpheroidRemnant           , &
         &                                                                         massSpheroidRemnant            , massSatellite                    , &
         &                                                                         radiusSatellite                , massSpheroidSatellite
    class           (nodeComponentDisk                 ), pointer               :: diskHost                       , diskSatelite
    class           (nodeComponentSpheroid             ), pointer               :: spheroidHost                   , spheroidSatellite
    type            (rootFinder                        ), save                  :: finder
    !$omp threadprivate(finder)
    double precision                                                            :: massComponent                  , factorDarkMatterDiskHost         , &
         &                                                                         radiusHalfMassDiskHost         , factorDarkMatterSpheroidHost     , &
         &                                                                         radiusHalfMassSpheroidHost     , factorDarkMatterDiskSatellite    , &
         &                                                                         radiusHalfMassDiskSatellite    , factorDarkMatterSpheroidSatellite, &
         &                                                                         radiusHalfMassSpheroidSatellite
    integer                                                                     :: destinationGasSatellite, destinationGasHost       , &
         &                                                                         destinationStarsHost   , destinationStarsSatellite
    logical                                                                     :: mergerIsMajor

    ! Find how mass is moved by the merger.
    call self%mergerMassMovements_%get(nodeSatellite,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    ! Initialize our root finder.
    if (.not.finder%isInitialized()) then
       call finder%rootFunction(cole2000HalfMassRadiusRoot                      )
       call finder%tolerance   (toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6)
       call finder%rangeExpand (                                                             &
            &                   rangeExpandUpward            =2.0d0                        , &
            &                   rangeExpandDownward          =0.5d0                        , &
            &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                   rangeExpandType              =rangeExpandMultiplicative      &
            &                  )
    end if
    ! Get the disk and spheroid components of host and satellite.
    diskHost          => nodeHost     %disk    ()
    spheroidHost      => nodeHost     %spheroid()
    diskSatelite      => nodeSatellite%disk    ()
    spheroidSatellite => nodeSatellite%spheroid()
    ! Find the baryonic masses of the two galaxies.
    massSatellite=Galactic_Structure_Enclosed_Mass(nodeSatellite,massType=massTypeGalactic)
    massHost     =Galactic_Structure_Enclosed_Mass(     nodeHost,massType=massTypeGalactic)
    ! Compute dark matter factors. These are the specific angular momenta of components divided by sqrt(G M r) where M is the
    ! component mass and r its half-mass radius. We use a weighted average of these factors to infer the specific angular momentum
    ! of the remnant from its mass and radius.
    massComponent                  =          spheroidHost%massStellar   () &
         &                          +         spheroidHost%massGas       ()
    radiusHalfMassSpheroidHost     =          spheroidHost%halfMassRadius()
    if (radiusHalfMassSpheroidHost > 0.0d0 .and. massComponent > 0.0d0) then
       factorDarkMatterSpheroidHost=           spheroidHost%angularMomentum()                           &
            &                       /(massComponent**1.5d0)                                                      &
            &                       /sqrt(gravitationalConstantGalacticus*          radiusHalfMassSpheroidHost)
    else
       factorDarkMatterSpheroidHost=0.0d0
    end if
    massComponent                  =          diskHost%massStellar   () &
         &                          +         diskHost%massGas       ()
    radiusHalfMassDiskHost         =          diskHost%halfMassRadius()
    if (radiusHalfMassDiskHost > 0.0d0 .and. massComponent > 0.0d0) then
       factorDarkMatterDiskHost         =          diskHost%angularMomentum()                           &
            &                            /(massComponent**1.5d0)                                                 &
            &                            /sqrt(gravitationalConstantGalacticus*         radiusHalfMassDiskHost)
    else
       factorDarkMatterDiskHost=0.0d0
    end if
    massComponent                  = spheroidSatellite%massStellar   () &
         &                          +spheroidSatellite%massGas       ()
    radiusHalfMassSpheroidSatellite= spheroidSatellite%halfMassRadius()
    if (radiusHalfMassSpheroidSatellite > 0.0d0 .and. massComponent > 0.0d0) then
       factorDarkMatterSpheroidSatellite= spheroidSatellite%angularMomentum()                           &
            &                            /(massComponent**1.5d0)                                                 &
            &                            /sqrt(gravitationalConstantGalacticus*radiusHalfMassSpheroidSatellite)
    else
       factorDarkMatterSpheroidSatellite=0.0d0
    end if
    massComponent                  =     diskSatelite%massStellar   () &
        &                           +    diskSatelite%massGas       ()
    radiusHalfMassDiskSatellite    =     diskSatelite%halfMassRadius()
    if (radiusHalfMassDiskSatellite > 0.0d0 .and. massComponent > 0.0d0) then
       factorDarkMatterDiskSatellite    =     diskSatelite%angularMomentum()                           &
            &                            /(massComponent**1.5d0)                                                 &
            &                            /sqrt(gravitationalConstantGalacticus*    radiusHalfMassDiskSatellite)
    else
       factorDarkMatterDiskSatellite=0.0d0
    end if
    ! Find the masses of material that will end up in the spheroid component of the remnant.
    select case (destinationGasHost)
    case (destinationMergerSpheroid)
       massSpheroidHost      =spheroidHost%massGas()                             +diskHost%massGas()
       factorAngularMomentum =spheroidHost%massGas()*factorDarkMatterSpheroidHost+diskHost%massGas()*factorDarkMatterDiskHost
       massGasSpheroidRemnant=spheroidHost%massGas()                             +diskHost%massGas()
       massSpheroidRemnant   =spheroidHost%massGas()                             +diskHost%massGas()
    case (destinationMergerDisk)
       massSpheroidHost      =0.0d0
       factorAngularMomentum =0.0d0
       massGasSpheroidRemnant=0.0d0
       massSpheroidRemnant   =0.0d0
    case (destinationMergerUnmoved)
       massSpheroidHost      =spheroidHost%massGas()
       factorAngularMomentum =spheroidHost%massGas()*factorDarkMatterSpheroidHost
       massGasSpheroidRemnant=spheroidHost%massGas()
       massSpheroidRemnant   =spheroidHost%massGas()
    case default
       call Galacticus_Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    select case (destinationStarsHost)
    case (destinationMergerSpheroid)
       massSpheroidHost     =massSpheroidHost     +spheroidHost%massStellar()                             +diskHost%massStellar()
       factorAngularMomentum=factorAngularMomentum+spheroidHost%massStellar()*factorDarkMatterSpheroidHost+diskHost%massStellar()*factorDarkMatterDiskHost
       massSpheroidRemnant  =massSpheroidRemnant  +spheroidHost%massStellar()                             +diskHost%massStellar()
    case (destinationMergerDisk)
       massSpheroidHost     =massSpheroidHost
       factorAngularMomentum=factorAngularMomentum
    case (destinationMergerUnmoved)
       massSpheroidHost     =massSpheroidHost     +spheroidHost%massStellar()
       factorAngularMomentum=factorAngularMomentum+spheroidHost%massStellar()*factorDarkMatterSpheroidHost
       massSpheroidRemnant  =massSpheroidRemnant  +spheroidHost%massStellar()
    case default
       call Galacticus_Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    select case (destinationGasSatellite)
    case (destinationMergerSpheroid)
       massSpheroidSatellite =                       spheroidSatellite%massGas()                                  +diskSatelite%massGas()
       factorAngularMomentum =factorAngularMomentum +spheroidSatellite%massGas()*factorDarkMatterSpheroidSatellite+diskSatelite%massGas()*factorDarkMatterDiskSatellite
       massGasSpheroidRemnant=massGasSpheroidRemnant+spheroidSatellite%massGas()                                  +diskSatelite%massGas()
       massSpheroidRemnant   =massSpheroidRemnant   +spheroidSatellite%massGas()                                  +diskSatelite%massGas()
    case (destinationMergerDisk)
       massSpheroidSatellite =0.0d0
       factorAngularMomentum =factorAngularMomentum
    case (destinationMergerUnmoved)
       massSpheroidSatellite=                        spheroidSatellite%massGas()
       factorAngularMomentum =factorAngularMomentum +spheroidSatellite%massGas()*factorDarkMatterSpheroidSatellite
       massGasSpheroidRemnant=massGasSpheroidRemnant+spheroidSatellite%massGas()
       massSpheroidRemnant   =massSpheroidRemnant   +spheroidSatellite%massGas()
    case default
       call Galacticus_Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    select case (destinationStarsSatellite)
    case (destinationMergerSpheroid)
       massSpheroidSatellite=massSpheroidSatellite+spheroidSatellite%massStellar()                                  +diskSatelite%massStellar()
       factorAngularMomentum=factorAngularMomentum+spheroidSatellite%massStellar()*factorDarkMatterSpheroidSatellite+diskSatelite%massStellar()*factorDarkMatterDiskSatellite
       massSpheroidRemnant  =massSpheroidRemnant  +spheroidSatellite%massStellar()                                  +diskSatelite%massStellar()
    case (destinationMergerDisk)
       massSpheroidSatellite=massSpheroidSatellite
       factorAngularMomentum=factorAngularMomentum
    case (destinationMergerUnmoved)
       massSpheroidSatellite=massSpheroidSatellite+spheroidSatellite%massStellar()
       massSpheroidRemnant  =massSpheroidRemnant  +spheroidSatellite%massStellar()
       factorAngularMomentum=factorAngularMomentum+spheroidSatellite%massStellar()*factorDarkMatterSpheroidSatellite
    case default
       call Galacticus_Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    ! Compute the half-mass radii of the material that will end up in the remnant spheroid.
    ! Host node.
    if (massSpheroidHost > 0.0d0) then
       cole2000Node       => nodeHost
       cole2000DestinationGas =  destinationGasHost
       cole2000DestinationStars=  destinationStarsHost
       cole2000HalfMass   =  0.0d0 ! Set to zero here so that cole2000HalfMassRadiusRoot() returns the actual half mass.
       cole2000HalfMass   =  0.5d0*cole2000HalfMassRadiusRoot(radiusLarge)
       if (cole2000HalfMassRadiusRoot(0.0d0) <= 0.0d0) then
          radiusHost=finder%find(rootGuess=Galactic_Structure_Radius_Enclosing_Mass(                                 &
               &                                                                                   cole2000Node    , &
               &                                                                    fractionalMass=0.50d0          , &
               &                                                                    massType      =massTypeGalactic  &
               &                                                                   )                                 &
               &                )
       else
          radiusHost      =0.0d0
          massSpheroidHost=0.0d0
       end if
    else
       radiusHost=0.0d0
    end if
    if (massSpheroidSatellite > 0.0d0) then
       cole2000Node       => nodeSatellite
       cole2000DestinationGas =  destinationGasSatellite
       cole2000DestinationStars=  destinationStarsSatellite
       cole2000HalfMass   =  0.0d0 ! Set to zero here so that cole2000HalfMassRadiusRoot() returns the actual half mass.
       cole2000HalfMass   =  0.50d0*cole2000HalfMassRadiusRoot(radiusLarge)
       if (cole2000HalfMassRadiusRoot(0.0d0) <= 0.0d0) then
          radiusSatellite=finder%find(rootGuess=Galactic_Structure_Radius_Enclosing_Mass(                                 &
               &                                                                                        cole2000Node    , &
               &                                                                         fractionalMass=0.50d0          , &
               &                                                                         massType      =massTypeGalactic  &
               &                                                                        )                                 &
               &                     )
       else
          radiusSatellite      =0.0d0
          massSpheroidSatellite=0.0d0
       end if
    else
       radiusSatellite=0.0d0
    end if
    ! Compute the angular momentum factor.
    if (massSpheroidSatellite+massSpheroidHost > 0.0d0) then
       factorAngularMomentum=factorAngularMomentum/(massSpheroidSatellite+massSpheroidHost)
    else
       factorAngularMomentum=1.0d0
    end if
    ! Compute the mass of the host spheroid before the merger.
    massSpheroidHostPreMerger=spheroidHost%massStellar()+spheroidHost%massGas()
    return
  end subroutine cole2000Get

  double precision function cole2000HalfMassRadiusRoot(radius)
    !% Function used in root finding for progenitor galaxy half-mass radii.
    use Satellite_Merging_Mass_Movements
    use Galactic_Structure_Options
    use Galactic_Structure_Enclosed_Masses
    implicit none
    double precision, intent(in   ) :: radius

    ! Initialize enclosed mass to negative of the half mass.
    cole2000HalfMassRadiusRoot=-cole2000HalfMass
    ! Account for gas mass.
    select case (cole2000DestinationGas)
    case (destinationMergerSpheroid)
       cole2000HalfMassRadiusRoot=+cole2000HalfMassRadiusRoot                                                                                         &
            &                     +Galactic_Structure_Enclosed_Mass(cole2000Node,radius,componentType=componentTypeSpheroid,massType=massTypeGaseous) &
            &                     +Galactic_Structure_Enclosed_Mass(cole2000Node,radius,componentType=componentTypeDisk    ,massType=massTypeGaseous)
    case (destinationMergerUnmoved    )
       cole2000HalfMassRadiusRoot=+cole2000HalfMassRadiusRoot                                                                                         &
            &                     +Galactic_Structure_Enclosed_Mass(cole2000Node,radius,componentType=componentTypeSpheroid,massType=massTypeGaseous)
    end select
    ! Account for stellar mass.
    select case (cole2000DestinationStars)
    case (destinationMergerSpheroid)
       cole2000HalfMassRadiusRoot=+cole2000HalfMassRadiusRoot                                                                                         &
            &                     +Galactic_Structure_Enclosed_Mass(cole2000Node,radius,componentType=componentTypeSpheroid,massType=massTypeStellar) &
            &                     +Galactic_Structure_Enclosed_Mass(cole2000Node,radius,componentType=componentTypeDisk    ,massType=massTypeStellar)
    case (destinationMergerUnmoved    )
       cole2000HalfMassRadiusRoot=+cole2000HalfMassRadiusRoot                                                                                         &
            &                     +Galactic_Structure_Enclosed_Mass(cole2000Node,radius,componentType=componentTypeSpheroid,massType=massTypeStellar)
    end select
    return
  end function cole2000HalfMassRadiusRoot
