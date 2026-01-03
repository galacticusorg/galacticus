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
  Implements a merger progenitor properties class which uses a standard calculation.
  !!}

  use :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass

  !![
  <mergerProgenitorProperties name="mergerProgenitorPropertiesStandard">
   <description>
    A merger progenitor properties class which implements a standard method to compute progenitor properties. Masses of
    progenitors are set to
    \begin{equation}
     M_\mathrm{host|satellite} = \sum_{i=\mathrm{disk|spheroid}} \sum_{j=\mathrm{stars|gas}} M_{i,j},
    \end{equation}
    where $M_{i,j}$ is the mass of mass type $j$ in \gls{component} $i$. Masses of progenitors that will end up in the remnant
    spheroid are set to
    \begin{equation}
     M_\mathrm{spheroid\,\,host|satellite} = \sum_{i=\mathrm{disk|spheroid}} \sum_{j=\mathrm{stars|gas}} M_{i,j} \delta_{i,j},
    \end{equation}
    where $\delta_{i,j}=0$ of mass type $j$ in \gls{component} $i$ will end up in the remnant spheroid and $0$ otherwise. Radii
    of material that will end up in the spheroid are set to
    \begin{equation}
     r_\mathrm{host|satellite} = {1 \over M_\mathrm{spheroid\,\,host|satellite}} \sum_{i=\mathrm{disk|spheroid}}
     \sum_{j=\mathrm{stars|gas}} M_{i,j} r_{1/2\,\,i,j} \delta_{i,j}.
    \end{equation}
    Finally, the angular momentum factor is set to
    \begin{equation}
     f_\mathrm{AM\,\,host|satellite} = {1 \over M_\mathrm{spheroid\,\,host|satellite}} \sum_{i=\mathrm{disk|spheroid}}
     \sum_{j=\mathrm{stars|gas}} M_{i,j} {J_{i,j} \over \mathrm{G} M^{3/2}_{i,j} r_{1/2\,\,i,j}} \delta_{i,j},
    \end{equation}
    where $J_{i,j}$ is the angular momentum or pseudo-angular momentum of mass type $j$ in \gls{component} $i$.
   </description>
  </mergerProgenitorProperties>
  !!]
  type, extends(mergerProgenitorPropertiesClass) :: mergerProgenitorPropertiesStandard
     !!{
     A merger progenitor properties class which uses a standard calculation.
     !!}
     private
     class(mergerMassMovementsClass), pointer :: mergerMassMovements_ => null()
   contains
     final     ::        standardDestructor
     procedure :: get => standardGet
  end type mergerProgenitorPropertiesStandard

  interface mergerProgenitorPropertiesStandard
     !!{
     Constructors for the \refClass{mergerProgenitorPropertiesStandard} merger progenitor properties class.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface mergerProgenitorPropertiesStandard

contains

  function standardConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerProgenitorPropertiesStandard} merger progenitor properties class which takes a parameter list as input.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Error_Report            , Component_List
    use :: Galacticus_Nodes, only : defaultDiskComponent    , defaultSpheroidComponent
    use :: Input_Parameters, only : inputParameter          , inputParameters
    implicit none
    type (mergerProgenitorPropertiesStandard)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(mergerMassMovementsClass          ), pointer       :: mergerMassMovements_

    if     (                                                                                                                                                           &
         &  .not.                                                                                                                                                      &
         &       (                                                                                                                                                     &
         &        defaultDiskComponent    %    massStellarIsGettable().and.                                                                                            &
         &        defaultDiskComponent    %        massGasIsGettable().and.                                                                                            &
         &        defaultDiskComponent    % halfMassRadiusIsGettable().and.                                                                                            &
         &        defaultDiskComponent    %angularMomentumIsGettable()                                                                                                 &
         &  )                                                                                                                                                          &
         & ) call Error_Report                                                                                                                                         &
         &        (                                                                                                                                                    &
         &         'this method requires that massStellar, massGas, halfMassRadius, and angularMomentum properties must all be gettable for the disk component.'    // &
         &         Component_List(                                                                                                                                     &
         &                        'disk'                                                                                                                            ,  &
         &                        defaultDiskComponent    %    massStellarAttributeMatch(requireGettable=.true.).intersection.                                         &
         &                        defaultDiskComponent    %        massGasAttributeMatch(requireGettable=.true.).intersection.                                         &
         &                        defaultDiskComponent    % halfMassRadiusAttributeMatch(requireGettable=.true.).intersection.                                         &
         &                        defaultDiskComponent    %angularMomentumAttributeMatch(requireGettable=.true.)                                                       &
         &                       )                                                                                                                                  // &
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
         & ) call Error_Report                                                                                                                                         &
         &        (                                                                                                                                                    &
         &         'this method requires that massStellar, massGas, halfMassRadius, and angularMomentum properties must all be gettable for the spheroid component.'// &
         &         Component_List(                                                                                                                                     &
         &                        'spheroid'                                                                                                                        ,  &
         &                        defaultSpheroidComponent%    massStellarAttributeMatch(requireGettable=.true.).intersection.                                         &
         &                        defaultSpheroidComponent%        massGasAttributeMatch(requireGettable=.true.).intersection.                                         &
         &                        defaultSpheroidComponent% halfMassRadiusAttributeMatch(requireGettable=.true.).intersection.                                         &
         &                        defaultSpheroidComponent%angularMomentumAttributeMatch(requireGettable=.true.)                                                       &
         &                       )                                                                                                                                  // &
         &         {introspection:location}                                                                                                                            &
         &        )
    !![
    <objectBuilder class="mergerMassMovements" name="mergerMassMovements_" source="parameters"/>
    !!]
    self=mergerProgenitorPropertiesStandard(mergerMassMovements_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerMassMovements_"/>
    !!]
    return
  end function standardConstructorParameters

 function standardConstructorInternal(mergerMassMovements_) result(self)
    !!{
    Internal constructor for the \refClass{mergerProgenitorPropertiesStandard} merger progenitor properties class.
    !!}
    implicit none
    type (mergerProgenitorPropertiesStandard)                        :: self
    class(mergerMassMovementsClass          ), intent(in   ), target :: mergerMassMovements_
    !![
    <constructorAssign variables="*mergerMassMovements_"/>
    !!]

    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
    !!{
    Destructor for the \refClass{mergerProgenitorPropertiesStandard} merger progenitor properties class.
    !!}
    implicit none
    type(mergerProgenitorPropertiesStandard), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerMassMovements_"/>
    !!]
    return
  end subroutine standardDestructor

  subroutine standardGet(self,nodeSatellite,nodeHost,massSatellite,massHost,massSpheroidSatellite,massSpheroidHost,massSpheroidHostPreMerger,radiusSatellite,radiusHost,factorAngularMomentum,massSpheroidRemnant,massGasSpheroidRemnant)
    !!{
    Computes various properties of the progenitor galaxies useful for calculations of merger remnant sizes.
    !!}
    use :: Galactic_Structure_Options      , only : massTypeGalactic
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentDisk             , nodeComponentSpheroid    , treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerDisk         , destinationMergerSpheroid, destinationMergerUnmoved, enumerationDestinationMergerType
    use :: Mass_Distributions              , only : massDistributionClass
    implicit none
    class           (mergerProgenitorPropertiesStandard), intent(inout), target :: self
    type            (treeNode                          ), intent(inout), target :: nodeSatellite                  , nodeHost
    double precision                                    , intent(  out)         :: factorAngularMomentum          , massHost                         , &
         &                                                                         radiusHost                     , massSpheroidHost                 , &
         &                                                                         massSpheroidHostPreMerger      , massGasSpheroidRemnant           , &
         &                                                                         massSpheroidRemnant            , massSatellite                    , &
         &                                                                         radiusSatellite                , massSpheroidSatellite
    class           (massDistributionClass             ), pointer               :: massDistributionHost           , massDistributionSatellite
    class           (nodeComponentDisk                 ), pointer               :: diskHost                       , diskSatellite
    class           (nodeComponentSpheroid             ), pointer               :: spheroidHost                   , spheroidSatellite
    double precision                                    , parameter             :: massComponentMinimum=1.0d-30
    double precision                                                            :: massComponent                  , factorDarkMatterDiskHost         , &
         &                                                                         radiusHalfMassDiskHost         , factorDarkMatterSpheroidHost     , &
         &                                                                         radiusHalfMassSpheroidHost     , factorDarkMatterDiskSatellite    , &
         &                                                                         radiusHalfMassDiskSatellite    , factorDarkMatterSpheroidSatellite, &
         &                                                                         radiusHalfMassSpheroidSatellite
    type            (enumerationDestinationMergerType  )                        :: destinationGasSatellite        , destinationGasHost               , &
         &                                                                         destinationStarsHost           , destinationStarsSatellite
    logical                                                                     :: mergerIsMajor

    ! Find how mass is moved by the merger.
    call self%mergerMassMovements_%get(nodeSatellite,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    ! Get the disk and spheroid components of host and satellite.
    diskHost          => nodeHost     %disk    ()
    spheroidHost      => nodeHost     %spheroid()
    diskSatellite     => nodeSatellite%disk    ()
    spheroidSatellite => nodeSatellite%spheroid()
    ! Find the baryonic masses of the two galaxies.
    massDistributionSatellite => nodeSatellite            %massDistribution(massType=massTypeGalactic)
    massDistributionHost      => nodeHost                 %massDistribution(massType=massTypeGalactic)
    massSatellite             =  massDistributionSatellite%massTotal       (                         )
    massHost                  =  massDistributionHost     %massTotal       (                         )
    ! Compute dark matter factors. These are the specific angular momenta of components divided by sqrt(G M r) where M is the
    ! component mass and r its half-mass radius. We use a weighted average of these factors to infer the specific angular momentum
    ! of the remnant from its mass and radius.
    massComponent                  =+spheroidHost     %massStellar   () &
         &                          +spheroidHost     %massGas       ()
    radiusHalfMassSpheroidHost     =+spheroidHost     %halfMassRadius()
    if (radiusHalfMassSpheroidHost > 0.0d0 .and. massComponent > massComponentMinimum) then
       factorDarkMatterSpheroidHost     =+spheroidHost%angularMomentum()                                       &
            &                            /massComponent**1.5d0                                                 &
            &                            /sqrt(gravitationalConstant_internal*radiusHalfMassSpheroidHost     )
    else
       factorDarkMatterSpheroidHost     =+0.0d0
    end if
    massComponent                  =+    diskHost     %massStellar   () &
        &                           +    diskHost     %massGas       ()
    radiusHalfMassDiskHost         =+    diskHost     %halfMassRadius()
    if (radiusHalfMassDiskHost > 0.0d0 .and. massComponent > massComponentMinimum) then
       factorDarkMatterDiskHost         =+    diskHost%angularMomentum()                                       &
            &                            /massComponent**1.5d0                                                 &
            &                            /sqrt(gravitationalConstant_internal*radiusHalfMassDiskHost         )
    else
       factorDarkMatterDiskHost         =+0.0d0
    end if
    massComponent                  =+spheroidSatellite%massStellar   () &
         &                          +spheroidSatellite%massGas       ()
    radiusHalfMassSpheroidSatellite=+spheroidSatellite%halfMassRadius()
    if (radiusHalfMassSpheroidSatellite > 0.0d0 .and. massComponent > massComponentMinimum) then
       factorDarkMatterSpheroidSatellite=+spheroidSatellite%angularMomentum()                                  &
            &                            /massComponent**1.5d0                                                 &
            &                            /sqrt(gravitationalConstant_internal*radiusHalfMassSpheroidSatellite)
    else
       factorDarkMatterSpheroidSatellite=+0.0d0
    end if
    massComponent                  =+    diskSatellite%massStellar   () &
         &                          +    diskSatellite%massGas       ()
    radiusHalfMassDiskSatellite    =+    diskSatellite%halfMassRadius()
    if (radiusHalfMassDiskSatellite > 0.0d0 .and. massComponent > massComponentMinimum) then
       factorDarkMatterDiskSatellite    =+    diskSatellite%angularMomentum()                                  &
            &                            /massComponent**1.5d0                                                 &
            &                            /sqrt(gravitationalConstant_internal*radiusHalfMassDiskSatellite    )
    else
       factorDarkMatterDiskSatellite    =+0.0d0
    end if
    ! Find the masses of material that will end up in the spheroid component of the remnant.
    select case (destinationGasHost%ID)
    case (destinationMergerSpheroid%ID)
       massSpheroidHost      =spheroidHost%massGas()                             +diskHost%massGas()
       radiusHost            =spheroidHost%massGas()*radiusHalfMassSpheroidHost  +diskHost%massGas()*radiusHalfMassDiskHost
       factorAngularMomentum =spheroidHost%massGas()*factorDarkMatterSpheroidHost+diskHost%massGas()*factorDarkMatterDiskHost
       massGasSpheroidRemnant=spheroidHost%massGas()                             +diskHost%massGas()
       massSpheroidRemnant   =spheroidHost%massGas()                             +diskHost%massGas()
    case (destinationMergerDisk   %ID)
       massSpheroidHost      =0.0d0
       radiusHost            =0.0d0
       factorAngularMomentum =0.0d0
       massGasSpheroidRemnant=0.0d0
       massSpheroidRemnant   =0.0d0
    case (destinationMergerUnmoved%ID)
       massSpheroidHost      =spheroidHost%massGas()
       radiusHost            =spheroidHost%massGas()*radiusHalfMassSpheroidHost
       factorAngularMomentum =spheroidHost%massGas()*factorDarkMatterSpheroidHost
       massGasSpheroidRemnant=spheroidHost%massGas()
       massSpheroidRemnant   =spheroidHost%massGas()
    case default
       call Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    select case (destinationStarsHost%ID)
    case (destinationMergerSpheroid%ID)
       massSpheroidHost     =massSpheroidHost     +spheroidHost%massStellar()                             +diskHost%massStellar()
       radiusHost           =radiusHost           +spheroidHost%massStellar()*radiusHalfMassSpheroidHost  +diskHost%massStellar()*radiusHalfMassDiskHost
       factorAngularMomentum=factorAngularMomentum+spheroidHost%massStellar()*factorDarkMatterSpheroidHost+diskHost%massStellar()*factorDarkMatterDiskHost
       massSpheroidRemnant  =massSpheroidRemnant  +spheroidHost%massStellar()                             +diskHost%massStellar()
    case (destinationMergerDisk   %ID)
       massSpheroidHost     =massSpheroidHost
       radiusHost           =radiusHost
       factorAngularMomentum=factorAngularMomentum
    case (destinationMergerUnmoved%ID)
       massSpheroidHost     =massSpheroidHost     +spheroidHost%massStellar()
       radiusHost           =radiusHost           +spheroidHost%massStellar()*radiusHalfMassSpheroidHost
       factorAngularMomentum=factorAngularMomentum+spheroidHost%massStellar()*factorDarkMatterSpheroidHost
       massSpheroidRemnant  =massSpheroidRemnant  +spheroidHost%massStellar()
    case default
       call Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    select case (destinationGasSatellite%ID)
    case (destinationMergerSpheroid%ID)
       massSpheroidSatellite =                       spheroidSatellite%massGas()                                  +diskSatellite%massGas()
       radiusSatellite       =                       spheroidSatellite%massGas()*radiusHalfMassSpheroidSatellite  +diskSatellite%massGas()*radiusHalfMassDiskSatellite
       factorAngularMomentum =factorAngularMomentum +spheroidSatellite%massGas()*factorDarkMatterSpheroidSatellite+diskSatellite%massGas()*factorDarkMatterDiskSatellite
       massGasSpheroidRemnant=massGasSpheroidRemnant+spheroidSatellite%massGas()                                  +diskSatellite%massGas()
       massSpheroidRemnant   =massSpheroidRemnant   +spheroidSatellite%massGas()                                  +diskSatellite%massGas()
    case (destinationMergerDisk   %ID)
       massSpheroidSatellite =0.0d0
       radiusSatellite       =0.0d0
       factorAngularMomentum =factorAngularMomentum
    case (destinationMergerUnmoved%ID)
       massSpheroidSatellite =                       spheroidSatellite%massGas()
       radiusSatellite       =                       spheroidSatellite%massGas()*radiusHalfMassSpheroidSatellite
       factorAngularMomentum =factorAngularMomentum +spheroidSatellite%massGas()*factorDarkMatterSpheroidSatellite
       massGasSpheroidRemnant=massGasSpheroidRemnant+spheroidSatellite%massGas()
       massSpheroidRemnant   =massSpheroidRemnant   +spheroidSatellite%massGas()
    case default
       call Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    select case (destinationStarsSatellite%ID)
    case (destinationMergerSpheroid%ID)
       massSpheroidSatellite=massSpheroidSatellite+spheroidSatellite%massStellar()                                  +diskSatellite%massStellar()
       radiusSatellite      =radiusSatellite      +spheroidSatellite%massStellar()*radiusHalfMassSpheroidSatellite  +diskSatellite%massStellar()*radiusHalfMassDiskSatellite
       factorAngularMomentum=factorAngularMomentum+spheroidSatellite%massStellar()*factorDarkMatterSpheroidSatellite+diskSatellite%massStellar()*factorDarkMatterDiskSatellite
       massSpheroidRemnant  =massSpheroidRemnant  +spheroidSatellite%massStellar()                                  +diskSatellite%massStellar()
    case (destinationMergerDisk   %ID)
       massSpheroidSatellite=massSpheroidSatellite
       radiusSatellite      =radiusSatellite
       factorAngularMomentum=factorAngularMomentum
    case (destinationMergerUnmoved%ID)
       massSpheroidSatellite=massSpheroidSatellite+spheroidSatellite%massStellar()
       radiusSatellite      =radiusSatellite      +spheroidSatellite%massStellar()*radiusHalfMassSpheroidSatellite
       factorAngularMomentum=factorAngularMomentum+spheroidSatellite%massStellar()*factorDarkMatterSpheroidSatellite
       massSpheroidRemnant  =massSpheroidRemnant  +spheroidSatellite%massStellar()
    case default
       call Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    ! Compute the angular momentum factor.
    if (massSpheroidSatellite+massSpheroidHost > massComponentMinimum) then
       factorAngularMomentum=factorAngularMomentum/(massSpheroidSatellite+massSpheroidHost)
    else
       factorAngularMomentum=1.0d0
    end if
    ! Trap cases where radius is zero, but mass is finite (due to numerical inaccuracies).
    if (radiusHost      <= 0.0d0) massSpheroidHost     =0.0d0
    if (radiusSatellite <= 0.0d0) massSpheroidSatellite=0.0d0
    ! Compute the radii of the spheroid components.
    if (massSpheroidHost      > 0.0d0) radiusHost     =radiusHost     /massSpheroidHost
    if (massSpheroidSatellite > 0.0d0) radiusSatellite=radiusSatellite/massSpheroidSatellite
    ! Compute the mass of the host spheroid before the merger.
    massSpheroidHostPreMerger=spheroidHost%massStellar()+spheroidHost%massGas()
    if (radiusHost <= 0.0d0) massSpheroidHostPreMerger=0.0d0
    ! Clean up.
    !![
    <objectDestructor name="massDistributionSatellite"/>
    <objectDestructor name="massDistributionHost     "/>
    !!]
    return
  end subroutine standardGet
