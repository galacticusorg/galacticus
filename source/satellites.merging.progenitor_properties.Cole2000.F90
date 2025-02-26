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
  Implements a merger progenitor properties class which uses the algorithm of \cite{cole_hierarchical_2000}.
  !!}

  use :: Root_Finder                     , only : rootFinder
  use :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass, enumerationDestinationMergerType
  use :: Mass_Distributions              , only : massDistributionClass

  !![
  <mergerProgenitorProperties name="mergerProgenitorPropertiesCole2000">
   <description>
    A merger progenitor properties class which uses the algorithms of \cite{cole_hierarchical_2000} to compute progenitor
    properties. Masses of progenitors are set to
    \begin{equation}
     M_\mathrm{host|satellite} = \sum_{i=\mathrm{disk|spheroid}} \sum_{j=\mathrm{stars|gas}} M_{i,j},
    \end{equation}
    where $M_{i,j}$ is the mass of mass type $j$ in \gls{component} $i$. Masses of progenitors that will end up in the remnant
    spheroid are set to
    \begin{equation}
     M_\mathrm{spheroid\,\,host|satellite} = \sum_{i=\mathrm{disk|spheroid}} \sum_{j=\mathrm{stars|gas}} M_{i,j} \delta_{i,j},
    \end{equation}
    where $\delta_{i,j}=0$ of mass type $j$ in \gls{component} $i$ will end up in the remnant spheroid and $0$ otherwise. Radii
    of material that will end up in the spheroid are set by finding the solution to:
    \begin{equation}
    \sum_{i=\mathrm{disk|spheroid}} \sum_{j=\mathrm{stars|gas}} M_{i,j}(r) \delta_{i,j} = {1 \over 2}
    \sum_{i=\mathrm{disk|spheroid}} \sum_{j=\mathrm{stars|gas}} M_{i,j} \delta_{i,j},
    \end{equation}
    such that the radii are the half-mass radii of the material that will end up in the remnant spheroid. Finally, the angular
    momentum factor is set to
    \begin{equation}
     f_\mathrm{AM\,\,host|satellite} = {1 \over M_\mathrm{spheroid\,\,host|satellite}} \sum_{i=\mathrm{disk|spheroid}}
     \sum_{j=\mathrm{stars|gas}} M_{i,j} {J_{i,j} \over \mathrm{G} M^{3/2}_{i,j} r_{1/2\,\,i,j}} \delta_{i,j},
    \end{equation}
    where $J_{i,j}$ is the angular momentum or pseudo-angular momentum of mass type $j$ in \gls{component} $i$\footnote{This is
    technically not quite what \protect\cite{cole_hierarchical_2000} do. Instead, when computing the masses of the material
    which ends up in the spheroid they include twice the mass of dark matter (accounting for the effects of adiabatic
    contraction) within the half-mass radius of each galaxy (as calculated above). The final angular momentum is then
    $j=\sqrt{\mathrm{G} M_\mathrm{remnant} r_\mathrm{remnant}/2}$ (where $M_\mathrm{remnant}$ includes the contribution from
    dark matter and the factor of $2$ appears to make this the half-mass). This approach is currently not used in \protect\glc\
    since there is no way to get the mass of dark matter enclosed accounting for adiabatic contraction in the general
    case. This is a solvable problem, and so this algorithm is expected to be modified to match that of
    \protect\cite{cole_hierarchical_2000} precisely in a future version of \protect\glc.}.
   </description>
  </mergerProgenitorProperties>
  !!]
  type, extends(mergerProgenitorPropertiesClass) :: mergerProgenitorPropertiesCole2000
     !!{
     A merger progenitor properties class which uses the algorithm of \cite{cole_hierarchical_2000}.
     !!}
     private
     class(mergerMassMovementsClass), pointer :: mergerMassMovements_ => null()
     type (rootFinder              )          :: finder
   contains
     final     ::        cole2000Destructor
     procedure :: get => cole2000Get
  end type mergerProgenitorPropertiesCole2000

  interface mergerProgenitorPropertiesCole2000
     !!{
     Constructors for the {\normalfont \ttfamily cole2000} merger progenitor properties class.
     !!}
     module procedure cole2000ConstructorParameters
     module procedure cole2000ConstructorInternal
  end interface mergerProgenitorPropertiesCole2000

  ! Module global variables used in root finding.
  class           (mergerProgenitorPropertiesCole2000), pointer :: self_
  type            (treeNode                          ), pointer :: node_
  class           (massDistributionClass             ), pointer :: massDistributionSpheroidStellar_, massDistributionDiskStellar_, &
       &                                                           massDistributionSpheroidGaseous_, massDistributionDiskGaseous_
  type            (enumerationDestinationMergerType  )          :: destinationGas_                 , destinationStars_
  double precision                                              :: massHalf_
  !$omp threadprivate(self_,node_,destinationGas_,destinationStars_,massHalf_,massDistributionSpheroidStellar_,massDistributionDiskStellar_,massDistributionSpheroidGaseous_,massDistributionDiskGaseous_)

contains

  function cole2000ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily cole2000} merger progenitor properties class which takes a parameter list as input.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Error_Report            , Component_List
    use :: Galacticus_Nodes, only : defaultDiskComponent    , defaultSpheroidComponent
    use :: Input_Parameters, only : inputParameter          , inputParameters
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
    self=mergerProgenitorPropertiesCole2000(mergerMassMovements_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerMassMovements_"/>
    !!]
    return
  end function cole2000ConstructorParameters

 function cole2000ConstructorInternal(mergerMassMovements_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily cole2000} merger progenitor properties class.
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    implicit none
    type (mergerProgenitorPropertiesCole2000)                        :: self
    class(mergerMassMovementsClass          ), intent(in   ), target :: mergerMassMovements_
    !![
    <constructorAssign variables="*mergerMassMovements_"/>
    !!]
    
    self%finder=rootFinder(                                                             &
         &                 rootFunction                 =cole2000HalfMassRadiusRoot   , &
         &                 toleranceAbsolute            =0.0d+0                       , &
         &                 toleranceRelative            =1.0d-6                       , &
         &                 rangeExpandUpward            =2.0d+0                       , &
         &                 rangeExpandDownward          =0.5d+0                       , &
         &                 rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
         &                 rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
         &                 rangeExpandType              =rangeExpandMultiplicative      &
         &                )
    return
  end function cole2000ConstructorInternal

  subroutine cole2000Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily cole2000} merger progenitor properties class.
    !!}
    implicit none
    type(mergerProgenitorPropertiesCole2000), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerMassMovements_"/>
    !!]
    return
  end subroutine cole2000Destructor

  subroutine cole2000Get(self,nodeSatellite,nodeHost,massSatellite,massHost,massSpheroidSatellite,massSpheroidHost,massSpheroidHostPreMerger,radiusSatellite,radiusHost,factorAngularMomentum,massSpheroidRemnant,massGasSpheroidRemnant)
    !!{
    Computes various properties of the progenitor galaxies useful for calculations of merger remnant sizes.
    !!}
    use :: Galactic_Structure_Options      , only : componentTypeDisk             , componentTypeSpheroid    , massTypeGaseous         , massTypeStellar, &
         &                                          massTypeGalactic              , radiusLarge
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentDisk             , nodeComponentSpheroid    , treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerDisk         , destinationMergerSpheroid, destinationMergerUnmoved
    implicit none
    class           (mergerProgenitorPropertiesCole2000), intent(inout), target :: self
    type            (treeNode                          ), intent(inout), target :: nodeSatellite                  , nodeHost
    double precision                                    , intent(  out)         :: factorAngularMomentum          , massHost                         , &
         &                                                                         radiusHost                     , massSpheroidHost                 , &
         &                                                                         massSpheroidHostPreMerger      , massGasSpheroidRemnant           , &
         &                                                                         massSpheroidRemnant            , massSatellite                    , &
         &                                                                         radiusSatellite                , massSpheroidSatellite
    class           (nodeComponentDisk                 ), pointer               :: diskHost                       , diskSatelite
    class           (nodeComponentSpheroid             ), pointer               :: spheroidHost                   , spheroidSatellite
    class           (massDistributionClass             ), pointer               :: massDistributionHost           , massDistributionSatellite
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
    diskSatelite      => nodeSatellite%disk    ()
    spheroidSatellite => nodeSatellite%spheroid()
    ! Find the baryonic masses of the two galaxies.
    massDistributionSatellite => nodeSatellite            %massDistribution(massType=massTypeGalactic)
    massDistributionHost      => nodeHost                 %massDistribution(massType=massTypeGalactic)
    massSatellite             =  massDistributionSatellite%massTotal       (                         )
    massHost                  =  massDistributionHost     %massTotal       (                         )
    ! Compute dark matter factors. These are the specific angular momenta of components divided by sqrt(G M r) where M is the
    ! component mass and r its half-mass radius. We use a weighted average of these factors to infer the specific angular momentum
    ! of the remnant from its mass and radius.
    massComponent                  =          spheroidHost%massStellar   () &
         &                          +         spheroidHost%massGas       ()
    radiusHalfMassSpheroidHost     =          spheroidHost%halfMassRadius()
    if (radiusHalfMassSpheroidHost > 0.0d0 .and. massComponent > 0.0d0) then
       factorDarkMatterSpheroidHost=           spheroidHost%angularMomentum()                                  &
            &                       /(massComponent**1.5d0)                                                    &
            &                       /sqrt(gravitationalConstant_internal*          radiusHalfMassSpheroidHost)
    else
       factorDarkMatterSpheroidHost=0.0d0
    end if
    massComponent                  =          diskHost%massStellar   () &
         &                          +         diskHost%massGas       ()
    radiusHalfMassDiskHost         =          diskHost%halfMassRadius()
    if (radiusHalfMassDiskHost > 0.0d0 .and. massComponent > 0.0d0) then
       factorDarkMatterDiskHost         =          diskHost%angularMomentum()                                  &
            &                            /(massComponent**1.5d0)                                               &
            &                            /sqrt(gravitationalConstant_internal*         radiusHalfMassDiskHost)
    else
       factorDarkMatterDiskHost=0.0d0
    end if
    massComponent                  = spheroidSatellite%massStellar   () &
         &                          +spheroidSatellite%massGas       ()
    radiusHalfMassSpheroidSatellite= spheroidSatellite%halfMassRadius()
    if (radiusHalfMassSpheroidSatellite > 0.0d0 .and. massComponent > 0.0d0) then
       factorDarkMatterSpheroidSatellite= spheroidSatellite%angularMomentum()                                  &
            &                            /(massComponent**1.5d0)                                               &
            &                            /sqrt(gravitationalConstant_internal*radiusHalfMassSpheroidSatellite)
    else
       factorDarkMatterSpheroidSatellite=0.0d0
    end if
    massComponent                  =     diskSatelite%massStellar   () &
        &                           +    diskSatelite%massGas       ()
    radiusHalfMassDiskSatellite    =     diskSatelite%halfMassRadius()
    if (radiusHalfMassDiskSatellite > 0.0d0 .and. massComponent > 0.0d0) then
       factorDarkMatterDiskSatellite    =     diskSatelite%angularMomentum()                                   &
            &                            /(massComponent**1.5d0)                                               &
            &                            /sqrt(gravitationalConstant_internal*    radiusHalfMassDiskSatellite)
    else
       factorDarkMatterDiskSatellite=0.0d0
    end if
    ! Find the masses of material that will end up in the spheroid component of the remnant.
    select case (destinationGasHost%ID)
    case (destinationMergerSpheroid%ID)
       massSpheroidHost      =spheroidHost%massGas()                             +diskHost%massGas()
       factorAngularMomentum =spheroidHost%massGas()*factorDarkMatterSpheroidHost+diskHost%massGas()*factorDarkMatterDiskHost
       massGasSpheroidRemnant=spheroidHost%massGas()                             +diskHost%massGas()
       massSpheroidRemnant   =spheroidHost%massGas()                             +diskHost%massGas()
    case (destinationMergerDisk   %ID)
       massSpheroidHost      =0.0d0
       factorAngularMomentum =0.0d0
       massGasSpheroidRemnant=0.0d0
       massSpheroidRemnant   =0.0d0
    case (destinationMergerUnmoved%ID)
       massSpheroidHost      =spheroidHost%massGas()
       factorAngularMomentum =spheroidHost%massGas()*factorDarkMatterSpheroidHost
       massGasSpheroidRemnant=spheroidHost%massGas()
       massSpheroidRemnant   =spheroidHost%massGas()
    case default
       call Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    select case (destinationStarsHost%ID)
    case (destinationMergerSpheroid%ID)
       massSpheroidHost     =massSpheroidHost     +spheroidHost%massStellar()                             +diskHost%massStellar()
       factorAngularMomentum=factorAngularMomentum+spheroidHost%massStellar()*factorDarkMatterSpheroidHost+diskHost%massStellar()*factorDarkMatterDiskHost
       massSpheroidRemnant  =massSpheroidRemnant  +spheroidHost%massStellar()                             +diskHost%massStellar()
    case (destinationMergerDisk   %ID)
       massSpheroidHost     =massSpheroidHost
       factorAngularMomentum=factorAngularMomentum
    case (destinationMergerUnmoved%ID)
       massSpheroidHost     =massSpheroidHost     +spheroidHost%massStellar()
       factorAngularMomentum=factorAngularMomentum+spheroidHost%massStellar()*factorDarkMatterSpheroidHost
       massSpheroidRemnant  =massSpheroidRemnant  +spheroidHost%massStellar()
    case default
       call Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    select case (destinationGasSatellite%ID)
    case (destinationMergerSpheroid%ID)
       massSpheroidSatellite =                       spheroidSatellite%massGas()                                  +diskSatelite%massGas()
       factorAngularMomentum =factorAngularMomentum +spheroidSatellite%massGas()*factorDarkMatterSpheroidSatellite+diskSatelite%massGas()*factorDarkMatterDiskSatellite
       massGasSpheroidRemnant=massGasSpheroidRemnant+spheroidSatellite%massGas()                                  +diskSatelite%massGas()
       massSpheroidRemnant   =massSpheroidRemnant   +spheroidSatellite%massGas()                                  +diskSatelite%massGas()
    case (destinationMergerDisk   %ID)
       massSpheroidSatellite =0.0d0
       factorAngularMomentum =factorAngularMomentum
    case (destinationMergerUnmoved%ID)
       massSpheroidSatellite=                        spheroidSatellite%massGas()
       factorAngularMomentum =factorAngularMomentum +spheroidSatellite%massGas()*factorDarkMatterSpheroidSatellite
       massGasSpheroidRemnant=massGasSpheroidRemnant+spheroidSatellite%massGas()
       massSpheroidRemnant   =massSpheroidRemnant   +spheroidSatellite%massGas()
    case default
       call Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    select case (destinationStarsSatellite%ID)
    case (destinationMergerSpheroid%ID)
       massSpheroidSatellite=massSpheroidSatellite+spheroidSatellite%massStellar()                                  +diskSatelite%massStellar()
       factorAngularMomentum=factorAngularMomentum+spheroidSatellite%massStellar()*factorDarkMatterSpheroidSatellite+diskSatelite%massStellar()*factorDarkMatterDiskSatellite
       massSpheroidRemnant  =massSpheroidRemnant  +spheroidSatellite%massStellar()                                  +diskSatelite%massStellar()
    case (destinationMergerDisk   %ID)
       massSpheroidSatellite=massSpheroidSatellite
       factorAngularMomentum=factorAngularMomentum
    case (destinationMergerUnmoved%ID)
       massSpheroidSatellite=massSpheroidSatellite+spheroidSatellite%massStellar()
       massSpheroidRemnant  =massSpheroidRemnant  +spheroidSatellite%massStellar()
       factorAngularMomentum=factorAngularMomentum+spheroidSatellite%massStellar()*factorDarkMatterSpheroidSatellite
    case default
       call Error_Report('unrecognized moveTo descriptor'//{introspection:location})
    end select
    ! Compute the half-mass radii of the material that will end up in the remnant spheroid.
    ! Host node.
    if (massSpheroidHost > 0.0d0) then
       self_                            => self
       node_                            => nodeHost
       massDistributionSpheroidStellar_ => nodeHost%massDistribution(componentType=componentTypeSpheroid,massType=massTypeStellar)
       massDistributionDiskStellar_     => nodeHost%massDistribution(componentType=componentTypeDisk    ,massType=massTypeStellar)
       massDistributionSpheroidGaseous_ => nodeHost%massDistribution(componentType=componentTypeSpheroid,massType=massTypeGaseous)
       massDistributionDiskGaseous_     => nodeHost%massDistribution(componentType=componentTypeDisk    ,massType=massTypeGaseous)
       destinationGas_                  =  destinationGasHost
       destinationStars_                =  destinationStarsHost
       massHalf_                        =  0.0d0 ! Set to zero here so that cole2000HalfMassRadiusRoot() returns the actual half mass.
       massHalf_                        =  0.5d0*cole2000HalfMassRadiusRoot(radiusLarge)
       if (cole2000HalfMassRadiusRoot(0.0d0) <= 0.0d0) then
          radiusHost=self%finder%find(rootGuess=massDistributionHost%radiusEnclosingMass(massFractional=0.50d0))
       else
          radiusHost      =0.0d0
          massSpheroidHost=0.0d0
       end if
       !![
       <objectDestructor name="massDistributionSpheroidStellar_"/>
       <objectDestructor name="massDistributionDiskStellar_"    />
       <objectDestructor name="massDistributionSpheroidGaseous_"/>
       <objectDestructor name="massDistributionDiskGaseous_"    />
       !!]
    else
       radiusHost=0.0d0
    end if
    if (massSpheroidSatellite > 0.0d0) then
       self_                            => self
       node_                            => nodeSatellite
       massDistributionSpheroidStellar_ => nodeSatellite%massDistribution(componentType=componentTypeSpheroid,massType=massTypeStellar)
       massDistributionDiskStellar_     => nodeSatellite%massDistribution(componentType=componentTypeDisk    ,massType=massTypeStellar)
       massDistributionSpheroidGaseous_ => nodeSatellite%massDistribution(componentType=componentTypeSpheroid,massType=massTypeGaseous)
       massDistributionDiskGaseous_     => nodeSatellite%massDistribution(componentType=componentTypeDisk    ,massType=massTypeGaseous)
       destinationGas_                  =  destinationGasSatellite
       destinationStars_                =  destinationStarsSatellite
       massHalf_                        =  0.0d0 ! Set to zero here so that cole2000HalfMassRadiusRoot() returns the actual half mass.
       massHalf_                        =  0.50d0*cole2000HalfMassRadiusRoot(radiusLarge)
       if (cole2000HalfMassRadiusRoot(0.0d0) <= 0.0d0) then
          radiusSatellite=self%finder%find(rootGuess=massDistributionSatellite%radiusEnclosingMass(massFractional=0.50d0))
       else
          radiusSatellite      =0.0d0
          massSpheroidSatellite=0.0d0
       end if
       !![
       <objectDestructor name="massDistributionSpheroidStellar_"/>
       <objectDestructor name="massDistributionDiskStellar_"    />
       <objectDestructor name="massDistributionSpheroidGaseous_"/>
       <objectDestructor name="massDistributionDiskGaseous_"    />
       !!]
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
    ! Clean up.
    !![
    <objectDestructor name="massDistributionSatellite"/>
    <objectDestructor name="massDistributionHost     "/>
    !!]
    return
  end subroutine cole2000Get

  double precision function cole2000HalfMassRadiusRoot(radius)
    !!{
    Function used in root finding for progenitor galaxy half-mass radii.
    !!}
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerSpheroid, destinationMergerUnmoved
    implicit none
    double precision, intent(in   ) :: radius

    ! Initialize enclosed mass to negative of the half mass.
    cole2000HalfMassRadiusRoot=-massHalf_
    ! Account for gas mass.
    select case (destinationGas_%ID)
    case (destinationMergerSpheroid%ID)
       cole2000HalfMassRadiusRoot=+cole2000HalfMassRadiusRoot                                    &
            &                     +massDistributionSpheroidGaseous_%massEnclosedBySphere(radius) &
            &                     +massDistributionDiskGaseous_    %massEnclosedBySphere(radius)
    case (destinationMergerUnmoved %ID)
       cole2000HalfMassRadiusRoot=+cole2000HalfMassRadiusRoot                                    &
            &                     +massDistributionSpheroidGaseous_%massEnclosedBySphere(radius)
    end select
    ! Account for stellar mass.
    select case (destinationStars_%ID)
    case (destinationMergerSpheroid%ID)
       cole2000HalfMassRadiusRoot=+cole2000HalfMassRadiusRoot                                    &
            &                     +massDistributionSpheroidStellar_%massEnclosedBySphere(radius) &
            &                     +massDistributionDiskStellar_    %massEnclosedBySphere(radius)
    case (destinationMergerUnmoved %ID)
       cole2000HalfMassRadiusRoot=+cole2000HalfMassRadiusRoot                                    &
            &                     +massDistributionSpheroidStellar_%massEnclosedBySphere(radius)
    end select
    return
  end function cole2000HalfMassRadiusRoot
