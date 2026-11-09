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
  Implements a node operator class that applies stochastic gravitational acceleration from fuzzy-dark-matter granules to orbiting satellite halos.
  !!}
  
  use :: Dark_Matter_Particles   , only : darkMatterParticleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <nodeOperator name="nodeOperatorSatelliteGranule">
   <description>A node operator class that applies stochastic gravitational acceleration from fuzzy-dark-matter granules to orbiting satellite halos.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteGranule
     private
     double precision                                     :: massParticleSI
     class           (darkMatterParticleClass  ), pointer :: darkMatterParticle_   => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     integer                                              :: axnCoreID      , aynCoreID      , &
                     &                                       aznCoreID      , axnOuterID     , &
                     &                                       aynOuterID     , aznOuterID     , &
                     &                                       massCoreID     , radiusCoreID   , &
                     &                                       radiusSolitonID, tauID          , &
                     &                                       maxTimeID      , frequencyCoreID, &
                     &                                       stepNumID, solitonStatusID
     double precision                                     :: alphaCore      , alphaOuter
     logical                                              :: applyFiniteSize
   contains
     final     :: satelliteGranuleDestructor
     procedure :: acceleration                => granuleAcceleration
     procedure :: nodesMerge                  => satelliteGranuleNodesMerge
     procedure :: nodeInitialize              => satelliteGranuleNodeInitialize
     procedure :: differentialEvolution       => satelliteGranuleDifferentialEvolution
     procedure :: differentialEvolutionScales => satelliteGranuleDifferentialEvolutionScales
  end type nodeOperatorSatelliteGranule

  interface nodeOperatorSatelliteGranule
     module procedure satelliteGranuleConstructorParameters
     module procedure satelliteGranuleConstructorInternal
  end interface nodeOperatorSatelliteGranule

contains

  function satelliteGranuleConstructorParameters(parameters) result(self)
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorSatelliteGranule)                       :: self
    type            (inputParameters             ), intent(inout)        :: parameters
    class           (darkMatterParticleClass     ),              pointer :: darkMatterParticle_
    class           (darkMatterProfileDMOClass   ),              pointer :: darkMatterProfileDMO_
    double precision                                                     :: alphaCore            , alphaOuter
    logical                                                              :: applyFiniteSize
  
    !![
    <inputParameter>
      <name>alphaCore</name>
      <defaultValue>0.3d0</defaultValue>
      <source>parameters</source>
      <description>Dimensionless scaling factor for the granule-induced acceleration in the solitonic core region.</description>
    </inputParameter>
    <inputParameter>
      <name>alphaOuter</name>
      <defaultValue>5.0d0</defaultValue>
      <source>parameters</source>
      <description>Dimensionless scaling factor for the granule-induced acceleration in the outer halo region.</description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>applyFiniteSize</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If true, suppress the outer stochastic acceleration according to the satellite half-mass radius relative to the effective granule size. If false, use the original point-mass model.</description>
    </inputParameter>
    <objectBuilder class="darkMatterParticle"   name="darkMatterParticle_"   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self =nodeOperatorSatelliteGranule(                      &
         &                            darkMatterParticle_  , &
         &                            darkMatterProfileDMO_, &
         &                            alphaCore            , &
         &                            alphaOuter           , &
         &                            applyFiniteSize        &
         &                           )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_"   />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
  end function satelliteGranuleConstructorParameters

  function satelliteGranuleConstructorInternal(darkMatterParticle_,darkMatterProfileDMO_,alphaCore,alphaOuter,applyFiniteSize) result(self)
    use :: Dark_Matter_Particles         , only : darkMatterParticleFuzzyDarkMatter
    use :: Numerical_Constants_Prefixes  , only : kilo
    use :: Numerical_Constants_Physical  , only : speedLight
    use :: Numerical_Constants_Units     , only : electronVolt
    implicit none
    type            (nodeOperatorSatelliteGranule)                     :: self
    class           (darkMatterParticleClass     ), intent(in), target :: darkMatterParticle_
    class           (darkMatterProfileDMOClass   ), intent(in), target :: darkMatterProfileDMO_
    double precision                              , intent(in)         :: alphaCore            , alphaOuter
    logical                                       , intent(in)         :: applyFiniteSize

    !![
    <constructorAssign variables="*darkMatterParticle_,*darkMatterProfileDMO_,alphaCore,alphaOuter,applyFiniteSize"/>
    <addMetaProperty component="basic"             name="tau"                  id="self%tauID"            isEvolvable="yes"   isCreator="yes"/>
    <addMetaProperty component="basic"             name="maxTime"              id="self%maxTimeID"        isEvolvable="no"    isCreator="yes"/>
    <addMetaProperty component="basic"             name="frequencyCore"        id="self%frequencyCoreID"  isEvolvable="no"    isCreator="yes"/>
    <addMetaProperty component="basic"             name="axnCore"              id="self%axnCoreID"        rank="1"            isCreator="yes"/>
    <addMetaProperty component="basic"             name="aynCore"              id="self%aynCoreID"        rank="1"            isCreator="yes"/>
    <addMetaProperty component="basic"             name="aznCore"              id="self%aznCoreID"        rank="1"            isCreator="yes"/>
    <addMetaProperty component="basic"             name="axnOuter"             id="self%axnOuterID"       rank="1"            isCreator="yes"/>
    <addMetaProperty component="basic"             name="aynOuter"             id="self%aynOuterID"       rank="1"            isCreator="yes"/>
    <addMetaProperty component="basic"             name="aznOuter"             id="self%aznOuterID"       rank="1"            isCreator="yes"/>
    <addMetaProperty component="basic"             name="stepNum"              type="integer"             id="self%stepNumID" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="solitonMassCore"      id="self%massCoreID"       isEvolvable="no"    isCreator="no" />
    <addMetaProperty component="darkMatterProfile" name="solitonRadiusCore"    id="self%radiusCoreID"     isEvolvable="no"    isCreator="no" />
    <addMetaProperty component="darkMatterProfile" name="solitonRadiusSoliton" id="self%radiusSolitonID"  isEvolvable="no"    isCreator="no" />
    <addMetaProperty component="darkMatterProfile" name="solitonStatus"        id="self%solitonStatusID"  type="integer"         isCreator="no" />
    !!]

    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleFuzzyDarkMatter)
       ! darkMatterParticle%mass() is in keV. Convert the rest-mass energy to a particle mass in kg for use with the SI Planck constant.
       self%massParticleSI=+darkMatterParticle__%mass() &
              &            *kilo                        &
              &            *electronVolt                &
              &            /speedLight**2
    class default
       call Error_Report('expected a `darkMatterParticleFuzzyDarkMatter` dark matter particle object'//{introspection:location})
    end select
  end function satelliteGranuleConstructorInternal

  subroutine satelliteGranuleDestructor(self)
    implicit none
    type(nodeOperatorSatelliteGranule), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterParticle_"   />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
  end subroutine satelliteGranuleDestructor

  subroutine satelliteGranuleDifferentialEvolutionScales(self,node)
    !!{
    Set the absolute ODE solver scale for the dimensionless granule phase.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSatellite, treeNode
    implicit none
    class           (nodeOperatorSatelliteGranule), intent(inout)          :: self
    type            (treeNode                    ), intent(inout)          :: node
    class           (nodeComponentBasic          )               , pointer :: basic
    class           (nodeComponentSatellite      )               , pointer :: satellite

    if (.not. node%isSatellite()) return
    
    satellite => node%satellite()
    basic     => node%basic    ()
    call basic%floatRank0MetaPropertyScale(            &
         &                                 self%tauID, &
         &                                 1.0d0       &
         &                                )
    return
  end subroutine satelliteGranuleDifferentialEvolutionScales

  subroutine satelliteGranuleDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    use :: Galacticus_Nodes                , only : nodeComponentBasic, nodeComponentSatellite, nodeComponentDarkMatterProfile, treeNode
    use :: Dark_Matter_Profiles_Soliton_Status, only : enumerationSolitonStatusType, solitonStatusSolitonOnly
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class    (nodeOperatorSatelliteGranule), intent(inout), target  :: self
    type     (treeNode                    ), intent(inout), target  :: node
    type     (treeNode                    ),                pointer :: nodeHost
    logical                                , intent(inout)          :: interrupt
    procedure(interruptTask)               , intent(inout), pointer :: functionInterrupt
    class    (nodeComponentBasic          ),                pointer :: basic
    class    (nodeComponentSatellite      )               , pointer :: satellite
    class    (nodeComponentDarkMatterProfile)             , pointer :: darkMatterProfileHost
    integer                                , intent(in   )          :: propertyType
    double precision                                                :: tauRate, maxTime
    type(enumerationSolitonStatusType)                              :: solitonStatus
    double precision                                 , dimension(3) :: acceleration
    
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    if (.not. node%isSatellite()) return
    
    satellite => node%satellite()
    basic     => node%basic    ()

    ! Use the halo state rather than radiusSoliton as the branch condition. In the NFW-only state radiusSoliton is -1, but outer-halo granules are still present.
    nodeHost              => node%parent
    darkMatterProfileHost => nodeHost%darkMatterProfile()
    solitonStatus         =  enumerationSolitonStatusType(darkMatterProfileHost%integerRank0MetaPropertyGet(self%solitonStatusID))

    if (solitonStatus /= solitonStatusSolitonOnly) then
        ! Soliton+NFW or NFW-only: evolve the phase and compute granule acceleration.
        maxTime      = basic%floatRank0MetaPropertyGet(self%maxTimeID)
        tauRate      = 2.0d0*Pi/maxTime
        acceleration = self%acceleration(node)
    else
        ! Soliton-only: the outer FDM halo, and therefore the granule perturbation, is absent.
        tauRate      = 0.0d0
        acceleration = [0.0d0,0.0d0,0.0d0]
    end if


    call basic%floatRank0MetaPropertyRate( &
          & self%tauID,                    &
          & tauRate                        &
          &                              )
    call satellite%velocityRate(acceleration)
    return
  end subroutine satelliteGranuleDifferentialEvolution

  subroutine satelliteGranuleNodesMerge(self,node)
    !!{RST
    Initialize granule state when a node is about to become a satellite.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    class(nodeOperatorSatelliteGranule), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    call satelliteGranuleNodeInitialize(self,node)
    return
  end subroutine satelliteGranuleNodesMerge

  subroutine satelliteGranuleNodeInitialize(self,node)
    !!{
    Initialize the stochastic granule state of a satellite. Reset the phase,
    store the integration time and core fluctuation frequency, and draw independent dimensionless Gaussian Fourier coefficients for the core and outer stochastic fields.    
    !!}
    use :: Coordinates                     , only : coordinateSpherical         , assignment(=)
    use :: Galacticus_Nodes                , only : nodeComponentBasic          , treeNode
    use :: Statistics_Distributions        , only : distributionFunction1DNormal
    use :: Mass_Distributions              , only : massDistributionClass
    implicit none
    class            (nodeOperatorSatelliteGranule), intent(inout), target     :: self
    type             (treeNode                    ), intent(inout), target     :: node
    type             (treeNode                    ),                pointer    :: nodeHost
    class            (nodeComponentBasic          ),                pointer    :: basic                , basicHost
    class            (massDistributionClass       )               , pointer    :: massDistributionHost_
    type             (coordinateSpherical         )                            :: coordinates
    type             (distributionFunction1DNormal)                            :: normalCore           , normalOuter
    double precision                               , dimension(:), allocatable :: axnCore              , aynCore       , &
                     &                                                            aznCore              , axnOuter      , &
                     &                                                            aynOuter             , aznOuter
    double precision                                                           :: maxTime              , densityCentral, &
                     &                                                            frequencyCore
    integer                                                                    :: i                    , stepNum

    ! Initialize current satellites and secondary progenitors that can become
    ! satellites after a merger; skip primary or parentless central nodes.
    if (.not. node%isSatellite() .and. (node%isPrimaryProgenitor() .or. .not. associated(node%parent))) return

    if (node%isSatellite()) then
       nodeHost => node%parent
    else
       nodeHost => node%parent%firstChild
    end if
    
    do while (associated(nodeHost%parent))
       nodeHost => nodeHost%parent
    end do
    
    basicHost => nodeHost%basic()
    maxTime   = basicHost%time ()

    if (node%isSatellite()) then
       nodeHost => node%parent
    else
       nodeHost => node%parent%firstChild
    end if
    
    massDistributionHost_ => self%darkMatterProfileDMO_%get(nodeHost)
    coordinates           = [0.0d0,0.0d0,0.0d0]
    densityCentral        =+massDistributionHost_%density(coordinates) !Msun/Mpc^3
    
    !![
    <objectDestructor name="massDistributionHost_"/>
    !!]

    ! Core oscillation frequency; reference density 10^9 Msun/kpc^3 = 10^18 Msun/Mpc^3 (Dutta Chowdhury et al. 2021).
    frequencyCore = 10.94d0*(sqrt(densityCentral/1.0d18)) !Gyr^-1

    ! Truncate the Paley-Wiener Fourier series at m = 4*f_core*t_max.
    ! The factor 4 ensures empirical convergence of acceleration statistics.
    stepNum  = ceiling(4.0d0*frequencyCore*maxTime)

    ! Adjust the integration time to match the integer number of modes.
    maxTime = stepNum/(4.0d0*frequencyCore)
    
    basic              => node%basic()
    call basic%floatRank0MetaPropertySet(self%tauID           ,0.0d0       )
    call basic%floatRank0MetaPropertySet(self%frequencyCoreID ,frequencyCore)
    call basic%floatRank0MetaPropertySet(self%maxTimeID ,maxTime)
    call basic%integerRank0MetaPropertySet(self%stepNumID ,stepNum)

    ! The core acceleration is isotropic and a_eff^core denotes its total
    ! three-dimensional RMS amplitude, so each Cartesian component carries
    ! a factor 1/sqrt(3).
    normalCore = distributionFunction1DNormal(0.0d0, self%alphaCore**2/3.0d0)
    normalOuter = distributionFunction1DNormal(0.0d0, self%alphaOuter**2)

    allocate(axnCore(0:stepNum))
    allocate(aynCore(0:stepNum))
    allocate(aznCore(0:stepNum))
    allocate(axnOuter(0:stepNum))
    allocate(aynOuter(0:stepNum))
    allocate(aznOuter(0:stepNum))
    
    ! Draw independent dimensionless Gaussian Fourier coefficients for each
    ! Cartesian component of the core and outer stochastic fields.
    do i = 0, stepNum
       axnCore(i) = normalCore%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)
       aynCore(i) = normalCore%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)
       aznCore(i) = normalCore%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)

       axnOuter(i) = normalOuter%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)
       aynOuter(i) = normalOuter%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)
       aznOuter(i) = normalOuter%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)
    end do
    
    call basic%floatRank1MetaPropertySet(self%axnCoreID ,axnCore)
    call basic%floatRank1MetaPropertySet(self%aynCoreID ,aynCore)
    call basic%floatRank1MetaPropertySet(self%aznCoreID ,aznCore)
    call basic%floatRank1MetaPropertySet(self%axnOuterID,axnOuter)
    call basic%floatRank1MetaPropertySet(self%aynOuterID,aynOuter)
    call basic%floatRank1MetaPropertySet(self%aznOuterID,aznOuter)
    return
  end subroutine satelliteGranuleNodeInitialize

  function granuleAcceleration(self,node)
    !!{
    Compute the stochastic acceleration generated by fuzzy-dark-matter granules
in the host halo. Core and outer stochastic fields are evaluated separately,
and the outer contribution includes the velocity-dependent anisotropy.
    Units: km/s/Gyr (same as other dynamical friction acceleration models).
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSatellite, nodeComponentDarkMatterProfile, treeNode
    use :: Dark_Matter_Profiles_Soliton_Status, only : enumerationSolitonStatusType, solitonStatusNfwOnly
    use :: Coordinates                     , only : coordinateSpherical  , assignment(=)
    use :: Mass_Distributions              , only : massDistributionClass, kinematicsDistributionClass
    use :: Numerical_Constants_Astronomical, only : gigaYear             , gravitationalConstant_internal, MpcPerKmPerSToGyr
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Constants_Physical    , only : plancksConstant
    use :: Vectors                         , only : Vector_Magnitude
    implicit none
    double precision                                               , dimension(3):: granuleAcceleration
    class           (nodeOperatorSatelliteGranule  ), intent(inout), target      :: self
    type            (treeNode                      ), intent(inout)              :: node
    type            (treeNode                      )               , pointer     :: nodeHost
    class           (nodeComponentBasic            )               , pointer     :: basic
    class           (nodeComponentSatellite        )               , pointer     :: satellite
    class           (massDistributionClass         )               , pointer     :: massDistributionHost_
    class           (kinematicsDistributionClass   )               , pointer     :: kinematicsHost_
    class           (nodeComponentDarkMatterProfile)               , pointer     :: darkMatterProfileHost
    type            (coordinateSpherical           )                             :: coordinates            , coordinatesSoliton
    type            (enumerationSolitonStatusType  )                             :: solitonStatus
    double precision                                                             :: velocityDispersionJeans, velocityDispersionJeansSoliton
    double precision                                                             :: densityLocal           , densitySoliton               , &
            &                                                                       massCore               , radiusSoliton                , &
            &                                                                       radiusCore
    double precision                                                             :: amplCore               , amplOuter
    double precision                                                             :: frequencyCore          , frequencyOuter
    double precision                                                             :: radius                 , velocityNorm                 , &
            &                                                                       tau                    , maxTime
    double precision                                                             :: velocityDispersionGranule, finiteSizeFactor
    double precision                                                             :: Atotal                 , factorPerpendicular          , &
            &                                                                       factorParallel
    double precision                                               , dimension(3):: positionSatellite      , velocitySatellite
    double precision                                ,  dimension(:), allocatable :: arrayN
    double precision                                ,  dimension(:), allocatable :: axnCore                , aynCore                      , &
            &                                                                       aznCore                , axnOuter                     , &
            &                                                                       aynOuter               , aznOuter
    double precision                                                             :: axCore                 , ayCore                       , &
            &                                                                       azCore                 , axOuter                      , &
            &                                                                       ayOuter                , azOuter
    double precision                                               , dimension(3):: aRawCore               , aRawOuter                    , &
            &                                                                       aRawPar                , aRawPerp                     , &
            &                                                                       velocityDirection
    double precision                                               , dimension(3):: aPar                   , aPerp
    double precision                                               , parameter   :: fractionSmall = 1.0d-2
    integer                                                                      :: stepNum
    !$GLC attributes initialized :: axnCore, aynCore, aznCore, axnOuter, aynOuter, aznOuter
    
    if (.not. node%isSatellite()) then
      granuleAcceleration = [0.0d0, 0.0d0, 0.0d0]
      return
    end if

    satellite             => node            %satellite(                 )
    basic                 => node            %basic    (                 )
    positionSatellite     =  satellite       %position (                 )
    radius                =  Vector_Magnitude          (positionSatellite)
    velocitySatellite     =  satellite       %velocity (                 )
    velocityNorm          =  Vector_Magnitude          (velocitySatellite)

    nodeHost              => node%parent
    darkMatterProfileHost => nodeHost%darkMatterProfile(                 )

    ! The granules belong to the host FDM halo, so both the density and the
    ! Jeans velocity dispersion must be evaluated from the host profile.
    massDistributionHost_ => self%darkMatterProfileDMO_%get(nodeHost)
    kinematicsHost_       => massDistributionHost_%kinematicsDistribution()
    
    axnCore = basic%floatRank1MetaPropertyGet(self%axnCoreID)
    aynCore = basic%floatRank1MetaPropertyGet(self%aynCoreID)
    aznCore = basic%floatRank1MetaPropertyGet(self%aznCoreID)
    axnOuter = basic%floatRank1MetaPropertyGet(self%axnOuterID)
    aynOuter = basic%floatRank1MetaPropertyGet(self%aynOuterID)
    aznOuter = basic%floatRank1MetaPropertyGet(self%aznOuterID)

    maxTime = basic%floatRank0MetaPropertyGet(self%maxTimeID)
    frequencyCore = basic%floatRank0MetaPropertyGet(self%frequencyCoreID)
    tau           = basic%floatRank0MetaPropertyGet(self%tauID)
    stepNum       = basic%integerRank0MetaPropertyGet(self%stepNumID)
    arrayN        = arrayCos(tau,stepNum)
    solitonStatus = enumerationSolitonStatusType(darkMatterProfileHost%integerRank0MetaPropertyGet(self%solitonStatusID))

    if (solitonStatus == solitonStatusNfwOnly) then
       ! NFW-only: radiusCore and radiusSoliton are sentinel values (-1), so
       ! use only the density and Jeans velocity dispersion at the satellite.
       coordinates                    = [radius,0.0d0,0.0d0]
       velocityDispersionJeans        =+kinematicsHost_%velocityDispersion1D(coordinates,massDistributionHost_,massDistributionHost_) ! km/s
       densityLocal                   =+massDistributionHost_%density(coordinates) ! Msun/Mpc^3
       amplCore                       = 0.0d0
       amplOuter                      = amplitudeOuter(self,velocityDispersionJeans,densityLocal)
       frequencyOuter                 = self%massParticleSI*(velocityDispersionJeans*kilo)**2.0d0/0.35d0/sqrt(2.0d0)/plancksConstant*gigaYear
       velocityDispersionGranule      = velocityDispersionJeans
    else
       massCore      = darkMatterProfileHost%floatRank0MetaPropertyGet(self%massCoreID)
       radiusCore    = darkMatterProfileHost%floatRank0MetaPropertyGet(self%radiusCoreID)
       radiusSoliton = darkMatterProfileHost%floatRank0MetaPropertyGet(self%radiusSolitonID )

       ! Impose a small radial floor relative to the soliton core radius to avoid evaluating the density and Jeans velocity dispersion at r = 0.
       radius             = max(radius,fractionSmall*radiusCore)
       coordinates        = [radius,0.0d0,0.0d0]
       coordinatesSoliton = [radiusSoliton,0.0d0,0.0d0]

       velocityDispersionJeans        =+kinematicsHost_%velocityDispersion1D(coordinates       ,massDistributionHost_,massDistributionHost_) ! km/s
       velocityDispersionJeansSoliton =+kinematicsHost_%velocityDispersion1D(coordinatesSoliton,massDistributionHost_,massDistributionHost_) ! km/s

       densityLocal                   =+massDistributionHost_%density(coordinates       ) ! Msun/Mpc^3
       densitySoliton                 =+massDistributionHost_%density(coordinatesSoliton) ! Msun/Mpc^3

       if (radius<=radiusCore) then
          amplCore  = gravitationalConstant_internal*massCore/(radiusCore**2.0d0)/MpcPerKmPerSToGyr                                     ! km/s/Gyr
          amplOuter = amplitudeOuter(self,velocityDispersionJeansSoliton,densitySoliton)                                                ! km/s/Gyr
          frequencyOuter = self%massParticleSI*(velocityDispersionJeansSoliton*kilo)**2.0d0/0.35d0/sqrt(2.0d0)/plancksConstant*gigaYear ! Gyr^-1
          velocityDispersionGranule = velocityDispersionJeansSoliton
       else if (radius<=radiusSoliton) then
          amplCore  = 0.0d0
          amplOuter = amplitudeOuter(self,velocityDispersionJeansSoliton,densitySoliton)
          frequencyOuter = self%massParticleSI*(velocityDispersionJeansSoliton*kilo)**2.0d0/0.35d0/sqrt(2.0d0)/plancksConstant*gigaYear
          velocityDispersionGranule = velocityDispersionJeansSoliton
       else
          amplCore  = 0.0d0
          amplOuter = amplitudeOuter(self,velocityDispersionJeans,densityLocal)
          frequencyOuter = self%massParticleSI*(velocityDispersionJeans*kilo)**2.0d0/0.35d0/sqrt(2.0d0)/plancksConstant*gigaYear
          velocityDispersionGranule = velocityDispersionJeans
       end if
    end if

    ! Keep the tested point-mass path unchanged. For a finite-size satellite,
    ! multiply only the outer stochastic amplitude by the relative suppression
    ! inferred from the finite-size calibration. The Fourier process and its
    ! correlation frequency are unchanged.
    if (self%applyFiniteSize) then
       finiteSizeFactor=finiteSizeAmplitudeFactor(self,node,velocityDispersionGranule)
       amplOuter       =finiteSizeFactor*amplOuter
    end if

    !![
    <objectDestructor name="kinematicsHost_"      />
    <objectDestructor name="massDistributionHost_"/>
    !!]

    ! The Fourier coefficients and arrayN are dimensionless. Since frequency
    ! is in Gyr^-1 and maxTime is in Gyr, sqrt(frequency*maxTime) is also
    ! dimensionless. Therefore aRawCore and aRawOuter are dimensionless.
    axCore = dot_product(axnCore, arrayN) / sqrt(frequencyCore * maxTime)
    ayCore = dot_product(aynCore, arrayN) / sqrt(frequencyCore * maxTime)
    azCore = dot_product(aznCore, arrayN) / sqrt(frequencyCore * maxTime)

    axOuter = dot_product(axnOuter, arrayN) / sqrt(frequencyOuter * maxTime)
    ayOuter = dot_product(aynOuter, arrayN) / sqrt(frequencyOuter * maxTime)
    azOuter = dot_product(aznOuter, arrayN) / sqrt(frequencyOuter * maxTime)

    aRawCore = [axCore, ayCore, azCore]
    aRawOuter = [axOuter, ayOuter, azOuter]

    if (velocityNorm > 0.0d0) then
        velocityDirection = velocitySatellite/velocityNorm
    else
        velocityDirection = [0.0d0, 0.0d0, 0.0d0]
    end if

    ! velocityNorm and velocityDispersionJeans are both in km/s,
    ! so Atotal, factorPerpendicular, factorParallel, aPar, and aPerp are dimensionless.
    Atotal = sqrt( 1.0d0 + (velocityNorm / velocityDispersionJeans)**2 )
    factorPerpendicular  = sqrt( 1.0d0 + 0.5d0 * (velocityNorm / velocityDispersionJeans)**2 )
    factorParallel   = 1.0d0 / (3.0d0 / Atotal - 2.0d0 / factorPerpendicular)

    ! The velocity-dependent factor is applied only to a_raw_outer.
    ! The acceleration from the core is set to 0 at d > r_c.
    aRawPar  = dot_product(aRawOuter,velocityDirection)*velocityDirection
    aRawPerp = aRawOuter-aRawPar
    aPar     = aRawPar /sqrt(factorParallel )
    aPerp    = aRawPerp/sqrt(factorPerpendicular)

    ! amplCore and amplOuter are in km/s/Gyr, while all stochastic and
    ! anisotropy factors are dimensionless. The returned acceleration is
    ! therefore in km/s/Gyr, as required by satellite%velocityRate().
    granuleAcceleration = amplCore*aRawCore + amplOuter*(aPar+aPerp)
    return
  end function granuleAcceleration

  function arrayCos(tau, m) result(arr)
    !!{
    Return the dimensionless Fourier basis. Both ``tau`` and ``arr`` are
    dimensionless; ``m`` is the number of non-zero Fourier modes. The first
    element is 1/sqrt(2), followed by cos(n*tau/2) for n = 1, 2, ..., m.
    !!}
    implicit none
    integer         , intent(in)                :: m
    double precision, intent(in)                :: tau
    double precision, allocatable, dimension(:) :: arr
    integer                                     :: n

    allocate(arr(0:m))
    arr(0) = 1.0d0 / sqrt(2.0d0)
    do n = 1, m
        arr(n) = cos(dble(n) * tau / 2.0d0)
    end do
    return
  end function arrayCos

  function amplitudeOuter(self,velocityDispersionJeans,densityBackground)
    !!{
    Compute the outer-halo effective granule acceleration amplitude.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal, megaParsec, MpcPerKmPerSToGyr
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Constants_Physical    , only : plancksConstant
    implicit none
    class           (nodeOperatorSatelliteGranule), intent(inout), target :: self
    double precision                              , intent(in   )         :: velocityDispersionJeans, densityBackground
    double precision                                                      :: amplitudeOuter         , LambdaDB         , &
            &                                                                deff                   , meff                                  
    
    ! plancksConstant [J s] / (massParticleSI [kg] * velocity [m/s])
    ! gives LambdaDB in metres; division by megaParsec converts it to Mpc.
    LambdaDB    = plancksConstant/self%massParticleSI/(velocityDispersionJeans*kilo)/megaParsec ! Mpc
    deff        = 0.35d0*LambdaDB                                                               ! Mpc
    meff        = densityBackground*(0.282d0*LambdaDB)**3.0d0                                   ! M_sun

    ! gravitationalConstant_internal*meff/deff^2 has units (km/s)^2/Mpc. Convert to the velocity evolution unit km/s/Gyr.
    amplitudeOuter = gravitationalConstant_internal*meff/(deff**2.0d0)/MpcPerKmPerSToGyr        ! km/s/Gyr
    return
  end function amplitudeOuter

  double precision function finiteSizeAmplitudeFactor(self,node,velocityDispersionJeans)
    !!{RST
    Return the finite-size suppression of the outer granule acceleration.

    The satellite half-mass radius is evaluated from its current bound dark
    matter mass. The effective granule size is
    :math:`d_\mathrm{eff}=0.35\lambda_\mathrm{dB}`. Satellites with
    :math:`r_\mathrm{half}/d_\mathrm{eff}\leq0.40` retain the point-mass
    amplitude. Above this transition, the acceleration amplitude is suppressed
    as :math:`(r_\mathrm{half}/0.40d_\mathrm{eff})^{-0.229}`. Because a
    velocity-diffusion coefficient is proportional to acceleration variance,
    the corresponding diffusion coefficient is multiplied by the square of
    this factor.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentSatellite, treeNode
    use :: Galactic_Structure_Options      , only : componentTypeAll      , massTypeDark, radiusLarge
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Constants_Physical    , only : plancksConstant
    implicit none
    class           (nodeOperatorSatelliteGranule), intent(inout), target   :: self
    type            (treeNode                    ), intent(inout), target   :: node
    double precision                              , intent(in   )           :: velocityDispersionJeans
    class           (nodeComponentSatellite      )               , pointer  :: satellite
    class           (massDistributionClass       )               , pointer  :: massDistributionSatellite_
    double precision                                                        :: deff                      , lambdaDB             , &
            &                                                                  massBound                 , massMaximum          , &
            &                                                                  massHalf                  , massAtOrigin         , &
            &                                                                  radiusHalf                , sizeRatio
    double precision                                            , parameter :: sizeRatioTransition=0.40d0, exponentSize=-0.229d0   !The transition ratio and suppression exponent are determined from numerical calibration of the finite-size correction.

    satellite                 => node% satellite      ()
    massDistributionSatellite_=> node%massDistribution(componentType=componentTypeAll,massType=massTypeDark)

    ! Limit the requested bound mass to the finite mass represented by the
    ! satellite profile, then locate the radius containing half of that mass.
    massMaximum =massDistributionSatellite_%massEnclosedBySphere(radius=radiusLarge)
    massBound   =max(0.0d0,min(satellite%boundMass(),massMaximum))
    massHalf    =0.5d0*massBound
    massAtOrigin=massDistributionSatellite_%massEnclosedBySphere(radius=0.0d0)
    if (massHalf > massAtOrigin) then
       radiusHalf=massDistributionSatellite_%radiusEnclosingMass(mass=massHalf)
    else
       radiusHalf=0.0d0
    end if
    !![
    <objectDestructor name="massDistributionSatellite_"/>
    !!]

    lambdaDB =plancksConstant/self%massParticleSI/(velocityDispersionJeans*kilo)/megaParsec ! Mpc
    deff     =0.35d0*lambdaDB                                                               ! Mpc
    sizeRatio=radiusHalf/deff

    if (sizeRatio <= sizeRatioTransition) then
       finiteSizeAmplitudeFactor=1.0d0
    else
       finiteSizeAmplitudeFactor=(sizeRatio/sizeRatioTransition)**exponentSize
    end if
    return
  end function finiteSizeAmplitudeFactor
