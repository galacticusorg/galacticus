!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
  Implements a black hole seed based on collapse of nuclear star clusters due to runaway stellar collisions.
  !!}
 
  use :: Mass_Distributions , only : massDistributionClass  , kinematicsDistributionClass
  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <blackHoleSeeds name="blackHoleSeedsVergara2023">
   <description>
    A model of black hole seeds in which seeds are formed due to the collapse of nuclear star clusters into a black hole. 
     Based on the model of \cite{Vergara_2023} and \cite{Escala_2021}.
   </description>
  </blackHoleSeeds>
  !!]

  type, extends(blackHoleSeedsClass) :: blackHoleSeedsVergara2023
     !!{
     A node operator class that handle the collapse of nuclear star clusters into a black hole.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     double precision                                   :: massSingleStar                       , radiusSingleStar                 , &
         &                                                 massEfficiency                       , radiusEfficiency                 , &
         &                                                 massThreshold
     integer                                            :: timeStellarMassFormedNSCID           , stellarMassFormedNSCID           , &
         &                                                 radiusNuclearStarClustersID          , blackHoleSeedMassID              , &
         &                                                 velocityNuclearStarClustersID        , ageNuclearStarClustersID         , &
         &                                                 gasMassNuclearStarClustersID         , criticalMassNuclearStarClustersID, &
         &                                                 redshiftBlackHoleSeedFormationID     , stellarMassNuclearStarClustersID , &                                                                    
         &                                                 mergerTreeWeightNuclearStarClustersID 
   contains
     final     ::                     Vergara2023Destructor              
     procedure :: mass             => Vergara2023Mass
     procedure :: spin             => Vergara2023Spin
     procedure :: formationChannel => Vergara2023FormationChannel
  end type blackHoleSeedsVergara2023
  
  interface blackHoleSeedsVergara2023
     !!{
     Constructors for the {\normalfont \ttfamily Vergara2023} black hole seeds class.
     !!}
     module procedure Vergara2023ConstructorParameters
     module procedure Vergara2023ConstructorInternal
  end interface blackHoleSeedsVergara2023

  class (massDistributionClass ), pointer :: massDistribution_, massDistributionStellarNuclearStarCluster_ 
  !$omp threadprivate(massDistribution_,massDistributionStellarNuclearStarCluster_)

contains

  function Vergara2023ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily Vergara2023} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (blackHoleSeedsVergara2023)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass  ), pointer       :: cosmologyFunctions_
    double precision                                           :: massSingleStar, radiusSingleStar, &
       &                                                          massEfficiency, radiusEfficiency, &
       &                                                          massThreshold 

    !![
    <inputParameter>
      <name>massSingleStar</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Specifies the mass of a single star in Solar units.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusSingleStar</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Specifies the radius of a single star in Solar units.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massEfficiency</name>
      <defaultValue>1.0d-1</defaultValue>
      <description>Specifies the efficiency of the mass converted into a black hole seed.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusEfficiency</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Specifies the efficiency of the radius used to compute the critical mass.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massThreshold</name>
      <defaultValue>1.0d3</defaultValue>
      <description>Specifies the minimum stellar mass to apply the seeding prescription.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=blackHoleSeedsVergara2023(massSingleStar, radiusSingleStar, massEfficiency, radiusEfficiency, massThreshold, cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function Vergara2023ConstructorParameters
  
  function Vergara2023ConstructorInternal(massSingleStar, radiusSingleStar, massEfficiency, radiusEfficiency, massThreshold,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily Vergara2023} node operator class.
    !!}
    implicit none
    type            (blackHoleSeedsVergara2023)                        :: self
    class           (cosmologyFunctionsClass  ), intent(in   ), target :: cosmologyFunctions_
    double precision                           , intent(in   )         :: massSingleStar
    double precision                           , intent(in   )         :: radiusSingleStar
    double precision                           , intent(in   )         :: massEfficiency
    double precision                           , intent(in   )         :: radiusEfficiency
    double precision                           , intent(in   )         :: massThreshold
    !![
    <constructorAssign variables="massSingleStar, radiusSingleStar, massEfficiency, radiusEfficiency, massThreshold, *cosmologyFunctions_"/>
    !!]
    !![
    <addMetaProperty   component="NSC"        name="agesStellarMassFormed"               id="self%stellarMassFormedNSCID"                isEvolvable="yes" isCreator="no" />
    <addMetaProperty   component="NSC"        name="agesTimeStellarMassFormed"           id="self%timeStellarMassFormedNSCID"            isEvolvable="yes" isCreator="no" />
    <addMetaProperty   component="NSC"        name="blackHoleSeedMassFormed"             id="self%blackHoleSeedMassID"                   isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty   component="NSC"        name="ageNuclearStarClusters"              id="self%ageNuclearStarClustersID"              isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty   component="NSC"        name="radiusNuclearStarClusters"           id="self%radiusNuclearStarClustersID"           isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty   component="NSC"        name="gasMassNuclearStarClusters"          id="self%gasMassNuclearStarClustersID"          isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty   component="NSC"        name="velocityNuclearStarClusters"         id="self%velocityNuclearStarClustersID"         isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty   component="NSC"        name="redshiftBlackHoleSeedFormation"      id="self%redshiftBlackHoleSeedFormationID"      isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty   component="NSC"        name="stellarMassNuclearStarClusters"      id="self%stellarMassNuclearStarClustersID"      isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty   component="NSC"        name="criticalMassNuclearStarClusters"     id="self%criticalMassNuclearStarClustersID"     isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty   component="NSC"        name="mergerTreeWeightNuclearStarClusters" id="self%mergerTreeWeightNuclearStarClustersID" isEvolvable="no"  isCreator="yes"/>
    !!]
    return
  end function Vergara2023ConstructorInternal

  subroutine Vergara2023Destructor(self)
      !!{
      Destructor for the {\normalfont \ttfamily Vergara2023} black hole seeds class.
      !!}
      implicit none 
      type(blackHoleSeedsVergara2023), intent(inout) :: self
      
      !![
      <objectDestructor name="self%cosmologyFunctions_"/>
      !!]
      return
  end subroutine Vergara2023Destructor

  double precision function Vergara2023Mass(self,node) result(mass)
      !!{
        Compute the nuclear star cluster collapse condition.
      !!}
    use :: Galacticus_Nodes                , only : nodeComponentNSC               , nodeComponentBasic, nodeComponentNSCStandard, treeNode  
    use :: Abundances_Structure            , only : operator(*)
    use :: Numerical_Constants_Math        , only : Pi
    use :: Galactic_Structure_Options      , only : componentTypenuclearStarCluster, massTypeStellar
    use :: Numerical_Constants_Prefixes    , only : mega                           , kilo
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal , radiusSolar       , megaParsec              , gigaYear, &
        &                                           parsec
    implicit none
    class(blackHoleSeedsVergara2023), intent(inout)          :: self
    type (treeNode                 ), intent(inout)          :: node
    class(nodeComponentNSC         )               , pointer :: nuclearStarCluster
    class(nodeComponentBasic       )               , pointer :: basic
    double precision                                         :: radiusNuclearStarCluster         , velocityNuclearStarCluster    , &
        &                                                       massStellarNuclearStarCluster    , massCriticalNuclearStarCluster, &
        &                                                       Theta                            , crossSectionNuclearStarCluster, &
        &                                                       massTimeStellarNuclearStarCluster, ageNuclearStarCluster         , &
        &                                                       time                             , redshift
    double precision                                         :: velocity        = 100.0d0 !km s¯¹
    
    ! Get the nuclear star cluster component.
    nuclearStarCluster => node%NSC()

    ! Detect the type of the nuclear star cluster component.
    select type (nuclearStarCluster)
      class default
          ! Generic type, do nothing.
          mass = 0.0d0
          return
      class is (nodeComponentNSCStandard)
          ! Standard class, get the properties of the nuclear star cluster component.
          radiusNuclearStarCluster                   =  self%radiusEfficiency*nuclearStarCluster%radius()
          
          ! Unphysical nuclear star cluster, do nothing.
          if (nuclearStarCluster%massStellar()<=0.0d0.or.nuclearStarCluster%radius()<=0.0d0) then
            mass=0.0d0
            return
          end if 
          velocityNuclearStarCluster                 = sqrt(gravitationalConstant_internal*nuclearStarCluster%massStellar()/radiusNuclearStarCluster)
          massStellarNuclearStarCluster     = nuclearStarCluster%floatRank0MetaPropertyGet(self%    stellarMassFormedNSCID)
          massTimeStellarNuclearStarCluster = nuclearStarCluster%floatRank0MetaPropertyGet(self%timeStellarMassFormedNSCID)
          
          basic=> node%basic()
          time =  basic%time()

          ! Get the age of the nuclear star cluster.
          if ( 0.0d0 < massStellarNuclearStarCluster) then 
            ageNuclearStarCluster= +time                              &
               &                   -massTimeStellarNuclearStarCluster &
               &                   /massStellarNuclearStarCluster
          else 
            ageNuclearStarCluster = 0.0d0
          end if 
          ! Do nothing if the nuclear star cluster has an unphysicall age or already formed a black hole seed.
          if (ageNuclearStarCluster<=0.0d0.or.nuclearStarCluster%isCollapsed()) then
            mass = 0.0d0
            return
          end if 
          ! Safronov number defined by Binney & Tremaine (2008, https://ui.adsabs.harvard.edu/abs/2008gady.book.....B/abstract)
          Theta=9.54d0*(self%massSingleStar/self%radiusSingleStar)*(velocity/velocityNuclearStarCluster)**2.0d0                   !Adimensional
          ! Probabilistic mean free path defined as in Landau & Lifshitz (1980, https://ui.adsabs.harvard.edu/abs/1981PhT....34a..74L/abstract) and Shu (1991, https://ui.adsabs.harvard.edu/abs/1991pav..book.....S/abstract)
          crossSectionNuclearStarCluster  = 16.0d0* sqrt(Pi)*(1+Theta)*(self%radiusSingleStar*radiusSolar/parsec)**2.0d0      !pc²
          ! Critical mass computation using equation (3) in the model of M.C. Vergara, A. Escala, D.R.G. Schleicher and B. Reinoso. (2023, https://ui.adsabs.harvard.edu/abs/2023MNRAS.522.4224V/abstract)
          massCriticalNuclearStarCluster  = (mega*radiusNuclearStarCluster)**(7.0d0/3.0d0)*((4.0d0*Pi*self%massSingleStar)/(3.0d0*crossSectionNuclearStarCluster*ageNuclearStarCluster*sqrt((gravitationalConstant_internal*megaParsec*(kilo*gigaYear)**2.0d0)*parsec**-3.0d0)))**(2.0d0/3.0d0)
          
          ! Generic type - interrupt and create a standard Black Hole if nuclear star cluster mass is greater than the critical mass.
          if (0.0d0<massCriticalNuclearStarCluster.and.massCriticalNuclearStarCluster<= nuclearStarCluster%massStellar() .and. self%massThreshold <= nuclearStarCluster%massStellar()) then
            
            redshift= self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
            
            call nuclearStarCluster%floatRank0MetaPropertySet(self%ageNuclearStarClustersID             , ageNuclearStarCluster           )
            call nuclearStarCluster%floatRank0MetaPropertySet(self%gasMassNuclearStarClustersID         , nuclearStarCluster%massGas    ())
            call nuclearStarCluster%floatRank0MetaPropertySet(self%stellarMassNuclearStarClustersID     , nuclearStarCluster%massStellar())
            call nuclearStarCluster%floatRank0MetaPropertySet(self%velocityNuclearStarClustersID        , velocityNuclearStarCluster      )
            call nuclearStarCluster%floatRank0MetaPropertySet(self%redshiftBlackHoleSeedFormationID     , redshift                        )
            call nuclearStarCluster%floatRank0MetaPropertySet(self%criticalMassNuclearStarClustersID    , massCriticalNuclearStarCluster  )
            call nuclearStarCluster%floatRank0MetaPropertySet(self%radiusNuclearStarClustersID          , radiusNuclearStarCluster        )
            call nuclearStarCluster%floatRank0MetaPropertySet(self%mergerTreeWeightNuclearStarClustersID, node%hostTree%volumeWeight      )
            ! Adjust stellar mass of the nuclear star cluster
            call nuclearStarCluster%massStellarSet      (                                          &
                &                                        ( 1.0d0-self%massEfficiency             ) &
                &                                         *nuclearStarCluster%massStellar      ()  &
                &                                       )
            ! Adjust stellar abundances of the nuclear star cluster 
            call nuclearStarCluster%abundancesStellarSet(                                          &
                &                                        ( 1.0d0-self%massEfficiency             ) &
                &                                         *nuclearStarCluster%abundancesStellar()  &
                &                                       )

            call nuclearStarCluster%isCollapsedSet(.true.)
            mass    = self%massEfficiency*nuclearStarCluster%massStellar()
            call nuclearStarCluster%floatRank0MetaPropertySet(self%blackHoleSeedMassID, mass)            
          else
            mass = 0.0d0
          end if 
    end select
    return
  end function Vergara2023Mass

  double precision function Vergara2023Spin(self,node) result(spin)
    !!{
    Compute the spin of the seed black hole.
    !!}
    implicit none
    class(blackHoleSeedsVergara2023), intent(inout) :: self
    type (treeNode                 ), intent(inout) :: node
    !$GLC attributes unused :: node
    ! Assume zero spin.
    spin=0.0d0
    return
  end function Vergara2023Spin

  integer function Vergara2023FormationChannel (self,node) result(enumeration)
    !!{
    Compute the spin of the seed black hole.
    !!}
    implicit none
    class(blackHoleSeedsVergara2023), intent(inout) :: self
    type (treeNode                 ), intent(inout) :: node
    !$GLC attributes unused :: node
    ! Assume zero spin.
    enumeration=1
    return
  end function Vergara2023FormationChannel
