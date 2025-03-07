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
  Implements a property extractor class the properties of nuclear star cluster when a black hole seed is formed using the model of \cite{vergara_global_2023}.
  !!}
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorBlackHoleSeedingVergara2023">
   <description>
    A property extractor class for the properties of the nuclear star cluster at the moment of the black hole formation.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorBlackHoleSeedingVergara2023
     !!{
     A property extractor class for the velocity dispersion at a set of radii.
     !!}
     private
     integer  :: radiusNuclearStarClustersID          , blackHoleSeedMassID              , &
         &       velocityNuclearStarClustersID        , ageNuclearStarClustersID         , &
         &       gasMassNuclearStarClustersID         , criticalMassNuclearStarClustersID, &
         &       redshiftBlackHoleSeedFormationID     , stellarMassNuclearStarClustersID , &                                                                    
         &       mergerTreeWeightNuclearStarClustersID   
   contains
     procedure :: elementCount       => blackHoleSeedingVergara2023ElementCount
     procedure :: extract            => blackHoleSeedingVergara2023Extract
     procedure :: names              => blackHoleSeedingVergara2023Names
     procedure :: descriptions       => blackHoleSeedingVergara2023Descriptions
     procedure :: unitsInSI          => blackHoleSeedingVergara2023UnitsInSI
  end type nodePropertyExtractorBlackHoleSeedingVergara2023

  interface nodePropertyExtractorBlackHoleSeedingVergara2023
     !!{
     Constructors for the ``BlackHoleSeedingVergara2023'' output analysis class.
     !!}
     module procedure blackHoleSeedingVergara2023ConstructorParameters
     module procedure blackHoleSeedingVergara2023ConstructorInternal
  end interface nodePropertyExtractorBlackHoleSeedingVergara2023

contains

  function blackHoleSeedingVergara2023ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily BlackHoleSeedingVergara2023} property extractor class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorBlackHoleSeedingVergara2023)                :: self
    type (inputParameters                                 ), intent(inout) :: parameters

    self=nodePropertyExtractorBlackHoleSeedingVergara2023()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function blackHoleSeedingVergara2023ConstructorParameters

  function blackHoleSeedingVergara2023ConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily BlackHoleSeedingVergara2023} property extractor class.
    !!}
    implicit none
    type          (nodePropertyExtractorBlackHoleSeedingVergara2023) :: self
    !![
    <addMetaProperty   component="NSC"  name="blackHoleSeedMassFormed"             id="self%blackHoleSeedMassID"                   isEvolvable="no"  isCreator="no"/>
    <addMetaProperty   component="NSC"  name="velocityNuclearStarClusters"         id="self%velocityNuclearStarClustersID"         isEvolvable="no"  isCreator="no"/>
    <addMetaProperty   component="NSC"  name="ageNuclearStarClusters"              id="self%ageNuclearStarClustersID"              isEvolvable="no"  isCreator="no"/>
    <addMetaProperty   component="NSC"  name="radiusNuclearStarClusters"           id="self%radiusNuclearStarClustersID"           isEvolvable="no"  isCreator="no"/>
    <addMetaProperty   component="NSC"  name="gasMassNuclearStarClusters"          id="self%gasMassNuclearStarClustersID"          isEvolvable="no"  isCreator="no"/>
    <addMetaProperty   component="NSC"  name="stellarMassNuclearStarClusters"      id="self%stellarMassNuclearStarClustersID"      isEvolvable="no"  isCreator="no"/>
    <addMetaProperty   component="NSC"  name="criticalMassNuclearStarClusters"     id="self%criticalMassNuclearStarClustersID"     isEvolvable="no"  isCreator="no"/>
    <addMetaProperty   component="NSC"  name="redshiftBlackHoleSeedFormation"      id="self%redshiftBlackHoleSeedFormationID"      isEvolvable="no"  isCreator="no"/>
    <addMetaProperty   component="NSC"  name="mergerTreeWeightNuclearStarClusters" id="self%mergerTreeWeightNuclearStarClustersID" isEvolvable="no"  isCreator="no"/>
    !!]
    return
  end function blackHoleSeedingVergara2023ConstructorInternal

  integer function blackHoleSeedingVergara2023ElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily BlackHoleSeedingVergara2023} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorBlackHoleSeedingVergara2023), intent(inout) :: self
    double precision                                                  , intent(in   ) :: time
    !$GLC attributes unused :: time

    blackHoleSeedingVergara2023ElementCount=8
    return
  end function blackHoleSeedingVergara2023ElementCount

  function blackHoleSeedingVergara2023Extract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily BlackHoleSeedingVergara2023} property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC
    implicit none
    double precision                                                  , dimension(:) , allocatable :: blackHoleSeedingVergara2023Extract
    class           (nodePropertyExtractorBlackHoleSeedingVergara2023), intent(inout), target      :: self
    type            (treeNode                                        ), intent(inout), target      :: node
    double precision                                                  , intent(in   )              :: time
    type            (multiCounter                                    ), intent(inout), optional    :: instance
    class           (nodeComponentNSC                                )               , pointer     :: nuclearStarCluster

    !$GLC attributes unused :: time, instance

    allocate(blackHoleSeedingVergara2023Extract(8))
    nuclearStarCluster => node%NSC()
    select type (nuclearStarCluster)
    type is (nodeComponentNSC)
      ! Nuclear star cluster does not yet exist.
      blackHoleSeedingVergara2023Extract=[       & 
        &                                 0.0d0, &
        &                                 0.0d0, &
        &                                 0.0d0, &
        &                                 0.0d0, &
        &                                 0.0d0, &
        &                                 0.0d0, &
        &                                 0.0d0, &
        &                                 0.0d0  &
        &                                ]
    class default
      blackHoleSeedingVergara2023Extract=[                                                                                      &
       &                                  nuclearStarCluster%floatRank0MetaPropertyGet(self%redshiftBlackHoleSeedFormationID ), &
       &                                  nuclearStarCluster%floatRank0MetaPropertyGet(self%blackHoleSeedMassID              ), &
       &                                  nuclearStarCluster%floatRank0MetaPropertyGet(self%ageNuclearStarClustersID         ), &
       &                                  nuclearStarCluster%floatRank0MetaPropertyGet(self%radiusNuclearStarClustersID      ), &
       &                                  nuclearStarCluster%floatRank0MetaPropertyGet(self%velocityNuclearStarClustersID    ), &
       &                                  nuclearStarCluster%floatRank0MetaPropertyGet(self%gasMassNuclearStarClustersID     ), &
       &                                  nuclearStarCluster%floatRank0MetaPropertyGet(self%stellarMassNuclearStarClustersID ), &
       &                                  nuclearStarCluster%floatRank0MetaPropertyGet(self%criticalMassNuclearStarClustersID)  &
       &                                 ]
    end select
    return
  end function blackHoleSeedingVergara2023Extract

  subroutine blackHoleSeedingVergara2023Names(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily BlackHoleSeedingVergara2023} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorBlackHoleSeedingVergara2023), intent(inout)                             :: self
    double precision                                                  , intent(in   )                             :: time
    type            (varying_string                                  ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time
    
    allocate(names(8))
    names(1)=var_str('blackHoleFormationRedshift'    )
    names(2)=var_str('blackHoleSeedMass'             )
    names(3)=var_str('nuclearStarClusterAge'         )
    names(4)=var_str('nuclearStarClusterRadius'      )
    names(5)=var_str('nuclearStarClusterVelocity'    )
    names(6)=var_str('nuclearStarClusterGasMass'     )
    names(7)=var_str('nuclearStarClusterStellarMass' )
    names(8)=var_str('nuclearStarClusterCriticalMass')
    return
  end subroutine blackHoleSeedingVergara2023Names

  subroutine blackHoleSeedingVergara2023Descriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily blackHoleSeedingVergara2023} property.
    !!}
    implicit none
    class           (nodePropertyExtractorBlackHoleSeedingVergara2023), intent(inout)                             :: self
    double precision                                                  , intent(in   )                             :: time
    type            (varying_string                                  ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(8))
    descriptions(1)=var_str('Redshift at the formation of the black hole seed in the Vergara et al. (2024) model.'                         )
    descriptions(2)=var_str('Black hole seed mass in the Vergara et al. (2024) model [M⊙].'                                                )
    descriptions(3)=var_str('Mass-weighted age of the nuclear star cluster at which black hole seed is formed [Gyr].'                      )
    descriptions(4)=var_str('Radius of the nuclear star cluster used to compute the critical mass (multiplied by radius efficiency) [Mpc].')
    descriptions(5)=var_str('Velocity of the nuclear star cluster used to compute the critical mass [km s⁻¹].'                             )
    descriptions(6)=var_str('Gaseous mass of the nuclear star cluster used to compute the critical mass [M⊙].'                             )
    descriptions(7)=var_str('Stellar mass of the nuclear star cluster used to compute the critical mass [M⊙].'                             )
    descriptions(8)=var_str('Critical mass of the nuclear star cluster [M⊙].'                                                              )
    return
  end subroutine blackHoleSeedingVergara2023Descriptions

  function blackHoleSeedingVergara2023UnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily BlackHoleSeedingVergara2023} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec, gigayear
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                                  , allocatable  , dimension(:) :: blackHoleSeedingVergara2023UnitsInSI
    class           (nodePropertyExtractorBlackHoleSeedingVergara2023), intent(inout)               :: self
    double precision                                                  , intent(in   )               :: time
    !$GLC attributes unused :: time

    allocate(blackHoleSeedingVergara2023UnitsInSI(8))
    blackHoleSeedingVergara2023UnitsInSI = [            &
     &                                      1.0d0     , &
     &                                      massSolar , &
     &                                      gigayear  , &
     &                                      megaParsec, &
     &                                      kilo      , &
     &                                      massSolar , &
     &                                      massSolar , &
     &                                      massSolar   &
     &                                     ]
    return
  end function blackHoleSeedingVergara2023UnitsInSI
