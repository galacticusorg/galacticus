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
  Implements a node operator class that records details of galaxy-galaxy mergers.
  !!}

  use, intrinsic :: ISO_C_Binding, only : c_size_t

  !![
  <nodeOperator name="nodeOperatorGalaxyMergers">
    <description>
      A node operator class that records details of galaxy-galaxy mergers.  Data recorded consists of:
      output are:
      \begin{description}
       \item The index of the satellite halo in the merger;
       \item The stellar/gas mass of te host/satellite galaxy;
       \item The time at which the merger occurred.
      \end{description}
      Two \refClass{galacticFilterClass}es are accepted, via parameters {\normalfont \ttfamily [galacticFilterSatellite]} and {\normalfont
      \ttfamily [galacticFilterCentral]} which can be used to control which galaxies are included in the output.
  </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorGalaxyMergers
     !!{
     A node operator class that records the times of gas-mass-based major mergers between galaxies.
     !!}
     private
     integer                               :: galaxyMergerHostMassGasID              , galaxyMergerHostMassStellarID               , &
          &                                   galaxyMergerSatelliteMassGasID         , galaxyMergerSatelliteMassStellarID          , &
          &                                   galaxyMergerTimeID                     , galaxyMergerSatelliteIndexID
     integer(c_size_t           )          :: countMergersMaximum
     class  (galacticFilterClass), pointer :: galacticFilterSatellite_      => null(), galacticFilterCentral_             => null()
  contains
     final     ::             galaxyMergersDestructor
     procedure :: autoHook => galaxyMergersAutoHook
  end type nodeOperatorGalaxyMergers
  
  interface nodeOperatorGalaxyMergers
     !!{
     Constructors for the \refClass{nodeOperatorGalaxyMergers} node operator class.
     !!}
     module procedure galaxyMergersConstructorParameters
     module procedure galaxyMergersConstructorInternal
  end interface nodeOperatorGalaxyMergers
  
contains

  function galaxyMergersConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorGalaxyMergers} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorGalaxyMergers)                :: self
    type   (inputParameters          ), intent(inout) :: parameters
    integer(c_size_t                 )                :: countMergersMaximum
    class  (galacticFilterClass      ), pointer       :: galacticFilterSatellite_, galacticFilterCentral_
    
    !![
    <inputParameter>
      <name>countMergersMaximum</name>
      <defaultValue>huge(0_c_size_t)</defaultValue>
      <description>The maximum number of mergers to accumulate for each node. Defaults to the maximum possible.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="galacticFilter" parameterName="galacticFilterSatellite" name="galacticFilterSatellite_" source="parameters">
     <default>
      <galacticFilterSatellite value="always"/>
     </default>
    </objectBuilder>
    <objectBuilder class="galacticFilter" parameterName="galacticFilterCentral"   name="galacticFilterCentral_"   source="parameters">
     <default>
      <galacticFilterCentral value="always"/>
     </default>
    </objectBuilder>
    !!]
    self=nodeOperatorGalaxyMergers(countMergersMaximum,galacticFilterSatellite_,galacticFilterCentral_)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function galaxyMergersConstructorParameters

  function galaxyMergersConstructorInternal(countMergersMaximum,galacticFilterSatellite_,galacticFilterCentral_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorGalaxyMergers} node operator class.
    !!}
    use :: Galacticus_Nodes, only : defaultBasicComponent
    implicit none
    type   (nodeOperatorGalaxyMergers)                        :: self
    integer(c_size_t                 ), intent(in   )         :: countMergersMaximum
    class  (galacticFilterClass      ), intent(in   ), target :: galacticFilterSatellite_, galacticFilterCentral_
    !![
    <constructorAssign variables="countMergersMaximum, *galacticFilterSatellite_, *galacticFilterCentral_"/>
    !!]
    
    !![
    <addMetaProperty component="basic" name="galaxyMergerSatelliteIndex"       id="self%galaxyMergerSatelliteIndexID"       type="longInteger" rank="1" isCreator="yes"/>
    <addMetaProperty component="basic" name="galaxyMergersTime"                id="self%galaxyMergerTimeID"                                    rank="1" isCreator="yes"/>
    <addMetaProperty component="basic" name="galaxyMergerHostMassGas"          id="self%galaxyMergerHostMassGasID"                             rank="1" isCreator="yes"/>
    <addMetaProperty component="basic" name="galaxyMergerSatelliteMassGas"     id="self%galaxyMergerSatelliteMassGasID"                        rank="1" isCreator="yes"/>
    <addMetaProperty component="basic" name="galaxyMergerHostMassStellar"      id="self%galaxyMergerHostMassStellarID"                         rank="1" isCreator="yes"/>
    <addMetaProperty component="basic" name="galaxyMergerSatelliteMassStellar" id="self%galaxyMergerSatelliteMassStellarID"                    rank="1" isCreator="yes"/>
    !!]
    return
  end function galaxyMergersConstructorInternal

  subroutine galaxyMergersAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent, openMPThreadBindingAtLevel, dependencyDirectionAfter, dependencyRegEx
    implicit none
    class(nodeOperatorGalaxyMergers), intent(inout) :: self
    type (dependencyRegEx          ), dimension(1)  :: dependenciesSatelliteMerger 

    dependenciesSatelliteMerger(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
    call satelliteMergerEvent%attach(self,satelliteMerger,openMPThreadBindingAtLevel,label='preAnalysis:galaxyMergers',dependencies=dependenciesSatelliteMerger)
    return
  end subroutine galaxyMergersAutoHook
  
  subroutine galaxyMergersDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorGalaxyMergers} node operator class.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent
    implicit none
    type(nodeOperatorGalaxyMergers), intent(inout) :: self

    if (satelliteMergerEvent%isAttached(self,satelliteMerger)) call satelliteMergerEvent%detach(self,satelliteMerger)
    !![
    <objectDestructor name="self%galacticFilterSatellite_"/>
    <objectDestructor name="self%galacticFilterCentral_" />
    !!]
    return
  end subroutine galaxyMergersDestructor

  subroutine satelliteMerger(self,node)
    !!{
    Record galaxy-galaxy mergers.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDisk, nodeComponentSpheroid
    implicit none
    class           (*                    ), intent(inout)              :: self
    type            (treeNode             ), intent(inout), target      :: node
    type            (treeNode             )               , pointer     :: nodeHost
    class           (nodeComponentBasic   )               , pointer     :: basicHost
    class           (nodeComponentDisk    )               , pointer     :: disk                       , diskHost
    class           (nodeComponentSpheroid)               , pointer     :: spheroid                   , spheroidHost
    double precision                       , dimension(:) , allocatable :: massGasSatelliteCurrent    , massGasSatelliteNew    , &
         &                                                                 massGasHostCurrent         , massGasHostNew         , &
         &                                                                 massStellarSatelliteCurrent, massStellarSatelliteNew, &
         &                                                                 massStellarHostCurrent     , massStellarHostNew     , &
         &                                                                 timeCurrent                , timeNew
    integer         (kind_int8            ), dimension(:), allocatable  :: indexSatelliteCurrent      , indexSatelliteNew
    double precision                                                    :: massGasSatellite           , massGasHost            , &
         &                                                                 massStellarSatellite       , massStellarHost        , &
         &                                                                 time
    integer         (kind_int8            )                             :: indexSatellite
    
    select type (self)
    class is (nodeOperatorGalaxyMergers)
       nodeHost      =>  node    %mergesWith()
       if (.not.(self%galacticFilterSatellite_%passes(node).and.self%galacticFilterCentral_%passes(nodeHost))) return
       disk                 =>  node     %disk       ()
       spheroid             =>  node     %spheroid   ()
       diskHost             =>  nodeHost %disk       ()
       spheroidHost         =>  nodeHost %spheroid   ()
       basicHost            =>  nodeHost %basic      ()
       indexSatellite       =   node     %index      ()
       time                 =   basicHost%time       ()
       massGasSatellite     =  +disk     %massGas    ()+spheroid    %massGas    ()
       massGasHost          =  +diskHost %massGas    ()+spheroidHost%massGas    ()
       massStellarSatellite =  +disk     %massStellar()+spheroid    %massStellar()
       massStellarHost      =  +diskHost %massStellar()+spheroidHost%massStellar()
       ! Append to the list of mergers.
       indexSatelliteCurrent      =basicHost%longIntegerRank1MetaPropertyGet(self%galaxyMergerSatelliteIndexID      )
       timeCurrent                =basicHost%      floatRank1MetaPropertyGet(self%galaxyMergerTimeID                )
       massGasSatelliteCurrent    =basicHost%      floatRank1MetaPropertyGet(self%galaxyMergerSatelliteMassGasID    )
       massGasHostCurrent         =basicHost%      floatRank1MetaPropertyGet(self%galaxyMergerHostMassGasID         )
       massStellarSatelliteCurrent=basicHost%      floatRank1MetaPropertyGet(self%galaxyMergerSatelliteMassStellarID)
       massStellarHostCurrent     =basicHost%      floatRank1MetaPropertyGet(self%galaxyMergerHostMassStellarID     )
       if (size(timeCurrent) < self%countMergersMaximum) then
          allocate(indexSatelliteNew      (size(indexSatelliteCurrent      )+1_c_size_t))
          allocate(timeNew                (size(timeCurrent                )+1_c_size_t))
          allocate(massGasSatelliteNew    (size(massGasSatelliteCurrent    )+1_c_size_t))
          allocate(massGasHostNew         (size(massGasHostCurrent         )+1_c_size_t))
          allocate(massStellarSatelliteNew(size(massStellarSatelliteCurrent)+1_c_size_t))
          allocate(massStellarHostNew     (size(massStellarHostCurrent     )+1_c_size_t))
          indexSatelliteNew      (1_c_size_t:size(indexSatelliteCurrent      )           )=indexSatelliteCurrent      (          :                                 )
          timeNew                (1_c_size_t:size(timeCurrent                )           )=timeCurrent                (          :                                 )
          massGasSatelliteNew    (1_c_size_t:size(massGasSatelliteCurrent    )           )=massGasSatelliteCurrent    (          :                                 )
          massGasHostNew         (1_c_size_t:size(massGasHostCurrent         )           )=massGasHostCurrent         (          :                                 )
          massStellarSatelliteNew(1_c_size_t:size(massStellarSatelliteCurrent)           )=massStellarSatelliteCurrent(          :                                 )
          massStellarHostNew     (1_c_size_t:size(massStellarHostCurrent     )           )=massStellarHostCurrent     (          :                                 )
       else
          allocate(indexSatelliteNew      (size(indexSatelliteCurrent      )           ))
          allocate(timeNew                (size(timeCurrent                )           ))
          allocate(massGasSatelliteNew    (size(massGasSatelliteCurrent    )           ))
          allocate(massGasHostNew         (size(massGasHostCurrent         )           ))
          allocate(massStellarSatelliteNew(size(massStellarSatelliteCurrent)           ))
          allocate(massStellarHostNew     (size(massStellarHostCurrent     )           ))
          indexSatelliteNew      (1_c_size_t:size(indexSatelliteCurrent      )-1_c_size_t)=indexSatelliteCurrent      (2_c_size_t:size(indexSatelliteCurrent      ))
          timeNew                (1_c_size_t:size(timeCurrent                )-1_c_size_t)=timeCurrent                (2_c_size_t:size(timeCurrent                ))
          massGasSatelliteNew    (1_c_size_t:size(massGasSatelliteCurrent    )-1_c_size_t)=massGasSatelliteCurrent    (2_c_size_t:size(massGasSatelliteCurrent    ))
          massGasHostNew         (1_c_size_t:size(massGasHostCurrent         )-1_c_size_t)=massGasHostCurrent         (2_c_size_t:size(massGasHostCurrent         ))
          massStellarSatelliteNew(1_c_size_t:size(massStellarSatelliteCurrent)-1_c_size_t)=massStellarSatelliteCurrent(2_c_size_t:size(massStellarSatelliteCurrent))
          massStellarHostNew     (1_c_size_t:size(massStellarHostCurrent     )-1_c_size_t)=massStellarHostCurrent     (2_c_size_t:size(massStellarHostCurrent     ))
       end if
       indexSatelliteNew      (size(indexSatelliteNew      ))=indexSatellite
       timeNew                (size(timeNew                ))=time
       massGasSatelliteNew    (size(massGasSatelliteNew    ))=massGasSatellite
       massGasHostNew         (size(massGasHostNew         ))=massGasHost
       massStellarSatelliteNew(size(massStellarSatelliteNew))=massStellarSatellite
       massStellarHostNew     (size(massStellarHostNew     ))=massStellarHost
       call basicHost%longIntegerRank1MetaPropertySet(self%galaxyMergerSatelliteIndexID      ,indexSatelliteNew      )
       call basicHost%      floatRank1MetaPropertySet(self%galaxyMergerTimeID                ,timeNew                )
       call basicHost%      floatRank1MetaPropertySet(self%galaxyMergerSatelliteMassGasID    ,massGasSatelliteNew    )
       call basicHost%      floatRank1MetaPropertySet(self%galaxyMergerHostMassGasID         ,massGasHostNew         )
       call basicHost%      floatRank1MetaPropertySet(self%galaxyMergerSatelliteMassStellarID,massStellarSatelliteNew)
       call basicHost%      floatRank1MetaPropertySet(self%galaxyMergerHostMassStellarID     ,massStellarHostNew     )
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine satelliteMerger
