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

  !+    Contributions to this file made by: Yu Zhao, Claude

  !!{RST
  Implements a node operator class that triggers merging of satellites based on their orbital radius.
  !!}

  !![
  <nodeOperator name="nodeOperatorSatelliteMergingSoliton" docformat="rst">
   <description>
   A node operator class for fuzzy dark matter (FDM/wave dark matter) models that triggers satellite merging when the orbital radius falls below the sum of the soliton core radii of the host and satellite halos. ``recordMergedSubhaloProperties`` optionally records Keplerian orbital elements of merged subhalos; ``recordFirstLevelOnly`` restricts recording to first-level subhalos.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorSatelliteMergingRadius) :: nodeOperatorSatelliteMergingSoliton
     !!{RST
     A node operator that triggers satellite merging in FDM models when the orbital radius falls below the sum of the soliton core radii of the host and satellite.
     !!}
     private
     integer :: radiusCoreID  , massCoreID      , &
                randomOffsetID, massCoreNormalID
   contains
     final     ::                satelliteMergingSolitonDestructor
     procedure :: radiusMerge => satelliteMergingSolitonRadiusMerge
     procedure :: autoHook    => satelliteMergingSolitonAutoHook
  end type nodeOperatorSatelliteMergingSoliton
  
  interface nodeOperatorSatelliteMergingSoliton
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorSatelliteMergingSoliton` node operator class.
     !!}
     module procedure satelliteMergingSolitonConstructorParameters
     module procedure satelliteMergingSolitonConstructorInternal
  end interface nodeOperatorSatelliteMergingSoliton

  
contains

  function satelliteMergingSolitonConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorSatelliteMergingSoliton` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorSatelliteMergingSoliton)                :: self
    type   (inputParameters                    ), intent(inout) :: parameters
    logical                                                     :: recordMergedSubhaloProperties, recordFirstLevelOnly

    !![
    <inputParameter docformat="rst">
      <name>recordMergedSubhaloProperties</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, record the orbital properties of subhalo that merge.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>recordFirstLevelOnly</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, record only mergers with first-level subhalos relative to the host.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodeOperatorSatelliteMergingSoliton(recordMergedSubhaloProperties,recordFirstLevelOnly)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function satelliteMergingSolitonConstructorParameters

  function satelliteMergingSolitonConstructorInternal(recordMergedSubhaloProperties,recordFirstLevelOnly) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorSatelliteMergingSoliton` node operator class.
    !!}
    implicit none
    type   (nodeOperatorSatelliteMergingSoliton)                :: self
    logical                                     , intent(in   ) :: recordMergedSubhaloProperties, recordFirstLevelOnly
    !![
    <constructorAssign variables="recordMergedSubhaloProperties, recordFirstLevelOnly"/>
    !!]
    
    !![
    <addMetaProperty component="darkMatterProfile" name="solitonRandomOffset"   id="self%randomOffsetID"   isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonRadiusCore"     id="self%radiusCoreID"     isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonMassCore"       id="self%massCoreID"       isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonMassCoreNormal" id="self%massCoreNormalID" isEvolvable="yes" isCreator="no"/>
    !!]
    
    call self%recordingInitialize()
    return
  end function satelliteMergingSolitonConstructorInternal

  subroutine satelliteMergingSolitonAutoHook(self)
    !!{RST
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent, openMPThreadBindingAtLevel, dependencyDirectionAfter, dependencyRegEx
    implicit none
    class(nodeOperatorSatelliteMergingSoliton), intent   (inout) :: self
    type (dependencyRegEx                    ), dimension(    1) :: dependenciesSatelliteMerger

    dependenciesSatelliteMerger(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
    call satelliteMergerEvent%attach(self,satelliteMerger,openMPThreadBindingAtLevel,label='satelliteMergingSoliton',dependencies=dependenciesSatelliteMerger)
    return
  end subroutine satelliteMergingSolitonAutoHook

  subroutine satelliteMergingSolitonDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorSatelliteMergingSoliton` node operator class.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent
    implicit none
    type(nodeOperatorSatelliteMergingSoliton), intent(inout) :: self

    if (satelliteMergerEvent%isAttached(self,satelliteMerger)) call satelliteMergerEvent%detach(self,satelliteMerger)
    return
  end subroutine satelliteMergingSolitonDestructor

  double precision function satelliteMergingSolitonRadiusMerge(self,node) result(radiusMerge)
    !!{RST
    Compute the merging radius for a node.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentSatellite, nodeComponentDarkMatterProfile, treeNode
    use :: Galactic_Structure_Options, only : massTypeGalactic
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (nodeOperatorSatelliteMergingSoliton), intent(inout) :: self
    class           (nodeComponentDarkMatterProfile     ), pointer       :: darkMatterProfile, darkMatterProfileHost
    type            (treeNode                           ), intent(inout) :: node
    type            (treeNode                           ), pointer       :: nodeHost
    double precision                                                     :: radiusCoreHost   , radiusCoreSatellite

    ! Find the host node profile.
    nodeHost              => node    %mergesWith       ()
    darkMatterProfileHost => nodeHost%darkMatterProfile()
    ! Get the satellite profile.
    darkMatterProfile     => node    %darkMatterProfile()
    ! Compute the merging radius.
    radiusCoreHost     =+darkMatterProfileHost%floatRank0MetaPropertyGet(self%radiusCoreID)
    radiusCoreSatellite=+darkMatterProfile    %floatRank0MetaPropertyGet(self%radiusCoreID)
    radiusMerge        =+radiusCoreHost      &
         &              +radiusCoreSatellite
    return
  end function satelliteMergingSolitonRadiusMerge

  subroutine satelliteMerger(self,node)
    !!{RST
    Merge the solitonic cores of the satellite and host halos.
    !!}
    use :: Error             , only : Error_Report
    use :: Function_Classes  , only : functionClass
    use :: Galacticus_Nodes  , only : nodeComponentSatellite, nodeComponentDarkMatterProfile, treeNode
    use :: ISO_Varying_String, only : char
    implicit none
    class           (*                             ), intent(inout)         :: self
    type            (treeNode                      ), intent(inout), target :: node
    type            (treeNode                      ), pointer               :: nodeHost
    class           (nodeComponentDarkMatterProfile), pointer               :: darkMatterProfile         , darkMatterProfileHost
    double precision                                , parameter             :: fractionMassRetained=0.7d0
    double precision                                                        :: massCoreHost              , massCoreSatellite      , &
            &                                                                  massCoreNormalHost        , massCoreNormalSatellite

    select type (self)
    class is (nodeOperatorSatelliteMergingSoliton)
       ! Find the host node profile.
       nodeHost              => node    %mergesWith       ()
       darkMatterProfileHost => nodeHost%darkMatterProfile()
       ! Get the satellite profile.
       darkMatterProfile     => node    %darkMatterProfile()
       ! Compute the new core mass.
       massCoreNormalHost     = darkMatterProfileHost%floatRank0MetaPropertyGet(self%massCoreNormalID)
       massCoreNormalSatellite= darkMatterProfile    %floatRank0MetaPropertyGet(self%massCoreNormalID)
       massCoreHost           = darkMatterProfileHost%floatRank0MetaPropertyGet(self%massCoreID      )
       massCoreSatellite      = darkMatterProfile    %floatRank0MetaPropertyGet(self%massCoreID      )
       call darkMatterProfileHost%floatRank0MetaPropertySet(                                                                    &
            &                                                self%massCoreNormalID                                            , &
            &                                               +fractionMassRetained*(massCoreNormalHost+massCoreNormalSatellite)  &
            &                                              )
       call darkMatterProfileHost%floatRank0MetaPropertySet(                                                                    &
            &                                                self%massCoreID                                                  , &
            &                                               +fractionMassRetained*(massCoreHost      +massCoreSatellite      )  &
            &                                              )
       call darkMatterProfileHost%floatRank0MetaPropertySet(                                                                    &
            &                                                self%randomOffsetID                                              , &
            &                                               +log10(                                                             &
            &                                                      +(massCoreHost      +massCoreSatellite      )                &
            &                                                      /(massCoreNormalHost+massCoreNormalSatellite)                &
            &                                                     )                                                             &
            &                                              )
    class is (functionClass)
       call Error_Report('object is not of [nodeOperatorSatelliteMergingSoliton] class, but of ['//char(self%objectType())//'] class'//{introspection:location})
    class default
       call Error_Report('object is not of [nodeOperatorSatelliteMergingSoliton] class'//{introspection:location})
    end select
    return
  end subroutine satelliteMerger


