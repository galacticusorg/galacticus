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

  !+    Contributions to this file made by: Claude.

  !!{RST
  Implements a node operator class that triggers merging of satellites based on their orbital radius.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <nodeOperator name="nodeOperatorSatelliteMergingRadiusTrigger" docformat="rst">
    <description>
    A node operator class that merges satellite halos with their central when the orbital radius falls below the larger of ``radiusVirialFraction`` times the host virial radius (default 0.01), and ``radiusHalfMassFraction`` times the sum of the central and satellite galactic half-mass radii (default 1.0). Optionally records Keplerian orbital elements of merged subhalos when ``recordMergedSubhaloProperties`` is true; ``recordFirstLevelOnly`` restricts recording to first-level subhalos relative to the final host.
    </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorSatelliteMergingRadius) :: nodeOperatorSatelliteMergingRadiusTrigger
     !!{RST
     A node operator class that triggers merging of satellites based on their orbital radius.
     !!}
     private
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     double precision                                    :: radiusVirialFraction          , radiusHalfMassFraction
   contains
     final     ::                satelliteMergingRadiusTriggerDestructor
     procedure :: radiusMerge => satelliteMergingRadiusTriggerRadiusMerge
  end type nodeOperatorSatelliteMergingRadiusTrigger
  
  interface nodeOperatorSatelliteMergingRadiusTrigger
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorSatelliteMergingRadiusTrigger` node operator class.
     !!}
     module procedure satelliteMergingRadiusTriggerConstructorParameters
     module procedure satelliteMergingRadiusTriggerConstructorInternal
  end interface nodeOperatorSatelliteMergingRadiusTrigger

  
contains

  function satelliteMergingRadiusTriggerConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorSatelliteMergingRadiusTrigger` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorSatelliteMergingRadiusTrigger)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass                 ), pointer       :: darkMatterHaloScale_
    double precision                                                           :: radiusVirialFraction         , radiusHalfMassFraction
    logical                                                                    :: recordMergedSubhaloProperties, recordFirstLevelOnly

    !![
    <inputParameter docformat="rst">
      <name>radiusHalfMassFraction</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The fraction of the sum of the central and satellite half-mass radii below which satellites are merged.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>radiusVirialFraction</name>
      <defaultValue>0.01d0</defaultValue>
      <description>
      The fraction of the virial radius below which satellites are merged.
      </description>
      <source>parameters</source>
    </inputParameter>
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
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodeOperatorSatelliteMergingRadiusTrigger(radiusHalfMassFraction,radiusVirialFraction,recordMergedSubhaloProperties,recordFirstLevelOnly,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function satelliteMergingRadiusTriggerConstructorParameters

  function satelliteMergingRadiusTriggerConstructorInternal(radiusHalfMassFraction,radiusVirialFraction,recordMergedSubhaloProperties,recordFirstLevelOnly,darkMatterHaloScale_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorSatelliteMergingRadiusTrigger` node operator class.
    !!}
    implicit none
    type            (nodeOperatorSatelliteMergingRadiusTrigger)                        :: self
    double precision                                           , intent(in   )         :: radiusVirialFraction         , radiusHalfMassFraction
    logical                                                    , intent(in   )         :: recordMergedSubhaloProperties, recordFirstLevelOnly
    class           (darkMatterHaloScaleClass                 ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="radiusHalfMassFraction, radiusVirialFraction, recordMergedSubhaloProperties, recordFirstLevelOnly, *darkMatterHaloScale_"/>
    !!]
    
    call self%recordingInitialize()
    return
  end function satelliteMergingRadiusTriggerConstructorInternal
  
  subroutine satelliteMergingRadiusTriggerDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorSatelliteMergingRadiusTrigger` node operator class.
    !!}
    implicit none
    type(nodeOperatorSatelliteMergingRadiusTrigger), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine satelliteMergingRadiusTriggerDestructor

  double precision function satelliteMergingRadiusTriggerRadiusMerge(self,node)
    !!{RST
    Compute the merging radius for a node.
    !!}
    use :: Galacticus_Nodes          , only : treeNode
    use :: Galactic_Structure_Options, only : massTypeGalactic
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (nodeOperatorSatelliteMergingRadiusTrigger), intent(inout) :: self
    type            (treeNode                                 ), intent(inout) :: node
    type            (treeNode                                 ), pointer       :: nodeHost
    class           (massDistributionClass                    ), pointer       :: massDistribution_    , massDistributionHost_
    double precision                                                           :: radiusHalfMassCentral, radiusHalfMassSatellite

    ! Find the host node.
    nodeHost => node%mergesWith()
    ! Get mass distributions.
    massDistribution_     => node    %massDistribution(massType=massTypeGalactic)
    massDistributionHost_ => nodeHost%massDistribution(massType=massTypeGalactic)
    ! Get half-mass radii of central and satellite galaxies. We first check that the total mass in the galactic component is
    ! non-zero as we do not want to attempt to find the half-mass radius of the galactic component, if no galactic component
    ! exists. To correctly handle the case that numerical errors lead to a zero-size galactic component (the enclosed mass
    ! within zero radius is non-zero and equals to the total mass of this component), we do a further check that the enclosed
    ! mass within zero radius is smaller than half of the total mass in the galactic component.
    if     (                                                                                    &
         &             massDistributionHost_%massTotal()                                        &
         &   >                                                                                  &
         &   max(                                                                               &
         &       0.0d0,                                                                         &
         &       2.0d0*massDistributionHost_%massEnclosedBySphere(radius=0.0d0)                 &
         &      )                                                                               &
         & ) then
       radiusHalfMassCentral  =massDistributionHost_%radiusEnclosingMass(massFractional=0.5d0)
    else
       radiusHalfMassCentral  =0.0d0
    end if
    if     (                                                                                    &
         &             massDistribution_    %massTotal()                                        &
         &   >                                                                                  &
         &   max(                                                                               &
         &       0.0d0,                                                                         &
         &       2.0d0*massDistribution_    %massEnclosedBySphere(radius=0.0d0)                 &
         &      )                                                                               &
         & ) then
       radiusHalfMassSatellite=massDistribution_    %radiusEnclosingMass(massFractional=0.5d0)
    else
       radiusHalfMassSatellite=0.0d0
    end if
    !![
    <objectDestructor name="massDistribution_"    />
    <objectDestructor name="massDistributionHost_"/>
    !!]
    satelliteMergingRadiusTriggerRadiusMerge=max(                                                                &
         &                                       +  self%                     radiusHalfMassFraction             &
         &                                       *(                                                              &
         &                                         +                          radiusHalfMassSatellite            &
         &                                         +                          radiusHalfMassCentral              &
         &                                        )                                                            , &
         &                                         +self%                     radiusVirialFraction               &
         &                                         *self%darkMatterHaloScale_%radiusVirial           (nodeHost)  &
         &                                      )
    return
  end function satelliteMergingRadiusTriggerRadiusMerge
