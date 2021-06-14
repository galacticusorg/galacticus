!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

  !% Implements a node operator class that triggers merging of satellites based on their orbital radius.

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !# <nodeOperator name="nodeOperatorSatelliteMergingRadiusTrigger">
  !#  <description>A node operator class that triggers merging of satellites based on their orbital radius.</description>
  !# </nodeOperator>
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteMergingRadiusTrigger
     !% A node operator class that triggers merging of satellites based on their orbital radius.
     private
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     double precision                                    :: radiusVirialFraction
   contains
     !# <methods>
     !#   <method description="Compute the radius at which the satellite will be merged." method="radiusMerge" />
     !# </methods>
     final     ::                          satelliteMergingRadiusTriggerDestructor
     procedure :: differentialEvolution => satelliteMergingRadiusTriggerDifferentialEvolution
     procedure :: radiusMerge           => satelliteMergingRadiusTriggerRadiusMerge
  end type nodeOperatorSatelliteMergingRadiusTrigger
  
  interface nodeOperatorSatelliteMergingRadiusTrigger
     !% Constructors for the {\normalfont \ttfamily satelliteMergingRadiusTrigger} node operator class.
     module procedure satelliteMergingRadiusTriggerConstructorParameters
     module procedure satelliteMergingRadiusTriggerConstructorInternal
  end interface nodeOperatorSatelliteMergingRadiusTrigger
  
contains

  function satelliteMergingRadiusTriggerConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily satelliteMergingRadiusTrigger} node operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorSatelliteMergingRadiusTrigger)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass                 ), pointer       :: darkMatterHaloScale_
    double precision                                                           :: radiusVirialFraction

    !# <inputParameter>
    !#   <name>radiusVirialFraction</name>
    !#   <defaultValue>0.01d0</defaultValue>
    !#   <description>The fraction of the virial radius below which satellites are merged.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=nodeOperatorSatelliteMergingRadiusTrigger(radiusVirialFraction,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterHaloScale_"/>
    return
  end function satelliteMergingRadiusTriggerConstructorParameters

  function satelliteMergingRadiusTriggerConstructorInternal(radiusVirialFraction,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily satelliteMergingRadiusTrigger} node operator class.
    implicit none
    type            (nodeOperatorSatelliteMergingRadiusTrigger)                        :: self
    double precision                                           , intent(in   )         :: radiusVirialFraction
    class           (darkMatterHaloScaleClass                 ), intent(in   ), target :: darkMatterHaloScale_
    !# <constructorAssign variables="radiusVirialFraction, *darkMatterHaloScale_"/>
    
    return
  end function satelliteMergingRadiusTriggerConstructorInternal
  
  subroutine satelliteMergingRadiusTriggerDestructor(self)
    !% Destructor for the {\normalfont \ttfamily satelliteMergingRadiusTrigger} node operator class.
    implicit none
    type(nodeOperatorSatelliteMergingRadiusTrigger), intent(inout) :: self
    
    !# <objectDestructor name="self%darkMatterHaloScale_"/>
    return
  end subroutine satelliteMergingRadiusTriggerDestructor
  
  subroutine satelliteMergingRadiusTriggerDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !% Trigger merging of a satellite halo based on its orbital radius.
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    use :: Vectors         , only : Vector_Magnitude
    implicit none
    class           (nodeOperatorSatelliteMergingRadiusTrigger), intent(inout), target  :: self
    type            (treeNode                                 ), intent(inout)          :: node
    logical                                                    , intent(inout)          :: interrupt
    procedure       (interruptTask                            ), intent(inout), pointer :: functionInterrupt
    integer                                                    , intent(in   )          :: propertyType
    class           (nodeComponentSatellite                   )               , pointer :: satellite
    double precision                                           , dimension(3)           :: position
    double precision                                                                    :: radius
    !$GLC attributes unused :: propertyType
    
    if (.not.node%isSatellite()) return
    satellite => node     %satellite(        )
    position  =  satellite%position (        )
    radius    =  Vector_Magnitude   (position)
    ! Test for merging.
    if     (                                 &
         &   radius > 0.0d0                  &
         &  .and.                            &
         &   radius < self%radiusMerge(node) &
         & ) then
       ! Merging criterion met - trigger an interrupt.
       interrupt         =  .true.
       functionInterrupt => mergerTrigger
    end if
    return
  end subroutine satelliteMergingRadiusTriggerDifferentialEvolution

  subroutine mergerTrigger(node)
    !% Trigger a merger of the satellite by setting the time until merging to zero.
    use :: Galacticus_Nodes, only : nodeComponentSatellite, treeNode
    implicit none
    type (treeNode              ), intent(inout), target  :: node
    class(nodeComponentSatellite)               , pointer :: satellite

    satellite => node%satellite()
    call satellite%mergeTimeSet(0.0d0)
    return
  end subroutine mergerTrigger

  double precision function satelliteMergingRadiusTriggerRadiusMerge(self,node)
    !% Compute the merging radius for a node.
    use :: Galacticus_Nodes                  , only : treeNode
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass, Galactic_Structure_Radius_Enclosing_Mass
    use :: Galactic_Structure_Options        , only : massTypeGalactic                , radiusLarge
    implicit none
    class           (nodeOperatorSatelliteMergingRadiusTrigger), intent(inout) :: self
    type            (treeNode                                 ), intent(inout) :: node
    type            (treeNode                                 ), pointer       :: nodeHost
    double precision                                                           :: radiusHalfMassCentral, radiusHalfMassSatellite

    ! Find the host node.
    nodeHost => node%mergesWith()
    ! Get half-mass radii of central and satellite galaxies. We first check that the total mass in the galactic component
    ! (found by setting the radius to "radiusLarge") is non-zero as we do not want to attempt to find the half-mass radius
    ! of the galactic component, if no galactic component exists. To correctly handle the case that numerical errors lead
    ! to a zero-size galactic component (the enclosed mass within zero radius is non-zero and equals to the total mass of
    ! this component), we do a further check that the enclosed mass within zero radius is smaller than half of the total
    ! mass in the galactic component.
    if     (                                                                                             &
         &   Galactic_Structure_Enclosed_Mass(nodeHost,massType=massTypeGalactic,radius=radiusLarge)     &
         &   >                                                                                           &
         &   max(                                                                                        &
         &       0.0d0,                                                                                  &
         &       2.0d0*Galactic_Structure_Enclosed_Mass(nodeHost,massType=massTypeGalactic,radius=0.0d0) &
         &      )                                                                                        &
         & ) then
       radiusHalfMassCentral  =Galactic_Structure_Radius_Enclosing_Mass(nodeHost,fractionalMass=0.5d0,massType=massTypeGalactic)
    else
       radiusHalfMassCentral  =0.0d0
    end if
    if     (                                                                                             &
         &   Galactic_Structure_Enclosed_Mass(node    ,massType=massTypeGalactic,radius=radiusLarge)     &
         &   >                                                                                           &
         &   max(                                                                                        &
         &       0.0d0,                                                                                  &
         &       2.0d0*Galactic_Structure_Enclosed_Mass(node    ,massType=massTypeGalactic,radius=0.0d0) &
         &      )                                                                                        &
         & ) then
       radiusHalfMassSatellite=Galactic_Structure_Radius_Enclosing_Mass(node    ,fractionalMass=0.5d0,massType=massTypeGalactic)
    else
       radiusHalfMassSatellite=0.0d0
    end if    
    satelliteMergingRadiusTriggerRadiusMerge=max(                                                              &
         &                                       +                          radiusHalfMassSatellite            &
         &                                       +                          radiusHalfMassCentral            , &
         &                                       +self%                     radiusVirialFraction               &
         &                                       *self%darkMatterHaloScale_%virialRadius           (nodeHost)  &
         &                                      )
    return
  end function satelliteMergingRadiusTriggerRadiusMerge
