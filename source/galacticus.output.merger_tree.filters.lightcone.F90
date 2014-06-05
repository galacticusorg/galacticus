!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which filters output for lightcone geometry.

module Galacticus_Merger_Tree_Output_Filter_Lightcones
  !% Filters output for lightcone geometry.
  implicit none
  private
  public :: Galacticus_Merger_Tree_Output_Filter_Lightcone,Galacticus_Output_Tree_Lightcone,&
       & Galacticus_Output_Tree_Lightcone_Property_Count, Galacticus_Output_Tree_Lightcone_Names,&
       & Galacticus_Merger_Tree_Output_Filter_Lightcone_Initialize

  ! Number of lightcone properties.
  integer                                      , parameter :: lightconePropertyCount    =8

  ! Flags indicating if the module has been initialized and if this filter is active.
  logical                                                  :: lightconeFilterInitialized=.false.
  logical                                                  :: lightconeFilterActive

  ! Variables used to describe the lightcone geometry.
  double precision             , dimension(3  )            :: lightconeOrigin
  double precision             , dimension(3,3)            :: lightconeUnitVector
  double precision, allocatable, dimension(:)              :: lightconeMaximumDistance          , lightconeMinimumDistance        , &
       &                                                      lightconeTime
  double precision                                         :: lightconeFieldOfViewLength        , lightconeFieldOfViewTanHalfAngle, &
       &                                                      lightconeReplicationPeriod
  double precision                                         :: angularWeight
  integer                                                  :: lightconeGeometry
  integer                                      , parameter :: lightconeGeometrySquare   =1

  ! Variables that store the position, velocity and redshift of a galaxy that is found to be in the lightcone. These are then used
  ! to output that position.
  double precision             , dimension(3  )            :: lightconePosition                 , lightconeVelocity
  double precision                                         :: lightconeRedshift
  !$omp threadprivate(lightconePosition,lightconeVelocity,lightconeRedshift)
contains

  !# <mergerTreeOutputFilterInitialize>
  !#   <unitName>Galacticus_Merger_Tree_Output_Filter_Lightcone_Initialize</unitName>
  !# </mergerTreeOutputFilterInitialize>
  subroutine Galacticus_Merger_Tree_Output_Filter_Lightcone_Initialize(filterNames)
    !% Initializes the lightcone filter module.
    use ISO_Varying_String
    use Input_Parameters
    use FoX_dom
    use IO_XML
    use Cosmology_Functions
    use Numerical_Constants_Astronomical
    use Cosmological_Parameters
    use Galacticus_Error
    use Memory_Management
    use Vectors
    implicit none
    type            (varying_string), dimension(:), intent(in   ) :: filterNames
    type            (Node          ), pointer                     :: doc                            , thisItem
    type            (NodeList      ), pointer                     :: itemList
    integer                                                       :: iAxis                          , iOutput                     , &
         &                                                           ioErr                          , lengthUnitsHubbleExponent
    double precision                                              :: lengthUnitsInSI                , lightconeMaximumDistanceTemp, &
         &                                                           lightconeMinimumDistanceTemp   , lightconeTimeTemp           , &
         &                                                           unitConversionLength
    type            (varying_string)                              :: filterLightconeGeometryFileName
    logical                                                       :: filterLightconeFixedTime
    character       (len=11        )                              :: tagName
    character       (len=10        )                              :: geometryLabel

    ! Initialize the filter if necessary.
    if (.not.lightconeFilterInitialized) then
       ! Determine if this filter has been selected.
       lightconeFilterActive=any(filterNames == "lightcone")
       ! If this filter is active, read the geometric constraints.
       if (lightconeFilterActive) then
          ! Get name of lightcone geometry specification file.
          !@ <inputParameter>
          !@   <name>filterLightconeGeometryFileName</name>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The name of an XML file from which to read details of lightcone geometry.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('filterLightconeGeometryFileName',filterLightconeGeometryFileName)
          ! Get option controlling whether output should occur at a fixed time.
          !@ <inputParameter>
          !@   <name>filterLightconeFixedTime</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Specifies if lightcone output should occur at a fixed time (as opposed to the usual case where time evolves along the lightcone). Intended for the construction of lightcones with no evolution.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('filterLightconeFixedTime',filterLightconeFixedTime,defaultValue=.false.)

          ! Extract data from geometry specification file.
          !$omp critical (FoX_DOM_Access)
          doc => parseFile(char(filterLightconeGeometryFileName),iostat=ioErr)
          if (ioErr /= 0) call Galacticus_Error_Report('Galacticus_Merger_Tree_Output_Filter_Lightcone','Unable to parse geometry file')
          ! Get the unit vectors of the lightcone.
          do iAxis=1,3
             write (tagName,'(a,i1)') "unitVector",iAxis
             thisItem => XML_Get_First_Element_By_Tag_Name(doc,trim(tagName))
             lightconeUnitVector(:,iAxis)=Filter_Lightcone_Get_Coordinates(thisItem)
             lightconeUnitVector(:,iAxis)=lightconeUnitVector(:,iAxis)/Vector_Magnitude(lightconeUnitVector(:,iAxis))
          end do
          ! Get units information.
          thisItem => XML_Get_First_Element_By_Tag_Name(doc     ,"units/length/unitsInSI"     )
          call extractDataContent(thisItem,lengthUnitsInSI          )
          thisItem => XML_Get_First_Element_By_Tag_Name(doc     ,"units/length/hubbleExponent")
          call extractDataContent(thisItem,lengthUnitsHubbleExponent)
          unitConversionLength=lengthUnitsInSI*(Little_H_0()**lengthUnitsHubbleExponent)/megaParsec
          ! Get the origin of the lightcone.
          thisItem => XML_Get_First_Element_By_Tag_Name(doc,"origin")
          lightconeOrigin=Filter_Lightcone_Get_Coordinates(thisItem)*unitConversionLength
          ! Get the replication period.
          thisItem => XML_Get_First_Element_By_Tag_Name(doc,"boxLength")
          call extractDataContent(thisItem,lightconeReplicationPeriod)
          lightconeReplicationPeriod=lightconeReplicationPeriod*unitConversionLength
          ! Get the geometry.
          thisItem => XML_Get_First_Element_By_Tag_Name(doc     ,"fieldOfView/geometry")
          call extractDataContent(thisItem,geometryLabel)
          select case (trim(geometryLabel))
          case ("square")
             lightconeGeometry=lightconeGeometrySquare
             thisItem => XML_Get_First_Element_By_Tag_Name(doc     ,"fieldOfView/length")
             call extractDataContent(thisItem,lightconeFieldOfViewLength)
             lightconeFieldOfViewTanHalfAngle=tan(0.5d0*lightconeFieldOfViewLength)
             angularWeight=1.0d0/(lightconeFieldOfViewLength*(180.0d0/Pi))**2
          case default
             call Galacticus_Error_Report('Galacticus_Merger_Tree_Output_Filter_Lightcone','unknown field of view geometry')
          end select
          ! Get lightcone radial intervals.
          thisItem => XML_Get_First_Element_By_Tag_Name(doc,"outputs")
          call XML_Array_Read(thisItem,"redshift"       ,lightconeTime           )
          call XML_Array_Read(thisItem,"minimumDistance",lightconeMinimumDistance)
          call XML_Array_Read(thisItem,"maximumDistance",lightconeMaximumDistance)
          call destroy(doc)
          if     (                                                       &
               &   size(lightconeMinimumDistance) /= size(lightconeTime) &
               &  .or.                                                   &
               &   size(lightconeMaximumDistance) /= size(lightconeTime) &
               & ) call Galacticus_Error_Report('Galacticus_Merger_Tree_Output_Filter_Lightcone','size mismatch in outputs arrays')
          do iOutput=1,size(lightconeTime)
             lightconeTime(iOutput)=Cosmology_Age(Expansion_Factor_From_Redshift(lightconeTime(iOutput)))
          end do
          lightconeMinimumDistance=lightconeMinimumDistance*unitConversionLength
          lightconeMaximumDistance=lightconeMaximumDistance*unitConversionLength
          ! If requested, reduce the list of lightcone output times to the final time, and have
          ! it incorporate the full radial length of the lightcone. This allows the construction
          ! of lightcones with no evolution along the cone.
          if (filterLightconeFixedTime) then
             lightconeTimeTemp           =lightConeTime           (getLength(itemList))
             lightconeMinimumDistanceTemp=lightconeMinimumDistance(getLength(itemList))
             lightconeMaximumDistanceTemp=lightconeMaximumDistance(                  1)
             call Dealloc_Array(lightconeTime           )
             call Dealloc_Array(lightconeMinimumDistance)
             call Dealloc_Array(lightconeMaximumDistance)
             call Alloc_Array(lightconeTime           ,[1])
             call Alloc_Array(lightconeMinimumDistance,[1])
             call Alloc_Array(lightconeMaximumDistance,[1])
             lightConeTime           (1)=lightconeTimeTemp
             lightConeMinimumDistance(1)=lightconeMinimumDistanceTemp
             lightConeMaximumDistance(1)=lightconeMaximumDistanceTemp
          end if
          !$omp end critical (FoX_DOM_Access)
       end if
       ! Flag that this filter is now initialized.
       lightconeFilterInitialized=.true.
    end if
    return
  end subroutine Galacticus_Merger_Tree_Output_Filter_Lightcone_Initialize

  !# <mergerTreeOutputFilter>
  !#   <unitName>Galacticus_Merger_Tree_Output_Filter_Lightcone</unitName>
  !# </mergerTreeOutputFilter>
  subroutine Galacticus_Merger_Tree_Output_Filter_Lightcone(thisNode,doOutput)
    !% Determines whether {\tt thisNode} lies within a lightcone and, therefore, should be output.
    use Galacticus_Nodes
    use Arrays_Search
    use Cosmology_Functions
    implicit none
    type            (treeNode             )                , intent(inout), pointer :: thisNode
    logical                                                , intent(inout)          :: doOutput
    class           (nodeComponentBasic   )                               , pointer :: thisBasicComponent
    class           (nodeComponentPosition)                               , pointer :: thisPositionComponent
    double precision                       , parameter                              :: timeTolerance        =1.0d-3
    integer                                , dimension(3,2)                         :: periodicRange
    double precision                       , dimension(3  )                         :: galaxyPosition              , galaxyVelocity
    logical                                                                         :: galaxyIsInFieldOfView       , galaxyIsInLightcone
    integer                                                                         :: i                           , iAxis              , iOutput, &
         &                                                                             j                           , k

    ! Return immediately if this filter is not active.
    if (.not.lightconeFilterActive) return

    ! Get components.
    thisBasicComponent => thisNode%basic()

    ! Determine to which output this galaxy corresponds.
    iOutput=Search_Array_For_Closest(lightconeTime,thisBasicCOmponent%time(),timeTolerance)

    ! Determine range of possible replicants of this galaxy which could be in the lightcone.
    periodicRange(:,1)=floor  ((lightconeOrigin+lightconeMinimumDistance(iOutput)*lightconeUnitVector(:,1))/lightconeReplicationPeriod)-1
    periodicRange(:,2)=ceiling((lightconeOrigin+lightconeMaximumDistance(iOutput)*lightconeUnitVector(:,1))/lightconeReplicationPeriod)+0

    ! Get position of galaxy in original coordinates.
    thisPositionComponent => thisNode             %position()
    galaxyPosition        =  thisPositionComponent%position()/Expansion_Factor(lightconeTime(iOutput))

    ! Loop over all replicants.
    galaxyIsInLightcone=.false.
    do i=periodicRange(1,1),periodicRange(1,2)
       do j=periodicRange(2,1),periodicRange(2,2)
          do k=periodicRange(3,1),periodicRange(3,2)
             if (.not.galaxyIsInLightcone) then

                ! Compute position of galaxy in lightcone coordinate system.
                forall(iAxis=1:3)
                   lightconePosition(iAxis)=Dot_Product(galaxyPosition-lightconeOrigin+lightconeReplicationPeriod*dble([i,j,k]) &
                        & ,lightconeUnitVector(:,iAxis))
                end forall

                ! Test if galaxy lies within the correct angular window.
                select case (lightconeGeometry)
                case (lightconeGeometrySquare)
                   galaxyIsInFieldOfView=                                                                    &
                        & abs(lightconePosition(2)/lightconePosition(1)) < lightconeFieldOfViewTanHalfAngle &
                        &  .and.                                                                             &
                        & abs(lightconePosition(3)/lightconePosition(1)) < lightconeFieldOfViewTanHalfAngle
                end select

                ! Test if galaxy lies within appropriate radial range.
                galaxyIsInLightcone=      galaxyIsInFieldOfView                                     &
                     &              .and. lightconePosition(1) >  lightconeMinimumDistance(iOutput) &
                     &              .and. lightconePosition(1) <= lightconeMaximumDistance(iOutput)

                ! If the galaxy is in the lightcone compute also its velocity in the lightcone coordinate system and the redshift.
                if (galaxyIsInLightcone) then
                   ! Compute velocity of galaxy in lightcone coordinate system.
                   galaxyVelocity=thisPositionComponent%velocity()
                   forall(iAxis=1:3)
                      lightconeVelocity(iAxis)=Dot_Product(galaxyVelocity,lightconeUnitVector(:,iAxis))
                   end forall
                   ! Get redshift of the galaxy.
                   lightconeRedshift=Redshift_from_Expansion_Factor(   &
                        &             Expansion_Factor              (  &
                        &              Time_From_Comoving_Distance   ( &
                        &               lightconePosition(1)           &
                        &                                            ) &
                        &                                           )  &
                        &                                          )
                end if

             end if
          end do
       end do
    end do

    ! If the galaxy is not in the lightcone, filter it out.
    if (.not.galaxyIsInLightcone) doOutput=.false.
    return
  end subroutine Galacticus_Merger_Tree_Output_Filter_Lightcone

  function Filter_Lightcone_Get_Coordinates(enclosingItem)
    !% Extract a vector of coordinates from an XML DOM element.
    use Galacticus_Error
    use FoX_DOM
    use IO_XML
    implicit none
    double precision          , dimension(3)           :: Filter_Lightcone_Get_Coordinates
    type            (Node    ), intent(in   ), pointer :: enclosingItem

    if (XML_Array_Length(enclosingItem,"coordinate") /= 3) call Galacticus_Error_Report('Filter_Lightcone_Get_Coordinates'&
         &,'three axes must be specified')
    call XML_Array_Read_Static(enclosingItem,"coordinate",Filter_Lightcone_Get_Coordinates)
    return
  end function Filter_Lightcone_Get_Coordinates

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Lightcone_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Lightcone</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Lightcone_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of link properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode)              , intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI

    if (lightconeFilterActive) then
       !@ <outputPropertyGroup>
       !@   <name>lightconePosition</name>
       !@   <description>Lightcone X-Y-Z coordinates</description>
       !@   <outputType>nodeData</outputType>
       !@ </outputPropertyGroup>
       !@ <outputPropertyGroup>
       !@   <name>lightconeVelocity</name>
       !@   <description>Lightcone X-Y-Z velocity components</description>
       !@   <outputType>nodeData</outputType>
       !@ </outputPropertyGroup>
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>lightconePositionX</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Position of galaxy in lightcone (radial axis) [Mpc].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>lightconePosition</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='lightconePositionX'
       doublePropertyComments(doubleProperty)='Position of galaxy in lightcone (radial axis) [Mpc].'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>lightconePositionY</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Position of galaxy in lightcone (1st angular axis) [Mpc].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>lightconePosition</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='lightconePositionY'
       doublePropertyComments(doubleProperty)='Position of galaxy in lightcone (1st angular axis) [Mpc].'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>lightconePositionZ</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Position of galaxy in lightcone (2nd angular axis) [Mpc].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>lightconePosition</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='lightconePositionZ'
       doublePropertyComments(doubleProperty)='Position of galaxy in lightcone (2nd angular axis) [Mpc].'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>lightconeVelocityX</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Velocity of galaxy in lightcone (radial axis) [km/s].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>lightconeVelocity</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='lightconeVelocityX'
       doublePropertyComments(doubleProperty)='Velocity of galaxy in lightcone (radial axis) [km/s].'
       doublePropertyUnitsSI (doubleProperty)=kilo
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>lightconeVelocityY</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Velocity of galaxy in lightcone (1st angular axis) [km/s].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>lightconeVelocity</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='lightconeVelocityY'
       doublePropertyComments(doubleProperty)='Velocity of galaxy in lightcone (1st angular axis) [km/s].'
       doublePropertyUnitsSI (doubleProperty)=kilo
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>lightconeVelocityZ</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Velocity of galaxy in lightcone (2nd angular axis) [km/s].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>lightconeVelocity</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='lightconeVelocityZ'
       doublePropertyComments(doubleProperty)='Velocity of galaxy in lightcone (2nd angular axis) [km/s].'
       doublePropertyUnitsSI (doubleProperty)=kilo
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>lightconeRedshift</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Reshift of galaxy in lightcone.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='lightconeRedshift'
       doublePropertyComments(doubleProperty)='Reshift of galaxy in lightcone.'
       doublePropertyUnitsSI (doubleProperty)=0.0d0
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>angularWeight</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Number of such galaxies per unit area [degrees^-2].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='angularWeight'
       doublePropertyComments(doubleProperty)='Number of such galaxies per unit area [degrees⁻²].'
       doublePropertyUnitsSI (doubleProperty)=(180.0d0/Pi)**2
    end if
    return
  end subroutine Galacticus_Output_Tree_Lightcone_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Lightcone_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Lightcone</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Lightcone_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of link properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

    if (lightconeFilterActive) doublePropertyCount=doublePropertyCount+lightconePropertyCount
    return
  end subroutine Galacticus_Output_Tree_Lightcone_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Lightcone</unitName>
  !#  <sortName>Galacticus_Output_Tree_Lightcone</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Lightcone(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store link properties in the \glc\ output file buffers.
    use Galacticus_Nodes
    use Kind_Numbers
    implicit none
    double precision                , intent(in   )          :: time
    type            (treeNode      ), intent(inout), pointer :: thisNode
    integer                         , intent(inout)          :: doubleBufferCount     , doubleProperty, integerBufferCount, &
         &                                                      integerProperty
    integer         (kind=kind_int8), intent(inout)          :: integerBuffer    (:,:)
    double precision                , intent(inout)          :: doubleBuffer     (:,:)

   if (lightconeFilterActive) then
       doubleBuffer(doubleBufferCount,doubleProperty+1:doubleProperty+3)=lightconePosition
       doubleProperty=doubleProperty+3
       doubleBuffer(doubleBufferCount,doubleProperty+1:doubleProperty+3)=lightconeVelocity
       doubleProperty=doubleProperty+3
       doubleBuffer(doubleBufferCount,doubleProperty+1                 )=lightconeRedshift
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty+1                 )=angularWeight
       doubleProperty=doubleProperty+1
    end if
    return
  end subroutine Galacticus_Output_Tree_Lightcone

end module Galacticus_Merger_Tree_Output_Filter_Lightcones
