!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which filters output for lightcone geometry.

module Galacticus_Merger_Tree_Output_Filter_Lightcones
  !% Filters output for lightcone geometry.
  implicit none
  private
  public :: Galacticus_Merger_Tree_Output_Filter_Lightcone,Galacticus_Output_Tree_Lightcone,&
       & Galacticus_Output_Tree_Lightcone_Property_Count, Galacticus_Output_Tree_Lightcone_Names,&
       & Galacticus_Merger_Tree_Output_Filter_Lightcone_Initialize

  ! Number of lightcone properties.
  integer, parameter   :: lightconePropertyCount=8

  ! Flags indicating if the module has been initialized and if this filter is active.
  logical                                      :: lightconeFilterInitialized=.false.
  logical                                      :: lightconeFilterActive

  ! Variables used to describe the lightcone geometry.
  double precision, dimension(3  )             :: lightconeOrigin
  double precision, dimension(3,3)             :: lightconeUnitVector
  double precision, dimension(:),  allocatable :: lightconeTime,lightconeMinimumDistance,lightconeMaximumDistance
  double precision                             :: lightconeReplicationPeriod,lightconeFieldOfViewLength,lightconeFieldOfViewTanHalfAngle
  double precision                             :: angularWeight
  integer                                      :: lightconeGeometry
  integer,          parameter                  :: lightconeGeometrySquare=1

  ! Variables that store the position, velocity and redshift of a galaxy that is found to be in the lightcone. These are then used
  ! to output that position.
  double precision, dimension(3  )             :: lightconePosition,lightconeVelocity
  double precision                             :: lightconeRedshift
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
    use Cosmology_Functions
    use Numerical_Constants_Astronomical
    use Cosmological_Parameters
    use Galacticus_Error
    use Memory_Management
    implicit none
    type(varying_string), intent(in),    dimension(:) :: filterNames
    type(Node),           pointer                     :: doc,thisItem
    type(NodeList),       pointer                     :: itemList
    integer                                           :: iAxis,iOutput,ioErr,lengthUnitsHubbleExponent
    double precision                                  :: lengthUnitsInSI,unitConversionLength
    type(varying_string)                              :: filterLightconeGeometryFileName
    character(len=11)                                 :: tagName
    character(len=10)                                 :: geometryLabel

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
          !@ </inputParameter>
          call Get_Input_Parameter('filterLightconeGeometryFileName',filterLightconeGeometryFileName)

          ! Extract data from geometry specification file.
          !$omp critical (FoX_DOM_Access)
          doc => parseFile(char(filterLightconeGeometryFileName),iostat=ioErr)
          if (ioErr /= 0) call Galacticus_Error_Report('Galacticus_Merger_Tree_Output_Filter_Lightcone','Unable to parse geometry file')
          ! Get the unit vectors of the lightcone.
          do iAxis=1,3
             write (tagName,'(a,i1)') "unitVector",iAxis
             itemList => getElementsByTagname(doc,tagName)
             thisItem => item(itemList,0)
             lightconeUnitVector(:,iAxis)=Filter_Lightcone_Get_Coordinates(thisItem)
          end do
          ! Get units information.
          itemList => getElementsByTagname(doc     ,"units"         )
          thisItem => item(itemList,0)
          itemList => getElementsByTagname(thisItem,"length"        )
          thisItem => item(itemList,0)
          itemList => getElementsByTagname(thisItem,"unitsInSI"     )
          thisItem => item(itemList,0)
          call extractDataContent(thisItem,lengthUnitsInSI          )
          itemList => getElementsByTagname(doc     ,"units"         )
          thisItem => item(itemList,0)
          itemList => getElementsByTagname(thisItem,"length"        )
          thisItem => item(itemList,0)
          itemList => getElementsByTagname(thisItem,"hubbleExponent")
          thisItem => item(itemList,0)
          call extractDataContent(thisItem,lengthUnitsHubbleExponent)
          unitConversionLength=lengthUnitsInSI*(Little_H_0()**lengthUnitsHubbleExponent)/megaParsec
          ! Get the origin of the lightcone.
          itemList => getElementsByTagname(doc,"origin")
          thisItem => item(itemList,0)
          lightconeOrigin=Filter_Lightcone_Get_Coordinates(thisItem)*unitConversionLength
          ! Get the replication period.
          itemList => getElementsByTagname(doc,"boxLength")
          thisItem => item(itemList,0)
          call extractDataContent(thisItem,lightconeReplicationPeriod)
          lightconeReplicationPeriod=lightconeReplicationPeriod*unitConversionLength
          ! Get the geometry.
          itemList => getElementsByTagname(doc,"fieldOfView")
          thisItem => item(itemList,0)
          itemList => getElementsByTagname(thisItem,"geometry")
          thisItem => item(itemList,0)
          call extractDataContent(thisItem,geometryLabel)
          select case (trim(geometryLabel))
          case ("square")
             lightconeGeometry=lightconeGeometrySquare
             itemList => getElementsByTagname(doc,"fieldOfView")
             thisItem => item(itemList,0)
             itemList => getElementsByTagname(thisItem,"length")
             thisItem => item(itemList,0)
             call extractDataContent(thisItem,lightconeFieldOfViewLength)
             lightconeFieldOfViewTanHalfAngle=dtan(0.5d0*lightconeFieldOfViewLength)
             angularWeight=1.0d0/(lightconeFieldOfViewLength*(180.0d0/Pi))**2
          case default
             call Galacticus_Error_Report('Galacticus_Merger_Tree_Output_Filter_Lightcone','unknown field of view geometry')
          end select
          ! Get lightcone radial intervals.
          itemList => getElementsByTagname(doc,"outputs")
          thisItem => item(itemList,0)
          itemList => getElementsByTagname(thisItem,"redshift")
          call Alloc_Array(lightconeTime,[getLength(itemList)])
          do iOutput=1,getLength(itemList)
             thisItem => item(itemList,iOutput-1)
             call extractDataContent(thisItem,lightconeTime(iOutput))
             lightconeTime(iOutput)=Cosmology_Age(Expansion_Factor_From_Redshift(lightconeTime(iOutput)))
          end do
          itemList => getElementsByTagname(doc,"outputs")
          thisItem => item(itemList,0)
          itemList => getElementsByTagname(thisItem,"minimumDistance")
          if (getLength(itemList) /= size(lightconeTime)) call Galacticus_Error_Report('Galacticus_Merger_Tree_Output_Filter_Lightcone' &
               & ,'size mismatch in outputs arrays')
          call Alloc_Array(lightconeMinimumDistance,[getLength(itemList)])
          do iOutput=1,getLength(itemList)
             thisItem => item(itemList,iOutput-1)
             call extractDataContent(thisItem,lightconeMinimumDistance(iOutput))
          end do
          lightconeMinimumDistance=lightconeMinimumDistance*unitConversionLength
          itemList => getElementsByTagname(doc,"outputs")
          thisItem => item(itemList,0)
          itemList => getElementsByTagname(thisItem,"maximumDistance")
          if (getLength(itemList) /= size(lightconeTime)) call Galacticus_Error_Report('Galacticus_Merger_Tree_Output_Filter_Lightcone' &
               & ,'size mismatch in outputs arrays')
          call Alloc_Array(lightconeMaximumDistance,[getLength(itemList)])
          do iOutput=1,getLength(itemList)
             thisItem => item(itemList,iOutput-1)
             call extractDataContent(thisItem,lightconeMaximumDistance(iOutput))
          end do
          lightconeMaximumDistance=lightconeMaximumDistance*unitConversionLength
          call destroy(doc)
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
    use Tree_Nodes
    use Arrays_Search
    use Numerical_Constants_Astronomical
    use Cosmology_Functions
    implicit none
    type(treeNode),       intent(inout), pointer      :: thisNode
    logical,              intent(inout)               :: doOutput
    double precision,     parameter                   :: timeTolerance=1.0d-3
    integer,              dimension(3,3)              :: periodicRange
    double precision,     dimension(3  )              :: galaxyPosition,galaxyVelocity
    logical                                           :: galaxyIsInLightcone,galaxyIsInFieldOfView
    integer                                           :: iAxis,iOutput,i,j,k
    
    ! Return immediately if this filter is not active.
    if (.not.lightconeFilterActive) return

    ! Determine to which output this galaxy corresponds.
    iOutput=Search_Array_For_Closest(lightconeTime,Tree_Node_Time(thisNode),timeTolerance)

    ! Determine range of possible replicants of this galaxy which could be in the lightcone.
    periodicRange(:,1)=floor  ((lightconeOrigin+lightconeMinimumDistance(iOutput)*lightconeUnitVector(:,1))/lightconeReplicationPeriod)-1
    periodicRange(:,2)=ceiling((lightconeOrigin+lightconeMaximumDistance(iOutput)*lightconeUnitVector(:,1))/lightconeReplicationPeriod)+0

    ! Get position of galaxy in original coordinates.
    call Tree_Node_Position(thisNode,galaxyPosition)
    galaxyPosition=galaxyPosition/Expansion_Factor(lightconeTime(iOutput))

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
                        & dabs(lightconePosition(2)/lightconePosition(1)) < lightconeFieldOfViewTanHalfAngle &
                        &  .and.                                                                             &
                        & dabs(lightconePosition(3)/lightconePosition(1)) < lightconeFieldOfViewTanHalfAngle
                end select
                
                ! Test if galaxy lies within appropriate radial range.
                galaxyIsInLightcone=      galaxyIsInFieldOfView                                     &
                     &              .and. lightconePosition(1) >  lightconeMinimumDistance(iOutput) &
                     &              .and. lightconePosition(1) <= lightconeMaximumDistance(iOutput)

                ! If the galaxy is in the lightcone compute also its velocity in the lightcone coordinate system and the redshift.
                if (galaxyIsInLightcone) then
                   ! Compute velocity of galaxy in lightcone coordinate system.
                   call Tree_Node_Velocity(thisNode,galaxyVelocity)
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
    implicit none
    double precision, dimension(3)            :: Filter_Lightcone_Get_Coordinates
    type(Node),       pointer,     intent(in) :: enclosingItem
    type(Node),       pointer                 :: thisItem
    type(NodeList),   pointer                 :: coordinates
    integer                                   :: iAxis

    coordinates => getElementsByTagname(enclosingItem,"coordinate")
    if (getLength(coordinates) /= 3) call Galacticus_Error_Report('Filter_Lightcone_Get_Coordinates','three axes must be specified')
    do iAxis=1,3
       thisItem => item(coordinates,iAxis-1)    
       call extractDataContent(thisItem,Filter_Lightcone_Get_Coordinates(iAxis))
    end do
    return
  end function Filter_Lightcone_Get_Coordinates

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Lightcone_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Lightcone</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Lightcone_Names(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of link properties to be written to the \glc\ output file.
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI

    if (lightconeFilterActive) then
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='lightconePositionX'
       doublePropertyComments(doubleProperty)='Position of galaxy in lightcone (radial axis) [Mpc].'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='lightconePositionY'
       doublePropertyComments(doubleProperty)='Position of galaxy in lightcone (1st angular axis) [Mpc].'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='lightconePositionZ'
       doublePropertyComments(doubleProperty)='Position of galaxy in lightcone (2nd angular axis) [Mpc].'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='lightconeVelocityX'
       doublePropertyComments(doubleProperty)='Velocity of galaxy in lightcone (radial axis) [km/s].'
       doublePropertyUnitsSI (doubleProperty)=kilo
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='lightconeVelocityY'
       doublePropertyComments(doubleProperty)='Velocity of galaxy in lightcone (1st angular axis) [km/s].'
       doublePropertyUnitsSI (doubleProperty)=kilo
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='lightconeVelocityZ'
       doublePropertyComments(doubleProperty)='Velocity of galaxy in lightcone (2nd angular axis) [km/s].'
       doublePropertyUnitsSI (doubleProperty)=kilo
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='lightconeRedshift'
       doublePropertyComments(doubleProperty)='Reshift of galaxy in lightcone.'
       doublePropertyUnitsSI (doubleProperty)=0.0d0
       doubleProperty=doubleProperty+1
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
  subroutine Galacticus_Output_Tree_Lightcone_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of link properties to be written to the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

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
    use Tree_Nodes
    use Kind_Numbers
    implicit none
    double precision,        intent(in)             :: time
    type(treeNode),          intent(inout), pointer :: thisNode
    integer,                 intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)          :: integerBuffer(:,:)
    double precision,        intent(inout)          :: doubleBuffer(:,:)

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
