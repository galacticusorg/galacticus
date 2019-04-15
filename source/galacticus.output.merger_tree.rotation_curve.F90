!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which handles outputting of rotation curve data to the \glc\ output file.

module Galacticus_Output_Trees_Rotation_Curve
  !% Handles outputting of rotation curve data to the \glc\ output file.
  use ISO_Varying_String
  use Galactic_Structure_Radii_Definitions
  implicit none
  private
  public :: Galacticus_Output_Tree_Rotation_Curve, Galacticus_Output_Tree_Rotation_Curve_Property_Count, Galacticus_Output_Tree_Rotation_Curve_Names

  ! Flag indicating whether or not rotation curve information is to be output.
  logical :: outputRotationCurveData

  ! Flag indicating whether or not this module has been initialized.
  logical :: outputRotationCurveDataInitialized=.false.

  ! Radii at which rotation curve is to be output.
  integer :: radiiCount
  type   (radiusSpecifier), allocatable, dimension(:) :: radii
  type   (varying_string ), allocatable, dimension(:) :: outputRotationCurveRadii
  logical                                             :: darkMatterScaleRadiusIsNeeded, diskIsNeeded        , &
       &                                                 spheroidIsNeeded             , virialRadiusIsNeeded

contains

  subroutine Galacticus_Output_Tree_Rotation_Curve_Initialize()
    !% Initializes the module by determining whether or not rotation curve data should be output.
    use Input_Parameters
    use Galacticus_Error
    implicit none

    if (.not.outputRotationCurveDataInitialized) then
       !$omp critical(Galacticus_Output_Tree_Rotation_Curve_Initialize)
       if (.not.outputRotationCurveDataInitialized) then
          !# <inputParameter>
          !#   <name>outputRotationCurveData</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.false.</defaultValue>
          !#   <description>Specifies whether or not rotation curve data should be incldued in the output file.</description>
          !#   <group>output</group>
          !#   <source>globalParameters</source>
          !#   <type>boolean</type>
          !# </inputParameter>
          ! Read and parse parameter controlling radii of output.
          if (outputRotationCurveData) then
             if (.not.globalParameters%isPresent('outputRotationCurveRadii'))                                             &
                  & call Galacticus_Error_Report(                                                                         &
                  &                              'outputRotationCurveRadii must be specified for rotation curve output'// &
                  &                              {introspection:location}                                                 &
                  &                             )
             radiiCount=globalParameters%count('outputRotationCurveRadii')
             allocate(outputRotationCurveRadii(radiiCount))
             !# <inputParameter>
             !#   <name>outputRotationCurveRadii</name>
             !#   <cardinality>1..*</cardinality>
             !#   <description>Specifies the radii at which the rotation curve will be output.</description>
             !#   <group>output</group>
             !#   <source>globalParameters</source>
             !#   <type>string</type>
             !# </inputParameter>
             ! Parse the radii definitions.
             call Galactic_Structure_Radii_Definition_Decode(outputRotationCurveRadii,radii,diskIsNeeded,spheroidIsNeeded,virialRadiusIsNeeded,darkMatterScaleRadiusIsNeeded)
             deallocate(outputRotationCurveRadii)
          end if
          ! Flag that module is now initialized.
          outputRotationCurveDataInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Rotation_Curve_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Rotation_Curve_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Rotation_Curve_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Rotation_Curve</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Rotation_Curve_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of rotation curve properties to be written to the \glc\ output file.
    use Galacticus_Nodes                , only : treeNode
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode)              , intent(inout) :: thisNode
    double precision                        , intent(in   ) :: time
    integer                                 , intent(inout) :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout) :: doublePropertyComments , doublePropertyNames   , &
         &                                                     integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout) :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    integer                                                 :: i
    !GCC$ attributes unused :: thisNode, time, integerProperty, integerPropertyNames, integerPropertyComments, integerPropertyUnitsSI
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Rotation_Curve_Initialize

    ! Return property names if we are outputting rotation curve data.
    if (outputRotationCurveData) then
       do i=1,radiiCount
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='rotationCurve:'//char(radii(i)%name)
          doublePropertyComments(doubleProperty)='Rotation curve at a given radius'
          doublePropertyUnitsSI (doubleProperty)=kilo
       end do
    end if
    return
  end subroutine Galacticus_Output_Tree_Rotation_Curve_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Rotation_Curve_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Rotation_Curve</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Rotation_Curve_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of rotation curve properties to be written to the \glc\ output file.
    use Galacticus_Nodes, only : treeNode
    implicit none
    type            (treeNode), intent(inout) :: thisNode
    double precision          , intent(in   ) :: time
    integer                   , intent(inout) :: doublePropertyCount, integerPropertyCount
    !GCC$ attributes unused :: thisNode, time, integerPropertyCount
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Rotation_Curve_Initialize

    ! Increment property count if we are outputting rotation curve data.
    if (outputRotationCurveData) doublePropertyCount=doublePropertyCount+radiiCount
    return
  end subroutine Galacticus_Output_Tree_Rotation_Curve_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Rotation_Curve</unitName>
  !#  <sortName>Galacticus_Output_Tree_Rotation_Curve</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Rotation_Curve(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time,instance)
    !% Store rotation curve properties in the \glc\ output file buffers.
    use Galacticus_Nodes, only : treeNode, nodeComponentDisk, nodeComponentSpheroid, nodeComponentDarkMatterProfile
    use Kind_Numbers
    use Galactic_Structure_Rotation_Curves
    use Dark_Matter_Halo_Scales
    use Galactic_Structure_Options
    use Galactic_Structure_Enclosed_Masses
    use Multi_Counters
    implicit none
    double precision                                , intent(in   ) :: time
    type            (treeNode                      ), intent(inout) :: thisNode
    integer                                         , intent(inout) :: doubleBufferCount         , doubleProperty, integerBufferCount, &
         &                                                             integerProperty
    integer         (kind=kind_int8                ), intent(inout) :: integerBuffer        (:,:)
    double precision                                , intent(inout) :: doubleBuffer         (:,:)
    type            (multiCounter                  ), intent(inout) :: instance
    class           (nodeComponentDisk             ), pointer       :: thisDisk
    class           (nodeComponentSpheroid         ), pointer       :: thisSpheroid
    class           (nodeComponentDarkMatterProfile), pointer       :: thisDarkMatterProfile
    class           (darkMatterHaloScaleClass      ), pointer       :: darkMatterHaloScale_
    integer                                                         :: i
    double precision                                                :: radius                    , radiusVirial
    !GCC$ attributes unused :: time, integerProperty, integerBufferCount, integerBuffer, instance
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Rotation_Curve_Initialize
    ! Store property data if we are outputting rotation curve data.
    if (outputRotationCurveData) then
       ! Compute required quantities.
       darkMatterHaloScale_ => darkMatterHaloScale()
       radiusVirial         =  0.0d0
       if (         virialRadiusIsNeeded) radiusVirial          =  darkMatterHaloScale_%virialRadius(thisNode                    )
       if (                 diskIsNeeded) thisDisk              =>                                thisNode%disk             ()
       if (             spheroidIsNeeded) thisSpheroid          =>                                thisNode%spheroid         ()
       if (darkMatterScaleRadiusIsNeeded) thisDarkMatterProfile =>                                thisNode%darkMatterProfile()
       do i=1,radiiCount
          ! Find the radius.
          radius=radii(i)%value
          select case (radii(i)%type)
          case   (radiusTypeRadius                )
             ! Nothing to do.
          case   (radiusTypeVirialRadius          )
             radius=radius*radiusVirial
          case   (radiusTypeDarkMatterScaleRadius )
             radius=radius*thisDarkMatterProfile%         scale()
          case   (radiusTypeDiskRadius            )
             radius=radius*thisDisk             %        radius()
          case   (radiusTypeSpheroidRadius        )
             radius=radius*thisSpheroid         %        radius()
          case   (radiusTypeDiskHalfMassRadius    )
             radius=radius*thisDisk             %halfMassRadius()
          case   (radiusTypeSpheroidHalfMassRadius)
             radius=radius*thisSpheroid         %halfMassRadius()
          case   (radiusTypeGalacticMassFraction ,  &
               &  radiusTypeGalacticLightFraction )
             radius= radius                                   &
                  & *Galactic_Structure_Radius_Enclosing_Mass &
                  &  (                                        &
                  &   thisNode                             ,  &
                  &   fractionalMass=radii(i)%fraction     ,  &
                  &   massType      =massTypeGalactic      ,  &
                  &   componentType =componentTypeAll      ,  &
                  &   weightBy      =radii(i)%weightBy     ,  &
                  &   weightIndex   =radii(i)%weightByIndex   &
                  &  )
          end select
          ! Store the rotation curve.
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=Galactic_Structure_Rotation_Curve(                                  &
               &                                                                           thisNode                        , &
               &                                                                           radius                          , &
               &                                                                           componentType=radii(i)%component, &
               &                                                                           massType     =radii(i)%mass       &
               &                                                                          )
       end do
    end if
    return
  end subroutine Galacticus_Output_Tree_Rotation_Curve

end module Galacticus_Output_Trees_Rotation_Curve
