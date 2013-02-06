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

!% Contains a module which handles outputting of node virial data to the \glc\ output file.

module Galacticus_Output_Trees_Rotation_Curve
  !% Handles outputting of node virial data to the \glc\ output file.
  use ISO_Varying_String
  implicit none
  private
  public :: Galacticus_Output_Tree_Rotation_Curve, Galacticus_Output_Tree_Rotation_Curve_Property_Count, Galacticus_Output_Tree_Rotation_Curve_Names

  ! Flag indicating whether or not virial information is to be output.
  logical                                            :: outputRotationCurveData

  ! Flag indicating whether or not this module has been initialized.
  logical                                            :: outputRotationCurveDataInitialized=.false.

  ! Radii at which rotation curve is to be output.
  integer                                            :: radiiCount
  type radiusSpecifier
     type            (varying_string) :: name
     integer                          :: type,mass,component
     logical                          :: loaded
     double precision                 :: value
  end type radiusSpecifier
  type   (radiusSpecifier), dimension(:), allocatable :: radii
  type   (varying_string ), dimension(:), allocatable :: outputRotationCurveRadii
  logical                                             :: diskIsNeeded,spheroidIsNeeded,virialRadiusIsNeeded&
       &,darkMatterScaleRadiusIsNeeded
  integer                 , parameter                 :: radiusTypeRadius                =0
  integer                 , parameter                 :: radiusTypeVirialRadius          =1
  integer                 , parameter                 :: radiusTypeDarkMatterScaleRadius =2
  integer                 , parameter                 :: radiusTypeDiskRadius            =3
  integer                 , parameter                 :: radiusTypeSpheroidRadius        =4
  integer                 , parameter                 :: radiusTypeDiskHalfMassRadius    =5
  integer                 , parameter                 :: radiusTypeSpheroidHalfMassRadius=6

contains

  subroutine Galacticus_Output_Tree_Rotation_Curve_Initialize
    !% Initializes the module by determining whether or not virial data should be output.
    use Input_Parameters
    use Galacticus_Error
    use String_Handling
    use Galacticus_Nodes
    use Memory_Management
    use Galactic_Structure_Options
    implicit none
    type     (varying_string), dimension(5) :: radiusDefinition
    character(len=20        )               :: radiusLabel
    integer                                 :: i

    if (.not.outputRotationCurveDataInitialized) then
       !$omp critical(Galacticus_Output_Tree_Rotation_Curve_Initialize)
       if (.not.outputRotationCurveDataInitialized) then
          !@ <inputParameter>
          !@   <name>outputRotationCurveData</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not rotation curve data should be incldued in the output file.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputRotationCurveData',outputRotationCurveData,defaultValue=.false.)
          

          ! Read and parse parameter controlling radii of output.
          if (outputRotationCurveData) then
             if (.not.Input_Parameter_Is_Present('outputRotationCurveRadii'))                                            &
                  & call Galacticus_Error_Report(                                                                        &
                  &                              'Galacticus_Output_Tree_Rotation_Curve_Initialize'                    , &
                  &                              'outputRotationCurveRadii must be specified for rotation curve output'  &
                  &                             )
             radiiCount=Get_Input_Parameter_Array_Size('outputRotationCurveRadii')
             allocate(outputRotationCurveRadii(radiiCount))
             allocate(radii                   (radiiCount))
             !@ <inputParameter>
             !@   <name>outputRotationCurveRadii</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Specifies the radii at which the rotation curve will be output.
             !@   </description>
             !@   <type>string</type>
             !@   <cardinality>1..*</cardinality>
             !@   <group>output</group>
             !@ </inputParameter>
             call Get_Input_Parameter('outputRotationCurveRadii',outputRotationCurveRadii)
             ! Parse the radii definitions.
             diskIsNeeded                  =.false.
             spheroidIsNeeded              =.false.
             virialRadiusIsNeeded          =.false.
             darkMatterScaleRadiusIsNeeded =.false.
             do i=1,radiiCount
                radii(i)%name=outputRotationCurveRadii(i)
                call String_Split_Words(radiusDefinition,char(outputRotationCurveRadii(i)),':')
                select case (char(radiusDefinition(1)))
                case ('radius'                )
                   radii(i)%type=radiusTypeRadius
                case ('virialRadius'          )
                   radii(i)%type=radiusTypeVirialRadius
                   virialRadiusIsNeeded         =.true.
                case ('darkMatterScaleRadius' )
                   radii(i)%type=radiusTypeDarkMatterScaleRadius
                   darkMatterScaleRadiusIsNeeded=.true.
                   if (.not.defaultDarkMatterProfileComponent%scaleIsGettable         ())                  &
                        & call Galacticus_Error_Report(                                                    &
                        &                              'Galacticus_Output_Tree_Rotation_Curve_Initialize', &
                        &                              'disk half-mass radius is not gettable'             &
                        &                             )
                case ('diskRadius'            )
                   radii(i)%type=radiusTypeDiskRadius
                   diskIsNeeded                 =.true.
                   if (.not.defaultDiskComponent             %radiusIsGettable        ())                  &
                        & call Galacticus_Error_Report(                                                    &
                        &                              'Galacticus_Output_Tree_Rotation_Curve_Initialize', &
                        &                              'disk radius is not gettable'                       &
                        &                             )
                case ('spheroidRadius'        )
                   radii(i)%type=radiusTypeSpheroidRadius
                   spheroidIsNeeded             =.true.
                   if (.not.defaultSpheroidComponent         %radiusIsGettable        ())                  &
                        & call Galacticus_Error_Report(                                                    &
                        &                              'Galacticus_Output_Tree_Rotation_Curve_Initialize', &
                        &                              'spheroid radius is not gettable'                   &
                        &                             )
                case ('diskHalfMassRadius'    )
                   radii(i)%type=radiusTypeDiskHalfMassRadius
                   diskIsNeeded                 =.true.
                   if (.not.defaultDiskComponent             %halfMassRadiusIsGettable())                  &
                        & call Galacticus_Error_Report(                                                    &
                        &                              'Galacticus_Output_Tree_Rotation_Curve_Initialize', &
                        &                              'disk half-mass radius is not gettable'             &
                        &                             )
                case ('spheroidHalfMassRadius')
                   radii(i)%type=radiusTypeSpheroidHalfMassRadius
                   spheroidIsNeeded   =.true.
                   if (.not.defaultSpheroidComponent         %halfMassRadiusIsGettable())                  &
                        & call Galacticus_Error_Report(                                                    &
                        &                              'Galacticus_Output_Tree_Rotation_Curve_Initialize', &
                        &                              'spheroid half-mass radius is not gettable'         &
                        &                             )
                case default
                   call Galacticus_Error_Report('Galacticus_Output_Tree_Rotation_Curve_Initialize','unrecognized radius specifier')
                end select
                radii(i)%component=Galactic_Structure_Component_Type_Decode(char(radiusDefinition(2)))
                radii(i)%mass     =Galactic_Structure_Mass_Type_Decode     (char(radiusDefinition(3)))
                select case (char(radiusDefinition(4)))
                case ('loaded'  )
                   radii(i)%loaded=.true.
                case ('unloaded')
                   radii(i)%loaded=.false.
                case default
                   call Galacticus_Error_Report('Galacticus_Output_Tree_Rotation_Curve_Initialize','unrecognized loading specifier')
                end select
                radiusLabel=radiusDefinition(5)
                read (radiusLabel,*) radii(i)%value
             end do
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
    use Galacticus_Nodes
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    implicit none
    type(treeNode),   intent(inout), pointer      :: thisNode
    double precision, intent(in   )               :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI
    integer                                       :: i

    ! Initialize the module.
    call Galacticus_Output_Tree_Rotation_Curve_Initialize

    ! Return property names if we are outputting rotation curve data.
    if (outputRotationCurveData) then
       !@ <outputProperty>
       !@   <name>rotationCurve</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Rotation curve at a given radius.</description>
       !@   <label>???</label>
       !@ </outputProperty>
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
    use Galacticus_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in   )          :: time
    integer,          intent(inout)          :: integerPropertyCount,doublePropertyCount
    
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
       &,doubleBufferCount,doubleBuffer,time)
    !% Store rotation curve properties in the \glc\ output file buffers.
    use Galacticus_Nodes
    use Kind_Numbers
    use Galactic_Structure_Rotation_Curves
    use Dark_Matter_Halo_Scales
    implicit none
    double precision                                , intent(in   )          :: time
    type            (treeNode                      ), intent(inout), pointer :: thisNode
    integer                                         , intent(inout)          :: integerProperty,integerBufferCount,doubleProperty&
         & ,doubleBufferCount
    integer         (kind=kind_int8                ), intent(inout)          :: integerBuffer(:,:)
    double precision                                , intent(inout)          :: doubleBuffer (:,:)
    class           (nodeComponentDisk             ),                pointer :: thisDisk
    class           (nodeComponentSpheroid         ),                pointer :: thisSpheroid
    class           (nodeComponentDarkMatterProfile),                pointer :: thisDarkMatterProfile
    integer                                                                  :: i
    double precision                                                         :: radius,radiusVirial
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Rotation_Curve_Initialize
    ! Store property data if we are outputting rotation curve data.
    if (outputRotationCurveData) then
       ! Compute required quantities.
       if (         virialRadiusIsNeeded) radiusVirial          =  Dark_Matter_Halo_Virial_Radius(thisNode                    )
       if (                 diskIsNeeded) thisDisk              =>                                thisNode%disk             ()
       if (             spheroidIsNeeded) thisSpheroid          =>                                thisNode%spheroid         ()
       if (darkMatterScaleRadiusIsNeeded) thisDarkMatterProfile =>                                thisNode%darkMatterProfile()
       do i=1,radiiCount
          ! Find the radius.
          radius=radii(i)%value
          select case (radii(i)%type)
          case (radiusTypeRadius                )
             ! Nothing to do.
          case (radiusTypeVirialRadius          )
             radius=radius*radiusVirial
          case (radiusTypeDarkMatterScaleRadius )
             radius=radius*thisDarkMatterProfile%         scale()
          case (radiusTypeDiskRadius            )
             radius=radius*thisDisk             %        radius()
          case (radiusTypeSpheroidRadius        )
             radius=radius*thisSpheroid         %        radius()
          case (radiusTypeDiskHalfMassRadius    )
             radius=radius*thisDisk             %halfMassRadius()
          case (radiusTypeSpheroidHalfMassRadius)
             radius=radius*thisSpheroid         %halfMassRadius()
          end select
          ! Store the rotation curve.
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=                                                                    &
               &                                         Galactic_Structure_Rotation_Curve(                                  &
               &                                                                           thisNode                        , &
               &                                                                           radius                          , &
               &                                                                           massType=radii(i)%mass          , &
               &                                                                           componentType=radii(i)%component, &
               &                                                                           haloLoaded=radii(i)%loaded        &
               &                                                                          )
       end do
    end if
    return
  end subroutine Galacticus_Output_Tree_Rotation_Curve

end module Galacticus_Output_Trees_Rotation_Curve
