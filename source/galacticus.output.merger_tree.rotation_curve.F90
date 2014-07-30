!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  implicit none
  private
  public :: Galacticus_Output_Tree_Rotation_Curve, Galacticus_Output_Tree_Rotation_Curve_Property_Count, Galacticus_Output_Tree_Rotation_Curve_Names

  ! Flag indicating whether or not rotation curve information is to be output.
  logical :: outputRotationCurveData

  ! Flag indicating whether or not this module has been initialized.
  logical :: outputRotationCurveDataInitialized=.false.

  ! Radii at which rotation curve is to be output.
  integer :: radiiCount
  type radiusSpecifier
     type            (varying_string) :: name
     integer                          :: component    , mass , type, weightBy, &
          &                              weightByIndex
     logical                          :: loaded
     double precision                 :: fraction     , value
  end type radiusSpecifier
  type   (radiusSpecifier)           , allocatable, dimension(:) :: radii
  type   (varying_string )           , allocatable, dimension(:) :: outputRotationCurveRadii
  logical                                                        :: darkMatterScaleRadiusIsNeeded     , diskIsNeeded        , &
       &                                                            spheroidIsNeeded                  , virialRadiusIsNeeded
  integer                 , parameter                            :: radiusTypeRadius                =0
  integer                 , parameter                            :: radiusTypeVirialRadius          =1
  integer                 , parameter                            :: radiusTypeDarkMatterScaleRadius =2
  integer                 , parameter                            :: radiusTypeDiskRadius            =3
  integer                 , parameter                            :: radiusTypeSpheroidRadius        =4
  integer                 , parameter                            :: radiusTypeDiskHalfMassRadius    =5
  integer                 , parameter                            :: radiusTypeSpheroidHalfMassRadius=6
  integer                 , parameter                            :: radiusTypeGalacticMassFraction  =7
  integer                 , parameter                            :: radiusTypeGalacticLightFraction =8

contains

  subroutine Galacticus_Output_Tree_Rotation_Curve_Initialize
    !% Initializes the module by determining whether or not rotation curve data should be output.
    use Input_Parameters
    use Galacticus_Error
    use String_Handling
    use Galacticus_Nodes
    use Galactic_Structure_Options
    use Stellar_Luminosities_Structure
    implicit none
    type     (varying_string), dimension(5) :: radiusDefinition
    type     (varying_string), dimension(3) :: fractionDefinition
    type     (varying_string)               :: valueDefinition
    character(len=20        )               :: fractionLabel     , radiusLabel
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
                call String_Split_Words(radiusDefinition,char(outputRotationCurveRadii(i)),':',bracketing="{}")
                ! Detect cases which specify radius via a mass or light fraction. In either case, extract the fraction (and filter
                ! for light fractions).
                valueDefinition=radiusDefinition(1)
                if (extract(valueDefinition,1,21) == 'galacticLightFraction') then
                   call String_Split_Words(fractionDefinition,char(valueDefinition),'{}')
                   radiusDefinition(1)='galacticLightFraction'
                end if
                if (extract(valueDefinition,1,20) == 'galacticMassFraction' ) then
                   call String_Split_Words(fractionDefinition,char(valueDefinition),'{}')
                   radiusDefinition(1)='galacticMassFraction'
                end if
                ! Parse the radius definition.
                select case (char(radiusDefinition(1)))
                case ('radius'                )
                   radii(i)%type=radiusTypeRadius
                case ('virialRadius'          )
                   radii(i)%type=radiusTypeVirialRadius
                   virialRadiusIsNeeded         =.true.
                case ('darkMatterScaleRadius' )
                   radii(i)%type=radiusTypeDarkMatterScaleRadius
                   darkMatterScaleRadiusIsNeeded=.true.
                   if (.not.defaultDarkMatterProfileComponent%scaleIsGettable         ())                                                  &
                        & call Galacticus_Error_Report                                                                                     &
                        &      (                                                                                                           &
                        &       'Galacticus_Output_Tree_Rotation_Curve_Initialize'                                                       , &
                        &       'dark matter profile scale radius is not gettable.'//                                                      &
                        &        Galacticus_Component_List(                                                                                &
                        &                                  'darkMatterProfile'                                                           , &
                        &                                   defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)  &
                        &                                  )                                                                               &
                        &      )
                case ('diskRadius'            )
                   radii(i)%type=radiusTypeDiskRadius
                   diskIsNeeded                 =.true.
                   if (.not.defaultDiskComponent             %radiusIsGettable        ())                                                  &
                        & call Galacticus_Error_Report                                                                                     &
                        &(                                                                                                                 &
                        &                              'Galacticus_Output_Tree_Rotation_Curve_Initialize'                                , &
                        &                              'disk radius is not gettable.'//                                                    &
                        &        Galacticus_Component_List(                                                                                &
                        &                                  'disk'                                                                        , &
                        &                                   defaultDiskComponent    %        radiusAttributeMatch(requireGettable=.true.)  &
                        &                                  )                                                                               &
                        &                             )
                case ('spheroidRadius'        )
                   radii(i)%type=radiusTypeSpheroidRadius
                   spheroidIsNeeded             =.true.
                   if (.not.defaultSpheroidComponent         %radiusIsGettable        ())                                                  &
                        & call Galacticus_Error_Report                                                                                     &
                        &(                                                                                                                 &
                        &                              'Galacticus_Output_Tree_Rotation_Curve_Initialize'                                , &
                        &                              'spheroid radius is not gettable.'//                                                &
                        &        Galacticus_Component_List(                                                                                &
                        &                                  'spheroid'                                                                    , &
                        &                                   defaultSpheroidComponent%        radiusAttributeMatch(requireGettable=.true.)  &
                        &                                  )                                                                               &
                        &                             )
                case ('diskHalfMassRadius'    )
                   radii(i)%type=radiusTypeDiskHalfMassRadius
                   diskIsNeeded                 =.true.
                   if (.not.defaultDiskComponent             %halfMassRadiusIsGettable())                                                  &
                        & call Galacticus_Error_Report                                                                                     &
                        &(                                                                                                                 &
                        &                              'Galacticus_Output_Tree_Rotation_Curve_Initialize'                                , &
                        &                              'disk half-mass radius is not gettable.'//                                          &
                        &        Galacticus_Component_List(                                                                                &
                        &                                  'disk'                                                                        , &
                        &                                   defaultDiskComponent    %halfMassRadiusAttributeMatch(requireGettable=.true.)  &
                        &                                  )                                                                               &
                        &                             )
                case ('spheroidHalfMassRadius')
                   radii(i)%type=radiusTypeSpheroidHalfMassRadius
                   spheroidIsNeeded   =.true.
                   if (.not.defaultSpheroidComponent         %halfMassRadiusIsGettable())                                                  &
                        & call Galacticus_Error_Report                                                                                     &
                        &(                                                                                                                 &
                        &                              'Galacticus_Output_Tree_Rotation_Curve_Initialize'                                , &
                        &                              'spheroid half-mass radius is not gettable.'//                                      &
                        &        Galacticus_Component_List(                                                                                &
                        &                                  'spheroid'                                                                    , &
                        &                                   defaultSpheroidComponent%halfMassRadiusAttributeMatch(requireGettable=.true.)  &
                        &                                  )                                                                               &
                        &                             )
                case ('galacticMassFraction'  )
                   radii(i)%type=radiusTypeGalacticMassFraction
                   fractionLabel=fractionDefinition(2)
                   read (fractionLabel,*) radii(i)%fraction
                   radii(i)%weightBy     =weightByMass
                   radii(i)%weightByIndex=weightIndexNull
                case ('galacticLightFraction' )
                   radii(i)%type=radiusTypeGalacticLightFraction
                   fractionLabel=fractionDefinition(2)
                   read (fractionLabel,*) radii(i)%fraction
                   radii(i)%weightBy      =weightByLuminosity
                   radii(i)%weightByIndex=unitStellarLuminosities%index(fractionDefinition(3))
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
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode)              , intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    integer                                                          :: i

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
       !@   <outputType>nodeData</outputType>
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
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

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
    use Galactic_Structure_Options
    use Galactic_Structure_Enclosed_Masses
    implicit none
    double precision                                , intent(in   )          :: time
    type            (treeNode                      ), intent(inout), pointer :: thisNode
    integer                                         , intent(inout)          :: doubleBufferCount         , doubleProperty, integerBufferCount, &
         &                                                                      integerProperty
    integer         (kind=kind_int8                ), intent(inout)          :: integerBuffer        (:,:)
    double precision                                , intent(inout)          :: doubleBuffer         (:,:)
    class           (nodeComponentDisk             )               , pointer :: thisDisk
    class           (nodeComponentSpheroid         )               , pointer :: thisSpheroid
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfile
    class           (darkMatterHaloScaleClass)               , pointer :: darkMatterHaloScale_
    integer                                                                  :: i
    double precision                                                         :: radius                    , radiusVirial

    ! Initialize the module.
    call Galacticus_Output_Tree_Rotation_Curve_Initialize
    ! Store property data if we are outputting rotation curve data.
    if (outputRotationCurveData) then
       ! Compute required quantities.
       darkMatterHaloScale_ => darkMatterHaloScale()
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
          doubleBuffer(doubleBufferCount,doubleProperty)=                                                                    &
               &                                         Galactic_Structure_Rotation_Curve(                                  &
               &                                                                           thisNode                        , &
               &                                                                           radius                          , &
               &                                                                           componentType=radii(i)%component, &
               &                                                                           massType=radii(i)%mass          , &
               &                                                                           haloLoaded   =radii(i)%loaded     &
               &                                                                          )
       end do
    end if
    return
  end subroutine Galacticus_Output_Tree_Rotation_Curve

end module Galacticus_Output_Trees_Rotation_Curve
