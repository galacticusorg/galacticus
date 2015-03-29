!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

module Galacticus_Output_Trees_Projected_Density
  !% Handles outputting of projected density data to the \glc\ output file.
  use ISO_Varying_String
  implicit none
  private
  public :: Galacticus_Output_Tree_Projected_Density, Galacticus_Output_Tree_Projected_Density_Property_Count, &
       &    Galacticus_Output_Tree_Projected_Density_Names

  ! Flag indicating whether or not projected density information is to be output.
  logical :: outputProjectedDensityData

  ! Flag indicating whether or not this module has been initialized.
  logical :: outputProjectedDensityDataInitialized=.false.

  ! Radii at which projected density is to be output.
  integer :: radiiCount
  type radiusSpecifier
     type            (varying_string) :: name
     integer                          :: component    , mass , type, weightBy, &
          &                              weightByIndex
     logical                          :: loaded
     double precision                 :: fraction     , value
  end type radiusSpecifier
  type   (radiusSpecifier)           , allocatable, dimension(:) :: radii
  type   (varying_string )           , allocatable, dimension(:) :: outputProjectedDensityRadii
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

  subroutine Galacticus_Output_Tree_Projected_Density_Initialize()
    !% Initializes the module by determining whether or not projected density data should be output.
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

    if (.not.outputProjectedDensityDataInitialized) then
       !$omp critical(Galacticus_Output_Tree_Projected_Density_Initialize)
       if (.not.outputProjectedDensityDataInitialized) then
          !@ <inputParameter>
          !@   <name>outputProjectedDensityData</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not projected density data should be incldued in the output file.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputProjectedDensityData',outputProjectedDensityData,defaultValue=.false.)
          ! Read and parse parameter controlling radii of output.
          if (outputProjectedDensityData) then
             if (.not.Input_Parameter_Is_Present('outputProjectedDensityRadii'))                                               &
                  & call Galacticus_Error_Report(                                                                              &
                  &                              'Galacticus_Output_Tree_Projected_Density_Initialize'                       , &
                  &                              'outputProjectedDensityRadii must be specified for projected density output'  &
                  &                             )
             radiiCount=Get_Input_Parameter_Array_Size('outputProjectedDensityRadii')
             allocate(outputProjectedDensityRadii(radiiCount))
             allocate(radii                      (radiiCount))
             !@ <inputParameter>
             !@   <name>outputProjectedDensityRadii</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Specifies the radii at which the projected density will be output.
             !@   </description>
             !@   <type>string</type>
             !@   <cardinality>1..*</cardinality>
             !@   <group>output</group>
             !@ </inputParameter>
             call Get_Input_Parameter('outputProjectedDensityRadii',outputProjectedDensityRadii)
             ! Parse the radii definitions.
             diskIsNeeded                  =.false.
             spheroidIsNeeded              =.false.
             virialRadiusIsNeeded          =.false.
             darkMatterScaleRadiusIsNeeded =.false.
             do i=1,radiiCount
                radii(i)%name=outputProjectedDensityRadii(i)
                call String_Split_Words(radiusDefinition,char(outputProjectedDensityRadii(i)),':',bracketing="{}")
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
                        &       'Galacticus_Output_Tree_Projected_Density_Initialize'                                                    , &
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
                        &                              'Galacticus_Output_Tree_Projected_Density_Initialize'                             , &
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
                        &                              'Galacticus_Output_Tree_Projected_Density_Initialize'                             , &
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
                        &                              'Galacticus_Output_Tree_Projected_Density_Initialize'                             , &
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
                        &                              'Galacticus_Output_Tree_Projected_Density_Initialize'                             , &
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
                   call Galacticus_Error_Report('Galacticus_Output_Tree_Projected_Density_Initialize','unrecognized radius specifier "'//char(radiusDefinition(1))//'"')
                end select
                radii(i)%component=Galactic_Structure_Component_Type_Decode(char(radiusDefinition(2)))
                radii(i)%mass     =Galactic_Structure_Mass_Type_Decode     (char(radiusDefinition(3)))
                select case (char(radiusDefinition(4)))
                case ('loaded'  )
                   radii(i)%loaded=.true.
                case ('unloaded')
                   radii(i)%loaded=.false.
                case default
                   call Galacticus_Error_Report('Galacticus_Output_Tree_Projected_Density_Initialize','unrecognized loading specifier')
                end select
                radiusLabel=radiusDefinition(5)
                read (radiusLabel,*) radii(i)%value
             end do
             deallocate(outputProjectedDensityRadii)
          end if
          ! Flag that module is now initialized.
          outputProjectedDensityDataInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Projected_Density_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Projected_Density_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Projected_Density_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Projected_Density</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Projected_Density_Names(node,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of projected density properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode)              , intent(inout), pointer :: node
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    integer                                                          :: i

    ! Initialize the module.
    call Galacticus_Output_Tree_Projected_Density_Initialize

    ! Return property names if we are outputting projected density data.
    if (outputProjectedDensityData) then
       !@ <outputProperty>
       !@   <name>projectedDensity</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>projected density at a given radius.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       do i=1,radiiCount
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='projectedDensity:'//char(radii(i)%name)
          doublePropertyComments(doubleProperty)='projected density at a given radius'
          doublePropertyUnitsSI (doubleProperty)=kilo
       end do
    end if
    return
  end subroutine Galacticus_Output_Tree_Projected_Density_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Projected_Density_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Projected_Density</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Projected_Density_Property_Count(node,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of projected density properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    implicit none
    type            (treeNode), intent(inout), pointer :: node
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

    ! Initialize the module.
    call Galacticus_Output_Tree_Projected_Density_Initialize

    ! Increment property count if we are outputting projected density data.
    if (outputProjectedDensityData) doublePropertyCount=doublePropertyCount+radiiCount
    return
  end subroutine Galacticus_Output_Tree_Projected_Density_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Projected_Density</unitName>
  !#  <sortName>Galacticus_Output_Tree_Projected_Density</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Projected_Density(node,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store projected density properties in the \glc\ output file buffers.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Numerical_Integration
    use Galacticus_Nodes
    use Kind_Numbers
    use Galactic_Structure_Densities
    use Dark_Matter_Halo_Scales
    use Galactic_Structure_Options
    use Galactic_Structure_Enclosed_Masses
    implicit none
    double precision                                , intent(in   )          :: time
    type            (treeNode                      ), intent(inout), pointer :: node
    integer                                         , intent(inout)          :: doubleBufferCount         , doubleProperty, integerBufferCount, &
         &                                                                      integerProperty
    integer         (kind=kind_int8                ), intent(inout)          :: integerBuffer        (:,:)
    double precision                                , intent(inout)          :: doubleBuffer         (:,:)
    class           (nodeComponentDisk             )               , pointer :: disk
    class           (nodeComponentSpheroid         )               , pointer :: spheroid
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile
    class           (darkMatterHaloScaleClass      )               , pointer :: darkMatterHaloScale_
    integer                                                                  :: i
    double precision                                                         :: radiusProjected           , radiusOuter   , radiusVirial
    type            (c_ptr                         )                         :: parameterPointer
    type            (fgsl_function                 )                         :: integrandFunction
    type            (fgsl_integration_workspace    )                         :: integrationWorkspace

    ! Initialize the module.
    call Galacticus_Output_Tree_Projected_Density_Initialize()
    ! Store property data if we are outputting projected density data.
    if (outputProjectedDensityData) then
       ! Compute required quantities.
       darkMatterHaloScale_ => darkMatterHaloScale()
       if (         virialRadiusIsNeeded) radiusVirial      =  darkMatterHaloScale_%virialRadius(node                    )
       if (                 diskIsNeeded) disk              =>                                   node%disk             ()
       if (             spheroidIsNeeded) spheroid          =>                                   node%spheroid         ()
       if (darkMatterScaleRadiusIsNeeded) darkMatterProfile =>                                   node%darkMatterProfile()
       do i=1,radiiCount
          ! Find the projected radius.
          radiusProjected=radii(i)%value
          select case (radii(i)%type)
          case   (radiusTypeRadius                )
             ! Nothing to do.
          case   (radiusTypeVirialRadius          )
             radiusProjected=+radiusProjected*radiusVirial
          case   (radiusTypeDarkMatterScaleRadius )
             radiusProjected=+radiusProjected*darkMatterProfile%         scale()
          case   (radiusTypeDiskRadius            )
             radiusProjected=+radiusProjected*disk             %        radius()
          case   (radiusTypeSpheroidRadius        )
             radiusProjected=+radiusProjected*spheroid         %        radius()
          case   (radiusTypeDiskHalfMassRadius    )
             radiusProjected=+radiusProjected*disk             %halfMassRadius()
          case   (radiusTypeSpheroidHalfMassRadius)
             radiusProjected=+radiusProjected*spheroid         %halfMassRadius()
          case   (radiusTypeGalacticMassFraction ,  &
               &  radiusTypeGalacticLightFraction )
             radiusProjected=+radiusProjected                          &
                  &          *Galactic_Structure_Radius_Enclosing_Mass &
                  &           (                                        &
                  &            node                                 ,  &
                  &            fractionalMass=radii(i)%fraction     ,  &
                  &            massType      =massTypeGalactic      ,  &
                  &            componentType =componentTypeAll      ,  &
                  &            weightBy      =radii(i)%weightBy     ,  &
                  &            weightIndex   =radii(i)%weightByIndex   &
                  &           )
          end select
          ! Find the outer radius.
          radiusOuter=darkMatterHaloScale_%virialRadius(node)
          ! Store the projected density.
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=Integrate(                                                           &
               &                                                   radiusProjected                                          , &
               &                                                   radiusOuter                                              , &
               &                                                   Galacticus_Output_Tree_Projected_Density_Integrand       , &
               &                                                   parameterPointer                                         , &
               &                                                   integrandFunction                                        , &
               &                                                   integrationWorkspace                                     , &
               &                                                   toleranceAbsolute                                 =0.0d+0, &
               &                                                   toleranceRelative                                 =1.0d-3  &
               &                                                  )
          call Integrate_Done(integrandFunction,integrationWorkspace)
       end do
    end if
    return
    
  contains
    
    function Galacticus_Output_Tree_Projected_Density_Integrand(radius,parameterPointer) bind(c)
      !% Integrand function used for computing projected densities.
      use, intrinsic :: ISO_C_Binding
      use Galactic_Structure_Densities
      implicit none
      real(kind=c_double)        :: Galacticus_Output_Tree_Projected_Density_Integrand
      real(kind=c_double), value :: radius
      type(c_ptr        ), value :: parameterPointer

      if (radius <= radiusProjected) then
         Galacticus_Output_Tree_Projected_Density_Integrand=+0.0d0
      else
         Galacticus_Output_Tree_Projected_Density_Integrand=+2.0d0                                                        &
              &                                             *radius                                                       &
              &                                             /sqrt(                                                        &
              &                                                   +radius         **2                                     &
              &                                                   -radiusProjected**2                                     &
              &                                             )                                                             &
              &                                             *Galactic_Structure_Density(                                  &
              &                                                                         node                            , &
              &                                                                         [                                 &
              &                                                                          radius                         , &
              &                                                                          0.0d0                          , &
              &                                                                          0.0d0                            &
              &                                                                         ]                               , &
              &                                                                         componentType=radii(i)%component, &
              &                                                                         massType     =radii(i)%mass     , &
              &                                                                         haloLoaded   =radii(i)%loaded     &
              &                                                                        )
      end if
      return
    end function Galacticus_Output_Tree_Projected_Density_Integrand
    
  end subroutine Galacticus_Output_Tree_Projected_Density

end module Galacticus_Output_Trees_Projected_Density
