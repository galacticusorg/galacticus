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

!% Contains a module which handles outputting of velocity dispersion data to the \glc\ output file.

module Galacticus_Output_Trees_Velocity_Dispersion
  !% Handles outputting of velocity dispersion data to the \glc\ output file.
  use ISO_Varying_String
    use Galacticus_Nodes
  implicit none
  private
  public :: Galacticus_Output_Tree_Velocity_Dispersion, Galacticus_Output_Tree_Velocity_Dispersion_Property_Count,&
       & Galacticus_Output_Tree_Velocity_Dispersion_Names

  ! Flag indicating whether or not velocity dispersion information is to be output.
  logical :: outputVelocityDispersionData

  ! Flag indicating whether or not this module has been initialized.
  logical :: outputVelocityDispersionDataInitialized=.false.

  ! Radii at which velocity dispersion is to be output.
  integer :: radiiCount
  type radiusSpecifier
     type            (varying_string) :: name
     integer                          :: component            , direction    , integralWeightBy, &
          &                              integralWeightByIndex, mass         , type            , &
          &                              weightBy             , weightByIndex
     logical                          :: loaded
     double precision                 :: fraction             , value
  end type radiusSpecifier
  type            (radiusSpecifier)           , allocatable, dimension(:) :: radii
  type            (varying_string )           , allocatable, dimension(:) :: outputVelocityDispersionRadii
  logical                                                                 :: darkMatterScaleRadiusIsNeeded        , diskIsNeeded        , &
       &                                                                     spheroidIsNeeded                     , virialRadiusIsNeeded
  integer                          , parameter                            :: radiusTypeRadius                   =0
  integer                          , parameter                            :: radiusTypeVirialRadius             =1
  integer                          , parameter                            :: radiusTypeDarkMatterScaleRadius    =2
  integer                          , parameter                            :: radiusTypeDiskRadius               =3
  integer                          , parameter                            :: radiusTypeSpheroidRadius           =4
  integer                          , parameter                            :: radiusTypeDiskHalfMassRadius       =5
  integer                          , parameter                            :: radiusTypeSpheroidHalfMassRadius   =6
  integer                          , parameter                            :: radiusTypeGalacticMassFraction     =7
  integer                          , parameter                            :: radiusTypeGalacticLightFraction    =8
  integer                          , parameter                            :: directionRadial                    =0
  integer                          , parameter                            :: directionLineOfSight               =1
  integer                          , parameter                            :: directionLineOfSightInteriorAverage=2

  ! Module scope variables used in integrations.
  type            (treeNode       ), pointer                              :: activeNode
  logical                                                                 :: haloLoaded
  integer                                                                 :: componentType                        , massType            , &
       &                                                                     weightBy                             , weightIndex
  double precision                                                        :: radiusImpact                         , radiusOuter
  !$omp threadprivate(activeNode,radiusOuter,haloLoaded,massType,componentType,weightBy,weightIndex,radiusImpact)
contains

  subroutine Galacticus_Output_Tree_Velocity_Dispersion_Initialize
    !% Initializes the module by determining whether or not velocity dispersion should be output.
    use Input_Parameters
    use Galacticus_Error
    use String_Handling
    use Galactic_Structure_Options
    use Stellar_Population_Properties_Luminosities
    implicit none
    type     (varying_string), dimension(6) :: radiusDefinition
    type     (varying_string), dimension(3) :: fractionDefinition
    type     (varying_string), dimension(2) :: weightingDefinition
    type     (varying_string)               :: valueDefinition
    character(len=20        )               :: fractionLabel      , radiusLabel
    integer                                 :: i

    if (.not.outputVelocityDispersionDataInitialized) then
       !$omp critical(Galacticus_Output_Tree_Velocity_Dispersion_Initialize)
       if (.not.outputVelocityDispersionDataInitialized) then
          !@ <inputParameter>
          !@   <name>outputVelocityDispersionData</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not velocity dispersion data should be incldued in the output file.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputVelocityDispersionData',outputVelocityDispersionData,defaultValue=.false.)
          ! Read and parse parameter controlling radii of output.
          if (outputVelocityDispersionData) then
             if (.not.Input_Parameter_Is_Present('outputVelocityDispersionRadii'))                                                 &
                  & call Galacticus_Error_Report(                                                                                  &
                  &                              'Galacticus_Output_Tree_Velocity_Dispersion_Initialize'                         , &
                  &                              'outputVelocityDispersionRadii must be specified for velocity dispersion output'  &
                  &                             )
             radiiCount=Get_Input_Parameter_Array_Size('outputVelocityDispersionRadii')
             allocate(outputVelocityDispersionRadii(radiiCount))
             allocate(radii                        (radiiCount))
             !@ <inputParameter>
             !@   <name>outputVelocityDispersionRadii</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Specifies the radii at which the velocity dispersion will be output.
             !@   </description>
             !@   <type>string</type>
             !@   <cardinality>1..*</cardinality>
             !@   <group>output</group>
             !@ </inputParameter>
             call Get_Input_Parameter('outputVelocityDispersionRadii',outputVelocityDispersionRadii)
             ! Parse the radii definitions.
             diskIsNeeded                  =.false.
             spheroidIsNeeded              =.false.
             virialRadiusIsNeeded          =.false.
             darkMatterScaleRadiusIsNeeded =.false.
             do i=1,radiiCount
                radii(i)%name=outputVelocityDispersionRadii(i)
                call String_Split_Words(radiusDefinition,char(outputVelocityDispersionRadii(i)),':',bracketing="{}")
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
                   if (.not.defaultDarkMatterProfileComponent%scaleIsGettable         ())                       &
                        & call Galacticus_Error_Report(                                                         &
                        &                              'Galacticus_Output_Tree_Velocity_Dispersion_Initialize', &
                        &                              'disk half-mass radius is not gettable'                  &
                        &                             )
                case ('diskRadius'            )
                   radii(i)%type=radiusTypeDiskRadius
                   diskIsNeeded                 =.true.
                   if (.not.defaultDiskComponent             %radiusIsGettable        ())                       &
                        & call Galacticus_Error_Report(                                                         &
                        &                              'Galacticus_Output_Tree_Velocity_Dispersion_Initialize', &
                        &                              'disk radius is not gettable'                            &
                        &                             )
                case ('spheroidRadius'        )
                   radii(i)%type=radiusTypeSpheroidRadius
                   spheroidIsNeeded             =.true.
                   if (.not.defaultSpheroidComponent         %radiusIsGettable        ())                       &
                        & call Galacticus_Error_Report(                                                         &
                        &                              'Galacticus_Output_Tree_Velocity_Dispersion_Initialize', &
                        &                              'spheroid radius is not gettable'                        &
                        &                             )
                case ('diskHalfMassRadius'    )
                   radii(i)%type=radiusTypeDiskHalfMassRadius
                   diskIsNeeded                 =.true.
                   if (.not.defaultDiskComponent             %halfMassRadiusIsGettable())                       &
                        & call Galacticus_Error_Report(                                                         &
                        &                              'Galacticus_Output_Tree_Velocity_Dispersion_Initialize', &
                        &                              'disk half-mass radius is not gettable'                  &
                        &                             )
                case ('spheroidHalfMassRadius')
                   radii(i)%type=radiusTypeSpheroidHalfMassRadius
                   spheroidIsNeeded   =.true.
                   if (.not.defaultSpheroidComponent         %halfMassRadiusIsGettable())                       &
                        & call Galacticus_Error_Report(                                                         &
                        &                              'Galacticus_Output_Tree_Velocity_Dispersion_Initialize', &
                        &                              'spheroid half-mass radius is not gettable'              &
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
                   radii(i)%weightByIndex=Stellar_Population_Luminosities_Index(fractionDefinition(3))
                case default
                   call Galacticus_Error_Report('Galacticus_Output_Tree_Velocity_Dispersion_Initialize','unrecognized radius specifier')
                end select
                radii(i)%component=Galactic_Structure_Component_Type_Decode(char(radiusDefinition(2)))
                radii(i)%mass     =Galactic_Structure_Mass_Type_Decode     (char(radiusDefinition(3)))
                select case (char(radiusDefinition(4)))
                case ('loaded'  )
                   radii(i)%loaded=.true.
                case ('unloaded')
                   radii(i)%loaded=.false.
                case default
                   call Galacticus_Error_Report('Galacticus_Output_Tree_Velocity_Dispersion_Initialize','unrecognized loading specifier')
                end select
                ! Detect cases which specify the weighting for integrals over the velocity dispersion.
                valueDefinition=radiusDefinition(5)
                if (extract(valueDefinition,1,11) == 'lineOfSight') then
                   call String_Split_Words(weightingDefinition,char(valueDefinition),'{}')
                   radiusDefinition(5)=weightingDefinition(1)
                   if (weightingDefinition(2) == "mass" .or. weightingDefinition(2) == "") then
                      radii(i)%integralWeightBy     =weightByMass
                      radii(i)%integralWeightByIndex=weightIndexNull
                   else
                      radii(i)%integralWeightBy     =weightByLuminosity
                      radii(i)%integralWeightByIndex=Stellar_Population_Luminosities_Index(weightingDefinition(2))
                   end if
                end if
                ! Parse the direction definition.
                select case (char(radiusDefinition(5)))
                case ('radial'                    )
                   radii(i)%direction=directionRadial
                case ('lineOfSight'               )
                   radii(i)%direction=directionLineOfSight
                case ('lineOfSightInteriorAverage')
                   radii(i)%direction=directionLineOfSightInteriorAverage
                case default
                   call Galacticus_Error_Report('Galacticus_Output_Tree_Velocity_Dispersion_Initialize','unrecognized direction specifier')
                end select
                radiusLabel=radiusDefinition(6)
                read (radiusLabel,*) radii(i)%value
             end do
             deallocate(outputVelocityDispersionRadii)
          end if
          ! Flag that module is now initialized.
          outputVelocityDispersionDataInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Velocity_Dispersion_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Velocity_Dispersion_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Velocity_Dispersion_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Velocity_Dispersion</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Velocity_Dispersion_Names(thisNode,integerProperty,integerPropertyNames&
       &,integerPropertyComments,integerPropertyUnitsSI,doubleProperty ,doublePropertyNames,doublePropertyComments&
       &,doublePropertyUnitsSI,time)
    !% Set the names of velocity dispersion properties to be written to the \glc\ output file.
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
    call Galacticus_Output_Tree_Velocity_Dispersion_Initialize

    ! Return property names if we are outputting velocity dispersion data.
    if (outputVelocityDispersionData) then
       !@ <outputProperty>
       !@   <name>velocityDispersion</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Velocity dispersion at a given radius.</description>
       !@   <label>???</label>
       !@ </outputProperty>
       do i=1,radiiCount
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='velocityDispersion:'//char(radii(i)%name)
          doublePropertyComments(doubleProperty)='Velocity dispersion at a given radius'
          doublePropertyUnitsSI (doubleProperty)=kilo
       end do
    end if
    return
  end subroutine Galacticus_Output_Tree_Velocity_Dispersion_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Velocity_Dispersion_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Velocity_Dispersion</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Velocity_Dispersion_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of velocity dispersion properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

    ! Initialize the module.
    call Galacticus_Output_Tree_Velocity_Dispersion_Initialize

    ! Increment property count if we are outputting velocity dispersion data.
    if (outputVelocityDispersionData) doublePropertyCount=doublePropertyCount+radiiCount
    return
  end subroutine Galacticus_Output_Tree_Velocity_Dispersion_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Velocity_Dispersion</unitName>
  !#  <sortName>Galacticus_Output_Tree_Velocity_Dispersion</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Velocity_Dispersion(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store velocity dispersoin properties in the \glc\ output file buffers.
    use Kind_Numbers
    use Galactic_Structure_Velocity_Dispersions
    use Dark_Matter_Halo_Scales
    use FGSL
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    use Galactic_Structure_Options
    use Galactic_Structure_Enclosed_Masses
    implicit none
    double precision                                , intent(in   )          :: time
    type            (treeNode                      ), intent(inout), pointer :: thisNode
    integer                                         , intent(inout)          :: doubleBufferCount                , doubleProperty          , &
         &                                                                      integerBufferCount               , integerProperty
    integer         (kind=kind_int8                ), intent(inout)          :: integerBuffer        (:,:)
    double precision                                , intent(inout)          :: doubleBuffer         (:,:)
    class           (nodeComponentDisk             )               , pointer :: thisDisk
    class           (nodeComponentSpheroid         )               , pointer :: thisSpheroid
    class           (nodeComponentDarkMatterProfile)               , pointer :: thisDarkMatterProfile
    double precision                                , parameter              :: outerRadiusMultiplier     =10.0d0
    type            (c_ptr                         )                         :: parameterPointer
    type            (fgsl_function                 )                         :: integrandFunction
    type            (fgsl_integration_workspace    )                         :: integrationWorkspace
    integer                                                                  :: i
    logical                                                                  :: scaleIsZero
    double precision                                                         :: densityIntegrand                 , radius                  , &
         &                                                                      radiusFromFraction               , radiusVirial            , &
         &                                                                      radiusZero                       , velocityDensityIntegrand

    ! Initialize the module.
    call Galacticus_Output_Tree_Velocity_Dispersion_Initialize
    ! Store property data if we are outputting velocity dispersion data.
    if (outputVelocityDispersionData) then
       ! Compute required quantities.
       if (         virialRadiusIsNeeded) radiusVirial          =  Dark_Matter_Halo_Virial_Radius(thisNode                    )
       if (                 diskIsNeeded) thisDisk              =>                                thisNode%disk             ()
       if (             spheroidIsNeeded) thisSpheroid          =>                                thisNode%spheroid         ()
       if (darkMatterScaleRadiusIsNeeded) thisDarkMatterProfile =>                                thisNode%darkMatterProfile()
       do i=1,radiiCount
          ! Find the radius.
          scaleIsZero=.false.
          radius     =radii(i)%value
          select case (radii(i)%type)
          case   (radiusTypeRadius                )
             radiusOuter=    radius                                        *outerRadiusMultiplier
          case   (radiusTypeVirialRadius          )
             radius     =    radius*radiusVirial
             radiusOuter=max(radius,radiusVirial                          )*outerRadiusMultiplier
          case   (radiusTypeDarkMatterScaleRadius )
             radius     =    radius*thisDarkMatterProfile%         scale()
             radiusOuter=max(radius,thisDarkMatterProfile%         scale())*outerRadiusMultiplier
          case   (radiusTypeDiskRadius            )
             radius     =    radius*thisDisk             %        radius()
             radiusOuter=max(radius,thisDisk             %        radius())*outerRadiusMultiplier
             scaleIsZero=(thisDisk                       %        radius() <= 0.0d0)
          case   (radiusTypeSpheroidRadius        )
             radius     =    radius*thisSpheroid         %        radius()
             radiusOuter=max(radius,thisSpheroid         %        radius())*outerRadiusMultiplier
             scaleIsZero=(thisSpheroid                   %        radius() <= 0.0d0)
          case   (radiusTypeDiskHalfMassRadius    )
             radius     =    radius*thisDisk             %halfMassRadius()
             radiusOuter=max(radius,thisDisk             %halfMassRadius())*outerRadiusMultiplier
             scaleIsZero=(thisDisk                       %halfMassRadius() <= 0.0d0)
          case   (radiusTypeSpheroidHalfMassRadius)
             radius     =    radius*thisSpheroid         %halfMassRadius()
             radiusOuter=max(radius,thisSpheroid         %halfMassRadius())*outerRadiusMultiplier
             scaleIsZero=(thisSpheroid                   %halfMassRadius() <= 0.0d0)
          case   (radiusTypeGalacticMassFraction ,  &
               &  radiusTypeGalacticLightFraction )
             radiusFromFraction=                              &
                  &  Galactic_Structure_Radius_Enclosing_Mass &
                  &  (                                        &
                  &   thisNode                             ,  &
                  &   fractionalMass=radii(i)%fraction     ,  &
                  &   massType      =massTypeGalactic      ,  &
                  &   componentType =componentTypeAll      ,  &
                  &   weightBy      =radii(i)%weightBy     ,  &
                  &   weightIndex   =radii(i)%weightByIndex   &
                  &  )
             radius     =    radius*radiusFromFraction
             radiusOuter=max(radius,radiusFromFraction)*outerRadiusMultiplier
          end select
          ! Store the rotation curve.
          doubleProperty=doubleProperty+1
          if (scaleIsZero) then
             ! Do not compute dispersions if the component scale is zero.
             doubleBuffer(doubleBufferCount,doubleProperty)=0.0d0
          else
             select case (radii(i)%direction)
             case (directionRadial                    )
                ! Radial velocity dispersion.
                doubleBuffer(doubleBufferCount,doubleProperty)=                                                                   &
                     &                                   Galactic_Structure_Velocity_Dispersion(                                  &
                     &                                                                          thisNode                        , &
                     &                                                                          radius                          , &
                     &                                                                          radiusOuter                     , &
                     &                                                                          massType     =radii(i)%mass     , &
                     &                                                                          componentType=radii(i)%component, &
                     &                                                                          haloLoaded   =radii(i)%loaded     &
                     &                                                                         )
             case (directionLineOfSight               )
                ! Line-of-sight velocity dispersion.
                activeNode   => thisNode
                massType     =  radii(i)%mass
                componentType=  radii(i)%component
                haloLoaded   =  radii(i)%loaded
                weightBy     =  radii(i)%integralWeightBy
                weightIndex  =  radii(i)%integralWeightByIndex
                radiusImpact =  radius
                doubleBuffer(doubleBufferCount,doubleProperty)=Galacticus_Output_Trees_Line_of_Sight_Velocity_Dispersion(radius)
             case (directionLineOfSightInteriorAverage)
                ! Average over the line-of-sight velocity dispersion within the radius.
                activeNode   => thisNode
                massType     =  radii(i)%mass
                componentType=  radii(i)%component
                haloLoaded   =  radii(i)%loaded
                weightBy     =  radii(i)%integralWeightBy
                weightIndex  =  radii(i)%integralWeightByIndex
                radiusZero   =  0.0d0
                radiusImpact =  radius
                velocityDensityIntegrand=Integrate(radiusZero,radiusOuter&
                     &,Galacticus_Output_Trees_Vlcty_Dsprsn_Vlcty_Dnsty_Srfc_Intgrnd ,parameterPointer ,integrandFunction&
                     &,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
                call Integrate_Done(integrandFunction,integrationWorkspace)
                densityIntegrand        =Integrate(radiusZero,radiusOuter&
                     &,Galacticus_Output_Trees_Vlcty_Dsprsn_Dnsty_Srfc_Intgrnd ,parameterPointer ,integrandFunction&
                     &,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
                call Integrate_Done(integrandFunction,integrationWorkspace)
                if (velocityDensityIntegrand <= 0.0d0) then
                   doubleBuffer(doubleBufferCount,doubleProperty)=0.0d0
                else
                   doubleBuffer(doubleBufferCount,doubleProperty)=sqrt(velocityDensityIntegrand/densityIntegrand)
                end if
             end select
          end if
       end do
    end if
    return
  end subroutine Galacticus_Output_Tree_Velocity_Dispersion

  function Galacticus_Output_Trees_Vlcty_Dsprsn_Vlcty_Dnsty_Srfc_Intgrnd(radius,parameterPointer) bind(c)
    !% Integrand function used for integrating line-of-sight velocity dispersion over surface density.
    use, intrinsic :: ISO_C_Binding
    use Galactic_Structure_Densities
    use Galactic_Structure_Velocity_Dispersions
    implicit none
    real(kind=c_double             )        :: Galacticus_Output_Trees_Vlcty_Dsprsn_Vlcty_Dnsty_Srfc_Intgrnd
    real(kind=c_double             ), value :: radius
    type(c_ptr                     ), value :: parameterPointer

    if (radius <= 0.0d0) then
       Galacticus_Output_Trees_Vlcty_Dsprsn_Vlcty_Dnsty_Srfc_Intgrnd=0.0d0
    else
       Galacticus_Output_Trees_Vlcty_Dsprsn_Vlcty_Dnsty_Srfc_Intgrnd=                                    &
            &                        Spherical_Shell_Solid_Angle_In_Cylcinder(radius)                    &
            &                       *                                         radius **2                 &
            &                       *Galactic_Structure_Density            (                             &
            &                                                               activeNode                 , &
            &                                                               [radius,0.0d0,0.0d0]       , &
            &                                                               massType     =massType     , &
            &                                                               componentType=componentType, &
            &                                                               weightBy     =weightBy     , &
            &                                                               weightIndex  =weightIndex  , &
            &                                                               haloLoaded   =haloLoaded     &
            &                                                              )                             &
            &                       *Galactic_Structure_Velocity_Dispersion(                             &
            &                                                               activeNode                 , &
            &                                                               radius                     , &
            &                                                               radiusOuter                , &
            &                                                               massType     =massType     , &
            &                                                               componentType=componentType, &
            &                                                               haloLoaded   =haloLoaded     &
            &                                                              )**2
    end if
   return
  end function Galacticus_Output_Trees_Vlcty_Dsprsn_Vlcty_Dnsty_Srfc_Intgrnd

  function Galacticus_Output_Trees_Vlcty_Dsprsn_Dnsty_Srfc_Intgrnd(radius,parameterPointer) bind(c)
    !% Integrand function used for integrating line-of-sight surface density dispersion over area.
    use, intrinsic :: ISO_C_Binding
    use Galactic_Structure_Densities
    implicit none
    real(kind=c_double             )        :: Galacticus_Output_Trees_Vlcty_Dsprsn_Dnsty_Srfc_Intgrnd
    real(kind=c_double             ), value :: radius
    type(c_ptr                     ), value :: parameterPointer

    if (radius <= 0.0d0) then
       Galacticus_Output_Trees_Vlcty_Dsprsn_Dnsty_Srfc_Intgrnd=0.0d0
    else
       Galacticus_Output_Trees_Vlcty_Dsprsn_Dnsty_Srfc_Intgrnd=                               &
            &                         Spherical_Shell_Solid_Angle_In_Cylcinder(radius)        &
            &                        *                                         radius **2     &
            &                        *Galactic_Structure_Density(                             &
            &                                                    activeNode                 , &
            &                                                    [radius,0.0d0,0.0d0]       , &
            &                                                    massType     =massType     , &
            &                                                    componentType=componentType, &
            &                                                    weightBy     =weightBy     , &
            &                                                    weightIndex  =weightIndex  , &
            &                                                    haloLoaded   =haloLoaded     &
            &                                                   )
    end if
    return
  end function Galacticus_Output_Trees_Vlcty_Dsprsn_Dnsty_Srfc_Intgrnd

  double precision function Spherical_Shell_Solid_Angle_In_Cylcinder(radius)
    !% Computes the solid angle of a spherical shelll of given {\tt radius} that lies within a cylinder of radius {\tt
    !% radiusImpact}.
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in   ) :: radius

    ! Test size of sphere relative to cylinder.
    if (radius <= radiusImpact) then
       ! Sphere is entirely within the cylinder.
       Spherical_Shell_Solid_Angle_In_Cylcinder=4.0d0*Pi
    else
       ! Some of sphere lies outside of cylinder.
       Spherical_Shell_Solid_Angle_In_Cylcinder=4.0d0*Pi*(1.0d0-cos(asin(radiusImpact/radius)))
    end if
    return
  end function Spherical_Shell_Solid_Angle_In_Cylcinder

  double precision function Galacticus_Output_Trees_Line_of_Sight_Velocity_Dispersion(radius)
    !% Compute the line-of-sight velocity dispersion at the given {\tt radius}.
    use FGSL
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    implicit none
    double precision                            , intent(in   ) :: radius
    type            (c_ptr                     )                :: parameterPointer
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    double precision                                            :: densityIntegral     , velocityDensityIntegral

    velocityDensityIntegral=Integrate(radius,radiusOuter&
         &,Galacticus_Output_Trees_Vlcty_Dispersion_Vlcty_Dnsty_Intgrnd ,parameterPointer ,integrandFunction&
         &,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    densityIntegral        =Integrate(radius,radiusOuter,Galacticus_Output_Trees_Velocity_Dispersion_Density_Integrand&
         &,parameterPointer ,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    if (velocityDensityIntegral <= 0.0d0) then
       Galacticus_Output_Trees_Line_of_Sight_Velocity_Dispersion=0.0d0
    else
       Galacticus_Output_Trees_Line_of_Sight_Velocity_Dispersion=sqrt(velocityDensityIntegral/densityIntegral)
    end if
    return
  end function Galacticus_Output_Trees_Line_of_Sight_Velocity_Dispersion

  function Galacticus_Output_Trees_Velocity_Dispersion_Density_Integrand(radius,parameterPointer) bind(c)
    !% Integrand function used for computing line-of-sight velocity dispersions.
    use, intrinsic :: ISO_C_Binding
    use Galactic_Structure_Densities
    implicit none
    real(kind=c_double)        :: Galacticus_Output_Trees_Velocity_Dispersion_Density_Integrand
    real(kind=c_double), value :: radius
    type(c_ptr        ), value :: parameterPointer

    if (radius == 0.0d0) then
       Galacticus_Output_Trees_Velocity_Dispersion_Density_Integrand=0.0d0
    else
       Galacticus_Output_Trees_Velocity_Dispersion_Density_Integrand=                                    &
            &                        Galactic_Structure_Density           (                              &
            &                                                               activeNode                 , &
            &                                                               [radius,0.0d0,0.0d0]       , &
            &                                                               massType     =massType     , &
            &                                                               componentType=componentType, &
            &                                                               weightBy     =weightBy     , &
            &                                                               weightIndex  =weightIndex  , &
            &                                                               haloLoaded   =haloLoaded     &
            &                                                              )                             &
            &                       *     radius                                                         &
            &                       /sqrt(radius**2-radiusImpact**2)
    end if
    return
  end function Galacticus_Output_Trees_Velocity_Dispersion_Density_Integrand

  function Galacticus_Output_Trees_Vlcty_Dispersion_Vlcty_Dnsty_Intgrnd(radius,parameterPointer) bind(c)
    !% Integrand function used for computing line-of-sight velocity dispersions.
    use, intrinsic :: ISO_C_Binding
    use Galactic_Structure_Densities
    use Galactic_Structure_Velocity_Dispersions
    implicit none
    real(kind=c_double)        :: Galacticus_Output_Trees_Vlcty_Dispersion_Vlcty_Dnsty_Intgrnd
    real(kind=c_double), value :: radius
    type(c_ptr        ), value :: parameterPointer

    if (radius <= radiusImpact) then
       Galacticus_Output_Trees_Vlcty_Dispersion_Vlcty_Dnsty_Intgrnd=0.0d0
    else
       Galacticus_Output_Trees_Vlcty_Dispersion_Vlcty_Dnsty_Intgrnd=                                     &
            &                        Galactic_Structure_Density           (                              &
            &                                                               activeNode                 , &
            &                                                               [radius,0.0d0,0.0d0]       , &
            &                                                               massType     =massType     , &
            &                                                               componentType=componentType, &
            &                                                               weightBy     =weightBy     , &
            &                                                               weightIndex  =weightIndex  , &
            &                                                               haloLoaded   =haloLoaded     &
            &                                                              )                             &
            &                       *Galactic_Structure_Velocity_Dispersion(                             &
            &                                                               activeNode                 , &
            &                                                               radius                     , &
            &                                                               radiusOuter                , &
            &                                                               massType     =massType     , &
            &                                                               componentType=componentType, &
            &                                                               haloLoaded   =haloLoaded     &
            &                                                              )**2                          &
            &                       *     radius                                                         &
            &                       /sqrt(radius**2-radiusImpact**2)
    end if
    return
  end function Galacticus_Output_Trees_Vlcty_Dispersion_Vlcty_Dnsty_Intgrnd

end module Galacticus_Output_Trees_Velocity_Dispersion
