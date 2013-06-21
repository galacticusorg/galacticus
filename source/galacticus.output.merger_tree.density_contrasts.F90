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

!% Contains a module which handles outputting of node density contrast properties (radii and masses).

module Galacticus_Output_Trees_Density_Contrasts
  !% Handles outputting of node density contrast properties (radii and masses).
  use Galacticus_Nodes
  use, intrinsic :: ISO_C_Binding
  implicit none
  private
  public :: Galacticus_Output_Tree_Density_Contrast, Galacticus_Output_Tree_Density_Contrast_Property_Count, Galacticus_Output_Tree_Density_Contrast_Names

  ! Number of overdensit properties.
  integer                                                         :: densityContrastPropertyCount                                                       
  
  ! Flag indicating whether or not density contrast information is to be output.
  logical                                                         :: outputDensityContrastData                   , outputDensityContrastDataDarkOnly, & 
       &                                                             outputDensityContrastHaloLoaded                                                    
  
  ! Array of density contrasts for which to output data.
  integer                                                         :: densityContrastCount                                                               
  double precision                    , allocatable, dimension(:) :: outputDensityContrastValues                                                        
  
  ! Flag indicating whether or not this module has been initialized.
  logical                                                         :: outputDensityContrastDataInitialized=.false.                                       
  
  ! Label used to specify for what type of mass (dark or all) density contrasts are required.
  integer                                                         :: massTypeSelected                                                                   
  
  ! Reference density from which contrasts are defined.
  double precision                                                :: referenceDensity                                                                   
  
  ! Global variables used in root finding.
  type            (treeNode          ), pointer                   :: activeNode                                                                         
  class           (nodeComponentBasic), pointer                   :: activeBasicComponent                                                               
  double precision                                                :: meanDensityContrastTarget                                                          
  !$omp threadprivate(activeNode,activeBasicComponent,meanDensityContrastTarget)
contains

  subroutine Galacticus_Output_Tree_Density_Contrast_Initialize
    !% Initializes the module by determining whether or not density contrast data should be output and getting a list of overdensities.
    use Input_Parameters
    use Memory_Management
    use Galactic_Structure_Options
    use Cosmological_Parameters
    use Cosmology_Functions
    implicit none

    if (.not.outputDensityContrastDataInitialized) then
       !$omp critical(Galacticus_Output_Tree_Density_Contrast_Initialize)
       if (.not.outputDensityContrastDataInitialized) then
          !@ <inputParameter>
          !@   <name>outputDensityContrastData</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not density contrast data (i.e. radius and mass at a given density contrast) should be included in the output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputDensityContrastData',outputDensityContrastData,defaultValue=.false.)
          !@ <inputParameter>
          !@   <name>outputDensityContrastDataDarkOnly</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not density contrast data should be computed using the dark matter component alone.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputDensityContrastDataDarkOnly',outputDensityContrastDataDarkOnly,defaultValue=.false.)
          select case (outputDensityContrastDataDarkOnly)
          case (.true.)
             massTypeSelected=massTypeDark
             referenceDensity=(Omega_Matter()-Omega_b())*Critical_Density()
          case (.false.)
             massTypeSelected=massTypeAll
             referenceDensity= Omega_Matter()           *Critical_Density()
          end select
          
          ! Read density contrast values if necessary.
          if (outputDensityContrastData) then
             densityContrastCount=Get_Input_Parameter_Array_Size('outputDensityContrastValues')
             densityContrastPropertyCount=2*densityContrastCount
             call Alloc_Array(outputDensityContrastValues,[densityContrastCount])
             !@ <inputParameter>
             !@   <name>outputDensityContrastValues</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     A list of density contrasts at which to output data.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1..*</cardinality>
             !@   <group>output</group>
             !@ </inputParameter>
             call Get_Input_Parameter('outputDensityContrastValues',outputDensityContrastValues)
             !@ <inputParameter>
             !@   <name>outputDensityContrastHaloLoaded</name>
             !@   <defaultValue></defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Specifies whether baryonic loading of the halo should be accounted for when outputting density contrast data.
             !@   </description>
             !@   <type>boolean</type>
             !@   <cardinality>1..*</cardinality>
             !@   <group>output</group>
             !@ </inputParameter>
             call Get_Input_Parameter('outputDensityContrastHaloLoaded',outputDensityContrastHaloLoaded,defaultValue=.true.)
          end if
          
          ! Flag that module is now initialized.
          outputDensityContrastDataInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Density_Contrast_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Density_Contrast_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Density_Contrast_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Density_Contrast</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Density_Contrast_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty ,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of density contrast properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode)              , intent(inout), pointer :: thisNode                                           
    double precision                        , intent(in   )          :: time                                               
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty           
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , & 
         &                                                              integerPropertyComments, integerPropertyNames      
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI    
    integer                                                          :: iDensity                                           
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Density_Contrast_Initialize

    ! Return property names if we are outputting density contrast data.
    if (outputDensityContrastData) then
       do iDensity=1,densityContrastCount
          doubleProperty=doubleProperty+1
          write (doublePropertyNames   (doubleProperty),'(a,f5.1  )') 'nodeRadius'                             ,outputDensityContrastValues(iDensity)
          write (doublePropertyComments(doubleProperty),'(a,f5.1,a)') 'Radius enclosing a density contrast of ',outputDensityContrastValues(iDensity),' [Mpc].'
          doublePropertyUnitsSI(doubleProperty)=megaParsec
          doubleProperty=doubleProperty+1
          write (doublePropertyNames   (doubleProperty),'(a,f5.1  )') 'nodeMass'                               ,outputDensityContrastValues(iDensity)
          write (doublePropertyComments(doubleProperty),'(a,f5.1,a)') 'Mass within a density contrast of '     ,outputDensityContrastValues(iDensity),' [Msolar].'
          doublePropertyUnitsSI(doubleProperty)=massSolar
       end do
    end if
    return
  end subroutine Galacticus_Output_Tree_Density_Contrast_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Density_Contrast_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Density_Contrast</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Density_Contrast_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of density contrast properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode                                  
    double precision          , intent(in   )          :: time                                      
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount 
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Density_Contrast_Initialize

    ! Increment property count if we are outputting density contrast data.
    if (outputDensityContrastData) doublePropertyCount=doublePropertyCount+densityContrastPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Density_Contrast_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Density_Contrast</unitName>
  !#  <sortName>Galacticus_Output_Tree_Density_Contrast</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Density_Contrast(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store density contrast properties in the \glc\ output file buffers.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Dark_Matter_Halo_Scales
    use Kind_Numbers
    use Root_Finder
    use FGSL
    implicit none
    double precision                , intent(in   )          :: time                                                       
    type            (treeNode      ), intent(inout), pointer :: thisNode                                                   
    integer                         , intent(inout)          :: doubleBufferCount            , doubleProperty          , & 
         &                                                      integerBufferCount           , integerProperty             
    integer         (kind=kind_int8), intent(inout)          :: integerBuffer     (:,:)                                    
    double precision                , intent(inout)          :: doubleBuffer      (:,:)                                    
    double precision                , parameter              :: toleranceAbsolute      =0.0d0, toleranceRelative=1.0d-3    
    type            (rootFinder    ), save                   :: finder                                                     
    !$omp threadprivate(finder)
    integer                                                  :: iDensity                                                   
    double precision                                         :: enclosedMass                 , radius                  , & 
         &                                                      radiusMaximum                , radiusMinimum               
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Density_Contrast_Initialize
    ! Store property data if we are outputting density contrast data.
    if (outputDensityContrastData) then
       ! Initialize our root finder.
       if (.not.finder%isInitialized()) then
          call finder%rangeExpand (                                                             &
               &                   rangeExpandDownward          =0.5d0                        , &
               &                   rangeExpandUpward            =2.0d0                        , &
               &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
               &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
               &                   rangeExpandType              =rangeExpandMultiplicative      &
               &                  )
          call finder%rootFunction(Mean_Density_Contrast_Root         )
          call finder%tolerance   (toleranceAbsolute,toleranceRelative)
       end if       
       ! Make the node and basic component available to the root finding routine.
       activeNode           => thisNode
       activeBasicComponent => thisNode%basic()
       ! Loop over all requested densities.
       do iDensity=1,densityContrastCount
          meanDensityContrastTarget=outputDensityContrastValues(iDensity)
          ! Locate the root.
          radius=finder%find(rootGuess=Dark_Matter_Halo_Virial_Radius(thisNode))
          ! Compute the mass enclosed in this radius.
          enclosedMass=Galactic_Structure_Enclosed_Mass(                                               &
               &                                        thisNode                                     , &
               &                                        radius                                       , &
               &                                        componentType=componentTypeAll               , &
               &                                        massType     =massTypeSelected               , &
               &                                        haloLoaded   =outputDensityContrastHaloLoaded  &
               &                                       )
          ! Store the resulting radius and mass.
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=radius
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=enclosedMass
       end do
    end if
    return
  end subroutine Galacticus_Output_Tree_Density_Contrast

  double precision function Mean_Density_Contrast_Root(radius)
    !% Root function used in finding the radius that encloses a given density contrast.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Cosmology_Functions
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in   ) :: radius       
    double precision                :: enclosedMass 
    
    ! Solve for the radius enclosing the specified density contrast.
    enclosedMass              =Galactic_Structure_Enclosed_Mass(activeNode,radius,componentType&
         &=componentTypeAll,massType=massTypeSelected,haloLoaded=outputDensityContrastHaloLoaded)
    Mean_Density_Contrast_Root=3.0d0*enclosedMass/4.0d0/Pi/radius**3/(referenceDensity &
         &/Expansion_Factor(activeBasicComponent%time())**3)-meanDensityContrastTarget
    return
  end function Mean_Density_Contrast_Root
  
end module Galacticus_Output_Trees_Density_Contrasts
