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


!% Contains a module which handles outputting of node density contrast properties (radii and masses).

module Galacticus_Output_Trees_Density_Contrasts
  !% Handles outputting of node density contrast properties (radii and masses).
  use Tree_Nodes
  use, intrinsic :: ISO_C_Binding
  implicit none
  private
  public :: Galacticus_Output_Tree_Density_Contrast, Galacticus_Output_Tree_Density_Contrast_Property_Count, Galacticus_Output_Tree_Density_Contrast_Names

  ! Number of overdensit properties.
  integer                                     :: densityContrastPropertyCount

  ! Flag indicating whether or not density contrast information is to be output.
  logical                                     :: outputDensityContrastData,outputDensityContrastDataDarkOnly

  ! Array of density contrasts for which to output data.
  integer                                     :: densityContrastCount
  double precision, allocatable, dimension(:) :: outputDensityContrastValues

  ! Flag indicating whether or not this module has been initialized.
  logical                                     :: outputDensityContrastDataInitialized=.false.

  ! Label used to specify for what type of mass (dark or all) density contrasts are required.
  integer                                     :: massTypeSelected

  ! Reference density from which contrasts are defined.
  double precision                            :: referenceDensity

  ! Global variables used in root finding.
  type(treeNode),   pointer                   :: activeNode
  double precision                            :: meanDensityContrastTarget
  !$omp threadprivate(activeNode,meanDensityContrastTarget)

contains

  subroutine Galacticus_Output_Tree_Density_Contrast_Initialize
    !% Initializes the module by determining whether or not density contrast data should be output and getting a list of overdensities.
    use Input_Parameters
    use Memory_Management
    use Galactic_Structure_Options
    use Cosmological_Parameters
    use Cosmology_Functions
    implicit none

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
       end if

       ! Flag that module is now initialized.
       outputDensityContrastDataInitialized=.true.
    end if
    !$omp end critical(Galacticus_Output_Tree_Density_Contrast_Initialize)
    return
  end subroutine Galacticus_Output_Tree_Density_Contrast_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Density_Contrast_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Density_Contrast</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Density_Contrast_Names(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of density contrast properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI
    integer                                       :: iDensity

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
  subroutine Galacticus_Output_Tree_Density_Contrast_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of density contrast properties to be written to the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount
    
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
    double precision,        intent(in)             :: time
    type(treeNode),          intent(inout), pointer :: thisNode
    integer,                 intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)          :: integerBuffer(:,:)
    double precision,        intent(inout)          :: doubleBuffer(:,:)
    type(fgsl_function),     save                   :: rootFunction
    type(fgsl_root_fsolver), save                   :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    double precision,        parameter              :: toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3
    type(c_ptr)                                     :: parameterPointer
    integer                                         :: iDensity
    double precision                                :: radius,radiusMinimum,radiusMaximum,enclosedMass

    ! Initialize the module.
    call Galacticus_Output_Tree_Density_Contrast_Initialize

    ! Store property data if we are outputting density contrast data.
    if (outputDensityContrastData) then

       ! Make the node available to the root finding routine.
       activeNode => thisNode

       ! Loop over all requested densities.
       do iDensity=1,densityContrastCount
          meanDensityContrastTarget=outputDensityContrastValues(iDensity)
          
          ! Find suitable minimum and maximum radii that bound the sought root.
          radiusMinimum=Dark_Matter_Halo_Virial_Radius(thisNode)
          radiusMaximum=Dark_Matter_Halo_Virial_Radius(thisNode)
          do while (Mean_Density_Contrast_Root(radiusMinimum,parameterPointer) < 0.0d0)
             radiusMinimum=0.5d0*radiusMinimum
          end do
          do while (Mean_Density_Contrast_Root(radiusMaximum,parameterPointer) > 0.0d0)
             radiusMaximum=2.0d0*radiusMaximum
          end do
          
          ! Locate the root.
          radius=Root_Find(radiusMinimum,radiusMaximum,Mean_Density_Contrast_Root,parameterPointer,rootFunction,rootFunctionSolver &
               &,toleranceAbsolute,toleranceRelative)
          
          ! Compute the mass enclosed in this radius.
          enclosedMass=Galactic_Structure_Enclosed_Mass(thisNode,radius,massType=massTypeSelected,componentType=componentTypeAll)
          
          ! Store the resulting radius and mass.
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=radius
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=enclosedMass
       end do
    end if
    return
  end subroutine Galacticus_Output_Tree_Density_Contrast

  function Mean_Density_Contrast_Root(radius,parameterPointer) bind(c)
    !% Root function used in finding the radius that encloses a given densit contrast.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Cosmology_Functions
    use Numerical_Constants_Math
    implicit none
    real(c_double)        :: Mean_Density_Contrast_Root
    real(c_double), value :: radius
    type(c_ptr),    value :: parameterPointer
    real(c_double)        :: enclosedMass
    
    ! Solve for the radius enclosing the specified density contrast.
    enclosedMass              =Galactic_Structure_Enclosed_Mass(activeNode,radius,massType=massTypeSelected,componentType&
         &=componentTypeAll)
    Mean_Density_Contrast_Root=3.0d0*enclosedMass/4.0d0/Pi/radius**3/(referenceDensity &
         &/Expansion_Factor(Tree_Node_Time(activeNode))**3)-meanDensityContrastTarget
    return
  end function Mean_Density_Contrast_Root
  
end module Galacticus_Output_Trees_Density_Contrasts
