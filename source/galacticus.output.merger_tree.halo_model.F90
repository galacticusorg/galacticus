!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which handles outputting of data required by the halo model of galaxy clustering to the \glc\ output file.

module Galacticus_Output_Halo_Models
  !% Handles outputting of data required by the halo model of galaxy clustering to the \glc\ output file.
  implicit none
  private
  public :: Galacticus_Output_Halo_Model, Galacticus_Output_Halo_Model_Property_Count, Galacticus_Output_Halo_Model_Names,&
       & Galacticus_Linear_Power_Spectrum_Output, Galacticus_Growth_Factor_Output, Galacticus_Extra_Output_Halo_Fourier_Profile

  ! Number of halo model properties.
  integer, parameter :: haloModelIntegerPropertyCount=1
  integer, parameter :: haloModelDoublePropertyCount =1

  ! Flag indicating whether or not halo model information is to be output.
  logical            :: outputHaloModelData

  ! Flag indicating whether or not this module has been initialized.
  logical            :: outputHaloModelDataInitialized=.false.

  ! Parameters controlling output data.
  integer            :: haloModelWavenumberPointsPerDecade
  double precision   :: haloModelWavenumberMinimum,haloModelWavenumberMaximum

  ! Wavenumbers used in power spectra.
  integer                                     :: wavenumberCount
  double precision, allocatable, dimension(:) :: wavenumber

contains

  subroutine Galacticus_Output_Halo_Model_Initialize
    !% Initializes the module by determining whether or not halo model data should be output.
    use Input_Parameters
    implicit none

    !$omp critical(Galacticus_Output_Halo_Model_Initialize)
    if (.not.outputHaloModelDataInitialized) then
       !@ <inputParameter>
       !@   <name>outputHaloModelData</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not halo model data (bias, power spectra, etc.) should be included in the output.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('outputHaloModelData',outputHaloModelData,defaultValue=.false.)

       ! Read additional parameters if required.
       if (outputHaloModelData) then
          !@ <inputParameter>
          !@   <name>haloModelWavenumberPointsPerDecade</name>
          !@   <defaultValue>10</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The number of points per decade in wavenumber at which to tabulate power spectra for the halo model.
          !@   </description>
          !@   <type>integer</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('haloModelWavenumberPointsPerDecade',haloModelWavenumberPointsPerDecade,defaultValue=10)
          !@ <inputParameter>
          !@   <name>haloModelWavenumberMinimum</name>
          !@   <defaultValue>$10^{-3}$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The minimum wavenumber (in Mpc${^-1}$) at which to tabulate power spectra for the halo model.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('haloModelWavenumberMinimum',haloModelWavenumberMinimum,defaultValue=1.0d-3)
          !@ <inputParameter>
          !@   <name>haloModelWavenumberMaximum</name>
          !@   <defaultValue>$10^4$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The maximum wavenumber (in Mpc${^-1}$) at which to tabulate power spectra for the halo model.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('haloModelWavenumberMaximum',haloModelWavenumberMaximum,defaultValue=1.0d4)
       end if

       ! Flag that module is now initialized.
       outputHaloModelDataInitialized=.true.
    end if
    !$omp end critical(Galacticus_Output_Halo_Model_Initialize)
    return
  end subroutine Galacticus_Output_Halo_Model_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Halo_Model_Names</unitName>
  !#  <sortName>Galacticus_Output_Halo_Model</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Halo_Model_Names(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of halo model properties to be written to the \glc\ output file.
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI

    ! Initialize the module.
    call Galacticus_Output_Halo_Model_Initialize

    ! Return property names if we are outputting halo model data.
    if (outputHaloModelData) then
       doubleProperty =doubleProperty +1
       !@ <outputProperty>
       !@   <name>nodeBias</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>The linear bias for this node.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames    (doubleProperty)='nodeBias'
       doublePropertyComments (doubleProperty)='The linear bias for this node.'
       doublePropertyUnitsSI  (doubleProperty)=0.0d0
       integerProperty=integerProperty+1
       !@ <outputProperty>
       !@   <name>isolatedHostIndex</name>
       !@   <datatype>integer</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>The index of the isolated node which hosts this node.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       integerPropertyNames   (integerProperty)='isolatedHostIndex'
       integerPropertyComments(integerProperty)='The index of the isolated node which hosts this node.'
       integerPropertyUnitsSI (integerProperty)=0.0d0
    end if
    return
  end subroutine Galacticus_Output_Halo_Model_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Halo_Model_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Halo_Model</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Halo_Model_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of halo model properties to be written to the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount
    
    ! Initialize the module.
    call Galacticus_Output_Halo_Model_Initialize

    ! Increment property count if we are outputting halo model data.
    if (outputHaloModelData) then
       integerPropertyCount=integerPropertyCount+haloModelIntegerPropertyCount
       doublePropertyCount =doublePropertyCount+haloModelDoublePropertyCount
    end if
    return
  end subroutine Galacticus_Output_Halo_Model_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Halo_Model</unitName>
  !#  <sortName>Galacticus_Output_Halo_Model</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Halo_Model(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store halo model properties in the \glc\ output file buffers.
    use Tree_Nodes
    use Dark_Matter_Halo_Biases
    use Kind_Numbers
    implicit none
    double precision,        intent(in)             :: time
    type(treeNode),          intent(inout), pointer :: thisNode
    integer,                 intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)          :: integerBuffer(:,:)
    double precision,        intent(inout)          :: doubleBuffer(:,:)
    type(treeNode),          pointer                :: isolatedNode
    
    ! Initialize the module.
    call Galacticus_Output_Halo_Model_Initialize

    ! Store property data if we are outputting halo model data.
    if (outputHaloModelData) then
       isolatedNode => thisNode
       do while (isolatedNode%isSatellite())
          isolatedNode => isolatedNode%parentNode
       end do
       doubleProperty =doubleProperty +1
       doubleBuffer(doubleBufferCount  ,doubleProperty )=Dark_Matter_Halo_Bias(isolatedNode)
       integerProperty=integerProperty+1
       integerBuffer(integerBufferCount,integerProperty)=isolatedNode%index()
    end if
    return
  end subroutine Galacticus_Output_Halo_Model

  !# <outputFileOpenTask>
  !#  <unitName>Galacticus_Linear_Power_Spectrum_Output</unitName>
  !# </outputFileOpenTask>
  subroutine Galacticus_Linear_Power_Spectrum_Output
    !% Output the linear theory power spectrum to the main output file.
    use Galacticus_HDF5
    use IO_HDF5
    use Numerical_Ranges
    use Memory_Management
    use CDM_Power_Spectrum
    use Numerical_Constants_Astronomical
    implicit none
    double precision, allocatable, dimension(:) :: powerSpectrum
    integer                                     :: iWavenumber
    type(hdf5Object)                            :: haloModelGroup,haloModelDataset

    ! Initialize the module.
    call Galacticus_Output_Halo_Model_Initialize

    ! Store power specturm if we are outputting halo model data.
    if (outputHaloModelData) then
       !@ <outputType>
       !@   <name>haloModel</name>
       !@   <description>A collection of data (including biases and halo profiles) that can be used in halo model calculations of galaxy clustering.</description>
       !@ </outputType>

       ! Create a group for halo model data..
       !$omp critical (HDF5_Access)
       haloModelGroup=IO_HDF5_Open_Group(galacticusOutputFile,'haloModel','Halo model data.')
       !$omp end critical (HDF5_Access)

       ! Determine how many wavenumbers to tabulate at.
       wavenumberCount=int(dlog10(haloModelWavenumberMaximum/haloModelWavenumberMinimum)*dble(haloModelWavenumberPointsPerDecade))+1

       ! Allocate arrays for power spectrum.
       call Alloc_Array(wavenumber   ,[wavenumberCount])
       call Alloc_Array(powerSpectrum,[wavenumberCount])

       ! Build a grid of wavenumbers.
       wavenumber=Make_Range(haloModelWavenumberMinimum,haloModelWavenumberMaximum,wavenumberCount,rangeType=rangeTypeLogarithmic)

       ! Compute power spectrum at each wavenumber.
       do iWavenumber=1,wavenumberCount
          powerSpectrum(iWavenumber)=Power_Spectrum_CDM(wavenumber(iWavenumber))
       end do

       ! Store the power spectrum
       !$omp critical (HDF5_Access)
       !@ <outputProperty>
       !@   <name>wavenumber</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Wavenumbers at which dark matter halo Fourier profiles are stored.</description>
       !@   <label>???</label>
       !@   <outputType>haloModel</outputType>
       !@ </outputProperty>
       call haloModelGroup%writeDataset(wavenumber   ,'wavenumber','Wavenumber at which power spectrum is tabulated [Mpc⁻¹].',datasetReturned=haloModelDataset)
       call haloModelDataset%writeAttribute(1.0d0/megaParsec,'unitsInSI')
       call haloModelDataset%close()
       !@ <outputProperty>
       !@   <name>powerSpectrum</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Linear theory power spectrum.</description>
       !@   <label>???</label>
       !@   <outputType>haloModel</outputType>
       !@ </outputProperty>
       call haloModelGroup%writeDataset(powerSpectrum,'powerSpectrum','Linear theory power spectrum [Mpc³].',datasetReturned=haloModelDataset)
       call haloModelDataset%writeAttribute(megaParsec**3   ,'unitsInSI')
       call haloModelDataset%close()
       !$omp end critical (HDF5_Access)
       
       ! Deallocate arrays.
       call Dealloc_Array(powerSpectrum)

       ! Close the halo model group.
       !$omp critical (HDF5_Access)
       call haloModelGroup%close()
       !$omp end critical (HDF5_Access)

    end if
    return
  end subroutine Galacticus_Linear_Power_Spectrum_Output

  !# <outputGroupOutputTask>
  !#  <unitName>Galacticus_Growth_Factor_Output</unitName>
  !# </outputGroupOutputTask>
  subroutine Galacticus_Growth_Factor_Output(outputGroup,time)
    !% Output the linear theory power spectrum to the main output file.
    use Linear_Growth
    use IO_HDF5
    type(hdf5Object), intent(inout) :: outputGroup
    double precision, intent(in)    :: time

    ! Initialize the module.
    call Galacticus_Output_Halo_Model_Initialize

    ! Store growth factor if we are outputting halo model data.
    if (outputHaloModelData) then
       call outputGroup%writeAttribute(Linear_Growth_Factor                       (time),'linearGrowthFactor'             )
       call outputGroup%writeAttribute(Linear_Growth_Factor_Logarithmic_Derivative(time),'linearGrowthFactorLogDerivative')
    end if
    
    return
  end subroutine Galacticus_Growth_Factor_Output

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Galacticus_Extra_Output_Halo_Fourier_Profile</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Galacticus_Extra_Output_Halo_Fourier_Profile(thisNode,iOutput,treeIndex,nodePassesFilter)
    !% Store Fourier-space halo profiles to the output file.
    use Tree_Nodes
    use IO_HDF5
    use Galacticus_HDF5
    use ISO_Varying_String
    use String_Handling
    use Memory_Management
    use Dark_Matter_Profiles
    use Cosmology_Functions
    use Kind_Numbers
    implicit none
    type(treeNode),          intent(inout), pointer      :: thisNode
    integer,                 intent(in)                  :: iOutput
    integer(kind=kind_int8), intent(in)                  :: treeIndex
    logical,                 intent(in)                  :: nodePassesFilter
    double precision,        allocatable,   dimension(:) :: fourierProfile
    integer                                              :: iWavenumber
    double precision                                     :: expansionFactor
    type(varying_string)                                 :: groupName
    type(hdf5Object)                                     :: profilesGroup,treeGroup,outputGroup
 
    ! If halo model output was requested, output the Fourier-space halo profiles.
    if (nodePassesFilter.and.outputHaloModelData.and..not.thisNode%isSatellite()) then
       ! Create a group for the profile datasets.
       !$omp critical (HDF5_Access)
       profilesGroup=IO_HDF5_Open_Group(galacticusOutputFile,"haloModel","Halo model data.")
       groupName="Output"
       groupName=groupName//iOutput
       outputGroup=IO_HDF5_Open_Group(profilesGroup,char(groupName),"Fourier space density profiles of halos for all trees at each output.")
       groupName="mergerTree"
       groupName=groupName//treeIndex
       treeGroup=IO_HDF5_Open_Group(outputGroup,char(groupName),"Fourier space density profiles of halos for each tree.")
       !$omp end critical (HDF5_Access)
       ! Allocate array to store profile.
       call Alloc_Array(fourierProfile,[wavenumberCount])
       ! Get the expansion factor.
       expansionFactor=Expansion_Factor(Tree_Node_Time(thisNode))
       ! Construct profile. (Our wavenumbers are comoving, so we must convert them to physics
       ! coordinates before passing them to the dark matter profile k-space routine.)
       do iWavenumber=1,waveNumberCount
          fourierProfile(iWavenumber)=Dark_Matter_Profile_kSpace(thisNode,waveNumber(iWavenumber)/expansionFactor)
       end do
       ! Write dataset to the group.
       groupName="fourierProfile"
       groupname=groupName//thisNode%index()
       !$omp critical (HDF5_Access)
       call treeGroup%writeDataset(fourierProfile,char(groupName),"The Fourier-space density profile.")
       !$omp end critical (HDF5_Access)
       ! Deallocate profile array.
       call Dealloc_Array(fourierProfile)
       ! Close the profile group.
       !$omp critical (HDF5_Access)
       call treeGroup    %close()
       call outputGroup  %close()
       call profilesGroup%close()
       !$omp end critical (HDF5_Access)
    end if
    return
  end subroutine Galacticus_Extra_Output_Halo_Fourier_Profile
  
end module Galacticus_Output_Halo_Models
