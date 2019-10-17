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

!% Contains a module which handles outputting of data required by the halo model of galaxy clustering to the \glc\ output file.

module Galacticus_Output_Halo_Models
  !% Handles outputting of data required by the halo model of galaxy clustering to the \glc\ output file.
  implicit none
  private
  public :: Galacticus_Extra_Output_Halo_Fourier_Profile

  ! Flag indicating whether or not halo model information is to be output.
  logical                                     :: outputHaloModelData

  ! Flag indicating whether or not this module has been initialized.
  logical                                     :: outputHaloModelDataInitialized    =.false.

  ! Parameters controlling output data.
  integer                                     :: haloModelWavenumberPointsPerDecade
  double precision                            :: haloModelWavenumberMaximum                , haloModelWavenumberMinimum

  ! Wavenumbers used in power spectra.
  integer                                     :: wavenumberCount
  double precision, allocatable, dimension(:) :: wavenumber

contains

  subroutine Galacticus_Output_Halo_Model_Initialize
    !% Initializes the module by determining whether or not halo model data should be output.
    use :: Galacticus_HDF5                 , only : galacticusOutputFile
    use :: IO_HDF5                         , only : hdf5Access          , hdf5Object
    use :: Input_Parameters                , only : globalParameters    , inputParameter
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Ranges                , only : Make_Range          , rangeTypeLogarithmic
    implicit none
    type(hdf5Object) :: haloModelDataset, haloModelGroup

    if (.not.outputHaloModelDataInitialized) then
       !$omp critical(Galacticus_Output_Halo_Model_Initialize)
       if (.not.outputHaloModelDataInitialized) then
          !# <inputParameter>
          !#   <name>outputHaloModelData</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.false.</defaultValue>
          !#   <description>Specifies whether or not halo model data (bias, power spectra, etc.) should be included in the output.</description>
          !#   <group>output</group>
          !#   <source>globalParameters</source>
          !#   <type>boolean</type>
          !# </inputParameter>

          ! Read additional parameters if required.
          if (outputHaloModelData) then
             !# <inputParameter>
             !#   <name>haloModelWavenumberPointsPerDecade</name>
             !#   <cardinality>1</cardinality>
             !#   <defaultValue>10</defaultValue>
             !#   <description>The number of points per decade in wavenumber at which to tabulate power spectra for the halo model.</description>
             !#   <source>globalParameters</source>
             !#   <type>integer</type>
             !# </inputParameter>
             !# <inputParameter>
             !#   <name>haloModelWavenumberMinimum</name>
             !#   <cardinality>1</cardinality>
             !#   <defaultValue>1.0d-3</defaultValue>
             !#   <description>The minimum wavenumber (in Mpc${^-1}$) at which to tabulate power spectra for the halo model.</description>
             !#   <source>globalParameters</source>
             !#   <type>real</type>
             !# </inputParameter>
             !# <inputParameter>
             !#   <name>haloModelWavenumberMaximum</name>
             !#   <cardinality>1</cardinality>
             !#   <defaultValue>1.0d4</defaultValue>
             !#   <description>The maximum wavenumber (in Mpc${^-1}$) at which to tabulate power spectra for the halo model.</description>
             !#   <source>globalParameters</source>
             !#   <type>real</type>
             !# </inputParameter>
             ! Determine how many wavenumbers to tabulate at.
             wavenumberCount=int(log10(haloModelWavenumberMaximum/haloModelWavenumberMinimum)*dble(haloModelWavenumberPointsPerDecade))+1
             ! Allocate arrays for power spectrum.
             allocate(wavenumber(wavenumberCount))
             ! Build a grid of wavenumbers.
             wavenumber=Make_Range(haloModelWavenumberMinimum,haloModelWavenumberMaximum,wavenumberCount,rangeType=rangeTypeLogarithmic)
             call hdf5Access%set()
             haloModelGroup=galacticusOutputFile%openGroup("haloModel","Halo model data.")
             call haloModelGroup  %writeDataset  (wavenumber      ,'wavenumber','Wavenumber at which power spectrum is tabulated [Mpc⁻¹].',datasetReturned=haloModelDataset)
             call haloModelDataset%writeAttribute(1.0d0/megaParsec,'unitsInSI'                                                                                             )
             call haloModelDataset%close()
             call haloModelGroup  %close()
             call hdf5Access      %unset()
          end if
          ! Flag that module is now initialized.
          outputHaloModelDataInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Halo_Model_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Halo_Model_Initialize

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Galacticus_Extra_Output_Halo_Fourier_Profile</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Galacticus_Extra_Output_Halo_Fourier_Profile(node,iOutput,treeIndex,nodePassesFilter)
    !% Store Fourier-space halo profiles to the output file.
    use            :: Cosmology_Functions     , only : cosmologyFunctions  , cosmologyFunctionsClass
    use            :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMO, darkMatterProfileDMOClass
    use            :: Galacticus_HDF5         , only : galacticusOutputFile
    use            :: Galacticus_Nodes        , only : mergerTree          , nodeComponentBasic       , treeNode
    use            :: IO_HDF5                 , only : hdf5Object          , hdf5Access
    use, intrinsic :: ISO_C_Binding
    use            :: ISO_Varying_String
    use            :: Kind_Numbers            , only : kind_int8
    use            :: Memory_Management       , only : allocateArray       , deallocateArray
    use            :: String_Handling         , only : operator(//)
    implicit none
    type            (treeNode                 ), intent(inout), pointer      :: node
    integer         (c_size_t                 ), intent(in   )               :: iOutput
    integer         (kind=kind_int8           ), intent(in   )               :: treeIndex
    logical                                    , intent(in   )               :: nodePassesFilter
    type            (treeNode                 )               , pointer      :: hostNode
    class           (nodeComponentBasic       )               , pointer      :: basic
    class           (cosmologyFunctionsClass  )               , pointer      :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass)               , pointer      :: darkMatterProfileDMO_
    double precision                           , allocatable  , dimension(:) :: fourierProfile
    logical                                                                  :: nodeExistsInOutput
    integer                                                                  :: iWavenumber
    double precision                                                         :: expansionFactor
    type            (varying_string           )                              :: groupName            , dataSetName
    type            (hdf5Object               )                              :: outputGroup          , profilesGroup, &
         &                                                                      treeGroup

    call Galacticus_Output_Halo_Model_Initialize()
    ! For any node that passes the filter, we want to ensure that the host halo profile is output.
    ! Return immediately if halo model is not to be output or this node is filtered out.
    if (.not.(outputHaloModelData.and.nodePassesFilter)) return
    ! Get required objects.
    darkMatterProfileDMO_ => darkMatterProfileDMO()
     ! Find the host halo.
    hostNode => node
    do while (hostNode%isSatellite())
       hostNode => hostNode%parent
    end do
    ! Create a group for the profile datasets.
    call hdf5Access%set()
    profilesGroup=galacticusOutputFile%openGroup("haloModel","Halo model data.")
    groupName="Output"
    groupName=groupName//iOutput
    outputGroup=profilesGroup%openGroup(char(groupName),"Fourier space density profiles of halos for all trees at each output.")
    groupName="mergerTree"
    groupName=groupName//treeIndex
    treeGroup=outputGroup%openGroup(char(groupName),"Fourier space density profiles of halos for each tree.")
    ! Check if this halo has already been output.
    dataSetName="fourierProfile"
    dataSetName=groupName//hostNode%index()
    nodeExistsInOutput=treeGroup%hasDataset(char(dataSetName))
    if (.not.nodeExistsInOutput) then
       ! Allocate array to store profile.
       call allocateArray(fourierProfile,[wavenumberCount])
       ! Get the basic component.
       basic => hostNode%basic()
       ! Get the default cosmology functions object.
       cosmologyFunctions_ => cosmologyFunctions()
       ! Get the expansion factor.
       expansionFactor=cosmologyFunctions_%expansionFactor(basic%time())
       ! Construct profile. (Our wavenumbers are comoving, so we must convert them to physics
       ! coordinates before passing them to the dark matter profile k-space routine.)
       do iWavenumber=1,waveNumberCount
          fourierProfile(iWavenumber)=darkMatterProfileDMO_%kSpace(hostNode,waveNumber(iWavenumber)/expansionFactor)
       end do
       ! Write dataset to the group.
       call treeGroup%writeDataset(fourierProfile,char(dataSetName),"The Fourier-space density profile.")
       ! Deallocate profile array.
       call deallocateArray(fourierProfile)
    end if
    ! Close the profile group.
    call treeGroup    %close()
    call outputGroup  %close()
    call profilesGroup%close()
    call hdf5Access%unset()
    return
  end subroutine Galacticus_Extra_Output_Halo_Fourier_Profile

end module Galacticus_Output_Halo_Models
