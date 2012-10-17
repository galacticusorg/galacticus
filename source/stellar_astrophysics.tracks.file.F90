!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculation of stellar tracks.

module Stellar_Astrophysics_Tracks_File
  !% Implements stellar tracks.
  use FGSL
  implicit none
  private
  public :: Stellar_Tracks_Initialize_File

  ! Variables to hold stellar track data.
  double precision, allocatable :: stellarTrackLogMetallicities(:),stellarTrackInitialMasses(:,:),stellarTrackLuminosities(:,:,:)&
       &,stellarTrackTemperatures(:,:,:),stellarTrackAges(:,:,:)
  integer,          allocatable :: stellarTrackInitialMassCount(:),stellarTrackAgesCount(:,:)
  integer                       :: stellarTrackMetallicityCount

  ! Interpolation objects.
  type(fgsl_interp_accel)       :: interpolationAcceleratorMetallicity
  logical                       :: interpolationResetMetallicity=.true.

  ! The current file format version.
  integer,          parameter   :: fileFormatVersionCurrent=1

contains

  !# <stellarTracksMethod>
  !#  <unitName>Stellar_Tracks_Initialize_File</unitName>
  !# </stellarTracksMethod>
  subroutine Stellar_Tracks_Initialize_File(stellarTracksMethod,Stellar_Luminosity_Get,Stellar_Effective_Temperature_Get)
    !% Initialize the stellar tracks module.
    use Galacticus_Error
    use Memory_Management
    use IO_HDF5
    use ISO_Varying_String
    use Galacticus_Input_Paths
    use Input_Parameters
    use String_Handling
    implicit none
    type     (varying_string  ),          intent(in   ) :: stellarTracksMethod
    procedure(double precision), pointer, intent(inout) :: Stellar_Luminosity_Get,Stellar_Effective_Temperature_Get
    type     (varying_string  )                         :: stellarTracksFile,groupName
    integer                                             :: initialMassCount,initialMassCountMaximum ,ageCountMaximum&
         &,metallicityCountMaximum,fileFormatVersion
    type     (hdf5Object      )                         :: stellarTracks,metallicityGroup,massGroup,ageDataset
    logical                                             :: foundMetallicityGroup,foundMassGroup

    ! Check if our method is selected.
    if (stellarTracksMethod == 'file') then
       ! Set up procedure pointers.
       Stellar_Luminosity_Get            => Stellar_Luminosity_File
       Stellar_Effective_Temperature_Get => Stellar_Effective_Temperature_File

       ! Get the name of the file from which to read stellar tracks.
       !@ <inputParameter>
       !@   <name>stellarTracksFile</name>
       !@   <defaultValue>data/stellarAstrophysics/Stellar\_Tracks\_Padova.hdf5</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the HDF5 file from which to read stellar tracks.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('stellarTracksFile',stellarTracksFile,defaultValue=char(Galacticus_Input_Path())//'data/stellarAstrophysics/Stellar_Tracks_Padova.hdf5')

       ! Open the HDF5 file.
       !$omp critical(HDF5_Access)
       call stellarTracks%openFile(char(stellarTracksFile),readOnly=.true.)

       ! Check that this file has the correct format.
       call stellarTracks%readAttribute('fileFormat',fileFormatVersion,allowPseudoScalar=.true.)
       if (fileFormatVersion /= fileFormatVersionCurrent) call Galacticus_Error_Report('Stellar_Tracks_Initialize_File','format of stellar tracks file is out of date')

       ! Count up number of metallicities present, the number of stellar masses tabulated and the number of ages tabulated.
       metallicityCountMaximum=0
       initialMassCountMaximum=0
       ageCountMaximum        =0
       ! Count metallicity groups.
       foundMetallicityGroup=.true.
       do while (foundMetallicityGroup)
          groupName="metallicity"
          groupName=groupName//(metallicityCountMaximum+1)
          foundMetallicityGroup=stellarTracks%hasGroup(char(groupName))
          if (foundMetallicityGroup) then
             metallicityCountMaximum=metallicityCountMaximum+1
             metallicityGroup=stellarTracks%openGroup(char(groupName))
             ! Find mass groups.
             foundMassGroup=.true.
             do while (foundMassGroup)
                groupName="mass"
                groupName=groupName//(initialMassCountMaximum+1)
                foundMassGroup=metallicityGroup%hasGroup(char(groupName))
                if (foundMassGroup) then
                   initialMassCountMaximum=initialMassCountMaximum+1
                   massGroup=metallicityGroup%openGroup(char(groupName))
                   ageDataset=massGroup%openDataset('age')
                   ageCountMaximum=max(ageCountMaximum,int(ageDataset%size(1)))
                   call ageDataset%close()
                   call massGroup%close()
                end if
             end do
             call metallicityGroup%close()
          end if
       end do

       ! Allocate storage space for data.
       call Alloc_Array(stellarTrackLogMetallicities,[                                        metallicityCountMaximum])
       call Alloc_Array(stellarTrackInitialMassCount,[                                        metallicityCountMaximum])
       call Alloc_Array(stellarTrackInitialMasses   ,[                initialMassCountMaximum,metallicityCountMaximum])
       call Alloc_Array(stellarTrackAgesCount       ,[                initialMassCountMaximum,metallicityCountMaximum])
       call Alloc_Array(stellarTrackAges            ,[ageCountMaximum,initialMassCountMaximum,metallicityCountMaximum])
       call Alloc_Array(stellarTrackLuminosities    ,[ageCountMaximum,initialMassCountMaximum,metallicityCountMaximum])
       call Alloc_Array(stellarTrackTemperatures    ,[ageCountMaximum,initialMassCountMaximum,metallicityCountMaximum])

       ! Read in all data.
       do stellarTrackMetallicityCount=1,metallicityCountMaximum
          ! Open the metallicity group.
          groupName="metallicity"
          groupName=groupName//stellarTrackMetallicityCount
          metallicityGroup=stellarTracks%openGroup(char(groupName))
          ! Get the metallicity.
          call metallicityGroup%readDatasetStatic('metallicity', stellarTrackLogMetallicities(stellarTrackMetallicityCount:stellarTrackMetallicityCount))
          ! Count how many masses are tabulated at this metallicity.
          initialMassCount=0
          foundMassGroup=.true.
          do while (foundMassGroup)
             groupName="mass"
             groupName=groupName//(initialMassCount+1)
             foundMassGroup=metallicityGroup%hasGroup(char(groupName))
             if (foundMassGroup) initialMassCount=initialMassCount+1                
          end do
          stellarTrackInitialMassCount(stellarTrackMetallicityCount)=initialMassCount
          ! Loop through all tabulated masses.
          do initialMassCount=1,stellarTrackInitialMassCount(stellarTrackMetallicityCount)
             ! Open the mass group.
             groupName="mass"
             groupName=groupName//initialMassCount
             massGroup=metallicityGroup%openGroup(char(groupName))
             ! Get initial mass.
             call massGroup%readDatasetStatic('mass',stellarTrackInitialMasses(initialMassCount:initialMassCount,stellarTrackMetallicityCount))
             ! Read tracks.
             ageDataset=massGroup%openDataset('age')
             stellarTrackAgesCount(initialMassCount,stellarTrackMetallicityCount)=int(ageDataset%size(1))
             call ageDataset%close()
             call massGroup%readDatasetStatic('age'                 ,stellarTrackAges        (1:stellarTrackAgesCount(initialMassCount,stellarTrackMetallicityCount),initialMassCount,stellarTrackMetallicityCount))
             call massGroup%readDatasetStatic('luminosity'          ,stellarTrackLuminosities(1:stellarTrackAgesCount(initialMassCount,stellarTrackMetallicityCount),initialMassCount,stellarTrackMetallicityCount))
             call massGroup%readDatasetStatic('effectiveTemperature',stellarTrackTemperatures(1:stellarTrackAgesCount(initialMassCount,stellarTrackMetallicityCount),initialMassCount,stellarTrackMetallicityCount))
             call massGroup%close()
          end do
          call metallicityGroup%close()
       end do
       ! Convert metallicities to logarithmic scale.
       stellarTrackLogMetallicities=dlog(stellarTrackLogMetallicities)
       stellarTrackMetallicityCount=metallicityCountMaximum

       ! Close the file.
       call stellarTracks%close()
       !$omp end critical(HDF5_Access)

    end if
    return
  end subroutine Stellar_Tracks_Initialize_File

  double precision function Stellar_Luminosity_File(initialMass,metallicity,age)
    !% Return the bolometric luminosity (in $L_\odot$) for a star of given {\tt initialMass}, {\tt metallicity} and {\tt age}.
    implicit none
    double precision, intent(in) :: initialMass,metallicity,age
    integer                      :: interpolationIndicesMetallicity(2),interpolationIndicesMass(2,2),interpolationIndicesAge(2,2,2)
    double precision             :: interpolationFactorsMetallicity(2),interpolationFactorsMass(2,2),interpolationFactorsAge(2,2,2)
    logical                      :: metallicityOutOfRange,massOutOfRange,ageOutOfRange

    ! Get the interpolating factors.
    call Stellar_Tracks_Interpolation_Get(initialMass,metallicity,age,interpolationIndicesMetallicity,interpolationIndicesMass &
         &,interpolationIndicesAge ,interpolationFactorsMetallicity,interpolationFactorsMass,interpolationFactorsAge&
         &,metallicityOutOfRange ,massOutOfRange,ageOutOfRange)

    ! Do the interpolation.
    if (massOutOfRange.or.ageOutOfRange) then
       Stellar_Luminosity_File=0.0d0
    else
       Stellar_Luminosity_File=Stellar_Tracks_Interpolation_Do(interpolationIndicesMetallicity,interpolationIndicesMass &
            &,interpolationIndicesAge ,interpolationFactorsMetallicity,interpolationFactorsMass,interpolationFactorsAge&
            &,stellarTrackLuminosities)
    end if

    return
  end function Stellar_Luminosity_File

  double precision function Stellar_Effective_Temperature_File(initialMass,metallicity,age)
    !% Return the effective temperature (in Kelvin) for a star of given {\tt initialMass}, {\tt metallicity} and {\tt age}.
    implicit none
    double precision, intent(in) :: initialMass,metallicity,age
    integer                      :: interpolationIndicesMetallicity(2),interpolationIndicesMass(2,2),interpolationIndicesAge(2,2,2)
    double precision             :: interpolationFactorsMetallicity(2),interpolationFactorsMass(2,2),interpolationFactorsAge(2,2,2)
    logical                      :: metallicityOutOfRange,massOutOfRange,ageOutOfRange

    ! Get the interpolating factors.
    call Stellar_Tracks_Interpolation_Get(initialMass,metallicity,age,interpolationIndicesMetallicity,interpolationIndicesMass &
         &,interpolationIndicesAge ,interpolationFactorsMetallicity,interpolationFactorsMass,interpolationFactorsAge&
         &,metallicityOutOfRange ,massOutOfRange,ageOutOfRange)

    ! Do the interpolation.
    if (massOutOfRange.or.ageOutOfRange) then
       Stellar_Effective_Temperature_File=0.0d0
    else
       Stellar_Effective_Temperature_File=Stellar_Tracks_Interpolation_Do(interpolationIndicesMetallicity,interpolationIndicesMass &
            &,interpolationIndicesAge ,interpolationFactorsMetallicity,interpolationFactorsMass,interpolationFactorsAge&
            &,stellarTrackTemperatures)
    end if

    return
  end function Stellar_Effective_Temperature_File

  double precision function Stellar_Tracks_Interpolation_Do(interpolationIndicesMetallicity,interpolationIndicesMass &
       &,interpolationIndicesAge,interpolationFactorsMetallicity,interpolationFactorsMass,interpolationFactorsAge&
       &,stellarTracks)
    !% Using precomputed factors, interpolate in metallicity, mass and age in the given {\tt stellarTracks}.
    implicit none
    integer,          intent(in) :: interpolationIndicesMetallicity(2),interpolationIndicesMass(2,2),interpolationIndicesAge(2,2,2)
    double precision, intent(in) :: interpolationFactorsMetallicity(2),interpolationFactorsMass(2,2),interpolationFactorsAge(2,2,2)
    double precision, intent(in) :: stellarTracks(:,:,:)
    integer                      :: iMetallicity,iMass,iAge,jMetallicity,jMass,jAge

    Stellar_Tracks_Interpolation_Do=0.0d0
    do iMetallicity=1,2
       jMetallicity=interpolationIndicesMetallicity(iMetallicity)
       do iMass=1,2
          jMass=interpolationIndicesMass(iMetallicity,iMass)
          do iAge=1,2
             jAge=interpolationIndicesAge(iMetallicity,iMass,iAge)
             Stellar_Tracks_Interpolation_Do=Stellar_Tracks_Interpolation_Do+stellarTracks(jAge,jMass,jMetallicity)&
                  &*interpolationFactorsMetallicity(iMetallicity)*interpolationFactorsMass(iMetallicity,iMass)&
                  &*interpolationFactorsAge(iMetallicity,iMass,iAge)
          end do
       end do
    end do
    return
  end function Stellar_Tracks_Interpolation_Do

  subroutine Stellar_Tracks_Interpolation_Get(initialMass,metallicity,age,interpolationIndicesMetallicity,interpolationIndicesMass&
       &,interpolationIndicesAge ,interpolationFactorsMetallicity,interpolationFactorsMass,interpolationFactorsAge&
       &,metallicityOutOfRange,massOutOfRange,ageOutOfRange)
    !% Get interpolating factors for stellar tracks.
    use Numerical_Interpolation
    implicit none
    double precision,        intent(in)  :: initialMass,metallicity,age
    integer,                 intent(out) :: interpolationIndicesMetallicity(2),interpolationIndicesMass(2,2),interpolationIndicesAge(2,2,2)
    double precision,        intent(out) :: interpolationFactorsMetallicity(2),interpolationFactorsMass(2,2),interpolationFactorsAge(2,2,2)
    logical,                 intent(out) :: metallicityOutOfRange,massOutOfRange,ageOutOfRange
    integer                              :: iMetallicity,jMetallicity,iMass,jMass
    double precision                     :: logMetallicity
    type(fgsl_interp_accel), save        :: interpolationAcceleratorMass,interpolationAcceleratorAge
    logical,                 save        :: interpolationResetMass=.true.,interpolationResetAge=.true.
    !$omp threadprivate(interpolationAcceleratorMass,interpolationAcceleratorAge,interpolationResetMass,interpolationResetAge)

    !$omp critical (Stellar_Tracks_Interpolate)

    ! Assume everything is in range initially.
    metallicityOutOfRange=.false.
    massOutOfRange       =.false.
    ageOutOfRange        =.false.

    ! Interpolate in metallicity.
    if (metallicity <= 0.0d0) then
       logMetallicity=stellarTrackLogMetallicities(1)
    else
       logMetallicity=dlog(metallicity)
    end if
    if (logMetallicity < stellarTrackLogMetallicities(1)) then
       interpolationIndicesMetallicity=[    1,    2]
       interpolationFactorsMetallicity=[1.0d0,0.0d0]
       metallicityOutOfRange=.true.
    else if (logMetallicity > stellarTrackLogMetallicities(stellarTrackMetallicityCount)) then
       interpolationIndicesMetallicity=[stellarTrackMetallicityCount-1,stellarTrackMetallicityCount]
       interpolationFactorsMetallicity=[0.0d0,1.0d0]
       metallicityOutOfRange=.true.
    else
       interpolationIndicesMetallicity(1)=Interpolate_Locate(stellarTrackMetallicityCount,stellarTrackLogMetallicities&
            &,interpolationAcceleratorMetallicity,logMetallicity,reset=interpolationResetMetallicity)
       interpolationIndicesMetallicity(2)=interpolationIndicesMetallicity(1)+1
       interpolationFactorsMetallicity=Interpolate_Linear_Generate_Factors(stellarTrackMetallicityCount&
            &,stellarTrackLogMetallicities,interpolationIndicesMetallicity(1),logMetallicity)
    end if

    ! Loop over metallicities.
    do iMetallicity=1,2
       jMetallicity=interpolationIndicesMetallicity(iMetallicity)

       ! Interpolate in mass at each metallicity.
       if (initialMass < stellarTrackInitialMasses(1,jMetallicity)) then
          interpolationIndicesMass(iMetallicity,:)=[    1,    2]
          interpolationFactorsMass(iMetallicity,:)=[1.0d0,0.0d0]
          massOutOfRange=.true.
       else if (initialMass > stellarTrackInitialMasses(stellarTrackInitialMassCount(jMetallicity),jMetallicity)) then
          interpolationIndicesMass(iMetallicity,:)=[stellarTrackInitialMassCount(jMetallicity)-1,stellarTrackInitialMassCount(jMetallicity)]
          interpolationFactorsMass(iMetallicity,:)=[0.0d0,1.0d0]         
          massOutOfRange=.true.
       else
          interpolationResetMass=.true.
          interpolationIndicesMass(iMetallicity,1)=Interpolate_Locate(stellarTrackInitialMassCount(jMetallicity)&
               &,stellarTrackInitialMasses(:,jMetallicity),interpolationAcceleratorMass,initialMass,reset=interpolationResetMass)
          call Interpolate_Done(interpolationAccelerator=interpolationAcceleratorMass,reset=interpolationResetMass)
          interpolationIndicesMass(iMetallicity,2)=interpolationIndicesMass(iMetallicity,1)+1
          interpolationFactorsMass(iMetallicity,:)=Interpolate_Linear_Generate_Factors(stellarTrackInitialMassCount(jMetallicity)&
               &,stellarTrackInitialMasses(:,jMetallicity),interpolationIndicesMass(iMetallicity,1),initialMass)
       end if

       ! Loop over masses.
       do iMass=1,2
          jMass=interpolationIndicesMass(iMetallicity,iMass)

          ! Interpolate in age at each mass.
          if (age < stellarTrackAges(1,jMass,jMetallicity)) then
             interpolationIndicesAge(iMetallicity,iMass,:)=[    1,    2]
             interpolationFactorsAge(iMetallicity,iMass,:)=[1.0d0,0.0d0]
             ageOutOfRange=.true.
          else if (age > stellarTrackAges(stellarTrackAgesCount(jMass,jMetallicity),jMass,jMetallicity)) then
             interpolationIndicesAge(iMetallicity,iMass,:)=[stellarTrackAgesCount(jMass,jMetallicity)-1,stellarTrackAgesCount(jMass,jMetallicity)]
             interpolationFactorsAge(iMetallicity,iMass,:)=[0.0d0,1.0d0]
             ageOutOfRange=.true.
          else
             interpolationResetAge=.true.
             interpolationIndicesAge(iMetallicity,iMass,1)=Interpolate_Locate(stellarTrackAgesCount(jMass,jMetallicity)&
                  &,stellarTrackAges(:,jMass,jMetallicity) ,interpolationAcceleratorAge,age,reset=interpolationResetAge)
             call Interpolate_Done(interpolationAccelerator=interpolationAcceleratorAge,reset=interpolationResetAge)
             interpolationIndicesAge(iMetallicity,iMass,2)=interpolationIndicesAge(iMetallicity,iMass,1)+1
             interpolationFactorsAge(iMetallicity,iMass,:)=Interpolate_Linear_Generate_Factors(stellarTrackAgesCount(jMass,jMetallicity)&
                  &,stellarTrackAges(:,jMass,jMetallicity),interpolationIndicesAge(iMetallicity,iMass,1),age)
          end if
       end do

    end do
    !$omp end critical (Stellar_Tracks_Interpolate)

    return
  end subroutine Stellar_Tracks_Interpolation_Get

end module Stellar_Astrophysics_Tracks_File
