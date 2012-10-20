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

!% Contains a module which implements an intergalatic background radiation component read from a file.

module Radiation_IGB_File
  !% Implements an intergalatic background radiation component read from a file.
  use FGSL
  implicit none
  private
  public :: Radiation_IGB_File_Initialize,Radiation_IGB_File_Format_Version

  ! Flag indicating whether the module has been initialized yet.
  logical :: moduleInitialized=.false.

  ! Arrays holding the radiation data.
  integer                                       :: spectraTimesCount,spectraWavelengthsCount
  double precision, allocatable, dimension(:  ) :: spectraTimes,spectraWavelengths
  double precision, allocatable, dimension(:,:) :: spectra

  ! Interpolation structures.
  logical                                       :: interpolationReset=.true., interpolationResetTimes=.true.
  type(fgsl_interp_accel)                       :: interpolationAccelerator , interpolationAcceleratorTimes
  type(fgsl_interp      )                       :: interpolationObject

  ! Current file format version for intergalactic background radiation files.
  integer         , parameter                   :: fileFormatVersionCurrent=1

contains

  integer function Radiation_IGB_File_Format_Version()
    !% Return the current file format version of intergalactic background radiation files.
    implicit none

    Radiation_IGB_File_Format_Version=fileFormatVersionCurrent
    return
  end function Radiation_IGB_File_Format_Version

  !# <radiationIntergalacticBackgroundMethod>
  !#  <unitName>Radiation_IGB_File_Initialize</unitName>
  !# </radiationIntergalacticBackgroundMethod>
  subroutine Radiation_IGB_File_Initialize(radiationIntergalacticBackgroundMethod,Radiation_Set_Intergalactic_Background_Do,Radiation_Flux_Intergalactic_Background_Do)
    !% Initialize the intergalactic background radiation component from file module by reading in the data.
    use Input_Parameters
    use ISO_Varying_String
    use FoX_dom
    use Galacticus_Error
    use Memory_Management
    use Cosmology_Functions
    use Array_Utilities
    use Galacticus_Input_Paths
    implicit none
    type(varying_string),          intent(in)    :: radiationIntergalacticBackgroundMethod
    procedure(),          pointer, intent(inout) :: Radiation_Set_Intergalactic_Background_Do,Radiation_Flux_Intergalactic_Background_Do
    type(Node),           pointer                :: doc,thisSpectrum,thisWavelength,thisDatum
    type(NodeList),       pointer                :: spectraList,datumList,wavelengthList
    integer                                      :: ioErr,iSpectrum,jSpectrum,iWavelength,fileFormatVersion
    logical                                      :: timesIncreasing
    type(varying_string)                         :: radiationIGBFileName

    if (radiationIntergalacticBackgroundMethod == 'file') then
       Radiation_Set_Intergalactic_Background_Do  => Radiation_IGB_File_Set
       Radiation_Flux_Intergalactic_Background_Do => Radiation_IGB_File_Flux
       ! Get the name of the file from which to read data.
       !@ <inputParameter>
       !@   <name>radiationIGBFileName</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the file containing a tabulation of the radiation field.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('radiationIGBFileName',radiationIGBFileName,defaultValue=char(Galacticus_Input_Path())//"data/radiation/Cosmic_Background_Radiation_Haardt_Madau_2005_Quasars_Galaxies.xml")

       !$omp critical (FoX_DOM_Access)
       ! Parse the XML file.
       doc => parseFile(char(radiationIGBFileName),iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Radiation_Initialize_File','Unable to find or parse data file')
       ! Check the file format version of the file.
       datumList => getElementsByTagname(doc,"fileFormat")
       thisDatum => item(datumList,0)
       call extractDataContent(thisDatum,fileFormatVersion)
       if (fileFormatVersion /= fileFormatVersionCurrent) call Galacticus_Error_Report('Radiation_IGB_File_Initialize','file format version is out of date')
       ! Get a list of all spectra.
       spectraList => getElementsByTagname(doc,"spectra")
       ! Get the times from the file.
       spectraTimesCount=getLength(spectraList)
       call Alloc_Array(spectraTimes,[spectraTimesCount])
       do iSpectrum=1,spectraTimesCount
          thisSpectrum => item(spectraList,iSpectrum-1)
          datumList    => getElementsByTagname(thisSpectrum,"redshift")
          thisDatum    => item(datumList,0)
          call extractDataContent(thisDatum,spectraTimes(iSpectrum))
          ! Convert redshift to a time.
          spectraTimes(iSpectrum)=Cosmology_Age(Expansion_Factor_from_Redshift(spectraTimes(iSpectrum)))
       end do
       ! Check if the times are monotonically ordered.
       if (.not.Array_Is_Monotonic(spectraTimes)) call Galacticus_Error_Report('Radiation_Initialize_File','spectra must be monotonically ordered in time')
       timesIncreasing=Array_Is_Monotonic(spectraTimes,direction=directionIncreasing)
       ! Reverse times if necessary.
       if (.not.timesIncreasing) spectraTimes=Array_Reverse(spectraTimes)
       ! Get the wavelengths.
       wavelengthList => getElementsByTagname(doc,"wavelengths")
       thisWavelength => item(wavelengthList,0)
       datumList => getElementsByTagname(thisWavelength,"datum")
       spectraWavelengthsCount=getLength(datumList)
       ! Allocate arrays for wavelengths and spectra.
       call Alloc_Array(spectraWavelengths,[spectraWavelengthsCount                  ])
       call Alloc_Array(spectra           ,[spectraWavelengthsCount,spectraTimesCount])
       ! Extract the wavelengths.
       do iWavelength=1,spectraWavelengthsCount
          thisDatum => item(datumList,iWavelength-1)
         call extractDataContent(thisDatum,spectraWavelengths(iWavelength))
       end do
       ! Read spectra into arrays.
       do iSpectrum=1,spectraTimesCount
          ! Determine where to store this spectrum, depending on whether the times were stored in increasing or decreasing order.
          if (timesIncreasing) then
             jSpectrum=iSpectrum
          else
             jSpectrum=spectraTimesCount+1-iSpectrum
          end if
          ! Get the data.
          thisSpectrum => item(spectraList,iSpectrum-1)
          datumList    => getElementsByTagname(thisSpectrum,"datum")
          ! Check that we have the correct number of data.
          if (getLength(datumList) /= spectraWavelengthsCount) call Galacticus_Error_Report('Radiation_Initialize_File','all spectra must contain the same number of wavelengths')
          ! Extract the data.
          do iWavelength=1,spectraWavelengthsCount
             thisDatum => item(datumList,iWavelength-1)
             call extractDataContent(thisDatum,spectra(iWavelength,jSpectrum))
          end do
       end do
       ! Destroy the document.
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)
       
       ! Flag that the module is now initialized.
       moduleInitialized=.true.
    end if

    return
  end subroutine Radiation_IGB_File_Initialize

  subroutine Radiation_IGB_File_Set(thisNode,radiationProperties)
    !% Property setting routine for the radiation component from file method.
    use Tree_Nodes
    use Memory_Management
    implicit none
    type(treeNode),   intent(inout), pointer                   :: thisNode
    double precision, intent(inout), allocatable, dimension(:) :: radiationProperties

    ! Ensure that the properties array is allocated.
    if (.not.allocated(radiationProperties)) call Alloc_Array(radiationProperties,[1])

    ! Store the time for the radiation field.
    radiationProperties(1)=Tree_Node_Time(thisNode)

    return
  end subroutine Radiation_IGB_File_Set

  subroutine Radiation_IGB_File_Flux(radiationProperties,wavelength,radiationFlux)
    !% Flux method for the radiation component from file method.
    use Numerical_Interpolation
    implicit none
    double precision, intent(in)                   :: wavelength
    double precision, intent(in),   dimension(:)   :: radiationProperties
    double precision, intent(inout)                :: radiationFlux
    double precision, save,         dimension(0:1) :: hSpectrum
    !$omp threadprivate(hSpectrum)
    double precision,               dimension(0:1) :: hWavelength
    integer,          save                         :: iSpectrum
    !$omp threadprivate(iSpectrum)
    double precision, save                         :: previousTime=-1.0d0
    integer                                        :: jSpectrum,iWavelength,jWavelength

    ! Return if out of range.
    if     (    radiationProperties(1) < spectraTimes      (1                      ) &
         & .or. radiationProperties(1) > spectraTimes      (spectraTimesCount      ) &
         & .or. wavelength             < spectraWavelengths(1                      ) &
         & .or. wavelength             > spectraWavelengths(spectraWavelengthsCount) &
         & ) return
   
    !$omp critical (Radiation_Interpolation)
    ! Find interpolate in the array of times (only necessary if the time differs from that on the previous call).
    if (radiationProperties(1) /= previousTime) then
       iSpectrum=Interpolate_Locate(spectraTimesCount,spectraTimes,interpolationAcceleratorTimes,radiationProperties(1),reset=interpolationResetTimes)
       hSpectrum=Interpolate_Linear_Generate_Factors(spectraTimesCount,spectraTimes,iSpectrum,radiationProperties(1))
       previousTime=radiationProperties(1)
    end if

    ! Find interpolation in the array of wavelengths.
    iWavelength=Interpolate_Locate(spectraWavelengthsCount,spectraWavelengths,interpolationAccelerator,wavelength,reset=interpolationReset)
    hWavelength=Interpolate_Linear_Generate_Factors(spectraWavelengthsCount,spectraWavelengths,iWavelength,wavelength)
    do jSpectrum=0,1
       do jWavelength=0,1
          radiationFlux=radiationFlux+hSpectrum(jSpectrum)*hWavelength(jWavelength)*spectra(jWavelength+iWavelength,jSpectrum&
               &+iSpectrum)
       end do
    end do
    !$omp end critical (Radiation_Interpolation)

    return
  end subroutine Radiation_IGB_File_Flux

end module Radiation_IGB_File
