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

!% Contains a module which implements calculations of power spectra and related properties for output.

module Power_Spectrum_Tasks
  !% Implements calculations of power spectra and related properties for output.
  use IO_HDF5
  implicit none
  private
  public :: Power_Spectrum_Compute, Power_Spectrum_Open_File, Power_Spectrum_Close_File, Power_Spectrum_Output
  
  ! HDF5 object for the output file.
  type            (hdf5Object), public                    :: powerSpectrumOutputFile                                  
  
  ! Arrays of power spectrum data.
  double precision            , allocatable, dimension(:) :: powerSpectrum_mass      , powerSpectrum_power        , & 
       &                                                     powerSpectrum_sigma     , powerSpectrum_sigmaGradient, & 
       &                                                     powerSpectrum_wavenumber                                 
  
contains
  
  subroutine Power_Spectrum_Open_File(outputFileName)
    !% Open the output file for power spectrum data.
    use ISO_Varying_String
    use HDF5
    implicit none
    type(varying_string), intent(in   ) :: outputFileName 
    
    ! Open the output file.
    call powerSpectrumOutputFile%openFile(char(outputFileName),overWrite=.true.,objectsOverwritable=.false.)
    
    ! Set default chunking and compression levels.
    call IO_HDF5_Set_Defaults(chunkSize=int(128,kind=hsize_t),compressionLevel=9)

    return
  end subroutine Power_Spectrum_Open_File

  subroutine Power_Spectrum_Close_File
    !% Close the output file for power spectrum data.
    implicit none
    
    call powerSpectrumOutputFile%close()
    return
  end subroutine Power_Spectrum_Close_File

  subroutine Power_Spectrum_Compute
    !% Computes power spectra and related properties for output.
    use Memory_Management
    use Numerical_Ranges
    use Input_Parameters
    use Power_Spectra
    use Numerical_Constants_Math
    use Cosmological_Parameters
    implicit none
    integer          :: iWavenumber                  , powerSpectraCount            , & 
         &              powerSpectraPointsPerDecade                                     
    double precision :: powerSpectraWavenumberMaximum, powerSpectraWavenumberMinimum    
    
    ! Find the wavenumber range and increment size.    !@ <inputParameter>    !@   <name>powerSpectraWavenumberMinimum</name>    !@   <defaultValue>$10^{-3}$ Mpc$^{-1}$</defaultValue>    !@   <attachedTo>module</attachedTo>    !@   <description>    !@     The minimum wavenumber at which to tabulate power spectra.    !@   </description>    !@   <type>real</type>    !@   <cardinality>1</cardinality>    !@ </inputParameter>
    call Get_Input_Parameter('powerSpectraWavenumberMinimum',powerSpectraWavenumberMinimum,defaultValue=1.0d-3)
    !@ <inputParameter>
    !@   <name>powerSpectraWavenumberMaximum</name>
    !@   <defaultValue>$10^{3}$ Mpc$^{-1}$</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The maximum wavenumber at which to tabulate power spectra.
    !@   </description>
    !@   <type>real</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    call Get_Input_Parameter('powerSpectraWavenumberMaximum',powerSpectraWavenumberMaximum,defaultValue=1.0d+3)
    !@ <inputParameter>
    !@   <name>powerSpectraPointsPerDecade</name>
    !@   <defaultValue>10</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The number of points per decade of wavenumber at which to tabulate power spectra.
    !@   </description>
    !@   <type>integer</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    call Get_Input_Parameter('powerSpectraPointsPerDecade',powerSpectraPointsPerDecade,defaultValue=10)

    ! Compute number of tabulation points.
    powerSpectraCount=int(log10(powerSpectraWavenumberMaximum/powerSpectraWavenumberMinimum)*dble(powerSpectraPointsPerDecade))+1

    ! Allocate arrays for power spectra.
    call Alloc_Array(powerSpectrum_Wavenumber   ,[powerSpectraCount])
    call Alloc_Array(powerSpectrum_Power        ,[powerSpectraCount])
    call Alloc_Array(powerSpectrum_Mass         ,[powerSpectraCount])
    call Alloc_Array(powerSpectrum_sigma        ,[powerSpectraCount])
    call Alloc_Array(powerSpectrum_sigmaGradient,[powerSpectraCount])
       
    ! Build a range of wavenumbers.
    powerSpectrum_Wavenumber(:)=Make_Range(powerSpectraWavenumberMinimum,powerSpectraWavenumberMaximum,powerSpectraCount,rangeTypeLogarithmic)
       
    ! Loop over all halo wavenumberes.
    do iWavenumber=1,powerSpectraCount
       ! Compute power spectrum.
       powerSpectrum_Power        (iWavenumber)=Power_Spectrum(powerSpectrum_Wavenumber(iWavenumber))
       ! Compute corresponding mass scale.
       powerSpectrum_Mass         (iWavenumber)=4.0d0*Pi*Omega_Matter()*Critical_Density()/3.0d0/powerSpectrum_Wavenumber(iWavenumber)**3
       ! Compute fluctuation on this mass scale.
       powerSpectrum_sigma        (iWavenumber)=Cosmological_Mass_Root_Variance                       (powerSpectrum_Mass(iWavenumber))
       ! Compute gradient of mass fluctuations.
       powerSpectrum_sigmaGradient(iWavenumber)=Cosmological_Mass_Root_Variance_Logarithmic_Derivative(powerSpectrum_Mass(iWavenumber))
    end do

    return
  end subroutine Power_Spectrum_Compute

  subroutine Power_Spectrum_Output
    !% Outputs power spectrum data.
    use Numerical_Constants_Astronomical
    implicit none
    type(hdf5Object) :: powerSpectrumGroup, thisDataset 
    
    ! Write power spectrum datasets.
    powerSpectrumGroup=powerSpectrumOutputFile%openGroup('powerSpectrum','Group containing datasets relating to&
         & the power spectrum.')

    ! Write the power spectrum data.
    call powerSpectrumGroup%writeDataset(powerSpectrum_Wavenumber,'wavenumber','The wavenumber.'  ,datasetReturned=thisDataset)
    call thisDataset%writeAttribute(1.0d0/megaParsec,'unitsInSI')
    call thisDataset%close()
    call powerSpectrumGroup%writeDataset(powerSpectrum_Power,'powerSpectrum','The power spectrum.',datasetReturned=thisDataset)
    call thisDataset%writeAttribute(megaParsec**3   ,'unitsInSI')
    call thisDataset%close()
    call powerSpectrumGroup%writeDataset(powerSpectrum_Mass,'mass','The corresponding mass scale.',datasetReturned=thisDataset)
    call thisDataset%writeAttribute(massSolar       ,'unitsInSI')
    call thisDataset%close()
    call powerSpectrumGroup%writeDataset(powerSpectrum_sigma        ,'sigma','The mass fluctuation on this scale.')
    call powerSpectrumGroup%writeDataset(powerSpectrum_sigmaGradient,'alpha','Logarithmic deriative of the mass flucation with respect to mass.')
 
    ! Close the datasets group.
    call powerSpectrumGroup%close()

    return
  end subroutine Power_Spectrum_Output

end module Power_Spectrum_Tasks
