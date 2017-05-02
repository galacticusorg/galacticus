!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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
  double precision            , allocatable, dimension(:  ) :: powerSpectrum_mass      , powerSpectrum_power        , &
       &                                                       powerSpectrum_sigma     , powerSpectrum_sigmaGradient, &
       &                                                       powerSpectrum_wavenumber, powerSpectrum_growthFactor , &
       &                                                       powerSpectrum_epochTime , powerSpectrum_epochRedshift
  double precision            , allocatable, dimension(:,:) :: powerSpectrum_nonLinear
  
  ! Options controlling output.
  logical                                                   :: powerSpectrumNonlinearOutput

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
    use, intrinsic :: ISO_C_Binding
    use               Memory_Management
    use               Numerical_Ranges
    use               Input_Parameters
    use               Power_Spectra
    use               Cosmological_Mass_Variance
    use               Numerical_Constants_Math
    use               Cosmology_Parameters
    use               Cosmology_Functions
    use               Galacticus_Output_Times
    use               Linear_Growth
    use               Power_Spectra_Nonlinear
    implicit none
    integer         (c_size_t                     )          :: iWavenumber                  , powerSpectraCount            , &
         &                                                      iOutput                      , outputCount
    integer                                                  :: powerSpectraPointsPerDecade
    double precision                                         :: powerSpectraWavenumberMaximum, powerSpectraWavenumberMinimum
    class           (cosmologyParametersClass     ), pointer :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_
    class           (powerSpectrumClass           ), pointer :: powerSpectrum_
    class           (linearGrowthClass            ), pointer :: linearGrowth_
    class           (powerSpectrumNonlinearClass  ), pointer :: powerSpectrumNonlinear_
    
    ! Find the wavenumber range and increment size.
    !@ <inputParameter>
    !@   <name>powerSpectraWavenumberMinimum</name>
    !@   <defaultValue>$10^{-3}$ Mpc$^{-1}$</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The minimum wavenumber at which to tabulate power spectra.
    !@   </description>
    !@   <type>real</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
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
    !@ <inputParameter>
    !@   <name>powerSpectrumNonlinearOutput</name>
    !@   <defaultValue>true</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@      If true the nonlinear power spectrum will be computed and output.
    !@   </description>
    !@   <type>boolean</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    call Get_Input_Parameter('powerSpectrumNonlinearOutput',powerSpectrumNonlinearOutput,defaultValue=.true.)

    ! Find number of output redshifts.
    outputCount      =Galacticus_Output_Time_Count()
    
    ! Compute number of tabulation points.
    powerSpectraCount=int(log10(powerSpectraWavenumberMaximum/powerSpectraWavenumberMinimum)*dble(powerSpectraPointsPerDecade))+1

    ! Allocate arrays for power spectra.
    call                                   allocateArray(powerSpectrum_Wavenumber   ,[powerSpectraCount            ])
    call                                   allocateArray(powerSpectrum_Power        ,[powerSpectraCount            ])
    call                                   allocateArray(powerSpectrum_Mass         ,[powerSpectraCount            ])
    call                                   allocateArray(powerSpectrum_sigma        ,[powerSpectraCount            ])
    call                                   allocateArray(powerSpectrum_sigmaGradient,[powerSpectraCount            ])
    call                                   allocateArray(powerSpectrum_growthFactor ,[                  outputCount])
    call                                   allocateArray(powerSpectrum_epochTime    ,[                  outputCount])
    call                                   allocateArray(powerSpectrum_epochRedshift,[                  outputCount])
    if (powerSpectrumNonlinearOutput) call allocateArray(powerSpectrum_nonLinear    ,[powerSpectraCount,outputCount])
    
    ! Build a range of wavenumbers.
    powerSpectrum_Wavenumber(:)=Make_Range(powerSpectraWavenumberMinimum,powerSpectraWavenumberMaximum,int(powerSpectraCount),rangeTypeLogarithmic)

    ! Get required objects.
    cosmologyParameters_      => cosmologyParameters     ()
    cosmologyFunctions_       => cosmologyFunctions      ()
    cosmologicalMassVariance_ => cosmologicalMassVariance()
    powerSpectrum_            => powerSpectrum           ()
    linearGrowth_             => linearGrowth            ()
    powerSpectrumNonLinear_   => powerSpectrumNonlinear  ()
    ! Iterate over all wavenumbers computing power spectrum and related quantities.
    do iWavenumber=1,powerSpectraCount
       ! Compute power spectrum.
       powerSpectrum_Power        (iWavenumber)=powerSpectrum_%power(powerSpectrum_Wavenumber(iWavenumber))
       ! Compute corresponding mass scale.
       powerSpectrum_Mass         (iWavenumber)=4.0d0*Pi*cosmologyParameters_%OmegaMatter()*cosmologyParameters_%densityCritical()/3.0d0/powerSpectrum_Wavenumber(iWavenumber)**3
       ! Compute fluctuation on this mass scale.
       powerSpectrum_sigma        (iWavenumber)=cosmologicalMassVariance_%rootVariance                   (powerSpectrum_Mass(iWavenumber))
       ! Compute gradient of mass fluctuations.
       powerSpectrum_sigmaGradient(iWavenumber)=cosmologicalMassVariance_%rootVarianceLogarithmicGradient(powerSpectrum_Mass(iWavenumber))
    end do
    ! Iterate over outputs.
    do iOutput=1,outputCount
       powerSpectrum_epochTime    (iOutput)=                                                       Galacticus_Output_Time(iOutput)
       powerSpectrum_epochRedshift(iOutput)=cosmologyFunctions_ %redshiftFromExpansionFactor(                                      &
            &                                cosmologyFunctions_%expansionFactor             (                                     &
            &                                                                                      Galacticus_Output_Time(iOutput) &
            &                                                                                )                                     &
            &                                                                               )
       powerSpectrum_growthFactor (iOutput)=linearGrowth_       %value                      ( time=Galacticus_Output_Time(iOutput))
       ! Iterate over all wavenumbers computing non-linear power spectrum.
       if (powerSpectrumNonlinearOutput) then
          do iWavenumber=1,powerSpectraCount
             powerSpectrum_nonLinear(iWavenumber,iOutput)=powerSpectrumNonlinear_%value(powerSpectrum_Wavenumber(iWavenumber),Galacticus_Output_Time(iOutput))
          end do
       end if
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
    call powerSpectrumGroup%writeDataset(powerSpectrum_sigma        ,'sigma'             ,'The mass fluctuation on this scale.'                                                          )
    call powerSpectrumGroup%writeDataset(powerSpectrum_sigmaGradient,'alpha'             ,'Logarithmic deriative of the mass flucation with respect to mass.'                            )
    call powerSpectrumGroup%writeDataset(powerSpectrum_growthFactor ,'linearGrowthFactor','The linear growth factor at each epoch.'                                                      )
    call powerSpectrumGroup%writeDataset(powerSpectrum_epochRedshift,'redshift'          ,'The redshift at each epoch.'                                                                  )
    call powerSpectrumGroup%writeDataset(powerSpectrum_epochTime    ,'time'              ,'The time at each epoch.'                                          ,datasetReturned=thisDataset)
    call thisDataset%writeAttribute(gigaYear       ,'unitsInSI')
    call thisDataset%close()
    if (powerSpectrumNonlinearOutput) then
       call powerSpectrumGroup%writeDataset(powerSpectrum_nonLinear,'powerSpectrumNonlinear','The non-linear power spectrum.',datasetReturned=thisDataset)
       call thisDataset%writeAttribute(megaParsec**3   ,'unitsInSI')
       call thisDataset%close()
    end if

    ! Close the datasets group.
    call powerSpectrumGroup%close()

    return
  end subroutine Power_Spectrum_Output

end module Power_Spectrum_Tasks
