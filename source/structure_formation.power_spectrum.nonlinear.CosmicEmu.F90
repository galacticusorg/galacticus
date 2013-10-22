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

!% Contains a module which implements the nonlinear power spectrum using the code of \cite{lawrence_coyote_2010}.

module Power_Spectra_Nonlinear_CosmicEmu
  !% Implements the nonlinear power spectrum using the code of \cite{lawrence_coyote_2010}.
  use FGSL
  implicit none
  private
  public :: Power_Spectrum_Nonlinear_CosmicEmu_Initialize

  ! Arrays to hold the power spectrum.
  integer                                                        :: wavenumberCount
  double precision                   , allocatable, dimension(:) :: powerSpectrumTable             , wavenumberTable

  ! Interpolators.
  type            (fgsl_interp      )                            :: interpolationObject
  type            (fgsl_interp_accel)                            :: interpolationAccelerator
  logical                                                        :: resetInterpolation      =.true.

contains

  !# <powerSpectrumNonlinearMethod>
  !#  <unitName>Power_Spectrum_Nonlinear_CosmicEmu_Initialize</unitName>
  !# </powerSpectrumNonlinearMethod>
  subroutine Power_Spectrum_Nonlinear_CosmicEmu_Initialize(powerSpectrumNonlinearMethod,Power_Spectrum_Nonlinear_Get)
    !% Initializes the ``CosmicEmu'' nonlinear power spectrum module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string  ), intent(in   )          :: powerSpectrumNonlinearMethod
    procedure(Power_Spectrum_Nonlinear_CosmicEmu), intent(inout), pointer :: Power_Spectrum_Nonlinear_Get

    if (powerSpectrumNonlinearMethod == 'CosmicEmu') Power_Spectrum_Nonlinear_Get => Power_Spectrum_Nonlinear_CosmicEmu
    return
  end subroutine Power_Spectrum_Nonlinear_CosmicEmu_Initialize

  double precision function Power_Spectrum_Nonlinear_CosmicEmu(waveNumber,time)
    !% Return a nonlinear power spectrum equal using the code of \cite{lawrence_coyote_2010}.
    use Numerical_Interpolation
    use Cosmology_Parameters
    use Cosmology_Functions
    use Galacticus_Error
    use FoX_wxml
    use Numerical_Comparison
    use Primordial_Power_Spectra
    use System_Command
    use ISO_Varying_String
    use Galacticus_Input_Paths
    use Power_Spectra
    use Input_Parameters
    use File_Utilities
    use Memory_Management
    implicit none
    double precision                , intent(in   ) :: time                    , waveNumber
    double precision                , save          :: timePrevious     =-1.0d0
    double precision                , parameter     :: wavenumberLong   =0.01d0, wavenumberShort  =1.0d0
    class(cosmologyFunctionsClass), pointer                    :: cosmologyFunctionsDefault
    class    (cosmologyParametersClass              )               , pointer :: thisCosmologyParameters
    double precision                                :: littleHubbleCMB         , redshift
    type            (varying_string)                :: parameterFile           , powerSpectrumFile
    type            (xmlf_t        )                :: parameterDoc
    character       (len=32        )                :: parameterLabel
    character       (len=128       )                :: powerSpectrumLine
    integer                                         :: iWavenumber             , powerSpectrumUnit

    ! If the time has changed, recompute the power spectrum.
    if (time /= timePrevious) then
       ! Get the default cosmology.
       thisCosmologyParameters => cosmologyParameters()
       ! Get the default cosmology functions object.
       cosmologyFunctionsDefault => cosmologyFunctions()     
 
       ! Store the new time and find the corresponding redshift.
       timePrevious=time
       redshift=cosmologyFunctionsDefault%redshiftFromExpansionFactor(cosmologyFunctionsDefault%expansionFactor(time))

       ! Check that this is a flat cosmology.
       if (Values_Differ(thisCosmologyParameters%OmegaMatter()+thisCosmologyParameters%OmegaDarkEnergy(),1.0d0,absTol=1.0d-3))                                       &
            & call Galacticus_Error_Report(                                                                    &
            &                              'Power_Spectrum_Nonlinear_CosmicEmu'                              , &
            &                              'this method is applicable only to flat matter+dark energy models'  &
            &                             )

       ! Check that the primordial power spectrum has no running of the spectral index.
       if     (                                                                                                              &
            &  Values_Differ(                                                                                                &
            &                Primordial_Power_Spectrum_Logarithmic_Derivative(waveNumberShort),                              &
            &                Primordial_Power_Spectrum_Logarithmic_Derivative(waveNumberLong ),                              &
            &                relTol=1.0d-3                                                                                   &
            &               )                                                                                                &
            & )                                                                                                              &
            & call Galacticus_Error_Report(                                                                                  &
            &                              'Power_Spectrum_Nonlinear_CosmicEmu'                                            , &
            &                              'this method is applicable only to models with no running of the spectral index'  &
            &                             )
       ! Generate a parameter file.
       powerSpectrumFile="powerSpectrum.txt"
       parameterFile    ="powerSpectrumParameters.xml"
       call xml_OpenFile(char(parameterFile),parameterDoc)
       call xml_NewElement(parameterDoc,"parameters")
       write (parameterLabel,'(f5.3)') thisCosmologyParameters%OmegaMatter   ()
       call Write_Parameter(parameterDoc,"Omega_Matter"             ,parameterLabel)
       write (parameterLabel,'(f6.4)') thisCosmologyParameters%OmegaBaryon   ()
       call Write_Parameter(parameterDoc,"Omega_b"                  ,parameterLabel)
       write (parameterLabel,'(f7.4)') thisCosmologyParameters%HubbleConstant()
       call Write_Parameter(parameterDoc,"H_0"                      ,parameterLabel)
       write (parameterLabel,'(f6.4)') sigma_8     ()
       call Write_Parameter(parameterDoc,"sigma_8"                  ,parameterLabel)
       write (parameterLabel,'(f6.4)') Primordial_Power_Spectrum_Logarithmic_Derivative(waveNumberShort)
       call Write_Parameter(parameterDoc,"powerSpectrumIndex"       ,parameterLabel)
       write (parameterLabel,'(f6.3)') -1.0d0
       call Write_Parameter(parameterDoc,"darkEnergyEquationOfState",parameterLabel)
       write (parameterLabel,'(f8.4)') redshift
       call Write_Parameter(parameterDoc,"redshift"                 ,parameterLabel)
       call xml_Close(parameterDoc)

       ! Generate the power spectrum.
       call System_Command_Do(Galacticus_Input_Path()//"scripts/aux/Cosmic_Emu_Driver.pl "//parameterFile//" "//powerSpectrumFile)

       ! Read the data file.
       wavenumberCount=Count_Lines_In_File(powerSpectrumFile,"#")
       if (allocated(wavenumberTable   )) call Dealloc_Array(wavenumberTable   )
       if (allocated(powerSpectrumTable)) call Dealloc_Array(powerSpectrumTable)
       call Alloc_Array(wavenumberTable   ,[wavenumberCount])
       call Alloc_Array(powerSpectrumTable,[wavenumberCount])
       open(newunit=powerSpectrumUnit,file=char(powerSpectrumFile),status='old',form='formatted')
       iWavenumber=0
       do while (iWavenumber < wavenumberCount)
          read (powerSpectrumUnit,'(a)') powerSpectrumLine
          if (powerSpectrumLine(1:1) == "#") then
             if (powerSpectrumLine(1:33) == "# dimensionless Hubble parameter") then
                read (powerSpectrumLine(index(powerSpectrumLine,":")+1:),*) littleHubbleCMB
                if (Values_Differ(littleHubbleCMB,thisCosmologyParameters%HubbleConstant(unitsLittleH),relTol=1.0d-2)) &
                     call Galacticus_Error_Report(                                                                         &
                     &                            'Power_Spectrum_Nonlinear_CosmicEmu'                                   , &
                     &                            'values of H_0 in Galacticus and CosmicEmu are significantly different'  &
                     &                           )
             end if
          else
             iWavenumber=iWavenumber+1
             read (powerSpectrumLine,*) wavenumberTable(iWavenumber),powerSpectrumTable(iWavenumber)
          end if
       end do
       close(powerSpectrumUnit)

       ! Convert to logarithmic values.
       wavenumberTable   =log(wavenumberTable   )
       powerSpectrumTable=log(powerSpectrumTable)

       ! Reset the interpolator.
       call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)
       resetInterpolation=.true.

       ! Destroy the parameter and power spectrum files.
       call System_Command_Do("rm -f "//parameterFile//" "//powerSpectrumFile)

    end if

    ! Interpolate in the tabulated data to get the power spectrum.
    Power_Spectrum_Nonlinear_CosmicEmu=exp(Interpolate(wavenumberCount,wavenumberTable,powerSpectrumTable ,interpolationObject&
         &,interpolationAccelerator,log(wavenumber),reset=resetInterpolation,extrapolationType=extrapolationTypeLinear))

    return
  end function Power_Spectrum_Nonlinear_CosmicEmu

end module Power_Spectra_Nonlinear_CosmicEmu
