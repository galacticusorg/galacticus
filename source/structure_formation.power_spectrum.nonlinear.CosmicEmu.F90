!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Contains a module which implements a nonlinear power spectrum class in which the nonlinear power spectrum is computed using the
!% code of \cite{lawrence_coyote_2010}.

  use FGSL
  use Cosmology_Parameters
  use Cosmology_Functions
  use Cosmological_Mass_Variance
  use Power_Spectra_Primordial

  !# <powerSpectrumNonlinear name="powerSpectrumNonlinearCosmicEmu">
  !#  <description>Provides a nonlinear power spectrum class in which the power spectrum is computed using the code of \cite{lawrence_coyote_2010}.</description>
  !# </powerSpectrumNonlinear>
  type, extends(powerSpectrumNonlinearClass) :: powerSpectrumNonlinearCosmicEmu
     !% A linear transfer function class.
     private
     integer                                                                      :: wavenumberCount
     double precision                                 , allocatable, dimension(:) :: powerSpectrumTable       , wavenumberTable
     type            (fgsl_interp                    )                            :: interpolationObject
     type            (fgsl_interp_accel              )                            :: interpolationAccelerator
     logical                                                                      :: resetInterpolation
     double precision                                                             :: timePrevious
     class           (cosmologyFunctionsClass        ), pointer                   :: cosmologyFunctions_
     class           (cosmologyParametersClass       ), pointer                   :: cosmologyParameters_
     class           (powerSpectrumPrimordialClass   ), pointer                   :: powerSpectrumPrimordial_
     class           (cosmologicalMassVarianceClass  ), pointer                   :: cosmologicalMassVariance_
   contains
     final     ::          cosmicEmuDestructor
     procedure :: value => cosmicEmuValue
  end type powerSpectrumNonlinearCosmicEmu

  interface powerSpectrumNonlinearCosmicEmu
     !% Constructors for the {\normalfont \ttfamily CosmicEmu} nonlinear power spectrum class.
     module procedure cosmicEmuConstructorParameters
     module procedure cosmicEmuConstructorInternal
  end interface powerSpectrumNonlinearCosmicEmu

  ! Wavenumber range used for testing shape of primordial power spectrum.
  double precision, parameter :: cosmicEmuWavenumberLong=0.01d0, cosmicEmuWavenumberShort=1.0d0

contains

  function cosmicEmuConstructorParameters(parameters)
    !% Constructor for the cosmicEmu nonlinear power spectrum class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type (powerSpectrumNonlinearCosmicEmu)                :: cosmicEmuConstructorParameters
    type (inputParameters                ), intent(inout) :: parameters
    class(cosmologyFunctionsClass        ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass       ), pointer       :: cosmologyParameters_
    class(powerSpectrumPrimordialClass   ), pointer       :: powerSpectrumPrimordial_
    class(cosmologicalMassVarianceClass  ), pointer       :: cosmologicalMassVariance_
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    ! Construct required objects.
    !# <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    !# <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    !# <objectBuilder class="powerSpectrumPrimordial"  name="powerSpectrumPrimordial_"  source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    ! Call the internal constructor.
    cosmicEmuConstructorParameters=cosmicEmuConstructorInternal(cosmologyFunctions_,cosmologyParameters_,powerSpectrumPrimordial_,cosmologicalMassVariance_)
    return
  end function cosmicEmuConstructorParameters

  function cosmicEmuConstructorInternal(cosmologyFunctions_,cosmologyParameters_,powerSpectrumPrimordial_,cosmologicalMassVariance_)
    !% Internal constructor for the {\normalfont \ttfamily CosmicEmu} nonlinear power spectrum class.
    use Galacticus_Error
    use Numerical_Comparison
    implicit none
    type (powerSpectrumNonlinearCosmicEmu)                        :: cosmicEmuConstructorInternal
    class(cosmologyFunctionsClass        ), intent(in   ), target :: cosmologyFunctions_
    class(cosmologyParametersClass       ), intent(in   ), target :: cosmologyParameters_
    class(powerSpectrumPrimordialClass   ), intent(in   ), target :: powerSpectrumPrimordial_
    class(cosmologicalMassVarianceClass  ), intent(in   ), target :: cosmologicalMassVariance_

    ! Initialize state.
    cosmicEmuConstructorInternal%timePrevious      =-1.0d0
    cosmicEmuConstructorInternal%resetInterpolation=.true.
    ! Store objects.
    cosmicEmuConstructorInternal%cosmologyFunctions_       => cosmologyFunctions_
    cosmicEmuConstructorInternal%cosmologyParameters_      => cosmologyParameters_
    cosmicEmuConstructorInternal%powerSpectrumPrimordial_  => powerSpectrumPrimordial_
    cosmicEmuConstructorInternal%cosmologicalMassVariance_ => cosmologicalMassVariance_
    ! Check that this is a flat cosmology.
    if     (                                                                                                &
         &  Values_Differ(                                                                                  &
         &                +cosmicEmuConstructorInternal%cosmologyParameters_%OmegaMatter    ()              &
         &                +cosmicEmuConstructorInternal%cosmologyParameters_%OmegaDarkEnergy()            , &
         &                       1.0d+0                                                                   , &
         &                absTol=1.0d-3                                                                     &
         &               )                                                                                  &
         & )                                                                                                &
         & call Galacticus_Error_Report(                                                                    &
         &                              'cosmicEmuConstructorInternal'                                    , &
         &                              'this method is applicable only to flat matter+dark energy models'  &
         &                             )
    ! Check that the primordial power spectrum has no running of the spectral index.
    if     (                                                                                                                     &
         &  Values_Differ(                                                                                                       &
         &                cosmicEmuConstructorInternal%powerSpectrumPrimordial_%logarithmicDerivative(cosmicEmuWavenumberShort), &
         &                cosmicEmuConstructorInternal%powerSpectrumPrimordial_%logarithmicDerivative(cosmicEmuWavenumberLong ), &
         &                relTol=1.0d-3                                                                                          &
         &               )                                                                                                       &
         & )                                                                                                                     &
         & call Galacticus_Error_Report(                                                                                         &
         &                              'cosmicEmuConstructorInternal'                                                         , &
         &                              'this method is applicable only to models with no running of the spectral index'         &
         &                             )
   return
  end function cosmicEmuConstructorInternal

  subroutine cosmicEmuDestructor(self)
    !% Destructor for the {\normalfont \ttfamily CosmicEmu} nonlinear power spectrum class.
    implicit none
    type(powerSpectrumNonlinearCosmicEmu), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"      />
    !# <objectDestructor name="self%cosmologyParameters_"     />
    !# <objectDestructor name="self%powerSpectrumPrimordial_" />
    !# <objectDestructor name="self%cosmologicalMassVariance_"/>
    return
  end subroutine cosmicEmuDestructor

  double precision function cosmicEmuValue(self,waveNumber,time)
    !% Return a nonlinear power spectrum equal using the code of \cite{lawrence_coyote_2010}.
    use Numerical_Interpolation
    use Galacticus_Error
    use FoX_wxml
    use Numerical_Comparison
    use System_Command
    use ISO_Varying_String
    use Galacticus_Input_Paths
    use Input_Parameters
    use Input_Parameters2
    use File_Utilities
    use Memory_Management
    use Table_Labels
    implicit none
    class           (powerSpectrumNonlinearCosmicEmu), intent(inout) :: self
    double precision                                 , intent(in   ) :: time                    , waveNumber
    double precision                                                 :: littleHubbleCMB         , redshift
    type            (varying_string                 )                :: parameterFile           , powerSpectrumFile
    type            (xmlf_t                         )                :: parameterDoc
    character       (len=32                         )                :: parameterLabel
    character       (len=128                        )                :: powerSpectrumLine
    integer                                                          :: iWavenumber             , powerSpectrumUnit
    type            (inputParameterList             )                :: parameters

    ! If the time has changed, recompute the power spectrum.
    !$omp critical(cosmicEmuValue)
    if (time /= self%timePrevious) then
       ! Store the new time and find the corresponding redshift.
       self%timePrevious=time
       redshift         =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
       ! Generate a parameter file.
       parameters=inputParameterList()
       write (parameterLabel,'(f5.3)') self%cosmologyParameters_     %OmegaMatter          (                        )
       call parameters%add("OmegaMatter"              ,parameterLabel)
       write (parameterLabel,'(f6.4)') self%cosmologyParameters_     %OmegaBaryon          (                        )
       call parameters%add("OmegaBaryon"              ,parameterLabel)
       write (parameterLabel,'(f7.4)') self%cosmologyParameters_     %HubbleConstant       (                        )
       call parameters%add("HubbleConstant"           ,parameterLabel)
       write (parameterLabel,'(f6.4)') self%cosmologicalMassVariance_%sigma8               (                        )
       call parameters%add("sigma_8"                  ,parameterLabel)
       write (parameterLabel,'(f6.4)') self%powerSpectrumPrimordial_ %logarithmicDerivative(cosmicEmuWavenumberShort)
       call parameters%add("powerSpectrumIndex"       ,parameterLabel)
       write (parameterLabel,'(f6.3)') -1.0d0
       call parameters%add("darkEnergyEquationOfState",parameterLabel)
       write (parameterLabel,'(f8.4)') redshift
       call parameters%add("redshift"                 ,parameterLabel)       
       powerSpectrumFile="powerSpectrum.txt"
       parameterFile    ="powerSpectrumParameters.xml"
       call xml_OpenFile(char(parameterFile),parameterDoc)
       call xml_NewElement(parameterDoc,"parameters")
       call parameters%serializeToXML(parameterDoc)
       call xml_Close(parameterDoc)
       ! Generate the power spectrum.
       call System_Command_Do(Galacticus_Input_Path()//"scripts/aux/Cosmic_Emu_Driver.pl "//parameterFile//" "//powerSpectrumFile)
       ! Read the data file.
       self%wavenumberCount=Count_Lines_In_File(powerSpectrumFile,"#")
       if (allocated(self%wavenumberTable   )) call Dealloc_Array(self%wavenumberTable   )
       if (allocated(self%powerSpectrumTable)) call Dealloc_Array(self%powerSpectrumTable)
       call Alloc_Array(self%wavenumberTable   ,[self%wavenumberCount])
       call Alloc_Array(self%powerSpectrumTable,[self%wavenumberCount])
       open(newunit=powerSpectrumUnit,file=char(powerSpectrumFile),status='old',form='formatted')
       iWavenumber=0
       do while (iWavenumber < self%wavenumberCount)
          read (powerSpectrumUnit,'(a)') powerSpectrumLine
          if (powerSpectrumLine(1:1) == "#") then
             if (powerSpectrumLine(1:33) == "# dimensionless Hubble parameter") then
                read (powerSpectrumLine(index(powerSpectrumLine,":")+1:),*) littleHubbleCMB
                if (Values_Differ(littleHubbleCMB,self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH),relTol=1.0d-2)) &
                     & call Galacticus_Error_Report(                                                                           &
                     &                              'cosmicEmuValue'                                                       ,   &
                     &                              'values of H_0 in Galacticus and CosmicEmu are significantly different'    &
                     &                             )
             end if
          else
             iWavenumber=iWavenumber+1
             read (powerSpectrumLine,*) self%wavenumberTable(iWavenumber),self%powerSpectrumTable(iWavenumber)
          end if
       end do
       close(powerSpectrumUnit)
       ! Convert to logarithmic values.
       self%wavenumberTable   =log(self%wavenumberTable   )
       self%powerSpectrumTable=log(self%powerSpectrumTable)
       ! Reset the interpolator.
       call Interpolate_Done(self%interpolationObject,self%interpolationAccelerator,self%resetInterpolation)
       self%resetInterpolation=.true.
       ! Destroy the parameter and power spectrum files.
       call System_Command_Do("rm -f "//parameterFile//" "//powerSpectrumFile)
    end if
    ! Interpolate in the tabulated data to get the power spectrum.
    cosmicEmuValue=exp(Interpolate(self%wavenumberTable,self%powerSpectrumTable,self%interpolationObject&
         &,self%interpolationAccelerator,log(wavenumber),reset=self%resetInterpolation,extrapolationType=extrapolationTypeExtrapolate))
    !$omp end critical(cosmicEmuValue)
    return
  end function cosmicEmuValue
