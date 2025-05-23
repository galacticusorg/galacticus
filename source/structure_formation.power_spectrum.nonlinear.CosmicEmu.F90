!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

!!{
Implements a nonlinear power spectrum class in which the nonlinear power spectrum is computed using the
code of \cite{lawrence_coyote_2010}.
!!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: File_Utilities            , only : lockDescriptor
  use :: Numerical_Interpolation   , only : interpolator
  use :: Power_Spectra_Primordial  , only : powerSpectrumPrimordialClass

  !![
  <powerSpectrumNonlinear name="powerSpectrumNonlinearCosmicEmu">
   <description>
    Provides a nonlinear power spectrum class in which the power spectrum is computed using the code of
    \cite{lawrence_coyote_2010}. The CosmicEmu code will be downloaded, compiled and run as necessary if this option is
    utilized.
   </description>
  </powerSpectrumNonlinear>
  !!]
  type, extends(powerSpectrumNonlinearClass) :: powerSpectrumNonlinearCosmicEmu
     !!{
     A linear transfer function class.
     !!}
     private
     integer                                                                    :: wavenumberCount
     double precision                               , allocatable, dimension(:) :: powerSpectrumTable                 , wavenumberTable
     type            (interpolator                 ), allocatable               :: interpolator_
     double precision                                                           :: timePrevious
     type            (lockDescriptor               )                            :: fileLock
     class           (cosmologyFunctionsClass      ), pointer                   :: cosmologyFunctions_       => null()
     class           (cosmologyParametersClass     ), pointer                   :: cosmologyParameters_      => null()
     class           (powerSpectrumPrimordialClass ), pointer                   :: powerSpectrumPrimordial_  => null()
     class           (cosmologicalMassVarianceClass), pointer                   :: cosmologicalMassVariance_ => null()
   contains
     final     ::          cosmicEmuDestructor
     procedure :: value => cosmicEmuValue
  end type powerSpectrumNonlinearCosmicEmu

  interface powerSpectrumNonlinearCosmicEmu
     !!{
     Constructors for the \refClass{powerSpectrumNonlinearCosmicEmu} nonlinear power spectrum class.
     !!}
     module procedure cosmicEmuConstructorParameters
     module procedure cosmicEmuConstructorInternal
  end interface powerSpectrumNonlinearCosmicEmu

  ! Wavenumber range used for testing shape of primordial power spectrum.
  double precision, parameter :: wavenumberLong=0.01d0, wavenumberShort=1.0d0

contains

  function cosmicEmuConstructorParameters(parameters) result(self)
    !!{
    Constructor for the cosmicEmu nonlinear power spectrum class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (powerSpectrumNonlinearCosmicEmu)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(cosmologyFunctionsClass        ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass       ), pointer       :: cosmologyParameters_
    class(powerSpectrumPrimordialClass   ), pointer       :: powerSpectrumPrimordial_
    class(cosmologicalMassVarianceClass  ), pointer       :: cosmologicalMassVariance_

    ! Check and read parameters.
    ! Construct required objects.
    !![
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="powerSpectrumPrimordial"  name="powerSpectrumPrimordial_"  source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    ! Call the internal constructor.
    self=powerSpectrumNonlinearCosmicEmu(cosmologyFunctions_,cosmologyParameters_,powerSpectrumPrimordial_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="powerSpectrumPrimordial_" />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function cosmicEmuConstructorParameters

  function cosmicEmuConstructorInternal(cosmologyFunctions_,cosmologyParameters_,powerSpectrumPrimordial_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the \refClass{powerSpectrumNonlinearCosmicEmu} nonlinear power spectrum class.
    !!}
    use :: Error               , only : Error_Report
    use :: Numerical_Comparison, only : Values_Differ
    implicit none
    type (powerSpectrumNonlinearCosmicEmu)                        :: self
    class(cosmologyFunctionsClass        ), intent(in   ), target :: cosmologyFunctions_
    class(cosmologyParametersClass       ), intent(in   ), target :: cosmologyParameters_
    class(powerSpectrumPrimordialClass   ), intent(in   ), target :: powerSpectrumPrimordial_
    class(cosmologicalMassVarianceClass  ), intent(in   ), target :: cosmologicalMassVariance_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *cosmologyParameters_, *powerSpectrumPrimordial_, *cosmologicalMassVariance_"/>
    !!]

    ! Initialize state.
    self%timePrevious=-1.0d0
    ! Check that this is a flat cosmology.
    if     (                                                                                      &
         &  Values_Differ(                                                                        &
         &                +self%cosmologyParameters_%OmegaMatter    ()                            &
         &                +self%cosmologyParameters_%OmegaDarkEnergy()                         ,  &
         &                       1.0d+0                                                        ,  &
         &                absTol=1.0d-3                                                           &
         &               )                                                                        &
         & )                                                                                      &
         & call Error_Report(                                                                     &
         &                   'this method is applicable only to flat matter+dark energy models'// &
         &                    {introspection:location}                                            &
         &                  )
    ! Check that the primordial power spectrum has no running of the spectral index.
    if     (                                                                                                    &
         &  Values_Differ(                                                                                      &
         &                self%powerSpectrumPrimordial_%logarithmicDerivative(wavenumberShort),                 &
         &                self%powerSpectrumPrimordial_%logarithmicDerivative(wavenumberLong ),                 &
         &                relTol=1.0d-3                                                                         &
         &               )                                                                                      &
         & )                                                                                                    &
         & call Error_Report(                                                                                   &
         &                   'this method is applicable only to models with no running of the spectral index'// &
         &                    {introspection:location}                                                          &
         &                  )
   return
  end function cosmicEmuConstructorInternal

  subroutine cosmicEmuDestructor(self)
    !!{
    Destructor for the \refClass{powerSpectrumNonlinearCosmicEmu} nonlinear power spectrum class.
    !!}
    implicit none
    type(powerSpectrumNonlinearCosmicEmu), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%powerSpectrumPrimordial_" />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine cosmicEmuDestructor

  double precision function cosmicEmuValue(self,waveNumber,time)
    !!{
    Return a nonlinear power spectrum equal using the code of \cite{lawrence_coyote_2010}.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsLittleH
    use :: Display             , only : displayMessage              , verbosityLevelWorking
    use :: File_Utilities      , only : Count_Lines_In_File         , Directory_Make       , File_Exists, File_Lock, &
          &                             File_Name_Temporary         , File_Remove          , File_Unlock
    use :: Error               , only : Error_Report
    use :: Input_Paths         , only : inputPath                   , pathTypeDataDynamic
    use :: ISO_Varying_String  , only : varying_string
    use :: Numerical_Comparison, only : Values_Differ
    use :: System_Command      , only : System_Command_Do
    use :: System_Download     , only : download
    use :: Table_Labels        , only : extrapolationTypeExtrapolate
    implicit none
    class           (powerSpectrumNonlinearCosmicEmu), intent(inout) :: self
    double precision                                 , intent(in   ) :: time             , waveNumber
    double precision                                                 :: littleHubbleCMB  , redshift
    type            (varying_string                 )                :: parameterFile    , powerSpectrumFile, &
         &                                                              parameters
    character       (len=32                         )                :: parameterLabel
    character       (len=128                        )                :: powerSpectrumLine
    integer                                                          :: iWavenumber      , powerSpectrumUnit, &
         &                                                              status

    ! If the time has changed, recompute the power spectrum.
    if (time /= self%timePrevious) then
       ! Store the new time and find the corresponding redshift.
       self%timePrevious=time
       redshift         =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
       ! Generate parameters and a file name for this power spectrum.
       call Directory_Make(inputPath(pathTypeDataDynamic)//"largeScaleStructure")
       powerSpectrumFile=inputPath(pathTypeDataDynamic)//"largeScaleStructure/powerSpectrumCosmicEmu"
       parameterFile    =File_Name_Temporary("cosmicEmuParameters")
       parameters       =''
       write (parameterLabel,'(f5.3)') +self%cosmologyParameters_     %OmegaMatter              (                                  )
       powerSpectrumFile=powerSpectrumFile//"_OmegaMatter"//trim(adjustl(parameterLabel))
       write (parameterLabel,'(f5.3)') +self%cosmologyParameters_     %OmegaMatter              (                                  )    &
            &                          *self%cosmologyParameters_     %HubbleConstant           (                hubbleUnitsLittleH)**2
       parameters=parameters//trim(adjustl(parameterLabel))//char(10)
       write (parameterLabel,'(f6.4)') +self%cosmologyParameters_     %OmegaBaryon              (                                  )
       powerSpectrumFile=powerSpectrumFile//"_OmegaBaryon"//trim(adjustl(parameterLabel))
       write (parameterLabel,'(f5.3)') +self%cosmologyParameters_     %OmegaBaryon              (                                  )    &
            &                          *self%cosmologyParameters_     %HubbleConstant           (                hubbleUnitsLittleH)**2
       parameters=parameters//trim(adjustl(parameterLabel))//char(10)
       write (parameterLabel,'(f7.4)') +self%cosmologyParameters_     %HubbleConstant           (                                  )
       powerSpectrumFile=powerSpectrumFile//"_HubbleConstant"//trim(adjustl(parameterLabel))
       write (parameterLabel,'(f6.4)') +self%powerSpectrumPrimordial_ %logarithmicDerivative    (                wavenumberShort   )
       powerSpectrumFile=powerSpectrumFile//"_powerSpectrumIndex"//trim(adjustl(parameterLabel))
       parameters=parameters//trim(adjustl(parameterLabel))//char(10)
       write (parameterLabel,'(f6.4)') +self%cosmologicalMassVariance_%sigma8                   (                                  )
       powerSpectrumFile=powerSpectrumFile//"_sigma8"//trim(adjustl(parameterLabel))
       parameters=parameters//trim(adjustl(parameterLabel))//char(10)
       write (parameterLabel,'(f6.3)') +self%cosmologyFunctions_      %equationOfStateDarkEnergy(expansionFactor=1.0d0             )
       powerSpectrumFile=powerSpectrumFile//"_w"//trim(adjustl(parameterLabel))
       parameters=parameters//trim(adjustl(parameterLabel))//char(10)
       write (parameterLabel,'(f7.4)') +redshift
       powerSpectrumFile=powerSpectrumFile//"_redshift"//trim(adjustl(parameterLabel))
       parameters=parameters//trim(adjustl(parameterLabel))//char(10)
       powerSpectrumFile=powerSpectrumFile//".txt"
       parameters=powerSpectrumFile//char(10)//parameters//'2'//char(10)
       ! Check for existence of the power spectrum, building it if necessary.
       call File_Lock(char(powerSpectrumFile),self%fileLock,lockIsShared=.true.)
       if (.not.File_Exists(char(powerSpectrumFile))) then
          call File_Unlock(self%fileLock)
          call File_Lock(char(powerSpectrumFile),self%fileLock,lockIsShared=.false.)
          open(newUnit=powerSpectrumUnit,file=char(parameterFile),status='unknown',form='formatted')
          write (powerSpectrumUnit,'(a)') char(parameters)
          close(powerSpectrumUnit)
          ! Check for presence of the executable.
          call Directory_Make(inputPath(pathTypeDataDynamic)//"CosmicEmu_v1.1")
          if (.not.File_Exists(inputPath(pathTypeDataDynamic)//"CosmicEmu_v1.1/emu.exe")) then
             ! Check for presence of the source code.
             if (.not.File_Exists(inputPath(pathTypeDataDynamic)//"CosmicEmu_v1.1/emu.c")) then
                ! Download the code.
                if (.not.File_Exists(inputPath(pathTypeDataDynamic)//"CosmicEmu_v1.1.tar.gz")) then
                   call displayMessage("downloading CosmicEmu code....",verbosityLevelWorking)
                   call download("http://www.hep.anl.gov/cosmology/CosmicEmu/CosmicEmu_v1.1.tar.gz",char(inputPath(pathTypeDataDynamic))//"CosmicEmu_v1.1.tar.gz",status=status)
                   if (status /= 0 .or. .not.File_Exists(inputPath(pathTypeDataDynamic)//"CosmicEmu_v1.1.tar.gz")) &
                        & call Error_Report("failed to download CosmicEmu code"//{introspection:location})
                end if
                ! Unpack the code.
                call displayMessage("unpacking CosmicEmu code....",verbosityLevelWorking)
                call System_Command_Do("tar -x -v -z -C "//inputPath(pathTypeDataDynamic)//"CosmicEmu_v1.1 -f "//inputPath(pathTypeDataDynamic)//"CosmicEmu_v1.1.tar.gz")
                if (.not.File_Exists(inputPath(pathTypeDataDynamic)//"CosmicEmu_v1.1/emu.c")) &
                     & call Error_Report("failed to unpack CosmicEmu code"//{introspection:location})
             end if
             ! Build the code.
             call displayMessage("compiling CosmicEmu code....",verbosityLevelWorking)
             call System_Command_Do("cd "//inputPath(pathTypeDataDynamic)//"CosmicEmu_v1.1; sed -i~ -r s/""^(\s*gcc.*\-lm)\s*$""/""\1 \-I\`gsl\-config \-\-prefix\`\n\n%.o: %.c\n\tgcc -c \$< -o \$\*\.o \-I\`gsl\-config \-\-prefix\`\n""/ makefile; sed -i~ -r s/""\-lgsl\s+\-lgslcblas\s+\-lm""/""\`gsl\-config \-\-libs\`""/ makefile; make");
             if (.not.File_Exists(inputPath(pathTypeDataDynamic)//"CosmicEmu_v1.1/emu.exe")) &
                  & call Error_Report("failed to build Cosmic_Emu code"//{introspection:location})
          end if
          ! Generate the power spectrum.
          call System_Command_Do(inputPath(pathTypeDataDynamic)//"CosmicEmu_v1.1/emu.exe < "//parameterFile)
          ! Destroy the parameter file.
          call File_Remove(char(parameterFile))
       end if
       ! Read the data file.
       self%wavenumberCount=Count_Lines_In_File(powerSpectrumFile,"#")
       if (allocated(self%wavenumberTable   )) deallocate(self%wavenumberTable   )
       if (allocated(self%powerSpectrumTable)) deallocate(self%powerSpectrumTable)
       allocate(self%wavenumberTable   (self%wavenumberCount))
       allocate(self%powerSpectrumTable(self%wavenumberCount))
       open(newunit=powerSpectrumUnit,file=char(powerSpectrumFile),status='old',form='formatted')
       iWavenumber=0
       do while (iWavenumber < self%wavenumberCount)
          read (powerSpectrumUnit,'(a)') powerSpectrumLine
          if (powerSpectrumLine(1:1) == "#") then
             if (powerSpectrumLine(1:33) == "# dimensionless Hubble parameter") then
                read (powerSpectrumLine(index(powerSpectrumLine,":")+1:),*) littleHubbleCMB
                if (Values_Differ(littleHubbleCMB,self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH),relTol=1.0d-2)) &
                     & call Error_Report(                                                                                      &
                     &                   'values of H₀ in Galacticus and CosmicEmu are significantly different'//              &
                     &                    {introspection:location}                                                             &
                     &                  )
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
       ! Build the interpolator.
       if (allocated(self%interpolator_)) deallocate(self%interpolator_)
       allocate(self%interpolator_)
       self%interpolator_=interpolator(self%wavenumberTable,self%powerSpectrumTable,extrapolationType=extrapolationTypeExtrapolate)
       call File_Unlock(self%fileLock)
    end if
    ! Interpolate in the tabulated data to get the power spectrum.
    cosmicEmuValue=exp(self%interpolator_%interpolate(log(wavenumber)))
    return
  end function cosmicEmuValue
