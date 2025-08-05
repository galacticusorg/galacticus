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
  Implements a chemical state class which utilizes the {\normalfont \scshape Cloudy} code to
  compute state in collisional ionization equilibrium.
  !!}

  !![
  <chemicalState name="chemicalStateAtomicCIECloudy">
   <description>
    A chemical state class that computes the chemical state using the {\normalfont \scshape Cloudy} code and under the
    assumption of collisional ionization equilibrium with no molecular contribution. Abundances are Solar, except for zero
    metallicity calculations which use {\normalfont \scshape Cloudy}'s ``primordial'' metallicity. The helium abundance for
    non-zero metallicity is scaled between primordial and Solar values linearly with metallicity. The {\normalfont \scshape
    Cloudy} code will be downloaded and run to compute the cooling function as needed, which will then be stored for future
    use. As this process is slow, a precomputed table is provided with \glc. If metallicities outside the range tabulated in
    this file are required it will be regenerated with an appropriate range.
   </description>
  </chemicalState>
  !!]
  type, extends(chemicalStateCIEFile) :: chemicalStateAtomicCIECloudy
     !!{
     A chemical state class which utilizes the {\normalfont \scshape Cloudy} code to compute
     state in collisional ionization equilibrium.
     !!}
     private
     logical :: initialized
   contains
     !![
     <methods>
       <method description="Run {\normalfont \scshape Cloudy} to tabulate chemical state as necessary." method="tabulate" />
     </methods>
     !!]
     procedure :: tabulate                           => atomicCIECloudyTabulate
     procedure :: electronDensity                    => atomicCIECloudyElectronDensity
     procedure :: electronDensityTemperatureLogSlope => atomicCIECloudyElectronDensityTemperatureLogSlope
     procedure :: electronDensityDensityLogSlope     => atomicCIECloudyElectronDensityDensityLogSlope
     procedure :: chemicalDensities                  => atomicCIECloudyChemicalDensities
  end type chemicalStateAtomicCIECloudy

  interface chemicalStateAtomicCIECloudy
     !!{
     Constructors for the \refClass{chemicalStateAtomicCIECloudy} chemical state class.
     !!}
     module procedure atomicCIECloudyConstructorParameters
     module procedure atomicCIECloudyConstructorInternal
  end interface chemicalStateAtomicCIECloudy

  ! File names for the cooling function and chemical state data.
  character       (len=56)           :: fileNameCoolingFunction  ='cooling/cooling_function_Atomic_CIE_Cloudy.hdf5'
  character       (len=56)           :: fileNameChemicalState    ='chemicalState/chemical_state_Atomic_CIE_Cloudy.hdf5'

  ! Maximum tabulated metallicity.
  double precision       , parameter :: metallicityMaximumDefault=30.0d0 ! Thirty times Solar.

  ! Maximum metallicity that we will ever tabulate. More than this should be physically implausible.
  double precision       , parameter :: metallicityMaximumLimit  =30.0d0 ! Thirty times Solar.

  ! Factor by which metallicity must exceed currently tabulated maximum before we retabulate.
  double precision       , parameter :: metallicityTolerance     = 0.1d0

contains

  function atomicCIECloudyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{chemicalStateAtomicCIECloudy} chemical state class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(chemicalStateAtomicCIECloudy)                :: self
    type(inputParameters             ), intent(inout) :: parameters

    self=chemicalStateAtomicCIECloudy()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function atomicCIECloudyConstructorParameters

  function atomicCIECloudyConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{chemicalStateAtomicCIECloudy} chemical state class.
    !!}
    use :: Chemical_Abundances_Structure, only : unitChemicalAbundances
    implicit none
    type(chemicalStateAtomicCIECloudy) :: self

    ! Initialize.
    self%electronDensityPrevious        =-1.0d0
    self%electronDensitySlopePrevious   =-1.0d0
    self%chemicalDensitiesPrevious      =-unitChemicalAbundances
    self%    metallicityPrevious        =-1.0d0
    self%    metallicitySlopePrevious   =-1.0d0
    self%    metallicityChemicalPrevious=-1.0d0
    self%    temperaturePrevious        =-1.0d0
    self%    temperatureSlopePrevious   =-1.0d0
    self%    temperatureChemicalPrevious=-1.0d0
    self%initialized                    =.false.
   return
  end function atomicCIECloudyConstructorInternal

  subroutine atomicCIECloudyTabulate(self,gasAbundances)
    !!{
    Create the chemical state.
    !!}
    use :: Abundances_Structure , only : Abundances_Get_Metallicity   , metallicityTypeLinearByMassSolar
    use :: File_Utilities       , only : File_Remove
    use :: Input_Paths          , only : inputPath                    , pathTypeDataStatic
    use :: Interfaces_Cloudy_CIE, only : Interface_Cloudy_CIE_Tabulate
    use :: String_Handling      , only : operator(//)
    implicit none
    class  (chemicalStateAtomicCIECloudy), intent(inout) :: self
    type   (abundances                  ), intent(in   ) :: gasAbundances
    logical                                              :: makeFile

    ! Determine if we need to retabulate the chemical state.
    if (.not.self%initialized) then
       makeFile=.true.
    else
       if (Abundances_Get_Metallicity(gasAbundances,metallicityTypeLinearByMassSolar) > 0.0d0) then
          makeFile=(                                                                                 &
               &     min(                                                                            &
               &         Abundances_Get_Metallicity(gasAbundances,metallicityTypeLinearByMassSolar), &
               &         metallicityMaximumLimit                                                     &
               &        )                                                                            &
               &    >                                                                                &
               &     self%metallicityMaximum*(1.0d0+metallicityTolerance)                            &
               &   )
       else
          makeFile=.false.
       end if
       ! Remove the chemical state file so that a new one will be created.
       if (makeFile) call File_Remove(inputPath(pathTypeDataStatic)//trim(fileNameChemicalState))
    end if
    ! Read the file if this module has not been initialized or if the metallicity is out of range.
    if (makeFile) then
       ! Determine maximum metallicity to which we should tabulate.
       if (Abundances_Get_Metallicity(gasAbundances,metallicityTypeLinearByMassSolar) > 0.0d0) then
          self%metallicityMaximum=min(                                                                                &
               &                      max(                                                                            &
               &                           metallicityMaximumDefault                                                , &
               &                          +3.0d0                                                                      &
               &                          *Abundances_Get_Metallicity(gasAbundances,metallicityTypeLinearByMassSolar) &
               &                         )                                                                          , &
               &                           metallicityMaximumLimit                                                    &
               &                     )
       else
          self%metallicityMaximum=metallicityMaximumDefault
       end if
       ! Generate the file.
       call Interface_Cloudy_CIE_Tabulate(                                                                    &
            &                             log10(self%metallicityMaximum                                    ), &
            &                                   inputPath(pathTypeDataStatic)//trim(fileNameCoolingFunction), &
            &                                   inputPath(pathTypeDataStatic)//trim(fileNameChemicalState  ), &
            &                                   fileFormatVersionCurrent                                      &
            &                            )
       ! Call routine to read in the tabulated data.
       call self%readFile(char(inputPath(pathTypeDataStatic)//trim(fileNameChemicalState)))
       ! Flag that chemical state is now initialized.
       self%initialized=.true.
    end if
    return
  end subroutine atomicCIECloudyTabulate

  double precision function atomicCIECloudyElectronDensity(self,numberDensityHydrogen,temperature,gasAbundances,radiation)
    !!{
    Return the electron density for collisional ionization equilibrium as computed by
    {\normalfont \scshape Cloudy}.
    !!}
    implicit none
    class           (chemicalStateAtomicCIECloudy), intent(inout) :: self
    double precision                              , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances                  ), intent(in   ) :: gasAbundances
    class           (radiationFieldClass         ), intent(inout) :: radiation

    call                           self%tabulate                            (                                  gasAbundances          )
    atomicCIECloudyElectronDensity=self%chemicalStateCIEFile%electronDensity(numberDensityHydrogen,temperature,gasAbundances,radiation)
    return
  end function atomicCIECloudyElectronDensity

  double precision function atomicCIECloudyElectronDensityTemperatureLogSlope(self,numberDensityHydrogen,temperature,gasAbundances,radiation)
    !!{
    Return the logarithmic slope of the electron density with respect to temperature for
    collisional ionization equilibrium as computed by {\normalfont \scshape Cloudy}.  read
    from a file.
    !!}
    implicit none
    class           (chemicalStateAtomicCIECloudy), intent(inout) :: self
    double precision                              , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances                  ), intent(in   ) :: gasAbundances
    class           (radiationFieldClass         ), intent(inout) :: radiation

    call                                              self%tabulate                                               (                                  gasAbundances          )
    atomicCIECloudyElectronDensityTemperatureLogSlope=self%chemicalStateCIEFile%electronDensityTemperatureLogSlope(numberDensityHydrogen,temperature,gasAbundances,radiation)
    return
  end function atomicCIECloudyElectronDensityTemperatureLogSlope

  double precision function atomicCIECloudyElectronDensityDensityLogSlope(self,numberDensityHydrogen,temperature,gasAbundances,radiation)
    !!{
    Return the logarithmic slope of the electron density with respect to density for
    collisional ionization equilibrium as computed by {\normalfont \scshape Cloudy}.
    !!}
    implicit none
    class           (chemicalStateAtomicCIECloudy), intent(inout) :: self
    double precision                              , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances                  ), intent(in   ) :: gasAbundances
    class           (radiationFieldClass         ), intent(inout) :: radiation

    call                                          self%tabulate                                           (                                  gasAbundances          )
    atomicCIECloudyElectronDensityDensityLogSlope=self%chemicalStateCIEFile%electronDensityDensityLogSlope(numberDensityHydrogen,temperature,gasAbundances,radiation)
    return
  end function atomicCIECloudyElectronDensityDensityLogSlope

  subroutine atomicCIECloudyChemicalDensities(self,chemicalDensities,numberDensityHydrogen,temperature,gasAbundances,radiation)
    !!{
    Return the densities of chemical species at the given temperature and hydrogen density for the specified set of abundances
    and radiation field. Units of the returned electron density are cm$^-3$.
    !!}
    implicit none
    class           (chemicalStateAtomicCIECloudy), intent(inout) :: self
    type            (chemicalAbundances          ), intent(inout) :: chemicalDensities
    double precision                              , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances                  ), intent(in   ) :: gasAbundances
    class           (radiationFieldClass         ), intent(inout) :: radiation

    call self%tabulate                              (                                                    gasAbundances          )
    call self%chemicalStateCIEFile%chemicalDensities(chemicalDensities,numberDensityHydrogen,temperature,gasAbundances,radiation)
    return
  end subroutine atomicCIECloudyChemicalDensities
