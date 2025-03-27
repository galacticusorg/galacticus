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
  Implements a cooling function class which utilizes the {\normalfont \scshape Cloudy} code to
  compute cooling in collisional ionization equilibrium.
  !!}

  !![
  <coolingFunction name="coolingFunctionAtomicCIECloudy">
   <description>
    A cooling function class that computes the cooling function using the {\normalfont \scshape Cloudy} code and under the
    assumption of collisional ionization equilibrium with no molecular contribution. Abundances are Solar, except for zero
    metallicity calculations which use {\normalfont \scshape Cloudy}'s {\normalfont \ttfamily primordial} metallicity. The helium abundance for
    non-zero metallicity is scaled between primordial and Solar values linearly with metallicity. The {\normalfont \scshape
    Cloudy} code will be downloaded and run to compute the cooling function as needed, which will then be stored for future
    use. As this process is slow, a precomputed table is provided with \glc. If metallicities outside the range tabulated in
    this file are required it will be regenerated with an appropriate range.
   </description>
  </coolingFunction>
  !!]
  type, extends(coolingFunctionCIEFile) :: coolingFunctionAtomicCIECloudy
     !!{
     A cooling function class which utilizes the {\normalfont \scshape Cloudy} code to compute
     the cooling function in collisional ionization equilibrium.
     !!}
     private
     logical :: initialized
   contains
     !![
     <methods>
       <method description="Run {\normalfont \scshape Cloudy} to tabulate the cooling function as necessary." method="tabulate" />
     </methods>
     !!]
     procedure :: tabulate                           => atomicCIECloudyTabulate
     procedure :: coolingFunction                    => atomicCIECloudyCoolingFunction
     procedure :: coolingFunctionFractionInBand      => atomicCIECloudyCoolingFunctionFractionInBand
     procedure :: coolingFunctionTemperatureLogSlope => atomicCIECloudyCoolingFunctionTemperatureLogSlope
     procedure :: coolingFunctionDensityLogSlope     => atomicCIECloudyCoolingFunctionDensityLogSlope
  end type coolingFunctionAtomicCIECloudy

  interface coolingFunctionAtomicCIECloudy
     !!{
     Constructors for the ``atomic CIE Cloudy'' cooling function class.
     !!}
     module procedure atomicCIECloudyConstructorParameters
     module procedure atomicCIECloudyConstructorInternal
  end interface coolingFunctionAtomicCIECloudy

  ! File names for the cooling function and chemical state data.
  character       (len=56)           :: coolingFunctionFileName='cooling/cooling_function_Atomic_CIE_Cloudy.hdf5'
  character       (len=56)           :: chemicalStateFileName  ='chemicalState/chemical_state_Atomic_CIE_Cloudy.hdf5'

  ! Maximum tabulated metallicity.
  double precision       , parameter :: metallicityMaximumDefault=30.0d0 ! Thirty times Solar.

  ! Maximum metallicity that we will ever tabulate. More than this should be physically implausible.
  double precision       , parameter :: metallicityMaximumLimit  =30.0d0 ! Thirty times Solar.

  ! Factor by which metallicity must exceed currently tabulated maximum before we retabulate.
  double precision       , parameter :: toleranceMetallicity     = 0.1d0

contains

  function atomicCIECloudyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``atomic CIE Cloudy'' cooling function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(coolingFunctionAtomicCIECloudy)                :: self
    type(inputParameters               ), intent(inout) :: parameters

    self=coolingFunctionAtomicCIECloudy()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function atomicCIECloudyConstructorParameters

  function atomicCIECloudyConstructorInternal() result(self)
    !!{
    Internal constructor for the ``atomic CIE Cloudy'' cooling function class.
    !!}
    implicit none
    type(coolingFunctionAtomicCIECloudy) :: self

    ! Initialize.
    self%temperaturePrevious     =-1.0d0
    self%metallicityPrevious     =-1.0d0
    self%temperatureSlopePrevious=-1.0d0
    self%metallicitySlopePrevious=-1.0d0
    self%initialized             =.false.
   return
  end function atomicCIECloudyConstructorInternal

  subroutine atomicCIECloudyTabulate(self,gasAbundances)
    !!{
    Create the cooling function.
    !!}
    use :: Abundances_Structure , only : Abundances_Get_Metallicity   , metallicityTypeLinearByMassSolar
    use :: File_Utilities       , only : File_Remove
    use :: Input_Paths          , only : inputPath                    , pathTypeDataStatic
    use :: Interfaces_Cloudy_CIE, only : Interface_Cloudy_CIE_Tabulate
    use :: String_Handling      , only : operator(//)
    implicit none
    class  (coolingFunctionAtomicCIECloudy), intent(inout) :: self
    type   (abundances                    ), intent(in   ) :: gasAbundances
    logical                                                :: makeFile

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
               &     self%metallicityMaximum*(1.0d0+toleranceMetallicity)                            &
               &   )
       else
          makeFile=.false.
       end if
       ! Remove the cooling function file so that a new one will be created.
       if (makeFile) call File_Remove(inputPath(pathTypeDataStatic)//trim(coolingFunctionFileName))
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
            &                                   inputPath(pathTypeDataStatic)//trim(coolingFunctionFileName), &
            &                                   inputPath(pathTypeDataStatic)//trim(chemicalStateFileName  ), &
            &                                   fileFormatVersionCurrent                                      &
            &                            )
       ! Call routine to read in the tabulated data.
       call self%readFile(inputPath(pathTypeDataStatic)//trim(coolingFunctionFileName))
       ! Flag that cooling function is now initialized.
       self%initialized=.true.
    end if
    return
  end subroutine atomicCIECloudyTabulate

  double precision function atomicCIECloudyCoolingFunction(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !!{
    Return the cooling function for collisional ionization equilibrium as computed by
    {\normalfont \scshape Cloudy}.
    !!}
    implicit none
    class           (coolingFunctionAtomicCIECloudy), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances                    ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances            ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass           ), intent(inout) :: radiation

    call                           self%tabulate                              (                                       gasAbundances                            )
    atomicCIECloudyCoolingFunction=self%coolingFunctionCIEFile%coolingFunction(node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    return
  end function atomicCIECloudyCoolingFunction

  double precision function atomicCIECloudyCoolingFunctionFractionInBand(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation,energyLow,energyHigh)
    !!{
    Return the fraction of the cooling luminosity due to emission in the given energy range as computed by
    {\normalfont \scshape Cloudy}.
    !!}
    implicit none
    class           (coolingFunctionAtomicCIECloudy), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: numberDensityHydrogen, temperature, &
         &                                                             energyLow            , energyHigh
    type            (abundances                    ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances            ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass           ), intent(inout) :: radiation

    call                                         self%tabulate                                            (                                       gasAbundances                                                 )
    atomicCIECloudyCoolingFunctionFractionInBand=self%coolingFunctionCIEFile%coolingFunctionFractionInBand(node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation,energyLow,energyHigh)
    return
  end function atomicCIECloudyCoolingFunctionFractionInBand

  double precision function atomicCIECloudyCoolingFunctionTemperatureLogSlope(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !!{
    Return the logarithmic slope of the cooling function with respect to temperature for
    collisional ionization equilibrium as computed by {\normalfont \scshape Cloudy}.  read
    from a file.
    !!}
    implicit none
    class           (coolingFunctionAtomicCIECloudy), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances                    ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances            ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass           ), intent(inout) :: radiation

    call                                              self%tabulate                                                 (                                       gasAbundances                            )
    atomicCIECloudyCoolingFunctionTemperatureLogSlope=self%coolingFunctionCIEFile%coolingFunctionTemperatureLogSlope(node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    return
  end function atomicCIECloudyCoolingFunctionTemperatureLogSlope

  double precision function atomicCIECloudyCoolingFunctionDensityLogSlope(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !!{
    Return the logarithmic slope of the cooling function with respect to density for
    collisional ionization equilibrium as computed by {\normalfont \scshape Cloudy}.
    !!}
    implicit none
    class           (coolingFunctionAtomicCIECloudy), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances                    ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances            ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass           ), intent(inout) :: radiation

    call                                          self%tabulate                                             (                                       gasAbundances                            )
    atomicCIECloudyCoolingFunctionDensityLogSlope=self%coolingFunctionCIEFile%coolingFunctionDensityLogSlope(node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    return
  end function atomicCIECloudyCoolingFunctionDensityLogSlope
