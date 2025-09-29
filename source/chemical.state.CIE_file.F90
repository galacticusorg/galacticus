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
  Implements a chemical state class which reads and interpolates a collisional ionization equilibrium chemical state from a file.
  !!}

  use :: Numerical_Interpolation, only : interpolator
  use :: Table_Labels           , only : enumerationExtrapolationTypeType

  !![
  <chemicalState name="chemicalStateCIEFile">
   <description>
    A chemical state class providing chemical state via interpolation of tabulated values read from file. The HDF5 file
    containing the table should have the following form:
    \begin{verbatim}
    HDF5 "chemicalState.hdf5" {
    GROUP "/" {
       ATTRIBUTE "fileFormat" {
          DATATYPE  H5T_STRING {
             STRSIZE 1;
             STRPAD H5T_STR_NULLTERM;
             CSET H5T_CSET_ASCII;
             CTYPE H5T_C_S1;
          }
          DATASPACE  SCALAR
          DATA {
          (0): "1"
          }
       }
       DATASET "electronDensity" {
          DATATYPE  H5T_IEEE_F64LE
          DATASPACE  SIMPLE { ( 7, 8 ) / ( 7, 8 ) }
       }
       DATASET "hiDensity" {
          DATATYPE  H5T_IEEE_F64LE
          DATASPACE  SIMPLE { ( 7, 8 ) / ( 7, 8 ) }
       }
       DATASET "hiiDensity" {
          DATATYPE  H5T_IEEE_F64LE
          DATASPACE  SIMPLE { ( 7, 8 ) / ( 7, 8 ) }
       }
       DATASET "metallicity" {
          DATATYPE  H5T_IEEE_F64LE
          DATASPACE  SIMPLE { ( 8 ) / ( 8 ) }
          ATTRIBUTE "extrapolateHigh" {
             DATATYPE  H5T_STRING {
                STRSIZE 3;
                STRPAD H5T_STR_NULLTERM;
                CSET H5T_CSET_ASCII;
                CTYPE H5T_C_S1;
             }
             DATASPACE  SCALAR
             DATA {
             (0): "fix"
             }
          }
          ATTRIBUTE "extrapolateLow" {
             DATATYPE  H5T_STRING {
                STRSIZE 3;
                STRPAD H5T_STR_NULLTERM;
                CSET H5T_CSET_ASCII;
                CTYPE H5T_C_S1;
             }
             DATASPACE  SCALAR
             DATA {
             (0): "fix"
             }
          }
       }
       DATASET "temperature" {
          DATATYPE  H5T_IEEE_F64LE
          DATASPACE  SIMPLE { ( 7 ) / ( 7 ) }
          ATTRIBUTE "extrapolateHigh" {
             DATATYPE  H5T_STRING {
                STRSIZE 3;
                STRPAD H5T_STR_NULLTERM;
                CSET H5T_CSET_ASCII;
                CTYPE H5T_C_S1;
             }
             DATASPACE  SCALAR
             DATA {
             (0): "fix"
             }
          }
          ATTRIBUTE "extrapolateLow" {
             DATATYPE  H5T_STRING {
                STRSIZE 3;
                STRPAD H5T_STR_NULLTERM;
                CSET H5T_CSET_ASCII;
                CTYPE H5T_C_S1;
             }
             DATASPACE  SCALAR
             DATA {
             (0): "fix"
             }
          }
       }
    }
    }
    \end{verbatim}
    The {\normalfont \ttfamily temperature} dataset should specify temperature (in Kelvin), while the {\normalfont \ttfamily
    metallicity} dataset should give the logarithmic metallicity relative to Solar (a value of -999 or less is taken to imply
    zero metallicity). The {\normalfont \ttfamily electronDensity} dataset should specify the number density of electrons
    relative to hydrogen at each temperature/metallicity pair. Optionally {\normalfont \ttfamily hiDensity} and {\normalfont
    \ttfamily hiiDensity} datasets may be added giving the number densities of H{\normalfont \scshape i} and H{\normalfont
    \scshape ii} relative to hydrogen respectively The {\normalfont \ttfamily extrapolateLow} and {\normalfont \ttfamily
    extrapolateHigh} attributes of the {\normalfont \ttfamily temperature} and {\normalfont \ttfamily metallicity} datasets
    specify how the cooling rate should be extrapolated in the low and high vale limits. Allowed options for these attributes
    are:
    \begin{description}
     \item[{\normalfont \ttfamily zero}] The electron density is set to zero beyond the relevant limit.
     \item[{\normalfont \ttfamily fixed}] The electron density is held fixed at the value at the relevant limit.
     \item[{\normalfont \ttfamily power law}] The electron density is extrapolated assuming a
     power-law dependence beyond the relevant limit. This option is only allowed if the
     electron density is everywhere positive.
    \end{description}
    If the electron density is everywhere positive the interpolation will be done in the
    logarithmic of temperature, metallicity\footnote{The exception is if the first electron
    density is tabulated for zero metallicity. In that case, a linear interpolation in
    metallicity is always used between zero and the first non-zero tabulated metallicity.}
    and electron density. Otherwise, interpolation is linear in these quantities. The
    electron density is scaled assuming a linear dependence on hydrogen density.
   </description>
   <runTimeFileDependencies paths="fileName"/>
  </chemicalState>
  !!]
  type, extends(chemicalStateClass) :: chemicalStateCIEFile
     !!{
     A chemical state class which interpolates state given in a file assuming collisional ionization equilibrium.
     !!}
     private
     type            (varying_string                  )                              :: fileName
     double precision                                                                :: metallicityMaximum               , metallicityMinimum         , &
          &                                                                             temperatureMaximum               , temperatureMinimum
     type            (enumerationExtrapolationTypeType)                              :: extrapolateMetallicityHigh       , extrapolateMetallicityLow  , &
          &                                                                             extrapolateTemperatureHigh       , extrapolateTemperatureLow
     logical                                                                         :: firstMetallicityIsZero           , gotHydrogenAtomic          , &
          &                                                                             gotHydrogenCation                , logarithmicTable
     integer                                                                         :: metallicityCount                 , temperatureCount
     double precision                                                                :: firstNonZeroMetallicity          , electronDensityPrevious    , &
          &                                                                             metallicityPrevious              , temperaturePrevious        , &
          &                                                                             electronDensitySlopePrevious     , metallicitySlopePrevious   , &
          &                                                                             temperatureSlopePrevious         , metallicityChemicalPrevious, &
          &                                                                             temperatureChemicalPrevious      , hMetallicityPrevious       , &
          &                                                                             hTemperaturePrevious
     integer         (c_size_t                        )                              :: iMetallicityPrevious             , iTemperaturePrevious
     type            (chemicalAbundances              )                              :: chemicalDensitiesPrevious
     double precision                                  , allocatable, dimension(:  ) :: metallicities                    , temperatures
     double precision                                  , allocatable, dimension(:,:) :: densityElectron                  , densityHydrogenAtomic      , &
          &                                                                             densityHydrogenCation
     type            (interpolator                    )                              :: interpolatorMetallicity          , interpolatorTemperature
     integer                                                                         :: atomicHydrogenCationChemicalIndex, atomicHydrogenChemicalIndex, &
          &                                                                             electronChemicalIndex
   contains
     !![
     <methods>
       <method description="Read the named chemical state file."                         method="readFile"            />
       <method description="Compute interpolating factors in a CIE chemical state file." method="interpolatingFactors"/>
       <method description="Interpolate in the given density table."                     method="interpolate"         />
     </methods>
     !!]
     procedure :: readFile                           => cieFileReadFile
     procedure :: interpolatingFactors               => cieFileInterpolatingFactors
     procedure :: interpolate                        => cieFileInterpolate
     procedure :: electronDensity                    => cieFileElectronDensity
     procedure :: electronDensityTemperatureLogSlope => cieFileElectronDensityTemperatureLogSlope
     procedure :: electronDensityDensityLogSlope     => cieFileElectronDensityDensityLogSlope
     procedure :: chemicalDensities                  => cieFileChemicalDensities
  end type chemicalStateCIEFile

  interface chemicalStateCIEFile
     !!{
     Constructors for the \refClass{chemicalStateCIEFile} chemical state class.
     !!}
     module procedure cieFileConstructorParameters
     module procedure cieFileConstructorInternal
  end interface chemicalStateCIEFile

  ! Current file format version for CIE chemical state files.
  integer, parameter :: fileFormatVersionCurrent=1

contains

  function cieFileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{chemicalStateCIEFile} chemical state class which takes a parameter set as input.
    !!}
    implicit none
    type(chemicalStateCIEFile)                :: self
    type(inputParameters     ), intent(inout) :: parameters
    type(varying_string      )                :: fileName

    !![
    <inputParameter>
      <name>fileName</name>
      <source>parameters</source>
      <description>The name of the file containing a tabulation of the collisional ionization equilibrium chemical state.</description>
    </inputParameter>
    !!]
    ! Construct the instance.
    self=cieFileConstructorInternal(char(fileName))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function cieFileConstructorParameters

  function cieFileConstructorInternal(fileName) result(self)
    !!{
    Internal constructor for the \refClass{chemicalStateCIEFile} chemical state class.
    !!}
    use :: Chemical_Abundances_Structure, only : unitChemicalAbundances
    implicit none
    type     (chemicalStateCIEFile)                :: self
    character(len=*               ), intent(in   ) :: fileName

    ! Read the file.
    self%fileName=fileName
    call self%readFile(fileName)
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
    return
  end function cieFileConstructorInternal

  double precision function cieFileElectronDensity(self,numberDensityHydrogen,temperature,gasAbundances,radiation)
    !!{
    Return the electron density by interpolating in tabulated CIE data read from a file.
    !!}
    use            :: Abundances_Structure, only : Abundances_Get_Metallicity, abundances                  , metallicityTypeLinearByMassSolar
    use, intrinsic :: ISO_C_Binding       , only : c_size_t
    use            :: Radiation_Fields    , only : radiationFieldClass
    use            :: Table_Labels        , only : extrapolationTypeFix      , extrapolationTypeExtrapolate, extrapolationTypeZero
    implicit none
    class           (chemicalStateCIEFile), intent(inout) :: self
    double precision                      , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances          ), intent(in   ) :: gasAbundances
    class           (radiationFieldClass ), intent(inout) :: radiation
    integer         (c_size_t            )                :: iMetallicity         , iTemperature
    double precision                                      :: hMetallicity         , hTemperature  , &
         &                                                   metallicityUse       , temperatureUse
    !$GLC attributes unused :: radiation

    ! Handle out of range temperatures.
    temperatureUse=temperature
    if (temperatureUse < self%temperatureMinimum) then
       select case (self%extrapolateTemperatureLow%ID)
       case (extrapolationTypeZero%ID)
          cieFileElectronDensity=0.0d0
          return
       case (extrapolationTypeFix%ID,extrapolationTypeExtrapolate%ID)
          temperatureUse=self%temperatureMinimum
       end select
    end if
    if (temperatureUse > self%temperatureMaximum) then
       select case (self%extrapolateTemperatureHigh%ID)
       case (extrapolationTypeZero%ID)
          cieFileElectronDensity=0.0d0
          return
       case (extrapolationTypeFix%ID,extrapolationTypeExtrapolate%ID)
          temperatureUse=self%temperatureMaximum
       end select
    end if
    ! Handle out of range metallicities.
    metallicityUse=Abundances_Get_Metallicity(gasAbundances,metallicityType=metallicityTypeLinearByMassSolar)
    if (metallicityUse < self%metallicityMinimum) then
       select case (self%extrapolateMetallicityLow%ID)
       case (extrapolationTypeZero%ID)
          cieFileElectronDensity=0.0d0
          return
       case (extrapolationTypeFix%ID)
          metallicityUse=self%metallicityMinimum
       end select
    end if
    if (metallicityUse > self%metallicityMaximum) then
       select case (self%extrapolateMetallicityHigh%ID)
       case (extrapolationTypeZero%ID)
          cieFileElectronDensity=0.0d0
          return
       case (extrapolationTypeFix%ID)
          metallicityUse=self%metallicityMaximum
       end select
    end if
    ! Check if we need to recompute the chemical state.
    if     (                                            &
         &   temperatureUse /= self%temperaturePrevious &
         &  .or.                                        &
         &   metallicityUse /= self%metallicityPrevious &
         & ) then
       ! Get the interpolation.
       call self%interpolatingFactors(temperatureUse,metallicityUse,iTemperature,hTemperature,iMetallicity,hMetallicity)
       ! Do the interpolation.
       self%electronDensityPrevious=self%interpolate(iTemperature,hTemperature,iMetallicity,hMetallicity,self%densityElectron)
       ! Store the temperature and metallicity for which calculation was performed.
       self%temperaturePrevious=temperatureUse
       self%metallicityPrevious=metallicityUse
    end if
    ! Scale to the specified density assuming two-body processes, in which case electron density scales with hydrogen density.
    cieFileElectronDensity=self%electronDensityPrevious*numberDensityHydrogen
    return
  end function cieFileElectronDensity

  double precision function cieFileElectronDensityTemperatureLogSlope(self,numberDensityHydrogen,temperature,gasAbundances,radiation)
    !!{
    Return the logarithmic slope of the electron density with respect to temperature by interpolating in tabulated CIE data
    read from a file.
    !!}
    use            :: Abundances_Structure, only : Abundances_Get_Metallicity, abundances                  , metallicityTypeLinearByMassSolar
    use, intrinsic :: ISO_C_Binding       , only : c_size_t
    use            :: Radiation_Fields    , only : radiationFieldClass
    use            :: Table_Labels        , only : extrapolationTypeFix      , extrapolationTypeExtrapolate, extrapolationTypeZero
    implicit none
    class           (chemicalStateCIEFile), intent(inout) :: self
    double precision                      , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances          ), intent(in   ) :: gasAbundances
    class           (radiationFieldClass ), intent(inout) :: radiation
    integer         (c_size_t            )                :: iMetallicity         , iTemperature
    double precision                                      :: hMetallicity         , hTemperature  , &
         &                                                   metallicityUse       , temperatureUse
    !$GLC attributes unused :: radiation, numberDensityHydrogen

    ! Handle out of range temperatures.
    temperatureUse=temperature
    if (temperatureUse < self%temperatureMinimum) then
       select case (self%extrapolateTemperatureLow%ID)
       case (extrapolationTypeZero%ID,extrapolationTypeFix%ID)
          cieFileElectronDensityTemperatureLogSlope=0.0d0
          return
       case (extrapolationTypeExtrapolate%ID)
          temperatureUse=self%temperatureMinimum
       end select
    end if
    if (temperatureUse > self%temperatureMaximum) then
       select case (self%extrapolateTemperatureHigh%ID)
       case (extrapolationTypeZero%ID,extrapolationTypeFix%ID)
          cieFileElectronDensityTemperatureLogSlope=0.0d0
          return
       case (extrapolationTypeExtrapolate%ID)
          temperatureUse=self%temperatureMaximum
       end select
    end if
    ! Handle out of range metallicities.
    metallicityUse=Abundances_Get_Metallicity(gasAbundances,metallicityType=metallicityTypeLinearByMassSolar)
    if (metallicityUse < self%metallicityMinimum) then
       select case (self%extrapolateMetallicityLow%ID)
       case (extrapolationTypeZero%ID)
          cieFileElectronDensityTemperatureLogSlope=0.0d0
          return
       case (extrapolationTypeFix%ID)
          metallicityUse=self%metallicityMinimum
       end select
    end if
    if (metallicityUse > self%metallicityMaximum) then
       select case (self%extrapolateMetallicityHigh%ID)
       case (extrapolationTypeZero%ID)
          cieFileElectronDensityTemperatureLogSlope=0.0d0
          return
       case (extrapolationTypeFix%ID)
          metallicityUse=self%metallicityMaximum
       end select
    end if
    ! Check if we need to recompute the chemical state.
    if     (                                                 &
         &   temperatureUse /= self%temperatureSlopePrevious &
         &  .or.                                             &
         &   metallicityUse /= self%metallicitySlopePrevious &
         & ) then
       ! Get the interpolation.
       call self%interpolatingFactors(temperatureUse,metallicityUse,iTemperature,hTemperature,iMetallicity,hMetallicity)
       ! Do the interpolation.
       self%electronDensitySlopePrevious=+(                                                       &
            &                              +(                                                     &
            &                                +self%densityElectron(iMetallicity  ,iTemperature+1) &
            &                                -self%densityElectron(iMetallicity  ,iTemperature  ) &
            &                               )                                                     &
            &                              *(1.0d0-hMetallicity)                                  &
            &                              +(                                                     &
            &                                +self%densityElectron(iMetallicity+1,iTemperature+1) &
            &                                -self%densityElectron(iMetallicity+1,iTemperature  ) &
            &                               )                                                     &
            &                              *(     +hMetallicity)                                  &
            &                             )                                                       &
            &                            /(                                                       &
            &                              +  self%temperatures   (iTemperature+1               ) &
            &                              -  self%temperatures   (iTemperature                 ) &
            &                             )
       ! Convert to logarithmic gradient if table was not stored logarithmically.
       if (.not.self%logarithmicTable)                                                                                      &
            & self  %electronDensitySlopePrevious=                                                                          &
            &  +self%electronDensitySlopePrevious                                                                           &
            &  /self%interpolate                 (iTemperature,hTemperature,iMetallicity,hMetallicity,self%densityElectron) &
            &  *temperature
       ! Store the temperature and metallicity for which calculation was performed.
       self%temperatureSlopePrevious=temperatureUse
       self%metallicitySlopePrevious=metallicityUse
    end if
    ! Return the stored value.
    cieFileElectronDensityTemperatureLogSlope=self%electronDensitySlopePrevious
    return
  end function cieFileElectronDensityTemperatureLogSlope

  double precision function cieFileElectronDensityDensityLogSlope(self,numberDensityHydrogen,temperature,gasAbundances,radiation)
    !!{
    Return the logarithmic slope of the electron density with respect to density assuming atomic CIE.
    !!}
    use :: Abundances_Structure, only : abundances
    use :: Radiation_Fields    , only : radiationFieldClass
    implicit none
    class           (chemicalStateCIEFile), intent(inout) :: self
    double precision                      , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances          ), intent(in   ) :: gasAbundances
    class           (radiationFieldClass ), intent(inout) :: radiation
    !$GLC attributes unused :: self, numberDensityHydrogen, temperature, gasAbundances, radiation

    ! Electron density always scales as total density under CIE conditions.
    cieFileElectronDensityDensityLogSlope=1.0d0
    return
  end function cieFileElectronDensityDensityLogSlope

  subroutine cieFileChemicalDensities(self,chemicalDensities,numberDensityHydrogen,temperature,gasAbundances,radiation)
    !!{
    Return the densities of chemical species at the given temperature and hydrogen density for the specified set of abundances
    and radiation field. Units of the returned electron density are cm$^-3$.
    !!}
    use            :: Abundances_Structure         , only : Abundances_Get_Metallicity, abundances                  , metallicityTypeLinearByMassSolar
    use            :: Chemical_Abundances_Structure, only : chemicalAbundances
    use, intrinsic :: ISO_C_Binding                , only : c_size_t
    use            :: Radiation_Fields             , only : radiationFieldClass
    use            :: Table_Labels                 , only : extrapolationTypeFix      , extrapolationTypeExtrapolate, extrapolationTypeZero
    implicit none
    class           (chemicalStateCIEFile), intent(inout) :: self
    type            (chemicalAbundances  ), intent(inout) :: chemicalDensities
    double precision                      , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances          ), intent(in   ) :: gasAbundances
    class           (radiationFieldClass ), intent(inout) :: radiation
    integer         (c_size_t            )                :: iMetallicity         , iTemperature
    double precision                                      :: hMetallicity         , hTemperature  , &
         &                                                   metallicityUse       , temperatureUse
    !$GLC attributes unused :: radiation

    ! Handle out of range temperatures.
    temperatureUse=temperature
    if (temperatureUse < self%temperatureMinimum) then
       select case (self%extrapolateTemperatureLow%ID)
       case (extrapolationTypeZero%ID)
          call chemicalDensities%reset()
          return
       case (extrapolationTypeFix%ID,extrapolationTypeExtrapolate%ID)
          temperatureUse=self%temperatureMinimum
       end select
    end if
    if (temperatureUse > self%temperatureMaximum) then
       select case (self%extrapolateTemperatureHigh%ID)
       case (extrapolationTypeZero%ID)
          call chemicalDensities%reset()
          return
       case (extrapolationTypeFix%ID,extrapolationTypeExtrapolate%ID)
          temperatureUse=self%temperatureMaximum
       end select
    end if
    ! Handle out of range metallicities.
    metallicityUse=Abundances_Get_Metallicity(gasAbundances,metallicityType=metallicityTypeLinearByMassSolar)
    if (metallicityUse < self%metallicityMinimum) then
       select case (self%extrapolateMetallicityLow%ID)
       case (extrapolationTypeZero%ID)
          call chemicalDensities%reset()
          return
       case (extrapolationTypeFix%ID)
          metallicityUse=self%metallicityMinimum
       end select
    end if
    if (metallicityUse > self%metallicityMaximum) then
       select case (self%extrapolateMetallicityHigh%ID)
       case (extrapolationTypeZero%ID)
          call chemicalDensities%reset()
          return
       case (extrapolationTypeFix%ID)
          metallicityUse=self%metallicityMaximum
       end select
    end if
    ! Check if we need to recompute the chemical state.
    if   (                                                    &
       &   temperatureUse /= self%temperatureChemicalPrevious &
       &  .or.                                                &
       &   metallicityUse /= self%metallicityChemicalPrevious &
       & ) then
       ! Get the interpolation.
       call self%interpolatingFactors(temperatureUse,metallicityUse,iTemperature,hTemperature,iMetallicity,hMetallicity)
       ! Reset densities to zero.
       call self%chemicalDensitiesPrevious%reset()
       ! Do the interpolation.
       if (self%electronChemicalIndex             > 0)                                                                                                    &
            & call self%                                                                                                                                  &
            &       chemicalDensitiesPrevious%                                                                                                            &
            &        abundanceSet(                                                                                                                        &
            &                     self%electronChemicalIndex                                                                                            , &
            &                     self%interpolate                      (iTemperature,hTemperature,iMetallicity,hMetallicity,self%densityElectron      )  &
            &                    )
       if (self%atomicHydrogenChemicalIndex       > 0)                                                                                                    &
            & call self%                                                                                                                                  &
            &       chemicalDensitiesPrevious%                                                                                                            &
            &        abundanceSet(                                                                                                                        &
            &                     self%atomicHydrogenChemicalIndex                                                                                      , &
            &                     self%interpolate                      (iTemperature,hTemperature,iMetallicity,hMetallicity,self%densityHydrogenAtomic)  &
            &                    )
       if (self%atomicHydrogenCationChemicalIndex > 0)                                                                                                    &
            & call self%                                                                                                                                  &
            &       chemicalDensitiesPrevious%                                                                                                            &
            &        abundanceSet(                                                                                                                        &
            &                     self%atomicHydrogenCationChemicalIndex                                                                                , &
            &                     self%interpolate                      (iTemperature,hTemperature,iMetallicity,hMetallicity,self%densityHydrogenCation)  &
            &                    )
       ! Store the temperature and metallicity for which calculation was performed.
       self%temperatureChemicalPrevious=temperatureUse
       self%metallicityChemicalPrevious=metallicityUse
    end if
    ! Scale to the specified density assuming two-body processes, in which case densities scale with hydrogen density.
    chemicalDensities=self%chemicalDensitiesPrevious*numberDensityHydrogen
    return
  end subroutine cieFileChemicalDensities

  subroutine cieFileInterpolatingFactors(self,temperature,metallicity,iTemperature,hTemperature,iMetallicity,hMetallicity)
    !!{
    Determine the interpolating parameters.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (chemicalStateCIEFile), intent(inout) :: self
    double precision                      , intent(in   ) :: metallicity   , temperature
    integer         (c_size_t            ), intent(  out) :: iMetallicity  , iTemperature
    double precision                      , intent(  out) :: hMetallicity  , hTemperature
    double precision                                      :: metallicityUse, temperatureUse

    ! Get interpolation in temperature.
    if (temperature /= self%temperaturePrevious) then
       temperatureUse=temperature
       if (self%logarithmicTable) temperatureUse=log(temperatureUse)
       self%iTemperaturePrevious=max(                                                         &
            &                        min(                                                     &
            &                            self%interpolatorTemperature%locate(temperatureUse), &
            &                            self%temperatureCount-1                              &
            &                           )                                                   , &
            &                            1                                                    &
            &                       )
       self%hTemperaturePrevious=+(     temperatureUse                             -self%temperatures(self%iTemperaturePrevious)) &
            &                    /(self%temperatures  (self%iTemperaturePrevious+1)-self%temperatures(self%iTemperaturePrevious))
    end if
    iTemperature=self%iTemperaturePrevious
    hTemperature=self%hTemperaturePrevious
    ! Get interpolation in metallicity.
    if (metallicity /= self%metallicityPrevious) then
       metallicityUse=max(metallicity,0.0d0)
       if     (                                               &
            &                    self%firstMetallicityIsZero  &
            &  .and.                                          &
            &   metallicityUse < self%firstNonZeroMetallicity &
            & ) then
          self%iMetallicityPrevious=+1
          self%hMetallicityPrevious=+metallicityUse               &
               &                    /self%firstNonZeroMetallicity
       else
          if (self%logarithmicTable) metallicityUse=log(metallicityUse)
          self%iMetallicityPrevious=max(                                                         &
               &                        min(                                                     &
               &                            self%interpolatorMetallicity%locate(metallicityUse), &
               &                            self%metallicityCount-1                              &
               &                           )                                                   , &
               &                            1                                                    &
               &                       )
          self%hMetallicityPrevious=+(     metallicityUse                             -self%metallicities(self%iMetallicityPrevious)) &
               &                    /(self%metallicities (self%iMetallicityPrevious+1)-self%metallicities(self%iMetallicityPrevious))
       end if
    end if
    iMetallicity=self%iMetallicityPrevious
    hMetallicity=self%hMetallicityPrevious
    return
  end subroutine cieFileInterpolatingFactors

  double precision function cieFileInterpolate(self,iTemperature,hTemperature,iMetallicity,hMetallicity,density)
    !!{
    Perform the interpolation.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (chemicalStateCIEFile)                , intent(inout) :: self
    integer         (c_size_t            )                , intent(in   ) :: iMetallicity, iTemperature
    double precision                                      , intent(in   ) :: hMetallicity, hTemperature
    double precision                      , dimension(:,:), intent(in   ) :: density

    ! Do the interpolation.
    cieFileInterpolate=+density(iMetallicity  ,iTemperature  )*(1.0d0-hTemperature)*(1.0d0-hMetallicity) &
         &             +density(iMetallicity+1,iTemperature  )*(1.0d0-hTemperature)*(      hMetallicity) &
         &             +density(iMetallicity  ,iTemperature+1)*(      hTemperature)*(1.0d0-hMetallicity) &
         &             +density(iMetallicity+1,iTemperature+1)*(      hTemperature)*(      hMetallicity)
    ! Exponentiate the result if the table was stored as the log.
    if (self%logarithmicTable) cieFileInterpolate=exp(cieFileInterpolate)
    return
  end function cieFileInterpolate

  subroutine cieFileReadFile(self,fileName)
    !!{
    Read in data from an chemical state file.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Display                      , only : displayIndent                       , displayUnindent     , verbosityLevelDebug         , displayGreen         , &
         &                                       displayReset
    use :: Error                        , only : Error_Report                        , errorStatusSuccess
    use :: HDF5_Access                  , only : hdf5Access
    use :: IO_HDF5                      , only : hdf5Object
    use :: ISO_Varying_String           , only : varying_string
    use :: Table_Labels                 , only : enumerationExtrapolationTypeEncode  , extrapolationTypeFix, extrapolationTypeExtrapolate, extrapolationTypeZero, &
         &                                       enumerationExtrapolationTypeDescribe
    implicit none
    class           (chemicalStateCIEFile), intent(inout) :: self
    character       (len=*               ), intent(in   ) :: fileName
    double precision                      , parameter     :: metallicityLogarithmicZero=-999.0d0
    type            (varying_string      )                :: limitType
    integer                                               :: fileFormatVersion                  , status
    type            (hdf5Object          )                :: chemicalStateFile                  , metallicityDataset, &
         &                                                   temperatureDataset
    
    !$ call hdf5Access%set()
    ! Parse the file.
    call displayIndent('Reading file: '//fileName,verbosityLevelDebug)
    chemicalStateFile=hdf5Object(fileName,readOnly=.true.)
    ! Check the file format version of the file.
    call chemicalStateFile%readAttribute('fileFormat',fileFormatVersion)
    if (fileFormatVersion /= fileFormatVersionCurrent) call Error_Report('file format version is out of date'//{introspection:location})
    ! Test for presence of hydrogen data.
    self%gotHydrogenAtomic=chemicalStateFile%hasDataset('hiDensity' )
    self%gotHydrogenCation=chemicalStateFile%hasDataset('hiiDensity')
    ! Read datasets.
    call                             chemicalStateFile%readDataset('temperature'    ,self%temperatures         )
    call                             chemicalStateFile%readDataset('metallicity'    ,self%metallicities        )
    call                             chemicalStateFile%readDataset('electronDensity',self%densityElectron      )
    if (self%gotHydrogenAtomic) call chemicalStateFile%readDataset('hiDensity'      ,self%densityHydrogenAtomic)
    if (self%gotHydrogenCation) call chemicalStateFile%readDataset('hiiDensity'     ,self%densityHydrogenCation)
    self%metallicityCount=size(self%metallicities)
    self%temperatureCount=size(self%temperatures )
    ! Unlog metallicities.
    where (self%metallicities > metallicityLogarithmicZero)
       self%metallicities=10.0d0**self%metallicities
    elsewhere
       self%metallicities= 0.0d0
    end where
    ! Extract extrapolation methods from the file.
    metallicityDataset=chemicalStateFile%openDataset('metallicity')
    call metallicityDataset%readAttribute('extrapolateLow' ,limitType,allowPseudoScalar=.true.)
    self%extrapolateMetallicityLow =enumerationExtrapolationTypeEncode(char(limitType),includesPrefix=.false.,status=status)
    if (status /= errorStatusSuccess) call Error_Report("low metallicity extrapolation type '" //char(limitType)//"' in file '"//trim(fileName)//"' is invalid"//char(10)//displayGreen()//"HELP:"//displayReset()//enumerationExtrapolationTypeDescribe()//{introspection:location})
    call metallicityDataset%readAttribute('extrapolateHigh',limitType,allowPseudoScalar=.true.)
    self%extrapolateMetallicityHigh=enumerationExtrapolationTypeEncode(char(limitType),includesPrefix=.false.,status=status)
    if (status /= errorStatusSuccess) call Error_Report("high metallicity extrapolation type '"//char(limitType)//"' in file '"//trim(fileName)//"' is invalid"//char(10)//displayGreen()//"HELP:"//displayReset()//enumerationExtrapolationTypeDescribe()//{introspection:location})
    temperatureDataset=chemicalStateFile%openDataset('temperature')
    call temperatureDataset%readAttribute('extrapolateLow' ,limitType,allowPseudoScalar=.true.)
    self%extrapolateTemperatureLow =enumerationExtrapolationTypeEncode(char(limitType),includesPrefix=.false.,status=status)
    if (status /= errorStatusSuccess) call Error_Report("low temperature extrapolation type '" //char(limitType)//"' in file '"//trim(fileName)//"' is invalid"//char(10)//displayGreen()//"HELP:"//displayReset()//enumerationExtrapolationTypeDescribe()//{introspection:location})
    call temperatureDataset%readAttribute('extrapolateHigh',limitType,allowPseudoScalar=.true.)
    self%extrapolateTemperatureHigh=enumerationExtrapolationTypeEncode(char(limitType),includesPrefix=.false.,status=status)
    if (status /= errorStatusSuccess) call Error_Report("high temperature extrapolation type '"//char(limitType)//"' in file '"//trim(fileName)//"' is invalid"//char(10)//displayGreen()//"HELP:"//displayReset()//enumerationExtrapolationTypeDescribe()//{introspection:location})
    ! Validate extrapolation methods.
    if     (                                                                 &
         &   self%extrapolateMetallicityLow  /= extrapolationTypeFix         &
         &  .and.                                                            &
         &   self%extrapolateMetallicityLow  /= extrapolationTypeZero        &
         &  .and.                                                            &
         &   self%extrapolateMetallicityLow  /= extrapolationTypeExtrapolate &
         & ) call Error_Report('extrapolation type not permitted'//{introspection:location})
    if     (                                                                 &
         &   self%extrapolateMetallicityHigh /= extrapolationTypeFix         &
         &  .and.                                                            &
         &   self%extrapolateMetallicityHigh /= extrapolationTypeZero        &
         &  .and.                                                            &
         &   self%extrapolateMetallicityHigh /= extrapolationTypeExtrapolate &
         & ) call Error_Report('extrapolation type not permitted'//{introspection:location})
    if     (                                                                 &
         &   self%extrapolateTemperatureLow  /= extrapolationTypeFix         &
         &  .and.                                                            &
         &   self%extrapolateTemperatureLow  /= extrapolationTypeZero        &
         &  .and.                                                            &
         &   self%extrapolateTemperatureLow  /= extrapolationTypeExtrapolate &
         & ) call Error_Report('extrapolation type not permitted'//{introspection:location})
    if     (                                                                 &
         &   self%extrapolateTemperatureHigh /= extrapolationTypeFix         &
         &  .and.                                                            &
         &   self%extrapolateTemperatureHigh /= extrapolationTypeZero        &
         &  .and.                                                            &
         &   self%extrapolateTemperatureHigh /= extrapolationTypeExtrapolate &
         & ) call Error_Report('extrapolation type not permitted'//{introspection:location})
    call displayUnindent('done',verbosityLevelDebug)
    !$ call hdf5Access%unset()
    ! Store table ranges for convenience.
    self%metallicityMinimum=self%metallicities(                    1)
    self%metallicityMaximum=self%metallicities(self%metallicityCount)
    self%temperatureMinimum=self%temperatures (                    1)
    self%temperatureMaximum=self%temperatures (self%temperatureCount)
    ! Decide whether or not to make the tables logarithmic.
    self%logarithmicTable=  all(self%densityElectron       > 0.0d0)                                   &
         &                .and.                                                                       &
         &                 (all(self%densityHydrogenAtomic > 0.0d0) .or. .not.self%gotHydrogenAtomic) &
         &                .and.                                                                       &
         &                 (all(self%densityHydrogenCation > 0.0d0) .or. .not.self%gotHydrogenCation)
    if (self%logarithmicTable) then
       self%firstMetallicityIsZero=(self%metallicities(1) == 0.0d0)
       if (self%firstMetallicityIsZero) self%firstNonZeroMetallicity=self%metallicities(2)
       where (self%metallicities > 0.0d0)
          self%metallicities=log(self%metallicities)
       elsewhere
          self%metallicities=metallicityLogarithmicZero
       end where
       self                            %         temperatures=log(self%         temperatures)
       self                            %      densityElectron=log(self%      densityElectron)
       if (self%gotHydrogenAtomic) self%densityHydrogenAtomic=log(self%densityHydrogenAtomic)
       if (self%gotHydrogenCation) self%densityHydrogenCation=log(self%densityHydrogenCation)
    else
       if     (                                                                 &
            &   self%extrapolateTemperatureLow  == extrapolationTypeExtrapolate &
            &  .or.                                                             &
            &   self%extrapolateTemperatureHigh == extrapolationTypeExtrapolate &
            &  .or.                                                             &
            &   self%extrapolateMetallicityLow  == extrapolationTypeExtrapolate &
            &  .or.                                                             &
            &   self%extrapolateMetallicityHigh == extrapolationTypeExtrapolate &
            & )                                                                 &
            & call Error_Report('extrapolation allowed only in loggable tables'//{introspection:location})
    end if
    ! Build interpolators.
    self%interpolatorTemperature=interpolator(self%temperatures )
    self%interpolatorMetallicity=interpolator(self%metallicities)
    ! Get chemical indices.
    self%electronChemicalIndex               =Chemicals_Index("Electron"            ,status)
    if (self%gotHydrogenAtomic) then
       self%atomicHydrogenChemicalIndex      =Chemicals_Index("AtomicHydrogen"      ,status)
    else
       self%atomicHydrogenChemicalIndex      =-1
    end if
    if (self%gotHydrogenCation) then
       self%atomicHydrogenCationChemicalIndex=Chemicals_Index("AtomicHydrogenCation",status)
    else
       self%atomicHydrogenCationChemicalIndex=-1
    end if
    return
  end subroutine cieFileReadFile
