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
  Implements a cooling function class which interpolates in a tabulated cooling function read from file.
  !!}

  use :: Numerical_Interpolation, only : interpolator
  use :: Table_Labels           , only : enumerationExtrapolationTypeType

  !![
  <coolingFunction name="coolingFunctionCIEFile">
   <description>
    A cooling function class providing a cooling function interpolated from a table read from file.  The HDF5 file containing
    the table should have the following form:
    \begin{verbatim}
    HDF5 "coolingFunction.hdf5" {
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
       DATASET "coolingRate" {
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
                STRSIZE 8;
                STRPAD H5T_STR_NULLTERM;
                CSET H5T_CSET_ASCII;
                CTYPE H5T_C_S1;
             }
             DATASPACE  SCALAR
             DATA {
             (0): "powerLaw"
             }
          }
          ATTRIBUTE "extrapolateLow" {
             DATATYPE  H5T_STRING {
                STRSIZE 8;
                STRPAD H5T_STR_NULLTERM;
                CSET H5T_CSET_ASCII;
                CTYPE H5T_C_S1;
             }
             DATASPACE  SCALAR
             DATA {
             (0): "powerLaw"
             }
          }
       }
       DATASET "energyContinuum" {
          DATATYPE  H5T_IEEE_F64LE
          DATASPACE  SIMPLE { ( 7 ) / ( 7 ) }
       }
       DATASET "powerEmittedFractionalCumulative" {
          DATATYPE  H5T_IEEE_F64LE
          DATASPACE  SIMPLE { ( 7, 8, 10 ) / ( 7, 8, 10 ) }
       }
    }
    }
    \end{verbatim}
    The {\normalfont \ttfamily temperature} dataset should specify temperature (in Kelvin), while the {\normalfont \ttfamily
    metallicity} dataset should give the logarithmic metallicity relative to Solar (a value of -999 or less is taken to imply
    zero metallicity). The {\normalfont \ttfamily coolingRate} dataset should specify the cooling function (in ergs cm$^3$
    s$^{-1}$ computed for a hydrogen density of 1 cm$^{-3}$) respectively at each temperature/metallicity pair. The
    {\normalfont \ttfamily extrapolateLow} and {\normalfont \ttfamily extrapolateHigh} attributes of the {\normalfont \ttfamily
    temperature} and {\normalfont \ttfamily metallicity} datasets specify how the cooling rate should be extrapolated in the
    low and high vale limits. Allowed options for these attributes are:
    \begin{description}
     \item[{\normalfont \ttfamily zero}] The cooling function is set to zero beyond the relevant limit.
     \item[{\normalfont \ttfamily fixed}] The cooling function is held fixed at the value at the relevant limit.
     \item[{\normalfont \ttfamily powerLaw}] The cooling function is extrapolated assuming a
     power-law dependence beyond the relevant limit. This option is only allowed if the
     cooling function is everywhere positive.
    \end{description}
    If the cooling function is everywhere positive the interpolation will be done in the
    logarithm of temperature, metallicity\footnote{The exception is if the first cooling
    function is tabulated for zero metallicity. In that case, a linear interpolation in
    metallicity is always used between zero and the first non-zero tabulated metallicity.}
    and cooling function. Otherwise, interpolation is linear in these quantities. The cooling
    function is scaled assuming a quadratic dependence on hydrogen density.
  
    The {\normalfont \ttfamily energyContinuum} and {\normalfont \ttfamily powerEmittedFractionalCumulative} are optional. If
    present, {\normalfont \ttfamily powerEmittedFractionalCumulative} gives the cumulative emitted power as a function of
    energy for each tabulated metallicity and temperature. The energies at which the emitted power is tabulated are given by
    {\normalfont \ttfamily energyContinuum}.
   </description>
   <runTimeFileDependencies paths="fileName"/>
  </coolingFunction>
  !!]
  type, extends(coolingFunctionClass) :: coolingFunctionCIEFile
     !!{
     A cooling function class which interpolates in a tabulated cooling function read from file.
     !!}
     private
     type            (varying_string                  )                                :: fileName
     double precision                                                                  :: metallicityMaximum              , metallicityMinimum          , &
          &                                                                               temperatureMaximum              , temperatureMinimum
     type            (enumerationExtrapolationTypeType)                                :: extrapolateMetallicityHigh      , extrapolateMetallicityLow   , &
          &                                                                               extrapolateTemperatureHigh      , extrapolateTemperatureLow
     logical                                                                           :: firstMetallicityIsZero          , logarithmicTable
     integer                                                                           :: metallicityCount                , temperatureCount
     double precision                                                                  :: firstNonZeroMetallicity         , hMetallicityPrevious       , &
          &                                                                               hTemperaturePrevious
     integer         (c_size_t                        )                                :: iMetallicityPrevious            , iTemperaturePrevious
     double precision                                  , allocatable, dimension(:    ) :: metallicities                   , temperatures                , &
          &                                                                               energyContinuum
     double precision                                  , allocatable, dimension(:,:  ) :: coolingFunctionTable
     double precision                                  , allocatable, dimension(:,:,:) :: powerEmittedFractionalCumulative
     type            (interpolator                    )                                :: interpolatorMetallicity         , interpolatorTemperature     , &
          &                                                                               interpolatorEnergy
     double precision                                                                  :: temperaturePrevious             , metallicityPrevious         , &
          &                                                                               temperatureSlopePrevious        , metallicitySlopePrevious    , &
          &                                                                               coolingFunctionPrevious         , coolingFunctionSlopePrevious
   contains
     !![
     <methods>
       <method description="Read the named cooling function file."                         method="readFile"            />
       <method description="Compute interpolating factors in a CIE cooling function file." method="interpolatingFactors"/>
       <method description="Interpolate in the cooling function."                          method="interpolate"         />
     </methods>
     !!]
     procedure :: readFile                           => cieFileReadFile
     procedure :: interpolatingFactors               => cieFileInterpolatingFactors
     procedure :: interpolate                        => cieFileInterpolate
     procedure :: coolingFunction                    => cieFileCoolingFunction
     procedure :: coolingFunctionTemperatureLogSlope => cieFileCoolingFunctionTemperatureLogSlope
     procedure :: coolingFunctionDensityLogSlope     => cieFileCoolingFunctionDensityLogSlope
     procedure :: coolingFunctionFractionInBand      => cieFileCoolingFunctionFractionInBand
  end type coolingFunctionCIEFile

  interface coolingFunctionCIEFile
     !!{
     Constructors for the \refClass{coolingFunctionCIEFile} cooling function class.
     !!}
     module procedure cieFileConstructorParameters
     module procedure cieFileConstructorInternal
  end interface coolingFunctionCIEFile

  ! Current file format version for CIE cooling function files.
  integer, parameter :: fileFormatVersionCurrent=1

contains

  function cieFileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{coolingFunctionCIEFile} cooling function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(coolingFunctionCIEFile)                :: self
    type(inputParameters       ), intent(inout) :: parameters
    type(varying_string        )                :: fileName
    
    !![
    <inputParameter>
      <name>fileName</name>
      <source>parameters</source>
      <description>The name of the file containing a tabulation of the collisional ionization equilibrium cooling function.</description>
    </inputParameter>
    !!]
    self=coolingFunctionCIEFile(fileName)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function cieFileConstructorParameters

  function cieFileConstructorInternal(fileName) result(self)
    !!{
    Internal constructor for the \refClass{coolingFunctionCIEFile} cooling function class.
    !!}
    implicit none
    type(coolingFunctionCIEFile)                :: self
    type(varying_string        ), intent(in   ) :: fileName
    !![
    <constructorAssign variables="fileName"/>
    !!]
    
    call self%readFile(fileName)
    self%temperaturePrevious     =-1.0d0
    self%metallicityPrevious     =-1.0d0
    self%temperatureSlopePrevious=-1.0d0
    self%metallicitySlopePrevious=-1.0d0
    return
  end function cieFileConstructorInternal

  double precision function cieFileCoolingFunction(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !!{
    Return the cooling function by interpolating in tabulated CIE data read from a file.
    !!}
    use            :: Abundances_Structure         , only : Abundances_Get_Metallicity, abundances                  , metallicityTypeLinearByMassSolar
    use            :: Chemical_Abundances_Structure, only : chemicalAbundances
    use, intrinsic :: ISO_C_Binding                , only : c_size_t
    use            :: Radiation_Fields             , only : radiationFieldClass
    use            :: Table_Labels                 , only : extrapolationTypeFix      , extrapolationTypeExtrapolate, extrapolationTypeZero
    implicit none
    class           (coolingFunctionCIEFile), intent(inout) :: self
    type            (treeNode              ), intent(inout) :: node
    double precision                        , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances            ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances    ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass   ), intent(inout) :: radiation
    integer         (c_size_t              )                :: iMetallicity         , iTemperature
    double precision                                        :: hMetallicity         , hTemperature  , &
         &                                                     metallicityUse       , temperatureUse
    !$GLC attributes unused :: node, chemicalDensities, radiation

    ! Handle out of range temperatures.
    temperatureUse=temperature
    if (temperatureUse < self%temperatureMinimum) then
       select case (self%extrapolateTemperatureLow%ID)
       case (extrapolationTypeZero%ID)
          cieFileCoolingFunction=0.0d0
          return
       case (extrapolationTypeFix%ID,extrapolationTypeExtrapolate%ID)
          temperatureUse=self%temperatureMinimum
       end select
    end if
    if (temperatureUse > self%temperatureMaximum) then
       select case (self%extrapolateTemperatureHigh%ID)
       case (extrapolationTypeZero%ID)
          cieFileCoolingFunction=0.0d0
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
          cieFileCoolingFunction=0.0d0
          return
       case (extrapolationTypeFix%ID)
          metallicityUse=self%metallicityMinimum
       end select
    end if
    if (metallicityUse > self%metallicityMaximum) then
       select case (self%extrapolateMetallicityHigh%ID)
       case (extrapolationTypeZero%ID)
          cieFileCoolingFunction=0.0d0
          return
       case (extrapolationTypeFix%ID)
          metallicityUse=self%metallicityMaximum
       end select
    end if
    ! Check if we need to recompute the cooling function.
    if     (                                            &
         &   temperatureUse /= self%temperaturePrevious &
         &  .or.                                        &
         &   metallicityUse /= self%metallicityPrevious &
         & ) then
       ! Get the interpolation.
       call self%interpolatingFactors(temperatureUse,metallicityUse,iTemperature,hTemperature,iMetallicity,hMetallicity)
       ! Do the interpolation.
       self%coolingFunctionPrevious=self%interpolate(iTemperature,hTemperature,iMetallicity,hMetallicity)
       ! Store the temperature and metallicity for which calculation was performed.
       self%temperaturePrevious=temperatureUse
       self%metallicityPrevious=metallicityUse
    end if
    ! Scale to the specified density assuming all processes are proportional to hydrogen density squared.
    cieFileCoolingFunction=+self%coolingFunctionPrevious    &
         &                 *numberDensityHydrogen       **2
    return
  end function cieFileCoolingFunction

  double precision function cieFileCoolingFunctionFractionInBand(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation,energyLow,energyHigh)
    !!{
    Return the fraction of the cooling function due to emission in the given energy band by interpolating in tabulated CIE data read from a file.
    !!}
    use            :: Abundances_Structure         , only : Abundances_Get_Metallicity, abundances                  , metallicityTypeLinearByMassSolar
    use            :: Chemical_Abundances_Structure, only : chemicalAbundances
    use            :: Error                        , only : Error_Report
    use, intrinsic :: ISO_C_Binding                , only : c_size_t
    use            :: Radiation_Fields             , only : radiationFieldClass
    use            :: Table_Labels                 , only : extrapolationTypeFix      , extrapolationTypeExtrapolate, extrapolationTypeZero
    implicit none
    class           (coolingFunctionCIEFile), intent(inout)  :: self 
    type            (treeNode              ), intent(inout)  :: node
    double precision                        , intent(in   )  :: numberDensityHydrogen, temperature        , &
         &                                                      energyLow            , energyHigh
    type            (abundances            ), intent(in   )  :: gasAbundances
    type            (chemicalAbundances    ), intent(in   )  :: chemicalDensities
    class           (radiationFieldClass   ), intent(inout)  :: radiation
    double precision                        , dimension(0:1) :: hEnergyLow           , hEnergyHigh
    integer         (c_size_t              )                 :: iMetallicity         , iTemperature       , &
         &                                                      iEnergyLow           , iEnergyHigh
    integer                                                  :: i
    double precision                                         :: hMetallicity         , hTemperature       , &
         &                                                      metallicityUse       , temperatureUse     , &
         &                                                      powerCumulativeLow   , powerCumulativeHigh
    !$GLC attributes unused :: node,chemicalDensities, radiation, numberDensityHydrogen

    ! Abort if cumulative power is not available.
    cieFileCoolingFunctionFractionInBand=0.0d0
    if (.not.allocated(self%powerEmittedFractionalCumulative)) call Error_Report('cumulative power data is not available'//{introspection:location})
    ! Handle out of range temperatures.
    temperatureUse=temperature
    if (temperatureUse < self%temperatureMinimum) then
       select case (self%extrapolateTemperatureLow%ID)
       case (extrapolationTypeZero%ID)
          return
       case (extrapolationTypeFix%ID,extrapolationTypeExtrapolate%ID)
          temperatureUse=self%temperatureMinimum
       end select
    end if
    if (temperatureUse > self%temperatureMaximum) then
       select case (self%extrapolateTemperatureHigh%ID)
       case (extrapolationTypeZero%ID)
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
          return
       case (extrapolationTypeFix%ID)
          metallicityUse=self%metallicityMinimum
       end select
    end if
    if (metallicityUse > self%metallicityMaximum) then
       select case (self%extrapolateMetallicityHigh%ID)
       case (extrapolationTypeZero%ID)
          return
       case (extrapolationTypeFix%ID)
          metallicityUse=self%metallicityMaximum
       end select
    end if
    ! Get the interpolation in temperature and metallicity.
    call self%interpolatingFactors(temperatureUse,metallicityUse,iTemperature,hTemperature,iMetallicity,hMetallicity)
    ! Get the interpolation in energy.
    call self%interpolatorEnergy%linearFactors(energyLow ,iEnergyLow ,hEnergyLow )
    call self%interpolatorEnergy%linearFactors(energyHigh,iEnergyHigh,hEnergyHigh)
    ! Do the interpolation.
    powerCumulativeLow =0.0d0
    powerCumulativeHigh=0.0d0
    do i=0,1
       powerCumulativeLow =+powerCumulativeLow                                                                                                                          &
            &              +self%powerEmittedFractionalCumulative(iMetallicity  ,iTemperature  ,iEnergyLow +i)*(1.0d0-hTemperature)*(1.0d0-hMetallicity)*hEnergyLow (i) &
            &              +self%powerEmittedFractionalCumulative(iMetallicity+1,iTemperature  ,iEnergyLow +i)*(1.0d0-hTemperature)*(      hMetallicity)*hEnergyLow (i) &
            &              +self%powerEmittedFractionalCumulative(iMetallicity  ,iTemperature+1,iEnergyLow +i)*(      hTemperature)*(1.0d0-hMetallicity)*hEnergyLow (i) &
            &              +self%powerEmittedFractionalCumulative(iMetallicity+1,iTemperature+1,iEnergyLow +i)*(      hTemperature)*(      hMetallicity)*hEnergyLow (i)
       powerCumulativeHigh=+powerCumulativeHigh                                                                                                                         &
            &              +self%powerEmittedFractionalCumulative(iMetallicity  ,iTemperature  ,iEnergyHigh+i)*(1.0d0-hTemperature)*(1.0d0-hMetallicity)*hEnergyHigh(i) &
            &              +self%powerEmittedFractionalCumulative(iMetallicity+1,iTemperature  ,iEnergyHigh+i)*(1.0d0-hTemperature)*(      hMetallicity)*hEnergyHigh(i) &
            &              +self%powerEmittedFractionalCumulative(iMetallicity  ,iTemperature+1,iEnergyHigh+i)*(      hTemperature)*(1.0d0-hMetallicity)*hEnergyHigh(i) &
            &              +self%powerEmittedFractionalCumulative(iMetallicity+1,iTemperature+1,iEnergyHigh+i)*(      hTemperature)*(      hMetallicity)*hEnergyHigh(i)
    end do
    cieFileCoolingFunctionFractionInBand=+powerCumulativeHigh &
         &                               -powerCumulativeLow
    return
  end function cieFileCoolingFunctionFractionInBand

  double precision function cieFileCoolingFunctionTemperatureLogSlope(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !!{
    Return the slope of the cooling function with respect to temperature by interpolating in tabulated CIE data
    read from a file.
    !!}
    use            :: Abundances_Structure         , only : Abundances_Get_Metallicity, abundances                  , metallicityTypeLinearByMassSolar
    use            :: Chemical_Abundances_Structure, only : chemicalAbundances
    use, intrinsic :: ISO_C_Binding                , only : c_size_t
    use            :: Radiation_Fields             , only : radiationFieldClass
    use            :: Table_Labels                 , only : extrapolationTypeFix      , extrapolationTypeExtrapolate, extrapolationTypeZero
    implicit none
    class           (coolingFunctionCIEFile), intent(inout) :: self
    type            (treeNode              ), intent(inout) :: node
    double precision                        , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances            ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances    ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass   ), intent(inout) :: radiation
    double precision                                        :: coolingFunction
    integer         (c_size_t              )                :: iMetallicity         , iTemperature
    double precision                                        :: hMetallicity         , hTemperature  , &
         &                                                     metallicityUse       , temperatureUse
    !$GLC attributes unused :: node
    
    ! Get the cooling function.
    coolingFunction=self%coolingFunction(node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    ! Handle out of range temperatures.
    temperatureUse=temperature
    if (temperatureUse < self%temperatureMinimum) then
       select case (self%extrapolateTemperatureLow%ID)
       case (extrapolationTypeZero%ID,extrapolationTypeFix%ID)
          cieFileCoolingFunctionTemperatureLogSlope=0.0d0
          return
       case (extrapolationTypeExtrapolate%ID)
          temperatureUse=self%temperatureMinimum
       end select
    end if
    if (temperatureUse > self%temperatureMaximum) then
       select case (self%extrapolateTemperatureHigh%ID)
       case (extrapolationTypeZero%ID,extrapolationTypeFix%ID)
          cieFileCoolingFunctionTemperatureLogSlope=0.0d0
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
          cieFileCoolingFunctionTemperatureLogSlope=0.0d0
          return
       case (extrapolationTypeFix%ID)
          metallicityUse=self%metallicityMinimum
       end select
    end if
    if (metallicityUse > self%metallicityMaximum) then
       select case (self%extrapolateMetallicityHigh%ID)
       case (extrapolationTypeZero%ID)
          cieFileCoolingFunctionTemperatureLogSlope=0.0d0
          return
       case (extrapolationTypeFix%ID)
          metallicityUse=self%metallicityMaximum
       end select
    end if
    ! Check if we need to recompute the cooling function.
    if     (                                                 &
         &   temperatureUse /= self%temperatureSlopePrevious &
         &  .or.                                             &
         &   metallicityUse /= self%metallicitySlopePrevious &
         & ) then

       ! Get the interpolation.
       call self%interpolatingFactors(temperatureUse,metallicityUse,iTemperature,hTemperature,iMetallicity,hMetallicity)
       ! Do the interpolation.
       self%coolingFunctionSlopePrevious=+(                                                            &
            &                              +(+self%coolingFunctionTable(iMetallicity  ,iTemperature+1) &
            &                                -self%coolingFunctionTable(iMetallicity  ,iTemperature  ) &
            &                               )                                                          &
            &                              *(1.0d0-hMetallicity)                                       &
            &                              +(+self%coolingFunctionTable(iMetallicity+1,iTemperature+1) &
            &                                -self%coolingFunctionTable(iMetallicity+1,iTemperature  ) &
            &                               )                                                          &
            &                              *(      hMetallicity)                                       &
            &                              )                                                           &
            &                              /(                                                          &
            &                                +self%temperatures   (iTemperature+1                    ) &
            &                                -self%temperatures   (iTemperature                      ) &
            &                               )

       ! Convert to logarithmic gradient if table was not stored logarithmically.
       if (.not.self%logarithmicTable)                                                &
            & self%coolingFunctionSlopePrevious=                                      &
            &  +self%coolingFunctionSlopePrevious                                     &
            &  *temperature                                                           &
            &  /self%interpolate(iTemperature,hTemperature,iMetallicity,hMetallicity)
       ! Store the temperature and metallicity for which calculation was performed.
       self%temperatureSlopePrevious=temperatureUse
       self%metallicitySlopePrevious=metallicityUse
    end if
    ! Return the stored value.
    cieFileCoolingFunctionTemperatureLogSlope=self%coolingFunctionSlopePrevious
    return
  end function cieFileCoolingFunctionTemperatureLogSlope

  double precision function cieFileCoolingFunctionDensityLogSlope(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !!{
    Return the logarithmic slope of the cooling function with respect to density.
    !!}
    use :: Abundances_Structure         , only : abundances
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (coolingFunctionCIEFile), intent(inout) :: self
    type            (treeNode              ), intent(inout) :: node
    double precision                        , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances            ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances    ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass   ), intent(inout) :: radiation
    !$GLC attributes unused :: self, node, numberDensityHydrogen, temperature, gasAbundances, chemicalDensities, radiation

    ! Logarithmic slope is always 2 for a CIE cooling function.
    cieFileCoolingFunctionDensityLogSlope=2.0d0
    return
  end function cieFileCoolingFunctionDensityLogSlope

  subroutine cieFileReadFile(self,fileName)
    !!{
    Read in data from a cooling function file.
    !!}
    use :: Display           , only : displayIndent                     , displayUnindent     , verbosityLevelWorking       , displayGreen         , &
         &                            displayReset
    use :: File_Utilities    , only : File_Name_Expand
    use :: Error             , only : Error_Report                      , errorStatusSuccess
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : varying_string
    use :: Table_Labels      , only : enumerationExtrapolationTypeEncode, extrapolationTypeFix, extrapolationTypeExtrapolate, extrapolationTypeZero, &
         &                            enumerationExtrapolationTypeDescribe
    implicit none
    class           (coolingFunctionCIEFile), intent(inout) :: self
    type            (varying_string        ), intent(in   ) :: fileName
    double precision                        , parameter     :: metallicityLogarithmicZero=-999.0d0
    type            (varying_string        )                :: limitType
    integer                                                 :: fileFormatVersion                  , status
    type            (hdf5Object            )                :: coolingFunctionFile                , metallicityDataset, &
         &                                                     temperatureDataset
    
    !$ call hdf5Access%set()
    ! Read the file.
    call displayIndent('Reading file: '//char(fileName),verbosityLevelWorking)
    coolingFunctionFile=hdf5Object(char(File_Name_Expand(char(fileName))),readOnly=.true.)
    ! Check the file format version of the file.
    call coolingFunctionFile%readAttribute('fileFormat',fileFormatVersion)
    if (fileFormatVersion /= fileFormatVersionCurrent) call Error_Report('file format version is out of date'//{introspection:location})
    ! Read datasets.
    call coolingFunctionFile%readDataset('temperature',self%temperatures        )
    call coolingFunctionFile%readDataset('metallicity',self%metallicities       )
    call coolingFunctionFile%readDataset('coolingRate',self%coolingFunctionTable)
    self%metallicityCount=size(self%metallicities)
    self%temperatureCount=size(self%temperatures )
    ! Unlog metallicities.
    where (self%metallicities > metallicityLogarithmicZero)
       self%metallicities=10.0d0**self%metallicities
    elsewhere
       self%metallicities=0.0d0
    end where
    ! Extract extrapolation methods from the file.
    metallicityDataset=coolingFunctionFile%openDataset('metallicity')
    call metallicityDataset%readAttribute('extrapolateLow' ,limitType,allowPseudoScalar=.true.)
    self%extrapolateMetallicityLow =enumerationExtrapolationTypeEncode(char(limitType),includesPrefix=.false.,status=status)
    if (status /= errorStatusSuccess) call Error_Report("low metallicity extrapolation type '" //char(limitType)//"' in file '"//trim(fileName)//"' is invalid"//char(10)//displayGreen()//"HELP:"//displayReset()//enumerationExtrapolationTypeDescribe()//{introspection:location})
    call metallicityDataset%readAttribute('extrapolateHigh',limitType,allowPseudoScalar=.true.)
    self%extrapolateMetallicityHigh=enumerationExtrapolationTypeEncode(char(limitType),includesPrefix=.false.,status=status)
    if (status /= errorStatusSuccess) call Error_Report("high metallicity extrapolation type '"//char(limitType)//"' in file '"//trim(fileName)//"' is invalid"//char(10)//displayGreen()//"HELP:"//displayReset()//enumerationExtrapolationTypeDescribe()//{introspection:location})
    temperatureDataset=coolingFunctionFile%openDataset('temperature')
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
    ! Read optional datasets.
    if (coolingFunctionFile%hasDataset('energyContinuum')) then
       call coolingFunctionFile%readDataset('energyContinuum'                 ,self%energyContinuum                 )
       call coolingFunctionFile%readDataset('powerEmittedFractionalCumulative',self%powerEmittedFractionalCumulative)
    end if
    call displayUnindent('done',verbosityLevelWorking)
    !$ call hdf5Access%unset()
    ! Store table ranges for convenience.
    self%metallicityMinimum=self%metallicities(                    1)
    self%metallicityMaximum=self%metallicities(self%metallicityCount)
    self%temperatureMinimum=self%temperatures (                    1)
    self%temperatureMaximum=self%temperatures (self%temperatureCount)
    ! Decide whether or not to make the tables logarithmic.
    self%logarithmicTable=all(self%coolingFunctionTable > 0.0d0)
    if (self%logarithmicTable) then
       self%firstMetallicityIsZero=(self%metallicities(1) == 0.0d0)
       if (self%firstMetallicityIsZero) self%firstNonZeroMetallicity=self%metallicities(2)
       where (self%metallicities > 0.0d0)
          self%metallicities=log(self%metallicities)
       elsewhere
          self%metallicities=metallicityLogarithmicZero
       end where
       self%temperatures        =log(self%temperatures        )
       self%coolingFunctionTable=log(self%coolingFunctionTable)
    else
       if     (                                                                &
            &  self%extrapolateTemperatureLow  == extrapolationTypeExtrapolate &
            &  .or.                                                            &
            &  self%extrapolateTemperatureHigh == extrapolationTypeExtrapolate &
            & )                                                                &
            & call Error_Report('extrapolation allowed only in loggable tables'//{introspection:location})
    end if
    if     (                                                                &
         &  self%extrapolateMetallicityLow  == extrapolationTypeExtrapolate &
         &   .or.                                                           &
         &  self%extrapolateMetallicityHigh == extrapolationTypeExtrapolate &
         & )                                                                &
         & call Error_Report('extrapolation not allowed in metallicity'//{introspection:location})
    ! Build interpolators.
    self       %interpolatorMetallicity=interpolator(self%metallicities  )
    self       %interpolatorTemperature=interpolator(self%temperatures   )
    if (allocated(self%energyContinuum)) &
         & self%interpolatorEnergy     =interpolator(self%energyContinuum)
    return
  end subroutine cieFileReadFile

  subroutine cieFileInterpolatingFactors(self,temperature,metallicity,iTemperature,hTemperature,iMetallicity,hMetallicity)
    !!{
    Determine the interpolating parameters.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (coolingFunctionCIEFile), intent(inout) :: self
    double precision                        , intent(in   ) :: metallicity   , temperature
    integer         (c_size_t              ), intent(  out) :: iMetallicity  , iTemperature
    double precision                        , intent(  out) :: hMetallicity  , hTemperature
    double precision                                        :: metallicityUse, temperatureUse

    ! Get interpolation in temperature.
    if (temperature /= self%temperaturePrevious) then
       temperatureUse=    temperature
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
       if (self%firstMetallicityIsZero .and. metallicityUse < self%firstNonZeroMetallicity) then
          self%iMetallicityPrevious=1
          self%hMetallicityPrevious=metallicityUse/self%firstNonZeroMetallicity
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

  double precision function cieFileInterpolate(self,iTemperature,hTemperature,iMetallicity,hMetallicity)
    !!{
    Perform the interpolation.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (coolingFunctionCIEFile), intent(inout) :: self
    integer         (c_size_t              ), intent(in   ) :: iMetallicity, iTemperature
    double precision                        , intent(in   ) :: hMetallicity, hTemperature

    ! Do the interpolation.
    cieFileInterpolate=+self%coolingFunctionTable(iMetallicity  ,iTemperature  )*(1.0d0-hTemperature)*(1.0d0-hMetallicity) &
         &             +self%coolingFunctionTable(iMetallicity+1,iTemperature  )*(1.0d0-hTemperature)*(      hMetallicity) &
         &             +self%coolingFunctionTable(iMetallicity  ,iTemperature+1)*(      hTemperature)*(1.0d0-hMetallicity) &
         &             +self%coolingFunctionTable(iMetallicity+1,iTemperature+1)*(      hTemperature)*(      hMetallicity)
    ! Exponentiate the result if the table was stored as the log.
    if (self%logarithmicTable) cieFileInterpolate=exp(cieFileInterpolate)
    return
  end function cieFileInterpolate
