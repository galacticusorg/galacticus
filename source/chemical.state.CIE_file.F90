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

  !% Implements a chemical state class which reads and interpolates a collisional ionization equilibrium chemical state from a file.

  use FGSL
  
  !# <chemicalState name="chemicalStateCIEFile" defaultThreadPrivate="yes">
  !#  <description>
  !#   Class providing chemical state via interpolation of tabulated values read from file. The XML file containing the table
  !#   should have the following form:
  !#   \begin{verbatim}
  !#    &lt;chemicalStates&gt;
  !#    &lt;chemicalState&gt;
  !#      &lt;temperature&gt;
  !#        &lt;datum&gt;10000.0&lt;/datum&gt;
  !#        &lt;datum&gt;15000.0&lt;/datum&gt;
  !#        .
  !#        .
  !#        .
  !#      &lt;/temperature&gt;
  !#      &lt;electronDensity&gt;
  !#        &lt;datum&gt;1.0e-23&lt;/datum&gt;
  !#        &lt;datum&gt;1.7e-23&lt;/datum&gt;
  !#        .
  !#        .
  !#        .
  !#      &lt;/electronDensity&gt;
  !#      &lt;hiDensity&gt;
  !#        &lt;datum&gt;0.966495864314214&lt;/datum&gt;
  !#        &lt;datum&gt;0.965828463162061&lt;/datum&gt;
  !#        .
  !#        .
  !#        .
  !#      &lt;/hiDensity&gt;
  !#      &lt;hiiDensity&gt;
  !#        &lt;datum&gt;0.033504135685786&lt;/datum&gt;
  !#        &lt;datum&gt;0.0341715368379391&lt;/datum&gt;
  !#        .
  !#        .
  !#        .
  !#      &lt;/hiiDensity&gt;
  !#      &lt;metallicity&gt;-4.0&lt;/metallicity&gt;
  !#    &lt;/chemicalState&gt;
  !#    &lt;chemicalState&gt;
  !#    .
  !#    .
  !#    .
  !#    &lt;/chemicalState&gt;
  !#    &lt;description&gt;Some description of what this chemical state is.&lt;/description&gt;
  !#    &lt;extrapolation&gt;
  !#      &lt;metallicity&gt;
  !#        &lt;limit&gt;low&lt;/limit&gt;
  !#        &lt;method&gt;fixed&lt;/method&gt;
  !#      &lt;/metallicity&gt;
  !#      &lt;metallicity&gt;
  !#        &lt;limit&gt;high&lt;/limit&gt;
  !#        &lt;method&gt;fixed&lt;/method&gt;
  !#      &lt;/metallicity&gt;
  !#      &lt;temperature&gt;
  !#        &lt;limit&gt;low&lt;/limit&gt;
  !#        &lt;method&gt;fixed&lt;/method&gt;
  !#      &lt;/temperature&gt;
  !#      &lt;temperature&gt;
  !#        &lt;limit&gt;high&lt;/limit&gt;
  !#        &lt;method&gt;fixed&lt;/method&gt;
  !#      &lt;/temperature&gt;
  !#    &lt;/extrapolation&gt;
  !#   &lt;/chemicalStates&gt;
  !#   \end{verbatim}
  !#   Each {\normalfont \ttfamily chemicalState} element should contain two lists (inside
  !#   {\normalfont \ttfamily temperature} and {\normalfont \ttfamily electronDensity} tags) of
  !#   {\normalfont \ttfamily datum} elements which specify temperature (in Kelvin) and electron
  !#   density (by number, relative to hydrogen) respectively, and a {\normalfont \ttfamily
  !#   metallicity} element which gives the logarithmic metallcity relative to Solar (a value of
  !#   -999 or less is taken to imply zero metallicity). Optionally, {\normalfont \ttfamily
  !#   hiDensity} and {\normalfont \ttfamily hiiDensity} elements may be added containing lists
  !#   of H{\normalfont \scshape i} and H{\normalfont \scshape ii} densities (by number,
  !#   relative to hydrogen) respectively. Any number of {\normalfont \ttfamily coolingFunction}
  !#   elements may appear, but they must be in order of increasing metallicity and must all
  !#   contain the same set of temperatures. The {\normalfont \ttfamily extrapolation} element
  !#   defines how the table is to be extrapolated in the {\normalfont \ttfamily low} and
  !#   {\normalfont \ttfamily high} limits of {\normalfont \ttfamily temperature} and
  !#   {\normalfont \ttfamily metallicity}. The {\normalfont \ttfamily method} elements can take
  !#   the following values:
  !#   \begin{description}
  !#    \item[{\normalfont \ttfamily zero}] The electron density is set to zero beyond the relevant limit.
  !#    \item[{\normalfont \ttfamily fixed}] The electron density is held fixed at the value at the relevant limit.
  !#    \item[{\normalfont \ttfamily power law}] The electron density is extrapolated assuming a
  !#    power-law dependence beyond the relevant limit. This option is only allowed if the
  !#    electron density is everywhere positive.
  !#   \end{description}
  !#   If the electron density is everywhere positive the interpolation will be done in the
  !#   logarithmic of temperature, metallicity\footnote{The exception is if the first electron
  !#   density is tabulated for zero metallicity. In that case, a linear interpolation in
  !#   metallicity is always used between zero and the first non-zero tabulated metallicity.}
  !#   and electron density. Otherwise, interpolation is linear in these quantities. The
  !#   electron density is scaled assuming a linear dependence on hydrogen density.
  !#  </description>
  !# </chemicalState>
  type, extends(chemicalStateClass) :: chemicalStateCIEFile
     !% A chemical state class which interpolates state given in a file assuming collisional ionization equilibrium.
     private
     type            (varying_string    )                              :: fileName
     double precision                                                  :: metallicityMaximum               , metallicityMinimum         , &
          &                                                               temperatureMaximum               , temperatureMinimum
     integer                                                           :: extrapolateMetallicityHigh       , extrapolateMetallicityLow  , &
          &                                                               extrapolateTemperatureHigh       , extrapolateTemperatureLow
     logical                                                           :: firstMetallicityIsZero           , gotHydrogenAtomic          , &
          &                                                               gotHydrogenCation                , logarithmicTable
     integer                                                           :: metallicityCount                 , temperatureCount
     double precision                                                  :: firstNonZeroMetallicity          , electronDensityPrevious    , &
          &                                                               metallicityPrevious              , temperaturePrevious        , &
          &                                                               electronDensitySlopePrevious     , metallicitySlopePrevious   , &
          &                                                               temperatureSlopePrevious         , metallicityChemicalPrevious, &
          &                                                               temperatureChemicalPrevious
     type            (chemicalAbundances)                              :: chemicalDensitiesPrevious
     double precision                    , allocatable, dimension(:  ) :: metallicities                    , temperatures
     double precision                    , allocatable, dimension(:,:) :: densityElectron                  , densityHydrogenAtomic      , &
          &                                                               densityHydrogenCation
     logical                                                           :: resetMetallicity                 , resetTemperature
     type            (fgsl_interp_accel )                              :: acceleratorMetallicity           , acceleratorTemperature
     integer                                                           :: atomicHydrogenCationChemicalIndex, atomicHydrogenChemicalIndex, &
          &                                                               electronChemicalIndex
   contains
     !@ <objectMethods>
     !@   <object>chemicalStateCIEFile</object>
     !@   <objectMethod>
     !@     <method>readFile</method>
     !@     <type>void</type>
     !@     <arguments>\textcolor{red}{\textless char(len=*)\textgreater} fileName\argin</arguments>
     !@     <description>Read the named chemical state file.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>interpolatingFactors</method>
     !@     <type>void</type>
     !@     <arguments>\doublezero\ temperature\argin, \doublezero\ metallicity\argin, \textcolor{red}{\textless integer(c\_size\_t)\textgreater} iTemperature\argout, \doublezero\ hTemperature\argout, \textcolor{red}{\textless integer(c\_size\_t)\textgreater} iMetallicity\argout, \doublezero\ hMetallicity\argout</arguments>
     !@     <description>Compute interpolating factors in a CIE chemical state file.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>interpolate</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless integer(c\_size\_t)\textgreater} iTemperature\argin, \doublezero\ hTemperature\argin, \textcolor{red}{\textless integer(c\_size\_t)\textgreater} iMetallicity\argin, \doublezero\ hMetallicity\argin, \doubletwo\ density\argin</arguments>
     !@     <description>Interpolate in the given density table.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                                       cieFileDestructor
     procedure :: readFile                           => cieFileReadFile
     procedure :: interpolatingFactors               => cieFileInterpolatingFactors
     procedure :: interpolate                        => cieFileInterpolate
     procedure :: electronDensity                    => cieFileElectronDensity
     procedure :: electronDensityTemperatureLogSlope => cieFileElectronDensityTemperatureLogSlope
     procedure :: electronDensityDensityLogSlope     => cieFileElectronDensityDensityLogSlope
     procedure :: chemicalDensities                  => cieFileChemicalDensities
     procedure :: descriptor                         => cieFileDescriptor
  end type chemicalStateCIEFile

  interface chemicalStateCIEFile
     !% Constructors for the ``CIE file'' chemical state class.
     module procedure cieFileConstructorParameters
     module procedure cieFileConstructorInternal
  end interface chemicalStateCIEFile

  ! Current file format version for CIE chemical state files.
  integer, parameter :: cieFileFormatVersionCurrent=1

contains

  function cieFileConstructorParameters(parameters)
    !% Constructor for the ``CIE file'' chemical state class which takes a parameter set as input.
    use Galacticus_Input_Paths
    implicit none
    type(chemicalStateCIEFile)                :: cieFileConstructorParameters
    type(inputParameters     ), intent(inout) :: parameters
    type(varying_string      )                :: fileName
    !# <inputParameterList label="allowedParameterNames" />

    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>fileName</name>
    !#   <source>parameters</source>
    !#   <description>The name of the file containing a tabulation of the collisional ionization equilibrium chemical state.</description>
    !#   <type>string</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    ! Construct the instance.    
    cieFileConstructorParameters=cieFileConstructorInternal(char(fileName))
    return
  end function cieFileConstructorParameters
  
  function cieFileConstructorInternal(fileName)
    !% Internal constructor for the ``CIE file'' chemical state class.
    implicit none
    type     (chemicalStateCIEFile)                :: cieFileConstructorInternal
    character(len=*               ), intent(in   ) :: fileName
    
    ! Read the file.
    cieFileConstructorInternal%fileName=fileName
    call cieFileConstructorInternal%readFile(fileName)
    ! Initialize.
    cieFileConstructorInternal%electronDensityPrevious        =-1.0d0
    cieFileConstructorInternal%electronDensitySlopePrevious   =-1.0d0
    cieFileConstructorInternal%chemicalDensitiesPrevious      =-unitChemicalAbundances
    cieFileConstructorInternal%    metallicityPrevious        =-1.0d0
    cieFileConstructorInternal%    metallicitySlopePrevious   =-1.0d0
    cieFileConstructorInternal%    metallicityChemicalPrevious=-1.0d0
    cieFileConstructorInternal%    temperaturePrevious        =-1.0d0
    cieFileConstructorInternal%    temperatureSlopePrevious   =-1.0d0
    cieFileConstructorInternal%    temperatureChemicalPrevious=-1.0d0
    cieFileConstructorInternal%resetMetallicity               =.true.
    cieFileConstructorInternal%resetTemperature               =.true.
    return
  end function cieFileConstructorInternal
  
  subroutine cieFileDestructor(self)
    !% Destructor for the ``CIE file'' chemical state class.
    use Numerical_Interpolation
    implicit none
    type(chemicalStateCIEFile), intent(inout) :: self

    ! Free all FGSL objects.
    call Interpolate_Done(                                                      &
         &                interpolationAccelerator=self%acceleratorMetallicity, &
         &                reset                   =self%      resetMetallicity  &
         &               )
    call Interpolate_Done(                                                      &
         &                interpolationAccelerator=self%acceleratorTemperature, &
         &                reset                   =self%      resetTemperature  &
         &               )
    return
  end subroutine cieFileDestructor

  double precision function cieFileElectronDensity(self,numberDensityHydrogen,temperature,gasAbundances,radiation)
    !% Return the electron density by interpolating in tabulated CIE data read from a file.
    use, intrinsic :: ISO_C_Binding
    use               Abundances_Structure
    use               Radiation_Structure
    use               Table_Labels
    implicit none
    class           (chemicalStateCIEFile), intent(inout) :: self
    double precision                      , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances          ), intent(in   ) :: gasAbundances
    type            (radiationStructure  ), intent(in   ) :: radiation
    integer         (c_size_t            )                :: iMetallicity         , iTemperature
    double precision                                      :: hMetallicity         , hTemperature  , &
         &                                                   metallicityUse       , temperatureUse
    !GCC$ attributes unused :: radiation
    
    ! Handle out of range temperatures.
    temperatureUse=temperature
    if (temperatureUse < self%temperatureMinimum) then
       select case (self%extrapolateTemperatureLow)
       case (extrapolationTypeZero)
          cieFileElectronDensity=0.0d0
          return
       case (extrapolationTypeFix,extrapolationTypePowerLaw)
          temperatureUse=self%temperatureMinimum
       end select
    end if
    if (temperatureUse > self%temperatureMaximum) then
       select case (self%extrapolateTemperatureHigh)
       case (extrapolationTypeZero)
          cieFileElectronDensity=0.0d0
          return
       case (extrapolationTypeFix,extrapolationTypePowerLaw)
          temperatureUse=self%temperatureMaximum
       end select
    end if
    ! Handle out of range metallicities.
    metallicityUse=Abundances_Get_Metallicity(gasAbundances,metallicityType=metallicityTypeLinearByMassSolar)
    if (metallicityUse < self%metallicityMinimum) then
       select case (self%extrapolateMetallicityLow)
       case (extrapolationTypeZero)
          cieFileElectronDensity=0.0d0
          return
       case (extrapolationTypeFix)
          metallicityUse=self%metallicityMinimum
       end select
    end if
    if (metallicityUse > self%metallicityMaximum) then
       select case (self%extrapolateMetallicityHigh)
       case (extrapolationTypeZero)
          cieFileElectronDensity=0.0d0
          return
       case (extrapolationTypeFix)
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
    !% Return the logarithmic slope of the electron density with respect to temperature by interpolating in tabulated CIE data
    !% read from a file.
    use, intrinsic :: ISO_C_Binding
    use               Abundances_Structure
    use               Radiation_Structure
    use               Table_Labels
    implicit none
    class           (chemicalStateCIEFile), intent(inout) :: self
    double precision                      , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances          ), intent(in   ) :: gasAbundances
    type            (radiationStructure  ), intent(in   ) :: radiation
    integer         (c_size_t            )                :: iMetallicity         , iTemperature
    double precision                                      :: hMetallicity         , hTemperature  , &
         &                                                   metallicityUse       , temperatureUse
    !GCC$ attributes unused :: radiation, numberDensityHydrogen

    ! Handle out of range temperatures.
    temperatureUse=temperature
    if (temperatureUse < self%temperatureMinimum) then
       select case (self%extrapolateTemperatureLow)
       case (extrapolationTypeZero,extrapolationTypeFix)
          cieFileElectronDensityTemperatureLogSlope=0.0d0
          return
       case (extrapolationTypePowerLaw)
          temperatureUse=self%temperatureMinimum
       end select
    end if
    if (temperatureUse > self%temperatureMaximum) then
       select case (self%extrapolateTemperatureHigh)
       case (extrapolationTypeZero,extrapolationTypeFix)
          cieFileElectronDensityTemperatureLogSlope=0.0d0
          return
       case (extrapolationTypePowerLaw)
          temperatureUse=self%temperatureMaximum
       end select
    end if
    ! Handle out of range metallicities.
    metallicityUse=Abundances_Get_Metallicity(gasAbundances,metallicityType=metallicityTypeLinearByMassSolar)
    if (metallicityUse < self%metallicityMinimum) then
       select case (self%extrapolateMetallicityLow)
       case (extrapolationTypeZero)
          cieFileElectronDensityTemperatureLogSlope=0.0d0
          return
       case (extrapolationTypeFix)
          metallicityUse=self%metallicityMinimum
       end select
    end if
    if (metallicityUse > self%metallicityMaximum) then
       select case (self%extrapolateMetallicityHigh)
       case (extrapolationTypeZero)
          cieFileElectronDensityTemperatureLogSlope=0.0d0
          return
       case (extrapolationTypeFix)
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
            &                                +self%densityElectron(iTemperature+1,iMetallicity  ) &
            &                                -self%densityElectron(iTemperature  ,iMetallicity  ) &
            &                               )                                                     &
            &                              *(1.0d0-hMetallicity)                                  &
            &                              +(                                                     &
            &                                +self%densityElectron(iTemperature+1,iMetallicity+1) &
            &                                -self%densityElectron(iTemperature  ,iMetallicity+1) &
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
    !% Return the logarithmic slope of the electron density with respect to density assuming atomic CIE.
    use Abundances_Structure
    use Radiation_Structure
    implicit none
    class           (chemicalStateCIEFile), intent(inout) :: self
    double precision                      , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances          ), intent(in   ) :: gasAbundances
    type            (radiationStructure  ), intent(in   ) :: radiation
    !GCC$ attributes unused :: self, numberDensityHydrogen, temperature, gasAbundances, radiation

    ! Electron density always scales as total density under CIE conditions.
    cieFileElectronDensityDensityLogSlope=1.0d0
    return
  end function cieFileElectronDensityDensityLogSlope

  subroutine cieFileChemicalDensities(self,chemicalDensities,numberDensityHydrogen,temperature,gasAbundances,radiation)
    !% Return the densities of chemical species at the given temperature and hydrogen density for the specified set of abundances
    !% and radiation field. Units of the returned electron density are cm$^-3$.
    use, intrinsic :: ISO_C_Binding
    use               Abundances_Structure
    use               Radiation_Structure
    use               Chemical_Abundances_Structure
    use               Table_Labels
    implicit none
    class           (chemicalStateCIEFile), intent(inout) :: self
    type            (chemicalAbundances  ), intent(inout) :: chemicalDensities
    double precision                      , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances          ), intent(in   ) :: gasAbundances
    type            (radiationStructure  ), intent(in   ) :: radiation
    integer         (c_size_t          )                  :: iMetallicity         , iTemperature
    double precision                                      :: hMetallicity         , hTemperature  , &
         &                                                   metallicityUse       , temperatureUse
    !GCC$ attributes unused :: radiation

    ! Handle out of range temperatures.
    temperatureUse=temperature
    if (temperatureUse < self%temperatureMinimum) then
       select case (self%extrapolateTemperatureLow)
       case (extrapolationTypeZero)
          call chemicalDensities%reset()
          return
       case (extrapolationTypeFix,extrapolationTypePowerLaw)
          temperatureUse=self%temperatureMinimum
       end select
    end if
    if (temperatureUse > self%temperatureMaximum) then
       select case (self%extrapolateTemperatureHigh)
       case (extrapolationTypeZero)
          call chemicalDensities%reset()
          return
       case (extrapolationTypeFix,extrapolationTypePowerLaw)
          temperatureUse=self%temperatureMaximum
       end select
    end if
    ! Handle out of range metallicities.
    metallicityUse=Abundances_Get_Metallicity(gasAbundances,metallicityType=metallicityTypeLinearByMassSolar)
    if (metallicityUse < self%metallicityMinimum) then
       select case (self%extrapolateMetallicityLow)
       case (extrapolationTypeZero)
          call chemicalDensities%reset()
          return
       case (extrapolationTypeFix)
          metallicityUse=self%metallicityMinimum
       end select
    end if
    if (metallicityUse > self%metallicityMaximum) then
       select case (self%extrapolateMetallicityHigh)
       case (extrapolationTypeZero)
          call chemicalDensities%reset()
          return
       case (extrapolationTypeFix)
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
    !% Determine the interpolating paramters.
    use, intrinsic :: ISO_C_Binding
    use               Numerical_Interpolation
    implicit none
    class           (chemicalStateCIEFile), intent(inout) :: self
    double precision                      , intent(in   ) :: metallicity   , temperature
    integer         (c_size_t            ), intent(  out) :: iMetallicity  , iTemperature
    double precision                      , intent(  out) :: hMetallicity  , hTemperature
    double precision                                      :: metallicityUse, temperatureUse

    ! Copy the input parameters.
    temperatureUse=    temperature
    metallicityUse=max(metallicity,0.0d0)
    ! Get interpolation in temperature.
    if (self%logarithmicTable) temperatureUse=log(temperatureUse)
    iTemperature=max(                                                    &
         &           min(                                                &
         &               Interpolate_Locate(                             &
         &                                  self%temperatures          , &
         &                                  self%acceleratorTemperature, &
         &                                  temperatureUse             , &
         &                                  self%resetTemperature        &
         &                                 )                           , &
         &               self%temperatureCount-1                         &
         &              )                                              , &
         &               1                                               &
         &          )
    hTemperature=+(     temperatureUse                -self%temperatures(iTemperature)) &
         &       /(self%temperatures  (iTemperature+1)-self%temperatures(iTemperature))
    ! Get interpolation in metallicity.
    if     (                                               &
         &                    self%firstMetallicityIsZero  &
         &  .and.                                          &
         &   metallicityUse < self%firstNonZeroMetallicity &
         & ) then
       iMetallicity=+1
       hMetallicity=+metallicityUse               &
            &       /self%firstNonZeroMetallicity
    else
       if (self%logarithmicTable) metallicityUse=log(metallicityUse)
       iMetallicity=max(                                                    &
            &           min(                                                &
            &               Interpolate_Locate(                             &
            &                                  self%metallicities         , &
            &                                  self%acceleratorMetallicity, &
            &                                  metallicityUse             , &
            &                                  self%resetMetallicity        &
            &                                 )                           , &
            &               self%metallicityCount-1                         &
            &              )                                              , &
            &               1                                               &
            &          )
       hMetallicity=+(     metallicityUse                -self%metallicities(iMetallicity)) &
            &       /(self%metallicities (iMetallicity+1)-self%metallicities(iMetallicity))
    end if
    return
  end subroutine cieFileInterpolatingFactors

  double precision function cieFileInterpolate(self,iTemperature,hTemperature,iMetallicity,hMetallicity,density)
    !% Perform the interpolation.
    use, intrinsic :: ISO_C_Binding
    implicit none
    class           (chemicalStateCIEFile)                , intent(inout) :: self
    integer         (c_size_t            )                , intent(in   ) :: iMetallicity, iTemperature
    double precision                                      , intent(in   ) :: hMetallicity, hTemperature
    double precision                      , dimension(:,:), intent(in   ) :: density

    ! Do the interpolation.
    cieFileInterpolate=+density(iTemperature  ,iMetallicity  )*(1.0d0-hTemperature)*(1.0d0-hMetallicity) &
         &             +density(iTemperature  ,iMetallicity+1)*(1.0d0-hTemperature)*(      hMetallicity) &
         &             +density(iTemperature+1,iMetallicity  )*(      hTemperature)*(1.0d0-hMetallicity) &
         &             +density(iTemperature+1,iMetallicity+1)*(      hTemperature)*(      hMetallicity)
    ! Exponentiate the result if the table was stored as the log.
    if (self%logarithmicTable) cieFileInterpolate=exp(cieFileInterpolate)
    return
  end function cieFileInterpolate

  subroutine cieFileReadFile(self,fileName)
    !% Read in data from an chemical state file.
    use Galacticus_Error
    use FoX_DOM
    use Memory_Management
    use Numerical_Comparison
    use Galacticus_Display
    use IO_XML
    use Table_Labels
    implicit none
    class           (chemicalStateCIEFile), intent(inout)             :: self
    character       (len=*               ), intent(in   )             :: fileName
    double precision                      , allocatable, dimension(:) :: temperaturesReference
    type            (node                ), pointer                   :: doc                                  , extrapolation               , &
         &                                                               extrapolationElement                 , metallicityElement          , &
         &                                                               thisChemicalState                    , thisElectronDensity         , &
         &                                                               thisHydrogenAtomicDensity            , thisHydrogenCationDensity   , &
         &                                                               thisTemperature                      , version
    type            (nodeList            ), pointer                   :: chemicalStateList                    , metallicityExtrapolationList, &
         &                                                               temperatureExtrapolationList
    double precision                      , parameter                 :: metallicityLogarithmicZero  =-999.0d0
    integer                                                           :: extrapolationMethod                  , fileFormatVersion           , &
         &                                                               iChemicalState                       , iExtrapolation              , &
         &                                                               ioErr
    character       (len=32        )                                  :: limitType

    !$omp critical (FoX_DOM_Access)
    ! Parse the XML file.
    call Galacticus_Display_Indent('Parsing file: '//fileName,verbosityDebug)
    doc => parseFile(fileName,iostat=ioErr)
    if (ioErr /= 0) call Galacticus_Error_Report('cieFileReadFile','unable to find chemical state file')
    ! Check the file format version of the file.
    version => XML_Get_First_Element_By_Tag_Name(doc,"fileFormat")
    call extractDataContent(version,fileFormatVersion)
    if (fileFormatVersion /= cieFileFormatVersionCurrent) call Galacticus_Error_Report('cieFileReadFile','file format version is out of date')
    ! Get a list of all <chemicalState> elements.
    chemicalStateList     => getElementsByTagname(doc,"chemicalState")
    self%metallicityCount =  getLength(chemicalStateList)
    ! Extract data from first chemical state and count number of temperatures present.
    thisChemicalState      => item(chemicalStateList,0)
    thisTemperature        => XML_Get_First_Element_By_Tag_Name(thisChemicalState,"temperature")
    self%temperatureCount  =  XML_Array_Length                 (thisTemperature  ,"datum"      )
    ! Allocate space for the table.
    if (allocated(self%metallicities        )) call deallocateArray(self%metallicities        )
    if (allocated(self%temperatures         )) call deallocateArray(self%temperatures         )
    if (allocated(self%densityElectron      )) call deallocateArray(self%densityElectron      )
    if (allocated(self%densityHydrogenAtomic)) call deallocateArray(self%densityHydrogenAtomic)
    if (allocated(self%densityHydrogenCation)) call deallocateArray(self%densityHydrogenCation)
    call allocateArray(self%metallicities  ,[                      self%metallicityCount])
    call allocateArray(self%temperatures   ,[self%temperatureCount                      ])
    call allocateArray(self%densityElectron,[self%temperatureCount,self%metallicityCount])
    ! Allocate space for atomic hydrogen density, if such data is included.
    self%gotHydrogenAtomic=(XML_Array_Length(doc,"hiDensity" ) > 0)
    if (self%gotHydrogenAtomic) call allocateArray(self%densityHydrogenAtomic,[self%temperatureCount,self%metallicityCount])
    ! Allocate space for ionized hydrogen density, if such data is included.
    self%gotHydrogenCation=(XML_Array_Length(doc,"hiiDensity") > 0)
    if (self%gotHydrogenCation) call allocateArray(self%densityHydrogenCation,[self%temperatureCount,self%metallicityCount])
    ! Extract data from the chemical states and populate metallicity and temperature arrays.
    allocate(temperaturesReference(0))
    do iChemicalState=0,self%metallicityCount-1
       ! Get required chemical state.
       thisChemicalState  => item(chemicalStateList,iChemicalState)
       ! Extract the metallicity from the <metallicity> element.
       metallicityElement => XML_Get_First_Element_By_Tag_Name(thisChemicalState,"metallicity")
       call extractDataContent(metallicityElement,self%metallicities(iChemicalState+1))
       ! Extract the data.
       thisTemperature                                       => XML_Get_First_Element_By_Tag_Name(thisChemicalState,"temperature"    )
       thisElectronDensity                                   => XML_Get_First_Element_By_Tag_Name(thisChemicalState,"electronDensity")
       if (self%gotHydrogenAtomic) thisHydrogenAtomicDensity => XML_Get_First_Element_By_Tag_Name(thisChemicalState,"hiDensity"      )
       if (self%gotHydrogenCation) thisHydrogenCationDensity => XML_Get_First_Element_By_Tag_Name(thisChemicalState,"hiiDensity"     )
       ! Check that number of temperatures is consistent.
       if (XML_Array_Length(thisTemperature,"datum") /= self%temperatureCount                        ) &
            & call Galacticus_Error_Report('cieFileReadFile','sizes of temperature grids must be the same for all metallicities')
       ! Check that number of chemical states matches number of temperatures.
       if (XML_Array_Length(thisTemperature,"datum") /= XML_Array_Length(thisElectronDensity,"datum")) &
            & call Galacticus_Error_Report('cieFileReadFile','sizes of temperature and electron density arrays must match'      )
       ! Store the chemical state.
       call                             XML_Array_Read_Static(thisTemperature          ,"datum",self%temperatures (:                 ))
       call                             XML_Array_Read_Static(thisElectronDensity      ,"datum",self%densityElectron      (:,iChemicalState+1))
       if (self%gotHydrogenAtomic) call XML_Array_Read_Static(thisHydrogenAtomicDensity,"datum",self%densityHydrogenAtomic(:,iChemicalState+1))
       if (self%gotHydrogenCation) call XML_Array_Read_Static(thisHydrogenCationDensity,"datum",self%densityHydrogenCation(:,iChemicalState+1))
       if (iChemicalState == 0) then
          ! Copy the temperatures so we can check subsequent temperature reads for consistency.
          deallocate(temperaturesReference)
          temperaturesReference=self%temperatures
       else
          ! Check that temperature grids are aligned.
          if (any(Values_Differ(self%temperatures,temperaturesReference,relTol=1.0d-6))) &
               & call Galacticus_Error_Report('cieFileReadFile','temperature grids mismatch')
       end if
    end do
    deallocate(temperaturesReference)
    where (self%metallicities > metallicityLogarithmicZero)
       self%metallicities=10.0d0**self%metallicities
    elsewhere
       self%metallicities= 0.0d0
    end where
    ! Extract extrapolation methods from the file.
    extrapolationElement         => XML_Get_First_Element_By_Tag_Name(doc,"extrapolation")
    metallicityExtrapolationList => getElementsByTagname(extrapolationElement,"metallicity")
    do iExtrapolation=0,getLength(metallicityExtrapolationList)-1
       extrapolation => item(metallicityExtrapolationList,iExtrapolation)
       call XML_Extrapolation_Element_Decode(extrapolation,limitType,extrapolationMethod,allowedMethods=[extrapolationTypeZero,extrapolationTypeFix,extrapolationTypePowerLaw])
       select case (trim(limitType))
       case ('low')
          self%extrapolateMetallicityLow =extrapolationMethod
       case ('high')
          self%extrapolateMetallicityHigh=extrapolationMethod
       case default
          call Galacticus_Error_Report('cieFileReadFile','unrecognized extrapolation limit')
       end select
    end do
    temperatureExtrapolationList => getElementsByTagname(extrapolationElement,"temperature")
    do iExtrapolation=0,getLength(temperatureExtrapolationList)-1
       extrapolation => item(temperatureExtrapolationList,iExtrapolation)
       call XML_Extrapolation_Element_Decode(extrapolation,limitType,extrapolationMethod,allowedMethods=[extrapolationTypeZero,extrapolationTypeFix,extrapolationTypePowerLaw])
       select case (trim(limitType))
       case ('low')
          self%extrapolateTemperatureLow =extrapolationMethod
       case ('high')
          self%extrapolateTemperatureHigh=extrapolationMethod
       case default
          call Galacticus_Error_Report('cieFileReadFile','unrecognized extrapolation limit')
       end select
    end do
    ! Destroy the document.
    call destroy(doc)
    call Galacticus_Display_Unindent('done',verbosityDebug)
    !$omp end critical (FoX_DOM_Access)
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
       if     (                                                              &
            &   self%extrapolateTemperatureLow  == extrapolationTypePowerLaw &
            &  .or.                                                          &
            &   self%extrapolateTemperatureHigh == extrapolationTypePowerLaw &
            &  .or.                                                          &
            &   self%extrapolateMetallicityLow  == extrapolationTypePowerLaw &
            &  .or.                                                          &
            &   self%extrapolateMetallicityHigh == extrapolationTypePowerLaw &
            & )                                                              &
            & call Galacticus_Error_Report('cieFileReadFile','power law extrapolation allowed only in loggable tables')
    end if
    ! Force interpolation accelerators to be reset.
    self%resetTemperature=.true.
    self%resetMetallicity=.true.
    ! Get chemical indices.
    self%electronChemicalIndex               =Chemicals_Index("Electron"            )
    if (self%gotHydrogenAtomic) then
       self%atomicHydrogenChemicalIndex      =Chemicals_Index("AtomicHydrogen"      )
    else
       self%atomicHydrogenChemicalIndex      =-1
    end if
    if (self%gotHydrogenCation) then
       self%atomicHydrogenCationChemicalIndex=Chemicals_Index("AtomicHydrogenCation")
    else
       self%atomicHydrogenCationChemicalIndex=-1
    end if
    return
  end subroutine cieFileReadFile

  subroutine cieFileDescriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    use FoX_DOM
    implicit none
    class(chemicalStateCIEFile), intent(inout) :: self
    type (inputParameters     ), intent(inout) :: descriptor
    type (inputParameters     )                :: subParameters

    call descriptor%addParameter("chemicalStateMethod","cieFile")
    subParameters=descriptor%subparameters("chemicalStateMethod")
    call subParameters%addParameter("fileName",char(self%fileName))
    return
  end subroutine cieFileDescriptor
