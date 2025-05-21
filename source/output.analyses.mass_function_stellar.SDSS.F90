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
Implements a stellar mass function output analysis class.
!!}


  !![
  <outputAnalysis name="outputAnalysisMassFunctionStellarSDSS">
   <description>
    An SDSS stellar mass function output analysis class, which utilizes the stellar mass function for $z\approx 0.07$ galaxies
    measured by \cite{li_distribution_2009} from the \gls{sdss}.
    
    Given a \glc\ model, total stellar masses of model galaxies are adjusted using:
    \begin{equation}
     M_\star \rightarrow \mathbf{C} \mathbf{L} \mathbf{G} \mathbf{S} M_\star,
    \end{equation}
    where the $\mathbf{S}$ operator is a multiplicative factor accounting for systematic errors in stellar mass determination and is
    equal to \citep{behroozi_comprehensive_2010}
    \begin{equation}
     \log_\mathrm{10} S = \sum_{i=0}^N s_i \log_\mathrm{10}^i \left({M_\star \over 10^{11.3}M_\odot}\right),
    \end{equation}
    where $s=${\normalfont \ttfamily [systematicErrorPolynomialCoefficient]}, the {\normalfont \bfseries G} operator is a
    multiplicative factor drawn from a log-normal distribution of width $\sigma(M)$~dex for each galaxy to mimic the effects of random
    errors on stellar masses (motivated by the discussion of \cite{behroozi_comprehensive_2010}), the {\normalfont \bfseries L}
    operator accounts for gravitational lensing, and the {\normalfont \bfseries C} operator accounts for the difference between model
    and observed cosmologies. The random error model is given by:
    \begin{equation}
     \sigma(M) = \hbox{min}\left[\sigma_\mathrm{max},\hbox{max}\left[\sigma_\mathrm{min},\sum_{i=0}^N r_i \log_\mathrm{10}^i \left({M_\star \over 10^{11.3}M_\odot}\right)\right]\right],
    \end{equation}
    where $r=${\normalfont \ttfamily [randomErrorPolynomialCoefficient]}, $\sigma_\mathrm{min}$={\normalfont \ttfamily
      [randomErrorMinimum]}, and $\sigma_\mathrm{max}$={\normalfont \ttfamily [randomErrorMaximum]}.
    
    The model masses are then used to construct a mass function by binning into a histogram using the masses reported by
    \cite{li_distribution_2009} (modified as described above) as the centers of the bins (with bin boundaries placed at the geometric
    means of consecutive bin centers).
   </description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisMassFunctionStellar) :: outputAnalysisMassFunctionStellarSDSS
     !!{
     An SDSS stellar mass function output analysis class.
     !!}
     private
     class           (gravitationalLensingClass), pointer                   :: gravitationalLensing_            => null()
     double precision                           , allocatable, dimension(:) :: randomErrorPolynomialCoefficient          , systematicErrorPolynomialCoefficient
     double precision                                                       :: randomErrorMinimum                        , randomErrorMaximum                  , &
          &                                                                    sizeSourceLensing
   contains
     final :: massFunctionStellarSDSSDestructor
  end type outputAnalysisMassFunctionStellarSDSS

  interface outputAnalysisMassFunctionStellarSDSS
     !!{
     Constructors for the \refClass{outputAnalysisMassFunctionStellarSDSS} output analysis class.
     !!}
     module procedure massFunctionStellarSDSSConstructorParameters
     module procedure massFunctionStellarSDSSConstructorInternal
  end interface outputAnalysisMassFunctionStellarSDSS

contains

  function massFunctionStellarSDSSConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisMassFunctionStellarSDSS} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisMassFunctionStellarSDSS)                              :: self
    type            (inputParameters                      ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass              ), pointer                     :: cosmologyFunctions_
    class           (outputTimesClass                     ), pointer                     :: outputTimes_
    class           (gravitationalLensingClass            ), pointer                     :: gravitationalLensing_
    double precision                                       , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient , systematicErrorPolynomialCoefficient
    integer                                                                              :: covarianceBinomialBinsPerDecade
    double precision                                                                     :: covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum   , &
         &                                                                                  randomErrorMinimum               , randomErrorMaximum                  , &
         &                                                                                  sizeSourceLensing

    ! Check and read parameters.
    if (parameters%isPresent(    'randomErrorPolynomialCoefficient')) then
       allocate(    randomErrorPolynomialCoefficient(parameters%count(    'randomErrorPolynomialCoefficient')))
    else
       allocate(    randomErrorPolynomialCoefficient(1                                                   ))
    end if
    if (parameters%isPresent('systematicErrorPolynomialCoefficient')) then
       allocate(systematicErrorPolynomialCoefficient(parameters%count('systematicErrorPolynomialCoefficient')))
    else
       allocate(systematicErrorPolynomialCoefficient(1                                                   ))
    end if
    !![
    <inputParameter>
      <name>randomErrorMinimum</name>
      <source>parameters</source>
      <variable>randomErrorMinimum</variable>
      <defaultValue>0.07d0</defaultValue>
      <description>The minimum random error for SDSS stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorMaximum</name>
      <source>parameters</source>
      <variable>randomErrorMaximum</variable>
      <defaultValue>0.07d0</defaultValue>
      <description>The minimum random error for SDSS stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>randomErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>randomErrorPolynomialCoefficient</variable>
      <defaultValue>[0.07d0]</defaultValue>
      <description>The coefficients of the random error polynomial for SDSS stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>systematicErrorPolynomialCoefficient</name>
      <source>parameters</source>
      <variable>systematicErrorPolynomialCoefficient</variable>
      <defaultValue>[0.0d0]</defaultValue>
      <description>The coefficients of the systematic error polynomial for SDSS stellar masses.</description>
    </inputParameter>
    <inputParameter>
      <name>sizeSourceLensing</name>
      <source>parameters</source>
      <variable>sizeSourceLensing</variable>
      <defaultValue>2.0d-3</defaultValue>
      <description>The characteristic source size for gravitational lensing calculations.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialBinsPerDecade</name>
      <source>parameters</source>
      <variable>covarianceBinomialBinsPerDecade</variable>
      <defaultValue>10</defaultValue>
      <description>The number of bins per decade of halo mass to use when constructing SDSS stellar mass function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMinimum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMinimum</variable>
      <defaultValue>1.0d8</defaultValue>
      <description>The minimum halo mass to consider when constructing SDSS stellar mass function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>covarianceBinomialMassHaloMaximum</name>
      <source>parameters</source>
      <variable>covarianceBinomialMassHaloMaximum</variable>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum halo mass to consider when constructing SDSS stellar mass function covariance matrices for main branch galaxies.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    <objectBuilder class="outputTimes"          name="outputTimes_"          source="parameters"/>
    <objectBuilder class="gravitationalLensing" name="gravitationalLensing_" source="parameters"/>
    !!]
    ! Build the object.
    self=outputAnalysisMassFunctionStellarSDSS(cosmologyFunctions_,gravitationalLensing_,outputTimes_,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"  />
    <objectDestructor name="outputTimes_"         />
    <objectDestructor name="gravitationalLensing_"/>
    !!]
    return
  end function massFunctionStellarSDSSConstructorParameters

  function massFunctionStellarSDSSConstructorInternal(cosmologyFunctions_,gravitationalLensing_,outputTimes_,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisMassFunctionStellarSDSS} output analysis class for internal use.
    !!}
    use :: Cosmology_Functions                   , only : cosmologyFunctionsClass                        , cosmologyFunctionsMatterLambda
    use :: Cosmology_Parameters                  , only : cosmologyParametersSimple
    use :: Galactic_Filters                      , only : galacticFilterStellarMass
    use :: Input_Paths                           , only : inputPath                                      , pathTypeDataStatic
    use :: Geometry_Surveys                      , only : surveyGeometryLiWhite2009SDSS
    use :: Gravitational_Lensing                 , only : gravitationalLensingClass
    use :: Output_Analysis_Distribution_Operators, only : distributionOperatorList                       , outputAnalysisDistributionOperatorGrvtnlLnsng, outputAnalysisDistributionOperatorRandomErrorPlynml, outputAnalysisDistributionOperatorSequence
    use :: Output_Analysis_Property_Operators    , only : outputAnalysisPropertyOperatorSystmtcPolynomial
    implicit none
    type            (outputAnalysisMassFunctionStellarSDSS              )                              :: self
    class           (cosmologyFunctionsClass                            ), intent(in   ), target       :: cosmologyFunctions_
    class           (outputTimesClass                                   ), intent(inout), target       :: outputTimes_
    class           (gravitationalLensingClass                          ), intent(in   ), target       :: gravitationalLensing_
    double precision                                                     , intent(in   )               :: randomErrorMinimum                                         , randomErrorMaximum                  , &
         &                                                                                                sizeSourceLensing
    double precision                                                     , intent(in   ), dimension(:) :: randomErrorPolynomialCoefficient                           , systematicErrorPolynomialCoefficient
    integer                                                              , intent(in   )               :: covarianceBinomialBinsPerDecade
    double precision                                                     , intent(in   )               :: covarianceBinomialMassHaloMinimum                          , covarianceBinomialMassHaloMaximum
    type            (galacticFilterStellarMass                          )               , pointer      :: galacticFilter_
    type            (surveyGeometryLiWhite2009SDSS                      )               , pointer      :: surveyGeometry_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    )               , pointer      :: outputAnalysisPropertyOperator_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml)               , pointer      :: outputAnalysisDistributionOperatorRandomErrorPlynml_
    type            (outputAnalysisDistributionOperatorGrvtnlLnsng      )               , pointer      :: outputAnalysisDistributionOperatorGrvtnlLnsng_
    type            (outputAnalysisDistributionOperatorSequence         )               , pointer      :: outputAnalysisDistributionOperator_
    type            (cosmologyParametersSimple                          )               , pointer      :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                     )               , pointer      :: cosmologyFunctionsData
    type            (distributionOperatorList                           )               , pointer      :: distributionOperatorSequence
    double precision                                                     , parameter                   :: errorPolynomialZeroPoint                            =11.3d+0
    !![
    <constructorAssign variables="randomErrorPolynomialCoefficient, systematicErrorPolynomialCoefficient, randomErrorMinimum, randomErrorMaximum, sizeSourceLensing, *gravitationalLensing_"/>
    !!]

    ! Build a filter which select galaxies with stellar mass 10⁶M☉ or greater.
    allocate(galacticFilter_)
    !![
    <referenceConstruct object="galacticFilter_" constructor="galacticFilterStellarMass(massThreshold=1.0d6)"/>
    !!]
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    !![
    <referenceConstruct object="cosmologyParametersData">
     <constructor>
      cosmologyParametersSimple(                            &amp;
        &amp;                   OmegaMatter    = 0.30000d0, &amp;
        &amp;                   OmegaDarkEnergy= 0.70000d0, &amp;
        &amp;                   HubbleConstant =70.00000d0, &amp;
        &amp;                   temperatureCMB = 2.72548d0, &amp;
        &amp;                   OmegaBaryon    = 0.04550d0  &amp;
        &amp;                  )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="cosmologyFunctionsData">
    <constructor>
      cosmologyFunctionsMatterLambda(                        &amp;
        &amp;                        cosmologyParametersData &amp;
        &amp;                       )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build the SDSS survey geometry of Li & White (2009) with their imposed redshift limits.
    allocate(surveyGeometry_)
    !![
    <referenceConstruct object="surveyGeometry_"                 constructor="surveyGeometryLiWhite2009SDSS(redshiftMinimum=1.0d-3,redshiftMaximum=5000.0d0,cosmologyFunctions_=cosmologyFunctions_)"/>
    !!]
    ! Create property operators.
    !! Systematic error model.
    allocate(outputAnalysisPropertyOperator_    )
    !![
    <referenceConstruct object="outputAnalysisPropertyOperator_" constructor="outputAnalysisPropertyOperatorSystmtcPolynomial(errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient)"/>
    !!]
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperatorRandomErrorPlynml_)
    !![
    <referenceConstruct object="outputAnalysisDistributionOperatorRandomErrorPlynml_">
     <constructor>
      outputAnalysisDistributionOperatorRandomErrorPlynml(                                  &amp;
        &amp;                                             randomErrorMinimum              , &amp;
        &amp;                                             randomErrorMaximum              , &amp;
        &amp;                                             errorPolynomialZeroPoint        , &amp;
        &amp;                                             randomErrorPolynomialCoefficient  &amp;
        &amp;                                            )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build a gravitational lensing distribution operator.
    allocate(outputAnalysisDistributionOperatorGrvtnlLnsng_)
    !![
    <referenceConstruct object="outputAnalysisDistributionOperatorGrvtnlLnsng_">
     <constructor>
      outputAnalysisDistributionOperatorGrvtnlLnsng(                       &amp;
        &amp;                                       gravitationalLensing_, &amp;
        &amp;                                       outputTimes_         , &amp;
        &amp;                                       sizeSourceLensing      &amp;
        &amp;                                      )
     </constructor>
    </referenceConstruct>
    !!]
    ! Construct sequence distribution operator.
    allocate(distributionOperatorSequence                )
    allocate(distributionOperatorSequence           %next)
    allocate(outputAnalysisDistributionOperator_     )
    distributionOperatorSequence            %operator_   => outputAnalysisDistributionOperatorRandomErrorPlynml_
    distributionOperatorSequence       %next%operator_   => outputAnalysisDistributionOperatorGrvtnlLnsng_
    !![
    <referenceConstruct object="outputAnalysisDistributionOperator_">
     <constructor>
      outputAnalysisDistributionOperatorSequence(                             &amp;
        &amp;                                    distributionOperatorSequence &amp;
        &amp;                                   )
     </constructor>
    </referenceConstruct>
    !!]
    ! Build the object.
    self%outputAnalysisMassFunctionStellar=                                                                                                                      &
         & outputAnalysisMassFunctionStellar(                                                                                                                    &
         &                                   var_str('LiWhite2009SDSS'                                               )                                         , &
         &                                   var_str('Stellar mass function for the Li \& White (2009) SDSS analysis')                                         , &
         &                                   char(inputPath(pathTypeDataStatic)//'/observations/massFunctionsStellar/Stellar_Mass_Function_Li_White_2009.hdf5'), &
         &                                   galacticFilter_                                                                                                   , &
         &                                   surveyGeometry_                                                                                                   , &
         &                                   cosmologyFunctions_                                                                                               , &
         &                                   cosmologyFunctionsData                                                                                            , &
         &                                   outputAnalysisPropertyOperator_                                                                                   , &
         &                                   outputAnalysisDistributionOperator_                                                                               , &
         &                                   outputTimes_                                                                                                      , &
         &                                   covarianceBinomialBinsPerDecade                                                                                   , &
         &                                   covarianceBinomialMassHaloMinimum                                                                                 , &
         &                                   covarianceBinomialMassHaloMaximum                                                                                   &
         &                                  )
    ! Clean up.
    !![
    <objectDestructor name="surveyGeometry_"                                     />
    <objectDestructor name="galacticFilter_"                                     />
    <objectDestructor name="cosmologyParametersData"                             />
    <objectDestructor name="cosmologyFunctionsData"                              />
    <objectDestructor name="outputAnalysisPropertyOperator_"                     />
    <objectDestructor name="outputAnalysisDistributionOperator_"                 />
    <objectDestructor name="outputAnalysisDistributionOperatorGrvtnlLnsng_"      />
    <objectDestructor name="outputAnalysisDistributionOperatorRandomErrorPlynml_"/>
    !!]
    nullify(distributionOperatorSequence)
    return
  end function massFunctionStellarSDSSConstructorInternal

  subroutine massFunctionStellarSDSSDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisMassFunctionStellarSDSS} output analysis class.
    !!}
    implicit none
    type(outputAnalysisMassFunctionStellarSDSS), intent(inout) :: self

    !![
    <objectDestructor name="self%gravitationalLensing_"/>
    !!]
    return
  end subroutine massFunctionStellarSDSSDestructor
