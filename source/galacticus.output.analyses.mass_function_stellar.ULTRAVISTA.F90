!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements a stellar mass function output analysis class for the ULTRAVISTA survey of \cite{muzzin_evolution_2013}.

  use Gravitational_Lensing
  
  !# <outputAnalysis name="outputAnalysisMassFunctionStellarULTRAVISTA" defaultThreadPrivate="yes">
  !#  <description>An ULTRAVISTA stellar mass function output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisMassFunctionStellar) :: outputAnalysisMassFunctionStellarULTRAVISTA
     !% An ULTRAVISTA stellar mass function output analysis class.
     private
  end type outputAnalysisMassFunctionStellarULTRAVISTA

  interface outputAnalysisMassFunctionStellarULTRAVISTA
     !% Constructors for the ``massFunctionStellarULTRAVISTA'' output analysis class.
     module procedure massFunctionStellarULTRAVISTAConstructorParameters
     module procedure massFunctionStellarULTRAVISTAConstructorInternal
  end interface outputAnalysisMassFunctionStellarULTRAVISTA

contains

  function massFunctionStellarULTRAVISTAConstructorParameters(parameters) result (self)
    !% Constructor for the ``massFunctionStellarULTRAVISTA'' output analysis class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (outputAnalysisMassFunctionStellarULTRAVISTA)                              :: self
    type            (inputParameters                            ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass                    ), pointer                     :: cosmologyFunctions_
    class           (gravitationalLensingClass                  ), pointer                     :: gravitationalLensing_
    double precision                                             , allocatable  , dimension(:) :: randomErrorPolynomialCoefficient , systematicErrorPolynomialCoefficient
    integer                                                                                    :: covarianceBinomialBinsPerDecade  , redshiftInterval
    double precision                                                                           :: covarianceBinomialMassHaloMinimum, covarianceBinomialMassHaloMaximum   , &
         &                                                                                        randomErrorMinimum               , randomErrorMaximum                  , &
         &                                                                                        sizeSourceLensing

    ! Check and read parameters.
    if (parameters%isPresent(    'randomErrorPolynomialCoefficient')) then
       allocate(    randomErrorPolynomialCoefficient(parameters%count(    'randomErrorPolynomialCoefficient')))
    else
       allocate(    randomErrorPolynomialCoefficient(1                                                       ))
    end if
    if (parameters%isPresent('systematicErrorPolynomialCoefficient')) then
       allocate(systematicErrorPolynomialCoefficient(parameters%count('systematicErrorPolynomialCoefficient')))
    else
       allocate(systematicErrorPolynomialCoefficient(1                                                       ))
    end if
    !# <inputParameter>
    !#   <name>redshiftInterval</name>
    !#   <source>parameters</source>
    !#   <variable>redshiftInterval</variable>
    !#   <description>The redshift interval (0-6) to use.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomErrorMinimum</name>
    !#   <source>parameters</source>
    !#   <variable>randomErrorMinimum</variable>
    !#   <defaultValue>0.1d0</defaultValue>
    !#   <description>The minimum random error for ULTRAVISTA stellar masses.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomErrorMaximum</name>
    !#   <source>parameters</source>
    !#   <variable>randomErrorMaximum</variable>
    !#   <defaultValue>0.1d0</defaultValue>
    !#   <description>The minimum random error for ULTRAVISTA stellar masses.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomErrorPolynomialCoefficient</name>
    !#   <source>parameters</source>
    !#   <variable>randomErrorPolynomialCoefficient</variable>
    !#   <defaultValue>[0.1d0]</defaultValue>
    !#   <description>The coefficients of the random error polynomial for ULTRAVISTA stellar masses.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>systematicErrorPolynomialCoefficient</name>
    !#   <source>parameters</source>
    !#   <variable>systematicErrorPolynomialCoefficient</variable>
    !#   <defaultValue>[0.0d0]</defaultValue>
    !#   <description>The coefficients of the systematic error polynomial for ULTRAVISTA stellar masses.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>sizeSourceLensing</name>
    !#   <source>parameters</source>
    !#   <variable>sizeSourceLensing</variable>
    !#   <defaultValue>2.0d-3</defaultValue>
    !#   <description>The characteristic source size for gravitational lensing calculations.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialBinsPerDecade</name>
    !#   <source>parameters</source>
    !#   <variable>covarianceBinomialBinsPerDecade</variable>
    !#   <defaultValue>10</defaultValue>
    !#   <description>The number of bins per decade of halo mass to use when constructing ULTRAVISTA stellar mass function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMinimum</name>
    !#   <source>parameters</source>
    !#   <variable>covarianceBinomialMassHaloMinimum</variable>
    !#   <defaultValue>1.0d8</defaultValue>
    !#   <description>The minimum halo mass to consider when constructing ULTRAVISTA stellar mass function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>covarianceBinomialMassHaloMaximum</name>
    !#   <source>parameters</source>
    !#   <variable>covarianceBinomialMassHaloMaximum</variable>
    !#   <defaultValue>1.0d16</defaultValue>
    !#   <description>The maximum halo mass to consider when constructing ULTRAVISTA stellar mass function covariance matrices for main branch galaxies.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    !# <objectBuilder class="gravitationalLensing" name="gravitationalLensing_" source="parameters"/>
    ! Build the object.
    self=outputAnalysisMassFunctionStellarULTRAVISTA(cosmologyFunctions_,gravitationalLensing_,redshiftInterval,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing)
    !# <inputParametersValidate source="parameters"/>
    return
  end function massFunctionStellarULTRAVISTAConstructorParameters

  function massFunctionStellarULTRAVISTAConstructorInternal(cosmologyFunctions_,gravitationalLensing_,redshiftInterval,randomErrorMinimum,randomErrorMaximum,randomErrorPolynomialCoefficient,systematicErrorPolynomialCoefficient,covarianceBinomialBinsPerDecade,covarianceBinomialMassHaloMinimum,covarianceBinomialMassHaloMaximum,sizeSourceLensing) result (self)
    !% Constructor for the ``massFunctionStellarULTRAVISTA'' output analysis class for internal use.
    use Input_Parameters
    use Galacticus_Input_Paths
    use Galacticus_Error
    use Output_Analysis_Distribution_Operators
    use Cosmology_Parameters
    use ISO_Varying_String
    use String_Handling
    implicit none
    type            (outputAnalysisMassFunctionStellarULTRAVISTA        )                              :: self
    class           (cosmologyFunctionsClass                            ), intent(in   ), target       :: cosmologyFunctions_
    class           (gravitationalLensingClass                          ), intent(in   ), target       :: gravitationalLensing_
    integer                                                              , intent(in   )               :: redshiftInterval
    double precision                                                     , intent(in   )               :: randomErrorMinimum                                         , randomErrorMaximum                  , &
         &                                                                                                sizeSourceLensing
    double precision                                                     , intent(in   ), dimension(:) :: randomErrorPolynomialCoefficient                           , systematicErrorPolynomialCoefficient
    integer                                                              , intent(in   )               :: covarianceBinomialBinsPerDecade
    double precision                                                     , intent(in   )               :: covarianceBinomialMassHaloMinimum                          , covarianceBinomialMassHaloMaximum
    type            (galacticFilterStellarMass                          )               , pointer      :: galacticFilter_
    type            (surveyGeometryMuzzin2013ULTRAVISTA                 )               , pointer      :: surveyGeometry_
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial    )               , pointer      :: outputAnalysisPropertyOperator_
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml)               , pointer      :: outputAnalysisDistributionOperatorRandomErrorPlynml_
    type            (outputAnalysisDistributionOperatorGrvtnlLnsng      )               , pointer      :: outputAnalysisDistributionOperatorGrvtnlLnsng_
    type            (outputAnalysisDistributionOperatorSequence         )               , pointer      :: outputAnalysisDistributionOperator_
    type            (cosmologyParametersSimple                          )               , pointer      :: cosmologyParametersData
    type            (cosmologyFunctionsMatterLambda                     )               , pointer      :: cosmologyFunctionsData
    type            (distributionOperatorList                           )               , pointer      :: distributionOperatorSequence
    double precision                                                     , parameter                   :: errorPolynomialZeroPoint                            =11.3d+0
    type            (varying_string                                     )                              :: fileName
    double precision                                                                                   :: massThreshold
    
    ! Create cosmological model in which data were analyzed.
    allocate(cosmologyParametersData)
    allocate(cosmologyFunctionsData )
    cosmologyParametersData=cosmologyParametersSimple     (                            &
         &                                                 OmegaMatter    = 0.30000d0, &
         &                                                 OmegaDarkEnergy= 0.70000d0, &
         &                                                 HubbleConstant =70.00000d0, &
         &                                                 temperatureCMB = 2.72548d0, &
         &                                                 OmegaBaryon    = 0.04550d0  &
         &                                                )
    cosmologyFunctionsData =cosmologyFunctionsMatterLambda(                            &
         &                                                 cosmologyParametersData     &
         &                                                )
    ! Determine the data file and mass threshold to use.
    select case (redshiftInterval)
    case (0)
       fileName     ='Stellar_Mass_Function_ULTRAVISTA_2013_z0.20_0.50.hdf5'
       massThreshold=3.0d6
    case (1)
       fileName     ='Stellar_Mass_Function_ULTRAVISTA_2013_z0.50_1.00.hdf5'
       massThreshold=3.0d7
    case (2)
       fileName     ='Stellar_Mass_Function_ULTRAVISTA_2013_z1.00_1.50.hdf5'
       massThreshold=3.0d7
    case (3)
       fileName     ='Stellar_Mass_Function_ULTRAVISTA_2013_z1.50_2.00.hdf5'
       massThreshold=3.0d8
    case (4)
       fileName     ='Stellar_Mass_Function_ULTRAVISTA_2013_z2.00_2.50.hdf5'
       massThreshold=3.0d8
    case (5)
       fileName     ='Stellar_Mass_Function_ULTRAVISTA_2013_z2.50_3.00.hdf5'
       massThreshold=1.0d9
    case (6)
       fileName     ='Stellar_Mass_Function_ULTRAVISTA_2013_z3.00_4.00.hdf5'
       massThreshold=1.0d9
    case default
       call Galacticus_Error_Report('0 ≤ redshiftInterval ≤ 6 is required'//{introspection:location})
    end select
    ! Build a filter which select galaxies with stellar mass 3⨉10⁶M☉ or greater.
    allocate(galacticFilter_)
    galacticFilter_=galacticFilterStellarMass(massThreshold=3.0d6)
    ! Build the ULTRAVISTA survey geometry of Muzzin et al. (2013) with their imposed redshift limits.
    allocate(surveyGeometry_)
    surveyGeometry_=surveyGeometryMuzzin2013ULTRAVISTA(redshiftBin=redshiftInterval,cosmologyFunctions_=cosmologyFunctions_)
    ! Create property operators.
    !! Systematic error model.
    allocate(outputAnalysisPropertyOperator_    )
    outputAnalysisPropertyOperator_    =outputAnalysisPropertyOperatorSystmtcPolynomial(errorPolynomialZeroPoint,systematicErrorPolynomialCoefficient)
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperatorRandomErrorPlynml_)
    outputAnalysisDistributionOperatorRandomErrorPlynml_ =  outputAnalysisDistributionOperatorRandomErrorPlynml (                                  &
         &                                                                                                       randomErrorMinimum              , &
         &                                                                                                       randomErrorMaximum              , &
         &                                                                                                       errorPolynomialZeroPoint        , &
         &                                                                                                       randomErrorPolynomialCoefficient  &
         &                                                                                                      )
    ! Build a gravitational lensing distribution operator.
    allocate(outputAnalysisDistributionOperatorGrvtnlLnsng_)
    outputAnalysisDistributionOperatorGrvtnlLnsng_       =  outputAnalysisDistributionOperatorGrvtnlLnsng       (                                  &
         &                                                                                                       gravitationalLensing_           , &
         &                                                                                                       sizeSourceLensing                 &
         &                                                                                                      )
    ! Construct sequence distribution operator.
    allocate(distributionOperatorSequence                )
    allocate(distributionOperatorSequence           %next)
    allocate(outputAnalysisDistributionOperator_     )
    distributionOperatorSequence            %operator_   => outputAnalysisDistributionOperatorRandomErrorPlynml_
    distributionOperatorSequence       %next%operator_   => outputAnalysisDistributionOperatorGrvtnlLnsng_
    outputAnalysisDistributionOperator_                  =  outputAnalysisDistributionOperatorSequence          (                                  &
         &                                                                                                       distributionOperatorSequence      &
         &                                                                                                      )
    ! Build the object.
    self%outputAnalysisMassFunctionStellar=                                                                                                                                  &
         & outputAnalysisMassFunctionStellar(                                                                                                                                &
         &                                   var_str('Muzzin2013ULTRAVISTAz')//redshiftInterval                                                                            , &
         &                                   var_str('Stellar mass function for the Muzzin et al. (2013) ULTRAVISTA analysis')                                             , &
         &                                   char(Galacticus_Input_Path()//'/data/observations/massFunctionsStellar/'//fileName), &
         &                                   galacticFilter_                                                                                                               , &
         &                                   surveyGeometry_                                                                                                               , &
         &                                   cosmologyFunctions_                                                                                                           , &
         &                                   cosmologyFunctionsData                                                                                                        , &
         &                                   outputAnalysisPropertyOperator_                                                                                               , &
         &                                   outputAnalysisDistributionOperator_                                                                                           , &
         &                                   covarianceBinomialBinsPerDecade                                                                                               , &
         &                                   covarianceBinomialMassHaloMinimum                                                                                             , &
         &                                   covarianceBinomialMassHaloMaximum                                                                                               &
         &                                  )
    ! Clean up.
    nullify(surveyGeometry_                                     )
    nullify(galacticFilter_                                     )
    nullify(cosmologyParametersData                             )
    nullify(cosmologyFunctionsData                              )
    nullify(outputAnalysisDistributionOperator_                 )
    nullify(outputAnalysisDistributionOperatorGrvtnlLnsng_      )
    nullify(outputAnalysisDistributionOperatorRandomErrorPlynml_)
    nullify(distributionOperatorSequence                        )
    return
  end function massFunctionStellarULTRAVISTAConstructorInternal
