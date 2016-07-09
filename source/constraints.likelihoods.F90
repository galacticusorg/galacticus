!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements likelihoods for use when constraining \glc.

module Constraints_Likelihoods
  !% Implements likelihoods for use when constraining \glc.
  use Linear_Algebra
  use Constraints_State
  use Constraints_Convergence
  use Constraints_Mappings
  use ISO_Varying_String
  use Pseudo_Random
  use FGSL
  use Nearest_Neighbors
  private
  public :: likelihoodNew

  ! Define the basic likelihood class.
  type, abstract, public :: likelihood
   contains
     !@ <objectMethods>
     !@   <object>likelihood</object>
     !@   <objectMethod>
     !@     <method>evaluate</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless class(state)\textgreater} simulationState\argin, \doublezero\ temperature\argin, \doublezero\ logLikelihoodCurrent\argin, \doublezero\ logLikelihoodCurrent\argin, \doublezero\ logPriorProposed\argin, \doublezero\ timeEvaluate\arginout, \doublezero\ [logLikelihoodVariance]\argout</arguments>
     !@     <description>Evaluate the model likelihood at the given {\normalfont \ttfamily simulationState} and return the log-likelihood.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>willEvaluate</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\textcolor{red}{\textless class(state)\textgreater} simulationState\argin, \doublezero\ temperature\argin, \doublezero\ logLikelihoodCurrent\argin, \doublezero\ logLikelihoodCurrent\argin, \doublezero\ logPriorProposed\argin</arguments>
     !@     <description>Return true if the likelihood will evaluate the model likelihood at the given {\normalfont \ttfamily simulationState}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>functionChanged</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Informs the likelihood object that the likelihood function may have changed.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure(likelihoodEvaluate       ), deferred :: evaluate
     procedure                                      :: willEvaluate    => likelihoodWillEvaluate
     procedure(likelihoodFunctionChanged), deferred :: functionChanged
  end type likelihood

  ! Interface for deferred functions.
  abstract interface
     double precision function likelihoodEvaluate(self,simulationState,parameterMappings,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance)
       import :: likelihood, state, convergence, mappingList
       class           (likelihood ), intent(inout)               :: self
       class           (state      ), intent(in   )               :: simulationState
       type            (mappingList), intent(in   ), dimension(:) :: parameterMappings
       class           (convergence), intent(inout)               :: simulationConvergence
       double precision             , intent(in   )               :: temperature          , logLikelihoodCurrent, &
            &                                                        logPriorCurrent      , logPriorProposed
       real                         , intent(inout)               :: timeEvaluate
       double precision             , intent(  out), optional     :: logLikelihoodVariance
     end function likelihoodEvaluate
  end interface
  abstract interface
     subroutine likelihoodFunctionChanged(self)
       import :: likelihood
       class(likelihood ), intent(inout) :: self
     end subroutine likelihoodFunctionChanged
  end interface

  ! Include all likelihood types.
  include 'constraints.likelihoods.multivariate_normal.type.inc'
  include 'constraints.likelihoods.multivariate_normal.stochastic.type.inc'
  include 'constraints.likelihoods.Galacticus.type.inc'
  include 'constraints.likelihoods.gaussian_regression.type.inc'
  include 'constraints.likelihoods.mass_function.type.inc'
  include 'constraints.likelihoods.halo_mass_function.type.inc'
  include 'constraints.likelihoods.projected_correlation_function.type.inc'
  include 'constraints.likelihoods.SED_fit.type.inc'
  include 'constraints.likelihoods.posterior_as_prior.type.inc'

contains

  function likelihoodNew(definition,configFileName) result (newLikelihood)
    !% Create a new likelihood from an XML definition.
    use FoX_DOM
    use IO_XML
    use ISO_Varying_String
    use Galacticus_Error
    use String_Handling
    implicit none
    class           (likelihood    ), pointer                     :: newLikelihood
    type            (node          ), pointer    , intent(in   )  :: definition
    type            (varying_string), optional   , intent(in   )  :: configFileName
    type            (node          ), pointer                     :: likelihoodMeanDefinition                       , likelihoodCovarianceDefinition                          , &
         &                                                           covarianceRow                                  , likelihoodRealizationCountDefinition                    , &
         &                                                           likelihoodRealizationCountMinimumDefinition    , likelihoodDefinition                                    , &
         &                                                           likelihoodEmulatorRebuildCountDefinition       , likelihoodPolynomialOrderDefinition                     , &
         &                                                           likelihoodSigmaBufferDefinition                , likelihoodLogLikelihoodBufferDefinition                 , &
         &                                                           likelihoodLogLikelihoodErrorToleranceDefinition, likelihoodReportCountDefinition                         , &
         &                                                           likelihoodEmulateOutliersDefinition            , likelihoodMassRangeMinimumDefinition                    , &
         &                                                           likelihoodHaloMassMinimumDefinition            , likelihoodHaloMassMaximumDefinition                     , &
         &                                                           likelihoodRedshiftMinimumDefinition            , likelihoodRedshiftMaximumDefinition                     , &
         &                                                           likelihoodUseSurveyLimitsDefinition            , likelihoodMassFunctionFileNameDefinition                , &
         &                                                           likelihoodModelSurfaceBrightnessDefinition     , likelihoodSurfaceBrightnessLimitDefinition              , &
         &                                                           likelihoodChainBaseNameDefinition              , likelihoodToleranceDefinition                           , &
         &                                                           likelihoodNeighborCountDefinition              , likelihoodProjectedCorrelationFunctionFileNameDefinition, &
         &                                                           likelihoodLineOfSightDepthDefinition           , likelihoodHalfIntegralDefinition                        , &
         &                                                           likelihoodExclusionsDefinition                 , likelihoodDumpEmulatorDefinition                        , &
         &                                                           likelihoodDelayIntervalDefinition              , likelihoodDummyEmulatorDefinition                       , &
         &                                                           likelihoodRedshiftDefinition                   , likelihoodMassParticleDefinition                        , &
         &                                                           likelihoodBinCountMinimumDefinition            , likelihoodMassFunctionTypeDefinition                    , &
         &                                                           likelihoodMassFunctionErrorModelDefinition
    type            (nodeList      ), pointer                     :: covarianceRows
    double precision                , allocatable, dimension(:  ) :: likelihoodMean
    double precision                , allocatable, dimension(:,:) :: likelihoodCovariance
    integer                         , allocatable, dimension(:  ) :: likelihoodExclusions
    integer                                                       :: i                                              , dimensionCount                            , &
         &                                                           likelihoodRealizationCount                     , likelihoodRealizationCountMinimum         , &
         &                                                           likelihoodEmulatorRebuildCount                 , likelihoodPolynomialOrder                 , &
         &                                                           likelihoodReportCount                          , likelihoodNeighborCount                   , &
         &                                                           likelihoodBinCountMinimum
    double precision                                              :: likelihoodSigmaBuffer                          , likelihoodLogLikelihoodBuffer             , &
         &                                                           likelihoodHaloMassMinimum                      , likelihoodHaloMassMaximum                 , &
         &                                                           likelihoodRedshiftMinimum                      , likelihoodRedshiftMaximum                 , &
         &                                                           likelihoodLogLikelihoodErrorTolerance          , likelihoodSurfaceBrightnessLimit          , &
         &                                                           likelihoodTolerance                            , likelihoodLineOfSightDepth                , &
         &                                                           likelihoodDelayInterval                        , likelihoodRedshift                        , &
         &                                                           likelihoodMassRangeMinimum                     , likelihoodMassParticle
    logical                                                       :: likelihoodEmulateOutliers                      , likelihoodUseSurveyLimits                 , &
         &                                                           likelihoodModelSurfaceBrightness               , likelihoodHalfIntegral                    , &
         &                                                           likelihoodDummyEmulator

    select case (char(XML_Extract_Text(XML_Get_First_Element_By_Tag_Name(definition,"type",directChildrenOnly=.true.))))
    case ("multivariateNormal")
       allocate(likelihoodMultivariateNormal :: newLikelihood)
       select type (newLikelihood)
       type is (likelihoodMultivariateNormal)
          likelihoodMeanDefinition       => XML_Get_First_Element_By_Tag_Name(definition,"mean"      )
          likelihoodCovarianceDefinition => XML_Get_First_Element_By_Tag_Name(definition,"covariance")
          covarianceRows                 => getElementsByTagName(likelihoodCovarianceDefinition,"row")
          dimensionCount=String_Count_Words(getTextContent(likelihoodMeanDefinition))
          allocate(likelihoodMean      (dimensionCount               ))
          allocate(likelihoodCovariance(dimensionCount,dimensionCount))
          call extractDataContent(likelihoodMeanDefinition,likelihoodMean) 
          do i=1,dimensionCount
             covarianceRow => item(covarianceRows,i-1)
             call extractDataContent(covarianceRow,likelihoodCovariance(i,:))
          end do
          newLikelihood=likelihoodMultivariateNormal(likelihoodMean,likelihoodCovariance)
          deallocate(likelihoodMean      )
          deallocate(likelihoodCovariance)
       end select
    case ("multivariateNormalStochastic")
       allocate(likelihoodMultivariateNormalStochastic :: newLikelihood)
       select type (newLikelihood)
       type is (likelihoodMultivariateNormalStochastic)
          likelihoodMeanDefinition                    => XML_Get_First_Element_By_Tag_Name(definition,"mean"                   )
          likelihoodCovarianceDefinition              => XML_Get_First_Element_By_Tag_Name(definition,"covariance"             )
          likelihoodRealizationCountDefinition        => XML_Get_First_Element_By_Tag_Name(definition,"realizationCount"       )
          likelihoodRealizationCountMinimumDefinition => XML_Get_First_Element_By_Tag_Name(definition,"realizationCountMinimum")
          covarianceRows                              => getElementsByTagName(likelihoodCovarianceDefinition,"row")
          call extractDataContent(likelihoodRealizationCountDefinition,likelihoodRealizationCount)
          dimensionCount=String_Count_Words(getTextContent(likelihoodMeanDefinition))
          allocate(likelihoodMean      (dimensionCount               ))
          allocate(likelihoodCovariance(dimensionCount,dimensionCount))
          call extractDataContent(likelihoodMeanDefinition,likelihoodMean) 
          do i=1,dimensionCount
             covarianceRow => item(covarianceRows,i-1)
             call extractDataContent(covarianceRow,likelihoodCovariance(i,:))
          end do
          newLikelihood=likelihoodMultivariateNormalStochastic(likelihoodMean,likelihoodCovariance,likelihoodRealizationCount,likelihoodRealizationCountMinimum)
          deallocate(likelihoodMean      )
          deallocate(likelihoodCovariance)
       end select
    case ("Galacticus")
       allocate(likelihoodGalacticus :: newLikelihood)
       select type (newLikelihood)
       type is (likelihoodGalacticus)
          likelihoodDelayIntervalDefinition => XML_Get_First_Element_By_Tag_Name(definition,"delayInterval")
          call extractDataContent(likelihoodDelayIntervalDefinition,likelihoodDelayInterval)
          newLikelihood=likelihoodGalacticus(configFileName,likelihoodDelayInterval)
       end select
    case ("gaussianRegression")
       allocate(likelihoodGaussianRegression :: newLikelihood)
       select type (newLikelihood)
       type is (likelihoodGaussianRegression)
          likelihoodDefinition                            => XML_Get_First_Element_By_Tag_Name(definition,"simulatorLikelihood"        )
          likelihoodEmulatorRebuildCountDefinition        => XML_Get_First_Element_By_Tag_Name(definition,"emulatorRebuildCount"       )
          likelihoodPolynomialOrderDefinition             => XML_Get_First_Element_By_Tag_Name(definition,"polynomialOrder"            )
          likelihoodSigmaBufferDefinition                 => XML_Get_First_Element_By_Tag_Name(definition,"sigmaBuffer"                )
          likelihoodLogLikelihoodBufferDefinition         => XML_Get_First_Element_By_Tag_Name(definition,"logLikelihoodBuffer"        )
          likelihoodLogLikelihoodErrorToleranceDefinition => XML_Get_First_Element_By_Tag_Name(definition,"logLikelihoodErrorTolerance")
          likelihoodReportCountDefinition                 => XML_Get_First_Element_By_Tag_Name(definition,"reportCount"                )
          likelihoodEmulateOutliersDefinition             => XML_Get_First_Element_By_Tag_Name(definition,"emulateOutliers"            )
          likelihoodDumpEmulatorDefinition                => XML_Get_First_Element_By_Tag_Name(definition,"dumpEmulator"               )
          likelihoodDummyEmulatorDefinition               => XML_Get_First_Element_By_Tag_Name(definition,"dummyEmulator"              )
          call extractDataContent(likelihoodEmulatorRebuildCountDefinition       ,likelihoodEmulatorRebuildCount       )
          call extractDataContent(likelihoodPolynomialOrderDefinition            ,likelihoodPolynomialOrder            )
          call extractDataContent(likelihoodSigmaBufferDefinition                ,likelihoodSigmaBuffer                )
          call extractDataContent(likelihoodLogLikelihoodBufferDefinition        ,likelihoodLogLikelihoodBuffer        )
          call extractDataContent(likelihoodLogLikelihoodErrorToleranceDefinition,likelihoodLogLikelihoodErrorTolerance)
          call extractDataContent(likelihoodReportCountDefinition                ,likelihoodReportCount                )
          call extractDataContent(likelihoodEmulateOutliersDefinition            ,likelihoodEmulateOutliers            )
          call extractDataContent(likelihoodDummyEmulatorDefinition              ,likelihoodDummyEmulator              )
          newLikelihood=likelihoodGaussianRegression(                                                  &
               &                                     likelihoodDefinition                            , &
               &                                     likelihoodEmulatorRebuildCount                  , &
               &                                     likelihoodPolynomialOrder                       , &
               &                                     likelihoodSigmaBuffer                           , &
               &                                     likelihoodLogLikelihoodBuffer                   , &
               &                                     likelihoodLogLikelihoodErrorTolerance           , &
               &                                     likelihoodReportCount                           , &
               &                                     likelihoodEmulateOutliers                       , &
               &                                     getTextContent(likelihoodDumpEmulatorDefinition), &
               &                                     likelihoodDummyEmulator                         , &
               &                                     configFileName                                    &
               &                                    )
       end select
    case ("posteriorPrior")
       allocate(likelihoodPosteriorPrior :: newLikelihood)
       select type (newLikelihood)
       type is (likelihoodPosteriorPrior)
          likelihoodDefinition              => XML_Get_First_Element_By_Tag_Name(definition,"wrappedLikelihood")
          likelihoodChainBaseNameDefinition => XML_Get_First_Element_By_Tag_Name(definition,"chainBaseName"    )
          likelihoodNeighborCountDefinition => XML_Get_First_Element_By_Tag_Name(definition,"neighborCount"    )
          likelihoodToleranceDefinition     => XML_Get_First_Element_By_Tag_Name(definition,"tolerance"        )
          likelihoodExclusionsDefinition    => XML_Get_First_Element_By_Tag_Name(definition,"exclusions"       )
          dimensionCount=String_Count_Words(getTextContent(likelihoodExclusionsDefinition))
          allocate(likelihoodExclusions(dimensionCount))
          call extractDataContent(likelihoodExclusionsDefinition   ,likelihoodExclusions   )
          call extractDataContent(likelihoodNeighborCountDefinition,likelihoodNeighborCount)
          call extractDataContent(likelihoodToleranceDefinition    ,likelihoodTolerance    )
          newLikelihood=likelihoodPosteriorPrior(                                                   &
               &                                 likelihoodDefinition                             , &
               &                                 getTextContent(likelihoodChainBaseNameDefinition), &
               &                                 likelihoodNeighborCount                          , &
               &                                 likelihoodTolerance                              , &
               &                                 configFileName                                   , &
               &                                 likelihoodExclusions                               &
               &                                )
       end select
    case ("massFunction")
       allocate(likelihoodMassFunction :: newLikelihood)
       select type (newLikelihood)
       type is (likelihoodMassFunction)
          likelihoodHaloMassMinimumDefinition        => XML_Get_First_Element_By_Tag_Name(definition,"haloMassMinimum"       )
          likelihoodHaloMassMaximumDefinition        => XML_Get_First_Element_By_Tag_Name(definition,"haloMassMaximum"       )
          likelihoodRedshiftMinimumDefinition        => XML_Get_First_Element_By_Tag_Name(definition,"redshiftMinimum"       )
          likelihoodRedshiftMaximumDefinition        => XML_Get_First_Element_By_Tag_Name(definition,"redshiftMaximum"       )
          likelihoodUseSurveyLimitsDefinition        => XML_Get_First_Element_By_Tag_Name(definition,"useSurveyLimits"       )
          likelihoodMassFunctionFileNameDefinition   => XML_Get_First_Element_By_Tag_Name(definition,"massFunctionFileName"  )
          call extractDataContent(likelihoodHaloMassMinimumDefinition       ,likelihoodHaloMassMinimum       )
          call extractDataContent(likelihoodHaloMassMaximumDefinition       ,likelihoodHaloMassMaximum       )
          call extractDataContent(likelihoodRedshiftMinimumDefinition       ,likelihoodRedshiftMinimum       )
          call extractDataContent(likelihoodRedshiftMaximumDefinition       ,likelihoodRedshiftMaximum       )
          call extractDataContent(likelihoodUseSurveyLimitsDefinition       ,likelihoodUseSurveyLimits       )
          if (XML_Path_Exists(definition,"modelSurfaceBrightness")) then
             likelihoodModelSurfaceBrightnessDefinition => XML_Get_First_Element_By_Tag_Name(definition,"modelSurfaceBrightness")
             likelihoodSurfaceBrightnessLimitDefinition => XML_Get_First_Element_By_Tag_Name(definition,"surfaceBrightnessLimit")
             call extractDataContent(likelihoodModelSurfaceBrightnessDefinition,likelihoodModelSurfaceBrightness)
             call extractDataContent(likelihoodSurfaceBrightnessLimitDefinition,likelihoodSurfaceBrightnessLimit)
          else
             likelihoodModelSurfaceBrightness=.false.
             likelihoodSurfaceBrightnessLimit=0.0d0
          end if
          newLikelihood=likelihoodMassFunction(                                                          &
               &                               likelihoodHaloMassMinimum                               , &
               &                               likelihoodHaloMassMaximum                               , &
               &                               likelihoodRedshiftMinimum                               , &
               &                               likelihoodRedshiftMaximum                               , &
               &                               likelihoodUseSurveyLimits                               , &
               &                               getTextContent(likelihoodMassFunctionFileNameDefinition), &
               &                               likelihoodModelSurfaceBrightness                        , &
               &                               likelihoodSurfaceBrightnessLimit                          &
               &                              )
       end select
    case ("haloMassFunction")
       allocate(likelihoodHaloMassFunction :: newLikelihood)
       select type (newLikelihood)
       type is (likelihoodHaloMassFunction)
          likelihoodRedshiftDefinition               => XML_Get_First_Element_By_Tag_Name(definition,"redshift"              )
          likelihoodMassRangeMinimumDefinition       => XML_Get_First_Element_By_Tag_Name(definition,"massRangeMinimum"      )
          likelihoodBinCountMinimumDefinition        => XML_Get_First_Element_By_Tag_Name(definition,"binCountMinimum"       )
          likelihoodMassFunctionTypeDefinition       => XML_Get_First_Element_By_Tag_Name(definition,"massFunctionType"      )
          likelihoodMassFunctionErrorModelDefinition => XML_Get_First_Element_By_Tag_Name(definition,"massFunctionErrorModel")
          likelihoodMassFunctionFileNameDefinition   => XML_Get_First_Element_By_Tag_Name(definition,"fileName"              )
          if (getTextContent(likelihoodMassFunctionErrorModelDefinition) == "sphericalOverdensity") then
             likelihoodMassParticleDefinition => XML_Get_First_Element_By_Tag_Name(definition,"massParticle")
             call extractDataContent(likelihoodMassParticleDefinition,likelihoodMassParticle)
          else
             likelihoodMassParticle=0.0d0
          end if
          call extractDataContent(likelihoodRedshiftDefinition        ,likelihoodRedshift        )
          call extractDataContent(likelihoodMassRangeMinimumDefinition,likelihoodMassRangeMinimum)
          call extractDataContent(likelihoodBinCountMinimumDefinition ,likelihoodBinCountMinimum )
          newLikelihood=likelihoodHaloMassFunction(                                                            &
               &                                   getTextContent(likelihoodMassFunctionFileNameDefinition  ), &
               &                                   likelihoodRedshift                                        , &
               &                                   likelihoodMassRangeMinimum                                , &
               &                                   likelihoodBinCountMinimum                                 , &
               &                                   getTextContent(likelihoodMassFunctionTypeDefinition      ), &
               &                                   getTextContent(likelihoodMassFunctionErrorModelDefinition), &
               &                                   likelihoodMassParticle                                      &
               &                                  )
       end select
    case ("projectedCorrelationFunction")
       allocate(likelihoodProjectedCorrelationFunction :: newLikelihood)
       select type (newLikelihood)
       type is (likelihoodProjectedCorrelationFunction)
          likelihoodHaloMassMinimumDefinition                      => XML_Get_First_Element_By_Tag_Name(definition,"haloMassMinimum"                     )
          likelihoodHaloMassMaximumDefinition                      => XML_Get_First_Element_By_Tag_Name(definition,"haloMassMaximum"                     )
          likelihoodLineOfSightDepthDefinition                     => XML_Get_First_Element_By_Tag_Name(definition,"lineOfSightDepth"                    )
          likelihoodHalfIntegralDefinition                         => XML_Get_First_Element_By_Tag_Name(definition,"halfIntegral"                        )
          likelihoodProjectedCorrelationFunctionFileNameDefinition => XML_Get_First_Element_By_Tag_Name(definition,"projectedCorrelationFunctionFileName")
          call extractDataContent(likelihoodHaloMassMinimumDefinition ,likelihoodHaloMassMinimum )
          call extractDataContent(likelihoodHaloMassMaximumDefinition ,likelihoodHaloMassMaximum )
          call extractDataContent(likelihoodLineOfSightDepthDefinition,likelihoodLineOfSightDepth)
          call extractDataContent(likelihoodHalfIntegralDefinition    ,likelihoodHalfIntegral    )
          newLikelihood=likelihoodProjectedCorrelationFunction(                                                                          &
               &                                               likelihoodHaloMassMinimum                                               , &
               &                                               likelihoodHaloMassMaximum                                               , &
               &                                               likelihoodLineOfSightDepth                                              , &
               &                                               likelihoodHalfIntegral                                                  , &
               &                                               getTextContent(likelihoodProjectedCorrelationFunctionFileNameDefinition)  &
               &                                              )
       end select
    case ("sedFit")
       allocate(likelihoodSEDFit :: newLikelihood)
       select type (newLikelihood)
       type is (likelihoodSEDFit)
          newLikelihood=likelihoodSEDFit(definition)
       end select
    case default
       call Galacticus_Error_Report('likelihoodNew','likelihood type is unrecognized')
    end select
    return
  end function likelihoodNew
  
  logical function likelihoodWillEvaluate(self,simulationState,parameterMappings,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
    !% Return true if the likelihood will be evaluated
    class           (likelihood ), intent(inout)               :: self
    class           (state      ), intent(in   )               :: simulationState
    type            (mappingList), intent(in   ), dimension(:) :: parameterMappings
    class           (convergence), intent(inout)               :: simulationConvergence
    double precision             , intent(in   )               :: temperature          , logLikelihoodCurrent, &
         &                                                        logPriorCurrent      , logPriorProposed
    !GCC$ attributes unused :: self, simulationState, parameterMappings, simulationConvergence, temperature, logLikelihoodCurrent, logPriorCurrent, logPriorProposed
    
    likelihoodWillEvaluate=.true.
    return
  end function likelihoodWillEvaluate
  
  ! Include all likelihood methods.
  include 'constraints.likelihoods.multivariate_normal.methods.inc'
  include 'constraints.likelihoods.multivariate_normal.stochastic.methods.inc'
  include 'constraints.likelihoods.Galacticus.methods.inc'
  include 'constraints.likelihoods.gaussian_regression.methods.inc'
  include 'constraints.likelihoods.mass_function.methods.inc'
  include 'constraints.likelihoods.halo_mass_function.methods.inc'
  include 'constraints.likelihoods.projected_correlation_function.methods.inc'
  include 'constraints.likelihoods.SED_fit.methods.inc'
  include 'constraints.likelihoods.posterior_as_prior.methods.inc'

end module Constraints_Likelihoods
