!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  use ISO_Varying_String
  use Pseudo_Random
  use FGSL
  private
  public :: likelihoodNew
  
  ! A very small log likelihood which is used as an approximation to zero likelihood.
  double precision, parameter, public :: logImpossible=-1.0d30

  ! Define the basic likelihood class.
  type, abstract, public :: likelihood
   contains
     !@ <objectMethods>
     !@   <object>likelihood</object>
     !@   <objectMethod>
     !@     <method>evaluate</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless class(state)\textgreater} simulationState\argin, \doublezero\ temperature\argin, \doublezero\ logLikelihoodCurrent\argin</arguments>
     !@     <description>Evaluate the model likelihood at the given {\tt simulationState} and return the log-likelihood.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure(likelihoodEvaluate), deferred :: evaluate
  end type likelihood

  ! Interface for deferred functions.
  abstract interface
     double precision function likelihoodEvaluate(self,simulationState,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
       import :: likelihood, state, convergence
       class           (likelihood ), intent(inout) :: self
       class           (state      ), intent(in   ) :: simulationState
       class           (convergence), intent(inout) :: simulationConvergence
       double precision             , intent(in   ) :: temperature          , logLikelihoodCurrent, &
            &                                          logPriorCurrent      , logPriorProposed
     end function likelihoodEvaluate
  end interface

  ! Include all likelihood types.
  include 'constraints.likelihoods.multivariate_normal.type.inc'
  include 'constraints.likelihoods.multivariate_normal.stochastic.type.inc'
  include 'constraints.likelihoods.Galacticus.type.inc'
  include 'constraints.likelihoods.gaussian_regression.type.inc'

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
    type            (node          ), pointer                     :: likelihoodMeanDefinition                       , likelihoodCovarianceDefinition         , &
         &                                                           covarianceRow                                  , likelihoodRealizationCountDefinition   , &
         &                                                           likelihoodRealizationCountMinimumDefinition    , likelihoodDefinition                   , &
         &                                                           likelihoodEmulatorRebuildCountDefinition       , likelihoodPolynomialOrderDefinition    , &
         &                                                           likelihoodSigmaBufferDefinition                , likelihoodLogLikelihoodBufferDefinition, &
         &                                                           likelihoodLogLikelihoodErrorToleranceDefinition, likelihoodReportCountDefinition        , &
         &                                                           likelihoodEmulateOutliersDefinition
    type            (nodeList      ), pointer                     :: covarianceRows
    double precision                , allocatable, dimension(:  ) :: likelihoodMean
    double precision                , allocatable, dimension(:,:) :: likelihoodCovariance
    integer                                                       :: i                                              , dimensionCount                         , &
         &                                                           likelihoodRealizationCount                     , likelihoodRealizationCountMinimum      , &
         &                                                           likelihoodEmulatorRebuildCount                 , likelihoodPolynomialOrder              , &
         &                                                           likelihoodReportCount
    double precision                                              :: likelihoodSigmaBuffer                          , likelihoodLogLikelihoodBuffer          , &
         &                                                           likelihoodLogLikelihoodErrorTolerance
    logical                                                       :: likelihoodEmulateOutliers

    select case (char(XML_Extract_Text(XML_Get_First_Element_By_Tag_Name(definition,"type"))))
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
          newLikelihood=likelihoodGalacticus(configFileName)
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
          call extractDataContent(likelihoodEmulatorRebuildCountDefinition       ,likelihoodEmulatorRebuildCount       )
          call extractDataContent(likelihoodPolynomialOrderDefinition            ,likelihoodPolynomialOrder            )
          call extractDataContent(likelihoodSigmaBufferDefinition                ,likelihoodSigmaBuffer                )
          call extractDataContent(likelihoodLogLikelihoodBufferDefinition        ,likelihoodLogLikelihoodBuffer        )
          call extractDataContent(likelihoodLogLikelihoodErrorToleranceDefinition,likelihoodLogLikelihoodErrorTolerance)
          call extractDataContent(likelihoodReportCountDefinition                ,likelihoodReportCount                )
          call extractDataContent(likelihoodEmulateOutliersDefinition            ,likelihoodEmulateOutliers            )
          newLikelihood=likelihoodGaussianRegression(                                       &
               &                                     likelihoodDefinition                 , &
               &                                     likelihoodEmulatorRebuildCount       , &
               &                                     likelihoodPolynomialOrder            , &
               &                                     likelihoodSigmaBuffer                , &
               &                                     likelihoodLogLikelihoodBuffer        , &
               &                                     likelihoodLogLikelihoodErrorTolerance, &
               &                                     likelihoodReportCount                , &
               &                                     likelihoodEmulateOutliers            , &
               &                                     configFileName                         &
               &                                    )
       end select
    case default
       call Galacticus_Error_Report('likelihoodNew','likelihood type is unrecognized')
    end select
    return
  end function likelihoodNew

  ! Include all likelihood methods.
  include 'constraints.likelihoods.multivariate_normal.methods.inc'
  include 'constraints.likelihoods.multivariate_normal.stochastic.methods.inc'
  include 'constraints.likelihoods.Galacticus.methods.inc'
  include 'constraints.likelihoods.gaussian_regression.methods.inc'

end module Constraints_Likelihoods
