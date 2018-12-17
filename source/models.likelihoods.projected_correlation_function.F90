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
  
  !% Implementation of a posterior sampling likelihood class which implements a likelihood for projected correlation functions.

  use Linear_Algebra
  
  !# <posteriorSampleLikelihood name="posteriorSampleLikelihoodPrjctdCorrelationFunction">
  !#  <description>A posterior sampling likelihood class which implements a likelihood for projected correlation functions.</description>
  !# </posteriorSampleLikelihood>
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodPrjctdCorrelationFunction
     !% Implementation of a posterior sampling likelihood class which implements a likelihood for projected correlation functions.
     private
     double precision                                              :: haloMassMinimum                     , haloMassMaximum             , &
          &                                                           lineOfSightDepth
     logical                                                       :: halfIntegral
     double precision                , dimension(:  ), allocatable :: separation                          , massMaximum                 , &
          &                                                           massMinimum
     double precision                , dimension(:,:), allocatable :: covarianceMatrix                    , projectedCorrelationFunction, &
          &                                                           projectedCorrelationFunctionObserved, integralConstraint
     type            (vector        )                              :: means
     type            (matrix        )                              :: covariance                          , inverseCovariance
     type            (varying_string)                              :: fileName
   contains
     procedure :: evaluate        => projectedCorrelationFunctionEvaluate
     procedure :: functionChanged => projectedCorrelationFunctionFunctionChanged
  end type posteriorSampleLikelihoodPrjctdCorrelationFunction

  interface posteriorSampleLikelihoodPrjctdCorrelationFunction
     !% Constructors for the {\normalfont \ttfamily projectedCorrelationFunction} posterior sampling convergence class.
     module procedure projectedCorrelationFunctionConstructorParameters
     module procedure projectedCorrelationFunctionConstructorInternal
  end interface posteriorSampleLikelihoodPrjctdCorrelationFunction

contains

  function projectedCorrelationFunctionConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily projectedCorrelationFunction} posterior sampling convergence class which builds the object from a
    !% parameter set.
    use Input_Parameters
    implicit none
    type            (posteriorSampleLikelihoodPrjctdCorrelationFunction)                :: self
    type            (inputParameters                                   ), intent(inout) :: parameters
    double precision                                                                    :: haloMassMinimum    , haloMassMaximum, &
         &                                                                                 lineOfSightDepth
    logical                                                                             :: halfIntegral
    type            (varying_string                                    )                :: fileName

    !# <inputParameter>
    !#   <name>haloMassMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The minimum halo mass over which to integrate.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>haloMassMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The maximum halo mass over which to integrate.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>lineOfSightDepth</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The line of sight depth over which to integrate.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>halfIntegral</name>
    !#   <cardinality>1</cardinality>
    !#   <description>If true, integrate only over positive line of sight depths.</description>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>fileName</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The name of the file containing the target projected correlation function.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    self=posteriorSampleLikelihoodPrjctdCorrelationFunction(haloMassMinimum,haloMassMaximum,lineOfSightDepth,halfIntegral,char(fileName))
    !# <inputParametersValidate source="parameters"/>
    return
  end function projectedCorrelationFunctionConstructorParameters

  function projectedCorrelationFunctionConstructorInternal(haloMassMinimum,haloMassMaximum,lineOfSightDepth,halfIntegral,fileName) result(self)
    !% Constructor for ``projectedCorrelationFunction'' posterior sampling likelihood class.
    use Galacticus_Paths
    use IO_HDF5
    use Memory_Management
    use Node_Component_Dark_Matter_Profile_Scale
    implicit none
    type            (posteriorSampleLikelihoodPrjctdCorrelationFunction)                :: self
    double precision                                                    , intent(in   ) :: haloMassMinimum    , haloMassMaximum, &
         &                                                                                 lineOfSightDepth
    logical                                                             , intent(in   ) :: halfIntegral
    character       (len=*                                             ), intent(in   ) :: fileName
    type            (hdf5Object                                        )                :: file
    !# <constructorAssign variables="haloMassMinimum, haloMassMaximum, lineOfSightDepth, halfIntegral, fileName"/>

    ! Read the projected correlation function file.
    !$ call hdf5Access%set()
    call file%openFile(char(galacticusPath(pathTypeDataStatic))//fileName,readOnly=.true.)
    call file%readDataset("separation"                          ,self%separation                          )
    call file%readDataset("projectedCorrelationFunctionObserved",self%projectedCorrelationFunctionObserved)
    call file%readDataset("covariance"                          ,self%covarianceMatrix                    )
    call file%readDataset("massMinimum"                         ,self%massMinimum                         )
    call file%readDataset("massMaximum"                         ,self%massMaximum                         )
    if (file%hasDataset("integralConstraint")) then
       call file%readDataset("integralConstraint"               ,self%integralConstraint                  )
    else
       call allocateArray(self%integralConstraint,shape(self%projectedCorrelationFunctionObserved))
       self%integralConstraint=1.0d0
    end if
    call file%close()
    !$ call hdf5Access%unset()
    ! Allocate storage for the model projected correlation function.
    call allocateArray(self%projectedCorrelationFunction,[size(self%separation),size(self%massMinimum)])
    ! Find the inverse covariance matrix.
    self%covariance       =self%covarianceMatrix
    self%inverseCovariance=self%covariance      %invert()
    call self%inverseCovariance%makeSemiPositiveDefinite()
    return
  end function projectedCorrelationFunctionConstructorInternal

  double precision function projectedCorrelationFunctionEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !% Return the log-likelihood for the projected correlation function likelihood function.
    use Posterior_Sampling_State
    use Posterior_Sampling_Convergence
    use Conditional_Mass_Functions
    use Cosmology_Functions
    use Galacticus_Error
    use Halo_Model_Projected_Correlations
    use Models_Likelihoods_Constants
    implicit none
    class           (posteriorSampleLikelihoodPrjctdCorrelationFunction), intent(inout)               :: self
    class           (posteriorSampleStateClass                         ), intent(inout)               :: simulationState
    type            (modelParameterList                                ), intent(in   ), dimension(:) :: modelParametersActive_  , modelParametersInactive_
    class           (posteriorSampleConvergenceClass                   ), intent(inout)               :: simulationConvergence
    double precision                                                    , intent(in   )               :: temperature             , logLikelihoodCurrent    , &
         &                                                                                               logPriorCurrent         , logPriorProposed
    real                                                                , intent(inout)               :: timeEvaluate
    double precision                                                    , intent(  out), optional     :: logLikelihoodVariance
    logical                                                             , intent(inout), optional     :: forceAcceptance
    double precision                                                    , allocatable  , dimension(:) :: stateVector
    type            (conditionalMassFunctionBehroozi2010               )                              :: conditionalMassFunction_
    type            (vector                                            )                              :: difference
    integer                                                                                           :: i
    !GCC$ attributes unused :: logLikelihoodCurrent, logPriorCurrent, simulationConvergence, temperature, timeEvaluate, modelParametersInactive_, forceAcceptance

    ! There is no variance in our likelihood estimate.
    if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
    ! Do not evaluate if the proposed prior is impossible.
    if (logPriorProposed <= logImpossible) then
       projectedCorrelationFunctionEvaluate=0.0d0
       return
    end if
    ! Construct the conditional mass function object.
    stateVector=simulationState%get()
    if (size(stateVector) /= 11) call Galacticus_Error_Report('11 parameters are required for this likelihood function'//{introspection:location})
    do i=1,size(stateVector)
       stateVector(i)=modelParametersActive_(i)%modelParameter_%unmap(stateVector(i))
    end do
    conditionalMassFunction_                                     &
         & =conditionalMassFunctionBehroozi2010(                 &
         &                                      stateVector( 1), &
         &                                      stateVector( 2), &
         &                                      stateVector( 3), &
         &                                      stateVector( 4), &
         &                                      stateVector( 5), &
         &                                      stateVector( 6), &
         &                                      stateVector( 7), &
         &                                      stateVector( 8), &
         &                                      stateVector( 9), &
         &                                      stateVector(10), &
         &                                      stateVector(11)  &
         &                                     )
    deallocate(stateVector)
    ! Compute the projected correlation function.
    do i=1,size(self%massMinimum)
       call Halo_Model_Projected_Correlation(                                        &
            &                                conditionalMassFunction_              , &
            &                                self%separation                       , &
            &                                self%massMinimum                 (  i), &
            &                                self%massMaximum                 (  i), &
            &                                self%haloMassMinimum                  , &
            &                                self%haloMassMaximum                  , &
            &                                self%lineOfSightDepth                 , &
            &                                self%halfIntegral                     , &
            &                                self%projectedCorrelationFunction(:,i)  &
            &                               )
       ! Apply the integral constraint.
       self%projectedCorrelationFunction(:,i)=self%projectedCorrelationFunction(:,i)/self%integralConstraint(:,i)
    end do
    ! Evaluate the log-likelihood.
    difference                          =reshape(                                                         &
         &                                        +     self%projectedCorrelationFunction                 &
         &                                        -     self%projectedCorrelationFunctionObserved       , &
         &                                       [                                                        &
         &                                         size(self%projectedCorrelationFunction        ,dim=1)  &
         &                                        *size(self%projectedCorrelationFunction        ,dim=2)  &
         &                                       ]                                                        &
         &                                      )
    projectedCorrelationFunctionEvaluate=-0.5d0*(difference*(self%inverseCovariance*difference))
    return
  end function projectedCorrelationFunctionEvaluate

  subroutine projectedCorrelationFunctionFunctionChanged(self)
    !% Respond to possible changes in the likelihood function.
    implicit none
    class(posteriorSampleLikelihoodPrjctdCorrelationFunction), intent(inout) :: self
    !GCC$ attributes unused :: self

    return
  end subroutine projectedCorrelationFunctionFunctionChanged
