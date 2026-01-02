!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  Implementation of a posterior sampling likelihood class which implements a likelihood for projected correlation functions.
  !!}

  use :: Cosmology_Functions       , only : cosmologyFunctions          , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Biases   , only : darkMatterHaloBias          , darkMatterHaloBiasClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScale         , darkMatterHaloScaleClass
  use :: Dark_Matter_Profile_Scales, only : darkMatterProfileScaleRadius, darkMatterProfileScaleRadiusClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMO        , darkMatterProfileDMOClass
  use :: Geometry_Surveys          , only : surveyGeometry              , surveyGeometryClass
  use :: Halo_Mass_Functions       , only : haloMassFunction            , haloMassFunctionClass
  use :: Linear_Algebra            , only : matrix                      , vector
  use :: Linear_Growth             , only : linearGrowth                , linearGrowthClass
  use :: Power_Spectra             , only : powerSpectrum               , powerSpectrumClass

  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodPrjctdCorrelationFunction">
   <description>
    The likelihood is computed as
    \begin{equation}
    \log \mathcal{L} = -{1\over2} \Delta \cdot \mathcal{C}^{-1} \cdot \Delta^\mathrm{T},
    \end{equation}
    where $\mathcal{C}$ is the covariance matrix, and $\Delta_i = w_i^\mathrm{model} - w_i^\mathrm{obs}$, $w_i^\mathrm{model}$
    is the computed projected correlation function at the $i^\mathrm{th}$ separation, and $w_i^\mathrm{obs}$ is the observed
    projected correlation function at the $i^\mathrm{th}$ separation. The projected correlation function is computed using the
    halo model and the parameterized conditional galaxy mass function of \cite[][see also \protect\cite{leauthaud_new_2011};
    \protect\refPhysics{conditionalMassFunctionBehroozi2010}]{behroozi_comprehensive_2010}.  The details of the projected
    correlation function calculation are specified by the following subparameters:
    \begin{description}
    \item[{\normalfont \ttfamily haloMass(Min|Max)imum}] The minimum/maximum halo mass over which to integrate in the halo model;
    \item[{\normalfont \ttfamily redshift(Min|Max)imum}] The minimum/maximum redshift over which to integrate in the halo model;
    \item[{\normalfont \ttfamily projectedCorrelationFunctionFileName}] The name of an HDF5 file containing the observed projected
      correlation function and its covariance matrix.
    \end{description}
    
    The HDF5 file specified by the {\normalfont \ttfamily projectedCorrelationFunctionFileName} element should contain a {\normalfont
    \ttfamily separation} dataset, giving the separations at which the projected correlation function is measured (in units of Mpc),
    a {\normalfont \ttfamily projectedCorrelationFunctionObserved} dataset giving the observed values of the projected correlation
    function at those separations (in units of Mpc), and a {\normalfont \ttfamily covariance} dataset, giving the covariance of the
    projected correlation function (in units of Mpc$^2$).
   </description>
   <runTimeFileDependencies paths="fileName"/>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodPrjctdCorrelationFunction
     !!{
     Implementation of a posterior sampling likelihood class which implements a likelihood for projected correlation functions.
     !!}
     private
     class           (powerSpectrumClass               ), pointer                     :: powerSpectrum_                       => null()
     class           (cosmologyFunctionsClass          ), pointer                     :: cosmologyFunctions_                  => null()
     class           (surveyGeometryClass              ), pointer                     :: surveyGeometry_                      => null()
     class           (darkMatterHaloScaleClass         ), pointer                     :: darkMatterHaloScale_                 => null()
     class           (haloMassFunctionClass            ), pointer                     :: haloMassFunction_                    => null()
     class           (darkMatterProfileDMOClass        ), pointer                     :: darkMatterProfileDMO_                => null()
     class           (darkMatterHaloBiasClass          ), pointer                     :: darkMatterHaloBias_                  => null()
     class           (darkMatterProfileScaleRadiusClass), pointer                     :: darkMatterProfileScaleRadius_        => null()
     double precision                                                                 :: haloMassMinimum                               , haloMassMaximum             , &
          &                                                                              lineOfSightDepth
     logical                                                                          :: halfIntegral
     double precision                                   , dimension(:  ), allocatable :: separation                                    , massMaximum                 , &
          &                                                                              massMinimum
     double precision                                   , dimension(:,:), allocatable :: covarianceMatrix                              , projectedCorrelationFunction, &
          &                                                                              projectedCorrelationFunctionObserved          , integralConstraint
     type            (vector                           )                              :: means
     type            (matrix                           )                              :: covariance
     type            (varying_string                   )                              :: fileName
   contains
     final     ::                    projectedCorrelationFunctionDestructor
     procedure :: evaluate        => projectedCorrelationFunctionEvaluate
     procedure :: functionChanged => projectedCorrelationFunctionFunctionChanged
  end type posteriorSampleLikelihoodPrjctdCorrelationFunction

  interface posteriorSampleLikelihoodPrjctdCorrelationFunction
     !!{
     Constructors for the \refClass{posteriorSampleLikelihoodPrjctdCorrelationFunction} posterior sampling convergence class.
     !!}
     module procedure projectedCorrelationFunctionConstructorParameters
     module procedure projectedCorrelationFunctionConstructorInternal
  end interface posteriorSampleLikelihoodPrjctdCorrelationFunction

contains

  function projectedCorrelationFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleLikelihoodPrjctdCorrelationFunction} posterior sampling convergence class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleLikelihoodPrjctdCorrelationFunction)                :: self
    type            (inputParameters                                   ), intent(inout) :: parameters
    class           (powerSpectrumClass                                ), pointer       :: powerSpectrum_
    class           (cosmologyFunctionsClass                           ), pointer       :: cosmologyFunctions_
    class           (surveyGeometryClass                               ), pointer       :: surveyGeometry_
    class           (darkMatterHaloScaleClass                          ), pointer       :: darkMatterHaloScale_
    class           (haloMassFunctionClass                             ), pointer       :: haloMassFunction_
    class           (darkMatterProfileDMOClass                         ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloBiasClass                           ), pointer       :: darkMatterHaloBias_
    class           (darkMatterProfileScaleRadiusClass                 ), pointer       :: darkMatterProfileScaleRadius_
    double precision                                                                    :: haloMassMinimum    , haloMassMaximum, &
         &                                                                                 lineOfSightDepth
    logical                                                                             :: halfIntegral
    type            (varying_string                                    )                :: fileName

    !![
    <inputParameter>
      <name>haloMassMinimum</name>
      <description>The minimum halo mass over which to integrate.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>haloMassMaximum</name>
      <description>The maximum halo mass over which to integrate.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>lineOfSightDepth</name>
      <description>The line of sight depth over which to integrate.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>halfIntegral</name>
      <description>If true, integrate only over positive line of sight depths.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>fileName</name>
      <description>The name of the file containing the target projected correlation function.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="powerSpectrum"                name="powerSpectrum_"                source="parameters"/>
    <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    <objectBuilder class="surveyGeometry"               name="surveyGeometry_"               source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    <objectBuilder class="haloMassFunction"             name="haloMassFunction_"             source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"         name="darkMatterProfileDMO_"         source="parameters"/>
    <objectBuilder class="darkMatterHaloBias"           name="darkMatterHaloBias_"           source="parameters"/>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    !!]
    self=posteriorSampleLikelihoodPrjctdCorrelationFunction(haloMassMinimum,haloMassMaximum,lineOfSightDepth,halfIntegral,char(fileName),powerSpectrum_,cosmologyFunctions_,surveyGeometry_,darkMatterHaloScale_,haloMassFunction_,darkMatterProfileDMO_,darkMatterHaloBias_,darkMatterProfileScaleRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="powerSpectrum_"               />
    <objectDestructor name="cosmologyFunctions_"          />
    <objectDestructor name="surveyGeometry_"              />
    <objectDestructor name="darkMatterHaloScale_"         />
    <objectDestructor name="haloMassFunction_"            />
    <objectDestructor name="darkMatterProfileDMO_"        />
    <objectDestructor name="darkMatterHaloBias_"          />
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    !!]
    return
  end function projectedCorrelationFunctionConstructorParameters

  function projectedCorrelationFunctionConstructorInternal(haloMassMinimum,haloMassMaximum,lineOfSightDepth,halfIntegral,fileName,powerSpectrum_,cosmologyFunctions_,surveyGeometry_,darkMatterHaloScale_,haloMassFunction_,darkMatterProfileDMO_,darkMatterHaloBias_,darkMatterProfileScaleRadius_) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleLikelihoodPrjctdCorrelationFunction} posterior sampling likelihood class.
    !!}
    use :: Input_Paths      , only : inputPath    , pathTypeDataStatic
    use :: HDF5_Access      , only : hdf5Access
    use :: IO_HDF5          , only : hdf5Object
    use :: Linear_Algebra   , only : assignment(=)
    implicit none
    type            (posteriorSampleLikelihoodPrjctdCorrelationFunction)                        :: self
    double precision                                                    , intent(in   )         :: haloMassMinimum    , haloMassMaximum, &
         &                                                                                         lineOfSightDepth
    logical                                                             , intent(in   )         :: halfIntegral
    character       (len=*                                             ), intent(in   )         :: fileName
    class           (powerSpectrumClass                                ), intent(in   ), target :: powerSpectrum_
    class           (cosmologyFunctionsClass                           ), intent(in   ), target :: cosmologyFunctions_
    class           (surveyGeometryClass                               ), intent(in   ), target :: surveyGeometry_
    class           (darkMatterHaloScaleClass                          ), intent(in   ), target :: darkMatterHaloScale_
    class           (haloMassFunctionClass                             ), intent(in   ), target :: haloMassFunction_
    class           (darkMatterProfileDMOClass                         ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloBiasClass                           ), intent(in   ), target :: darkMatterHaloBias_
    class           (darkMatterProfileScaleRadiusClass                 ), intent(in   ), target :: darkMatterProfileScaleRadius_
    type            (hdf5Object                                        )                        :: file
    !![
    <constructorAssign variables="haloMassMinimum, haloMassMaximum, lineOfSightDepth, halfIntegral, fileName, *powerSpectrum_, *cosmologyFunctions_, *surveyGeometry_, *darkMatterHaloScale_, *haloMassFunction_, *darkMatterProfileDMO_, *darkMatterHaloBias_, *darkMatterProfileScaleRadius_"/>
    !!]

    ! Read the projected correlation function file.
    !$ call hdf5Access%set()
    call file%openFile(char(inputPath(pathTypeDataStatic))//fileName,readOnly=.true.)
    call file%readDataset("separation"                          ,self%separation                          )
    call file%readDataset("projectedCorrelationFunctionObserved",self%projectedCorrelationFunctionObserved)
    call file%readDataset("covariance"                          ,self%covarianceMatrix                    )
    call file%readDataset("massMinimum"                         ,self%massMinimum                         )
    call file%readDataset("massMaximum"                         ,self%massMaximum                         )
    if (file%hasDataset("integralConstraint")) then
       call file%readDataset("integralConstraint"               ,self%integralConstraint                  )
    else
       allocate(self%integralConstraint,mold=self%projectedCorrelationFunctionObserved)
       self%integralConstraint=1.0d0
    end if
    call file%close()
    !$ call hdf5Access%unset()
    ! Allocate storage for the model projected correlation function.
    allocate(self%projectedCorrelationFunction(size(self%separation),size(self%massMinimum)))
    ! Build the covariance matrix.
    self%covariance=self%covarianceMatrix
    return
  end function projectedCorrelationFunctionConstructorInternal

  subroutine projectedCorrelationFunctionDestructor(self)
    !!{
    Destructor for the \refClass{posteriorSampleLikelihoodPrjctdCorrelationFunction} class.
    !!}
    implicit none
    type(posteriorSampleLikelihoodPrjctdCorrelationFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%powerSpectrum_"               />
    <objectDestructor name="self%cosmologyFunctions_"          />
    <objectDestructor name="self%surveyGeometry_"              />
    <objectDestructor name="self%darkMatterHaloScale_"         />
    <objectDestructor name="self%haloMassFunction_"            />
    <objectDestructor name="self%darkMatterProfileDMO_"        />
    <objectDestructor name="self%darkMatterHaloBias_"          />
    <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    !!]
    return
  end subroutine projectedCorrelationFunctionDestructor

  double precision function projectedCorrelationFunctionEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !!{
    Return the log-likelihood for the projected correlation function likelihood function.
    !!}
    use :: Conditional_Mass_Functions       , only : conditionalMassFunctionBehroozi2010
    use :: Error                            , only : Error_Report
    use :: Halo_Model_Projected_Correlations, only : Halo_Model_Projected_Correlation
    use :: Linear_Algebra                   , only : assignment(=)                      , operator(*)
    use :: Models_Likelihoods_Constants     , only : logImpossible
    use :: Posterior_Sampling_Convergence   , only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State         , only : posteriorSampleStateClass
    implicit none
    class           (posteriorSampleLikelihoodPrjctdCorrelationFunction), intent(inout), target       :: self
    class           (posteriorSampleStateClass                         ), intent(inout)               :: simulationState
    type            (modelParameterList                                ), intent(inout), dimension(:) :: modelParametersActive_  , modelParametersInactive_
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
    !$GLC attributes unused :: logLikelihoodCurrent, logPriorCurrent, simulationConvergence, temperature, timeEvaluate, modelParametersInactive_, forceAcceptance

    ! There is no variance in our likelihood estimate.
    if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
    ! Do not evaluate if the proposed prior is impossible.
    if (logPriorProposed <= logImpossible) then
       projectedCorrelationFunctionEvaluate=0.0d0
       return
    end if
    ! Construct the conditional mass function object.
    stateVector=simulationState%get()
    if (size(stateVector) /= 11) call Error_Report('11 parameters are required for this likelihood function'//{introspection:location})
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
       call Halo_Model_Projected_Correlation(                                         &
            &                                conditionalMassFunction_               , &
            &                                self%powerSpectrum_                    , &
            &                                self%cosmologyFunctions_               , &
            &                                self%surveyGeometry_                   , &
            &                                self%darkMatterHaloScale_              , &
            &                                self%haloMassFunction_                 , &
            &                                self%darkMatterProfileDMO_             , &
            &                                self%darkMatterHaloBias_               , &
            &                                self%darkMatterProfileScaleRadius_     , &
            &                                self%separation                        , &
            &                                self%massMinimum                  (  i), &
            &                                self%massMaximum                  (  i), &
            &                                self%haloMassMinimum                   , &
            &                                self%haloMassMaximum                   , &
            &                                self%lineOfSightDepth                  , &
            &                                self%halfIntegral                      , &
            &                                self%projectedCorrelationFunction (:,i)  &
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
    projectedCorrelationFunctionEvaluate=-0.5d0*self%covariance%covarianceProduct(difference)
    return
  end function projectedCorrelationFunctionEvaluate

  subroutine projectedCorrelationFunctionFunctionChanged(self)
    !!{
    Respond to possible changes in the likelihood function.
    !!}
    implicit none
    class(posteriorSampleLikelihoodPrjctdCorrelationFunction), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine projectedCorrelationFunctionFunctionChanged
