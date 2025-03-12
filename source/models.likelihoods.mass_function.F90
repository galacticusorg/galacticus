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
  Implementation of a posterior sampling likelihood class which implements a likelihood for mass functions.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  use :: Geometry_Surveys   , only : surveyGeometryClass
  use :: Halo_Mass_Functions, only : haloMassFunctionClass

  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodMassFunction">
   <description>
    The likelihood is computed as
    \begin{equation}
    \log \mathcal{L} = -{1\over2} \Delta \cdot \mathcal{C}^{-1} \cdot \Delta^\mathrm{T},
    \end{equation}
    where $\mathcal{C}$ is the covariance matrix, and $\Delta_i = w_i^\mathrm{model} - w_i^\mathrm{obs}$, $w_i^\mathrm{model}$
    is the computed mass function at the $i^\mathrm{th}$ separation, and $w_i^\mathrm{obs}$ is the observed mass function at
    the $i^\mathrm{th}$ separation. The mass function is computed using the halo model and the parameterized conditional galaxy
    mass function of \cite[][see also \protect\cite{leauthaud_new_2011};
    \protect\refPhysics{conditionalMassFunctionBehroozi2010}]{behroozi_comprehensive_2010}.  The details of the mass function calculation are specified by
    the following subparameters:
    \begin{description}
    \item[{\normalfont \ttfamily haloMass(Min|Max)imum}] The minimum/maximum halo mass over which to integrate in the halo model;
    \item[{\normalfont \ttfamily redshift(Min|Max)imum}] The minimum/maximum redshift over which to integrate in the halo model;
    \item[{\normalfont \ttfamily massFunctionFileName}] The name of an HDF5 file containing the observed mass function and its
      covariance matrix.
    \end{description}
    
    The HDF5 file specified by the {\normalfont \ttfamily massFunctionFileName} element should contain a {\normalfont \ttfamily mass}
    dataset, giving the masses at which the mass function is measured (in units of $M_\odot$), a {\normalfont \ttfamily
    massFunctionObserved} dataset giving the observed values of the mass function at those masses (in units of Mpc$^{-3}$ per
    $\log M$), and a {\normalfont \ttfamily covariance} dataset, giving the covariance of the mass function (in units of Mpc$^{-6}$).
   </description>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodMassFunction
     !!{
     Implementation of a posterior sampling likelihood class which implements a likelihood for mass functions.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer                     :: cosmologyFunctions_   => null()
     class           (haloMassFunctionClass  ), pointer                     :: haloMassFunction_     => null()
     class           (surveyGeometryClass    ), pointer                     :: surveyGeometry_       => null()
     logical                                                                :: useSurveyLimits                , modelSurfaceBrightness
     double precision                                                       :: haloMassMinimum                , haloMassMaximum       , &
          &                                                                    redshiftMinimum                , redshiftMaximum       , &
          &                                                                    logHaloMassMinimum             , logHaloMassMaximum    , &
          &                                                                    surfaceBrightnessLimit
     double precision                         , dimension(:  ), allocatable :: mass                           , massFunctionObserved  , &
          &                                                                    massMinimum                    , massMaximum           , &
          &                                                                    massFunction
     double precision                         , dimension(:,:), allocatable :: covarianceMatrix
     type            (vector                 )                              :: means
     type            (matrix                 )                              :: covariance
     type            (varying_string         )                              :: massFunctionFileName
   contains
     final     ::                    massFunctionDestructor
     procedure :: evaluate        => massFunctionEvaluate
     procedure :: functionChanged => massFunctionFunctionChanged
  end type posteriorSampleLikelihoodMassFunction

  interface posteriorSampleLikelihoodMassFunction
     !!{
     Constructors for the {\normalfont \ttfamily massFunction} posterior sampling convergence class.
     !!}
     module procedure massFunctionConstructorParameters
     module procedure massFunctionConstructorInternal
  end interface posteriorSampleLikelihoodMassFunction

contains

  function massFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily massFunction} posterior sampling convergence class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleLikelihoodMassFunction)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    double precision                                                       :: redshiftMinimum       , redshiftMaximum       , &
         &                                                                    haloMassMinimum       , haloMassMaximum       , &
         &                                                                    surfaceBrightnessLimit
    logical                                                                :: useSurveyLimits       , modelSurfaceBrightness
    type            (varying_string                       )                :: massFunctionFileName
    class           (cosmologyFunctionsClass              ), pointer       :: cosmologyFunctions_
    class           (haloMassFunctionClass                ), pointer       :: haloMassFunction_
    class           (surveyGeometryClass                  ), pointer       :: surveyGeometry_

    !![
    <inputParameter>
      <name>massFunctionFileName</name>
      <description>The name of the file containing the mass function.</description>
      <source>parameters</source>
    </inputParameter>
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
      <name>redshiftMinimum</name>
      <description>The minimum redshift over which to integrate.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshiftMaximum</name>
      <description>The maximum redshift over which to integrate.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>useSurveyLimits</name>
      <description>If true, limit redshift integration range based on survey geometry limits.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>modelSurfaceBrightness</name>
      <description>If true, model the effects of surface brightness incompleteness on the mass function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>surfaceBrightnessLimit</name>
      <description>The limiting surface brightness to which to integrate.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="haloMassFunction"   name="haloMassFunction_"   source="parameters"/>
    <objectBuilder class="surveyGeometry"     name="surveyGeometry_"     source="parameters"/>
    !!]
    self=posteriorSampleLikelihoodMassFunction(haloMassMinimum,haloMassMaximum,redshiftMinimum,redshiftMaximum,useSurveyLimits,char(massFunctionFileName),modelSurfaceBrightness,surfaceBrightnessLimit,cosmologyFunctions_,haloMassFunction_,surveyGeometry_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="haloMassFunction_"  />
    <objectDestructor name="surveyGeometry_"    />
    !!]
    return
  end function massFunctionConstructorParameters

  function massFunctionConstructorInternal(haloMassMinimum,haloMassMaximum,redshiftMinimum,redshiftMaximum,useSurveyLimits,massFunctionFileName,modelSurfaceBrightness,surfaceBrightnessLimit,cosmologyFunctions_,haloMassFunction_,surveyGeometry_) result(self)
    !!{
    Constructor for ``massFunction'' posterior sampling likelihood class.
    !!}
    use :: Display          , only : displayMessage, displayMagenta    , displayReset
    use :: Input_Paths      , only : inputPath     , pathTypeDataStatic
    use :: HDF5_Access      , only : hdf5Access
    use :: IO_HDF5          , only : hdf5Object
    use :: Linear_Algebra   , only : assignment(=)
    type            (posteriorSampleLikelihoodMassFunction)                              :: self
    double precision                                       , intent(in   )               :: redshiftMinimum        , redshiftMaximum       , &
         &                                                                                  haloMassMinimum        , haloMassMaximum       , &
         &                                                                                  surfaceBrightnessLimit
    logical                                                , intent(in   )               :: useSurveyLimits        , modelSurfaceBrightness
    character       (len=*                                ), intent(in   )               :: massFunctionFileName
    class           (cosmologyFunctionsClass              ), intent(in   ), target       :: cosmologyFunctions_
    class           (haloMassFunctionClass                ), intent(in   ), target       :: haloMassFunction_
    class           (surveyGeometryClass                  ), intent(in   ), target       :: surveyGeometry_
    double precision                                       , allocatable  , dimension(:) :: massBinWidth           , eigenValueArray
    type            (hdf5Object                           )                              :: massFunctionFile
    integer                                                                              :: i
    type            (matrix                               )                              :: eigenVectors
    type            (vector                               )                              :: eigenValues
    !![
    <constructorAssign variables="haloMassMinimum, haloMassMaximum, redshiftMinimum, redshiftMaximum, modelSurfaceBrightness, surfaceBrightnessLimit, useSurveyLimits, massFunctionFileName, *cosmologyFunctions_, *haloMassFunction_, *surveyGeometry_"/>
    !!]

    self%logHaloMassMinimum=log10(haloMassMinimum)
    self%logHaloMassMaximum=log10(haloMassMaximum)
    ! Read the mass function file.
    !$ call hdf5Access%set()
    massFunctionFile=hdf5Object(char(inputPath(pathTypeDataStatic))//massFunctionFileName,readOnly=.true.)
    call massFunctionFile%readDataset("mass"                ,self%mass                )
    call massFunctionFile%readDataset("massFunctionObserved",self%massFunctionObserved)
    call massFunctionFile%readDataset("covariance"          ,self%covarianceMatrix    )
    ! Compute mass bin limits.
    allocate(self%massFunction,mold=self%mass)
    allocate(self%massMinimum,mold=self%mass)
    allocate(self%massMaximum,mold=self%mass)
    if (massFunctionFile%hasDataset("massWidthObserved")) then
       call massFunctionFile%readDataset("massWidthObserved",massBinWidth)
       do i=1,size(self%mass)
          self%massMinimum(i)=self%mass(i)/sqrt(massBinWidth(i))
          self%massMaximum(i)=self%mass(i)*sqrt(massBinWidth(i))
       end do
       deallocate(massBinWidth)
    else
       do i=1,size(self%mass)
          if (i == 1) then
             self              %massMinimum(i  ) &
                  & =      self%mass       (i  ) &
                  & *sqrt(                       &
                  &        self%mass       (i  ) &
                  &       /self%mass       (i+1) &
                  &      )
          else
             self              %massMinimum(i  ) &
                  & =sqrt(                       &
                  &        self%mass       (i-1) &
                  &       *self%mass       (i  ) &
                  &      )
          end if
          if (i == size(self%mass)) then
             self              %massMaximum(i  ) &
                  & =      self%mass       (i  ) &
                  & *sqrt(                       &
                  &        self%mass       (i  ) &
                  &       /self%mass       (i-1) &
                  &      )
          else
             self              %massMaximum(i  ) &
                  & =sqrt(                       &
                  &        self%mass       (i  ) &
                  &       *self%mass       (i+1) &
                  &      )
          end if
       end do
    end if
    !$ call hdf5Access%unset()
    ! Find the inverse covariance matrix.
    self%covariance=self%covarianceMatrix
    ! Get eigenvalues and vectors of the covariance matrix.
    allocate(eigenValueArray(size(self%mass)))
    call self%covariance%eigenSystem(eigenVectors,eigenValues)
    eigenValueArray=eigenValues
    if (any(eigenValueArray < 0.0d0)) call displayMessage(displayMagenta()//'WARNING:'//displayReset()//' inverse covariance matrix is not semi-positive definite')
    deallocate(eigenValueArray)
    return
  end function massFunctionConstructorInternal

  subroutine massFunctionDestructor(self)
    !!{
    Destructor for ``massFunction'' posterior sampling likelihood class.
    !!}
    implicit none
    type(posteriorSampleLikelihoodMassFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    <objectDestructor name="self%haloMassFunction_"  />
    <objectDestructor name="self%surveyGeometry_"    />
    !!]
    return
  end subroutine massFunctionDestructor

  double precision function massFunctionEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !!{
    Return the log-likelihood for the mass function likelihood function.
    !!}
    use :: Conditional_Mass_Functions    , only : conditionalMassFunctionBehroozi2010
    use :: Error                         , only : Error_Report
    use :: Linear_Algebra                , only : assignment(=)                              , operator(*)
    use :: Mass_Function_Incompletenesses, only : massFunctionIncompletenessSurfaceBrightness
    use :: Models_Likelihoods_Constants  , only : logImpossible
    use :: Numerical_Integration         , only : integrator
    use :: Posterior_Sampling_Convergence, only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State      , only : posteriorSampleStateClass
    implicit none
    class           (posteriorSampleLikelihoodMassFunction      ), intent(inout), target       :: self
    class           (posteriorSampleStateClass                  ), intent(inout)               :: simulationState
    type            (modelParameterList                         ), intent(inout), dimension(:) :: modelParametersActive_         , modelParametersInactive_
    class           (posteriorSampleConvergenceClass            ), intent(inout)               :: simulationConvergence
    double precision                                             , intent(in   )               :: temperature                    , logLikelihoodCurrent    , &
         &                                                                                        logPriorCurrent                , logPriorProposed
    real                                                         , intent(inout)               :: timeEvaluate
    double precision                                             , intent(  out), optional     :: logLikelihoodVariance
    logical                                                      , intent(inout), optional     :: forceAcceptance
    double precision                                             , allocatable  , dimension(:) :: stateVector
    type            (conditionalMassFunctionBehroozi2010        )                              :: conditionalMassFunction_
    type            (massFunctionIncompletenessSurfaceBrightness)                              :: massFunctionIncompletenessModel
    type            (integrator                                 )                              :: integratorMassFunction         , integratorNormalization
    type            (vector                                     )                              :: difference
    integer                                                                                    :: i                              , j
    double precision                                                                           :: binTimeMinimum                 , binTimeMaximum          , &
         &                                                                                        timeMinimum                    , timeMaximum             , &
         &                                                                                        distanceMaximum                , time                    , &
         &                                                                                        massFunction                   , normalization
    !$GLC attributes unused :: logLikelihoodCurrent, logPriorCurrent, simulationConvergence, temperature, timeEvaluate, modelParametersInactive_, forceAcceptance

    ! There is no variance in our likelihood estimate.
    if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
    ! Do not evaluate if the proposed prior is impossible.
    if (logPriorProposed <= logImpossible) then
       massFunctionEvaluate=0.0d0
       return
    end if
    ! Construct the conditional mass function object.
    stateVector=simulationState%get()
    if     (                                                                  &
         &   (.not.self%modelSurfaceBrightness .and. size(stateVector) /= 11) &
         &  .or.                                                              &
         &   (     self%modelSurfaceBrightness .and. size(stateVector) /= 14) &
         & )                                                                  &
         & call Error_Report('11 or 14 parameters are required for this likelihood function'//{introspection:location})
    do i=1,size(stateVector)
       stateVector(i)=modelParametersActive_(i)%modelParameter_%unmap(stateVector(i))
    end do
    conditionalMassFunction_=conditionalMassFunctionBehroozi2010(                 &
         &                                                       stateVector( 1), &
         &                                                       stateVector( 2), &
         &                                                       stateVector( 3), &
         &                                                       stateVector( 4), &
         &                                                       stateVector( 5), &
         &                                                       stateVector( 6), &
         &                                                       stateVector( 7), &
         &                                                       stateVector( 8), &
         &                                                       stateVector( 9), &
         &                                                       stateVector(10), &
         &                                                       stateVector(11)  &
         &                                                      )
    ! Extract surface brightness parameters and get an incompleteness model.
    if (self%modelSurfaceBrightness)                                                  &
         & massFunctionIncompletenessModel                                            &
         &  =massFunctionIncompletenessSurfaceBrightness(                             &
         &                                               self%surfaceBrightnessLimit, &
         &                                               1.0d0                      , &
         &                                               stateVector(12)            , &
         &                                               stateVector(13)            , &
         &                                               stateVector(14)              &
         &                                               )
    ! Compute the time corresponding to the specified redshifts.
    timeMinimum=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%redshiftMaximum))
    timeMaximum=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%redshiftMinimum))
    ! Compute the mass function.
    integratorMassFunction =integrator(likelihoodMassFunctionTimeIntegrand             ,toleranceRelative=1.0d-3)
    integratorNormalization=integrator(likelihoodMassFunctionTimeNormalizationIntegrand,toleranceRelative=1.0d-3)
    do i=1,size(self%mass)
       ! A survey geometry is imposed - sum over fields.
       massFunction =0.0d0
       normalization=0.0d0
       do j=1,self%surveyGeometry_%fieldCount()
          ! Find the maximum distance at which a galaxy of the present mass can be detected in this survey.
          distanceMaximum=self%surveyGeometry_%distanceMaximum(sqrt(self%massMinimum(i)*self%massMaximum(i)),field=j)
          ! Set integration limits appropriately.
          binTimeMaximum=                                                            timeMaximum
          binTimeMinimum=min(                                                                      &
               &                                                                     timeMaximum , &
               &             max(                                                                  &
               &                                                                     timeMinimum , &
               &                 self%cosmologyFunctions_%timeAtDistanceComoving(distanceMaximum)  &
               &                )                                                                  &
               &            )
          ! Integrate the mass function over this time interval.
          massFunction =+massFunction                                                      &
               &        +self%surveyGeometry_   %solidAngle(j                            ) &
               &        *integratorMassFunction %integrate (binTimeMinimum,binTimeMaximum)
          normalization=+normalization                                                     &
               &        +self%surveyGeometry_   %solidAngle(j                            ) &
               &        *integratorNormalization%integrate (binTimeMinimum,binTimeMaximum)
       end do
       self%massFunction(i)=                          &
            &               +massFunction             &
            &               /normalization            &
            &               /log(                     &
            &                     self%massMaximum(i) &
            &                    /self%massMinimum(i) &
            &                   )
       ! Account for surface brightness incompleteness if required.
       if (self%modelSurfaceBrightness) self%massFunction(i)=self%massFunction(i)*massFunctionIncompletenessModel%completeness(self%mass(i))
    end do
    deallocate(stateVector)
    ! Evaluate the log-likelihood.
    difference          =self%massFunction-self%massFunctionObserved
    massFunctionEvaluate=-0.5d0*self%covariance%covarianceProduct(difference)
    return

  contains

    double precision function likelihoodMassFunctionTimeIntegrand(timePrime)
      !!{
      Integral over time.
      !!}
      use :: Error          , only : Error_Report, errorStatusSuccess
      use :: String_Handling, only : operator(//)
      implicit none
      double precision                , intent(in   ) :: timePrime
      type            (integrator    )                :: integrator_
      integer                                         :: errorStatus
      type            (varying_string)                :: message
      character       (len=14        )                :: label
      double precision                                :: haloMassMinimum, haloMassMaximum

      ! Check for zero contribution from the ends of our halo mass range.
      haloMassMinimum=10.0d0**self%logHaloMassMinimum
      haloMassMaximum=10.0d0**self%logHaloMassMaximum
      if          (                                                                                  &
           &        conditionalMassFunction_%massFunction(      haloMassMinimum,self%massMinimum(i)) &
           &         ==                                                                              &
           &        conditionalMassFunction_%massFunction(      haloMassMinimum,self%massMaximum(i)) &
           &       .and.                                                                             &
           &        conditionalMassFunction_%massFunction(      haloMassMinimum,self%massMinimum(i)) &
           &         ==                                                                              &
           &        0.0d0                                                                            &
           &       .and.                                                                             &
           &        conditionalMassFunction_%massFunction(      haloMassMaximum,self%massMinimum(i)) &
           &         ==                                                                              &
           &        conditionalMassFunction_%massFunction(      haloMassMaximum,self%massMaximum(i)) &
           &       .and.                                                                             &
           &        conditionalMassFunction_%massFunction(      haloMassMaximum,self%massMinimum(i)) &
           &         >                                                                               &
           &        0.0d0                                                                            &
           &      ) then
         do while (                                                                                  &
              &     conditionalMassFunction_%massFunction(2.0d0*haloMassMinimum,self%massMinimum(i)) &
              &      ==                                                                              &
              &     conditionalMassFunction_%massFunction(2.0d0*haloMassMinimum,self%massMaximum(i)) &
              &    .and.                                                                             &
              &     conditionalMassFunction_%massFunction(2.0d0*haloMassMinimum,self%massMinimum(i)) &
              &      ==                                                                              &
              &     0.0d0                                                                            &
              &   )
            haloMassMinimum=2.0d0*haloMassMinimum
         end do
         do while (                                                                                  &
              &     conditionalMassFunction_%massFunction(0.5d0*haloMassMaximum,self%massMinimum(i)) &
              &      ==                                                                              &
              &     conditionalMassFunction_%massFunction(0.5d0*haloMassMaximum,self%massMaximum(i)) &
              &    .and.                                                                             &
              &     conditionalMassFunction_%massFunction(0.5d0*haloMassMaximum,self%massMinimum(i)) &
              &      >                                                                               &
              &     0.0d0                                                                            &
              &   )
            haloMassMaximum=0.5d0*haloMassMaximum
         end do
      end if
      time                =timePrime
      integrator_                        = integrator                                               (likelihoodMassFunctionHaloMassIntegrand,toleranceRelative=1.0d-3                ,toleranceAbsolute=1.0d-9     )
      likelihoodMassFunctionTimeIntegrand=+integrator_                    %integrate                (log10(haloMassMinimum)                 ,                  log10(haloMassMaximum),status           =errorStatus) &
           &                              *self       %cosmologyFunctions_%comovingVolumeElementTime(timePrime                                                                                                     )
      if (errorStatus /= errorStatusSuccess) then
         message='integration failed - state vector follows'
         do i=1,size(stateVector)
            write (label,'(e12.6)') stateVector(i)
            message=message//char(10)//" state ["//i//"] = "//trim(label)
         end do
         call Error_Report(message//{introspection:location})
      end if
      return
    end function likelihoodMassFunctionTimeIntegrand

    double precision function likelihoodMassFunctionTimeNormalizationIntegrand(timePrime)
      !!{
      Normalization integral over time.
      !!}
      implicit none
      double precision, intent(in   ) :: timePrime

      likelihoodMassFunctionTimeNormalizationIntegrand=self%cosmologyFunctions_%comovingVolumeElementTime(timePrime)
      return
    end function likelihoodMassFunctionTimeNormalizationIntegrand

    double precision function likelihoodMassFunctionHaloMassIntegrand(logMass)
      !!{
      Integral over halo mass function.
      !!}
      implicit none
      double precision, intent(in   ) :: logMass
      double precision                :: mass

      mass=10.0d0**logMass
      likelihoodMassFunctionHaloMassIntegrand=+self%haloMassFunction_%differential(time,mass)                        &
           &                                  *                                         mass                         &
           &                                  *log(10.0d0)                                                           &
           &                                  *max(                                                                  &
           &                                       +0.0d0                                                          , &
           &                                       +conditionalMassFunction_%massFunction(mass,self%massMinimum(i))  &
           &                                       -conditionalMassFunction_%massFunction(mass,self%massMaximum(i))  &
           &                                      )
      return
    end function likelihoodMassFunctionHaloMassIntegrand

  end function massFunctionEvaluate

  subroutine massFunctionFunctionChanged(self)
    !!{
    Respond to possible changes in the likelihood function.
    !!}
    implicit none
    class(posteriorSampleLikelihoodMassFunction), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine massFunctionFunctionChanged
