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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implementation of a posterior sampling likelihood class which evaluates the likelihood using emulators of model predictions.
  !!}

  use :: Statistics_Emulators, only : emulatorClass

  type, public :: emulatedEmulatorList
     class(emulatorClass       ), pointer :: emulator_ => null()
     type (emulatedEmulatorList), pointer :: next      => null()
  end type emulatedEmulatorList

  type :: emulatedEmulatorState
     !!{RST
     The state of one emulator used by the emulated likelihood: the index of the active parameter supplying each of its inputs,
     its target data, and which bins of the target data are included in the likelihood.
     !!}
     integer         , allocatable, dimension(:  ) :: indexParameter
     double precision, allocatable, dimension(:  ) :: yTarget
     double precision, allocatable, dimension(:,:) :: covarianceTarget
     logical         , allocatable, dimension(:  ) :: included
  end type emulatedEmulatorState

  !![
  <enumeration docformat="rst">
   <name>emulatedLikelihoodForm</name>
   <description>
   Used to specify the form of the likelihood of the data compared to an emulated observable.
   </description>
   <encodeFunction>yes</encodeFunction>
   <entry label="gaussianDiagonal"  />
   <entry label="gaussianCovariance"/>
  </enumeration>
  !!]

  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodEmulated" docformat="rst">
   <description>
   A likelihood evaluated using emulators of the model's predictions (see :ref:`manual-sec-Emulation`), in place of running
   the model. Each ``emulator`` predicts one observable, in the (possibly transformed) space in which it was trained, as a
   function of the prior quantiles of the parameters on which it depends. The likelihood is the product of the likelihoods of
   the target data of each observable, which are therefore treated as independent.

   Each emulator input is supplied by the active parameter of the same name, and every active parameter must supply an input
   to at least one emulator. On first use, the priors of the active parameters are checked against those under which the
   emulators were trained: the prior cumulative probability of the parameter value at each training point must equal the
   quantile at which the emulator was trained.

   The form of the likelihood for each emulator is given by the corresponding entry in ``[likelihoodForms]`` (a single entry
   applies to every emulator):

   ``gaussianDiagonal``
     A Gaussian likelihood neglecting covariance between bins,

     .. math::

        \log \mathcal{L} = -{1 \over 2} \sum_i \left[ {(y_i - \mu_i)^2 \over s_i^2} + \log (2 \pi s_i^2) \right],

     where :math:`y_i` is the target datum in bin :math:`i`, :math:`\mu_i` the emulator's mean prediction, and :math:`s_i^2
     = \mathsf{C}_{ii} + \sigma_i^2` the sum of the variance of the target datum, and the variance of the emulator's
     prediction.

   ``gaussianCovariance``
     A Gaussian likelihood including the full covariance, :math:`\mathsf{C}`, of the target data,

     .. math::

        \log \mathcal{L} = -{1 \over 2} \left[ \Delta^\mathrm{T} \mathsf{S}^{-1} \Delta + \log |\mathsf{S}| + n \log 2 \pi
        \right],

     where :math:`\Delta = \mathbf{y} - \boldsymbol{\mu}`, :math:`\mathsf{S} = \mathsf{C} + \mathrm{diag}(\sigma^2)`, and
     :math:`n` is the number of bins.

   If ``[includeEmulatorVariance]`` is false, the variance of the emulator's prediction, :math:`\sigma_i^2`, is omitted from
   :math:`s_i^2` and :math:`\mathsf{S}`. Bins whose target datum or variance is not finite (for example, bins of a mass
   function in which no galaxies were observed, after a logarithmic transformation), or whose total variance is not positive,
   are excluded.

   The variance of the log-likelihood due to the uncertainty of the emulation is estimated to first order as
   :math:`\sum_i g_i^2 \sigma_i^2`, where :math:`\mathbf{g} = \mathsf{S}^{-1} \Delta` is the derivative of the log-likelihood
   with respect to the predicted means.
   </description>
   <linkedList type="emulatedEmulatorList" variable="emulators" next="next" object="emulator_" objectType="emulatorClass"/>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodEmulated
     !!{RST
     Implementation of a posterior sampling likelihood class which evaluates the likelihood using emulators of model predictions.
     !!}
     private
     type   (emulatedEmulatorList                 ), pointer                   :: emulators               => null()
     type   (enumerationEmulatedLikelihoodFormType), allocatable, dimension(:) :: likelihoodForms
     type   (emulatedEmulatorState                ), allocatable, dimension(:) :: states
     logical                                                                   :: includeEmulatorVariance          , initialized
   contains
     !![
     <methods docformat="rst">
       <method method="initialize" description="Map the inputs of each emulator to the active parameters, and check the priors of those parameters."/>
     </methods>
     !!]
     final     ::                    emulatedDestructor
     procedure :: evaluate        => emulatedEvaluate
     procedure :: functionChanged => emulatedFunctionChanged
     procedure :: initialize      => emulatedInitialize
  end type posteriorSampleLikelihoodEmulated

  interface posteriorSampleLikelihoodEmulated
     !!{RST
     Constructors for the :galacticus-class:`posteriorSampleLikelihoodEmulated` posterior sampling likelihood class.
     !!}
     module procedure emulatedConstructorParameters
     module procedure emulatedConstructorInternal
  end interface posteriorSampleLikelihoodEmulated

  ! The tolerance in prior quantiles used when checking that the priors of the active parameters are those under which the
  ! emulators were trained.
  double precision, parameter :: emulatedToleranceQuantile=1.0d-6

contains

  function emulatedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`posteriorSampleLikelihoodEmulated` posterior sampling likelihood class which builds
    the object from a parameter set.
    !!}
    use :: Error             , only : Error_Report
    use :: Input_Parameters  , only : inputParameter, inputParameters
    use :: ISO_Varying_String, only : char          , var_str
    implicit none
    type   (posteriorSampleLikelihoodEmulated    )                              :: self
    type   (inputParameters                      ), intent(inout)               :: parameters
    type   (emulatedEmulatorList                 ), pointer                     :: emulators              , emulator_
    type   (varying_string                       ), allocatable  , dimension(:) :: likelihoodForms
    type   (enumerationEmulatedLikelihoodFormType), allocatable  , dimension(:) :: likelihoodForms_
    logical                                                                     :: includeEmulatorVariance
    integer                                                                     :: countEmulators         , countForms, &
         &                                                                         i

    countEmulators=parameters%copiesCount('emulator',zeroIfNotPresent=.true.)
    if (countEmulators < 1) call Error_Report('at least one emulator must be given'//{introspection:location})
    emulators => null()
    emulator_ => null()
    do i=1,countEmulators
       if (associated(emulator_)) then
          allocate(emulator_%next)
          emulator_ => emulator_%next
       else
          allocate(emulators)
          emulator_ => emulators
       end if
       !![
       <objectBuilder class="emulator" name="emulator_%emulator_" source="parameters" copy="i" />
       !!]
    end do
    countForms=parameters%count('likelihoodForms',zeroIfNotPresent=.true.)
    if (countForms == 0) then
       allocate(likelihoodForms(1))
       likelihoodForms(1)=var_str('gaussianDiagonal')
    else
       if (countForms /= 1 .and. countForms /= countEmulators) call Error_Report('the number of `likelihoodForms` must be 1 or equal to the number of emulators'//{introspection:location})
       allocate(likelihoodForms(countForms))
       !![
       <inputParameter docformat="rst">
         <name>likelihoodForms</name>
         <description>
         The form of the likelihood (``gaussianDiagonal`` or ``gaussianCovariance``) for each emulator, in order. A single value
         applies to every emulator. If not given, ``gaussianDiagonal`` is used for every emulator.
         </description>
         <source>parameters</source>
       </inputParameter>
       !!]
    end if
    allocate(likelihoodForms_(countEmulators))
    do i=1,countEmulators
       likelihoodForms_(i)=enumerationEmulatedLikelihoodFormEncode(char(likelihoodForms(min(i,size(likelihoodForms)))),includesPrefix=.false.)
    end do
    !![
    <inputParameter docformat="rst">
      <name>includeEmulatorVariance</name>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, the variance of each emulator's predictions is added to the variance of the target data in the likelihood.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=posteriorSampleLikelihoodEmulated(emulators,likelihoodForms_,includeEmulatorVariance)
    !![
    <inputParametersValidate source="parameters" multiParameters="emulator"/>
    !!]
    return
  end function emulatedConstructorParameters

  function emulatedConstructorInternal(emulators,likelihoodForms,includeEmulatorVariance) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`posteriorSampleLikelihoodEmulated` posterior sampling likelihood class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (posteriorSampleLikelihoodEmulated    )                              :: self
    type   (emulatedEmulatorList                 ), intent(in   ), target       :: emulators
    type   (enumerationEmulatedLikelihoodFormType), intent(in   ), dimension(:) :: likelihoodForms
    logical                                       , intent(in   )               :: includeEmulatorVariance
    type   (emulatedEmulatorList                 ), pointer                     :: emulator_
    integer                                                                     :: countEmulators
    !![
    <constructorAssign variables="likelihoodForms, includeEmulatorVariance"/>
    !!]

    self     %emulators => emulators
    emulator_           => emulators
    countEmulators      =  0
    do while (associated(emulator_))
       countEmulators=countEmulators+1
       !![
       <referenceCountIncrement owner="emulator_" object="emulator_"/>
       !!]
       emulator_ => emulator_%next
    end do
    if (size(likelihoodForms) /= countEmulators) call Error_Report('a likelihood form must be given for each emulator'//{introspection:location})
    self%initialized=.false.
    return
  end function emulatedConstructorInternal

  subroutine emulatedDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`posteriorSampleLikelihoodEmulated` posterior sampling likelihood class.
    !!}
    implicit none
    type(posteriorSampleLikelihoodEmulated), intent(inout) :: self
    type(emulatedEmulatorList             ), pointer       :: emulator_, emulatorNext

    emulator_ => self%emulators
    do while (associated(emulator_))
       emulatorNext => emulator_%next
       !![
       <objectDestructor name="emulator_%emulator_"/>
       !!]
       deallocate(emulator_)
       emulator_ => emulatorNext
    end do
    return
  end subroutine emulatedDestructor

  subroutine emulatedInitialize(self,modelParametersActive_)
    !!{RST
    Map the inputs of each emulator to the active parameters, check that the priors of those parameters are those under which
    the emulators were trained, and find the bins of each target dataset to include in the likelihood.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : operator(//)  , operator(==), var_str
    use :: String_Handling   , only : operator(//)
    implicit none
    class           (posteriorSampleLikelihoodEmulated), intent(inout)                 :: self
    type            (modelParameterList               ), intent(in   ), dimension(:  ) :: modelParametersActive_
    type            (emulatedEmulatorList             ), pointer                       :: emulator_
    type            (varying_string                   ), allocatable  , dimension(:  ) :: names
    double precision                                   , allocatable  , dimension(:,:) :: quantiles             , values
    logical                                            , allocatable  , dimension(:  ) :: used
    character       (len=12                           )                                :: label
    double precision                                                                   :: differenceMaximum
    integer                                                                            :: countEmulators        , k     , &
         &                                                                                i                     , j     , &
         &                                                                                n

    ! Count emulators.
    countEmulators=0
    emulator_ => self%emulators
    do while (associated(emulator_))
       countEmulators=countEmulators+1
       emulator_ => emulator_%next
    end do
    if (allocated(self%states)) deallocate(self%states)
    allocate(self%states(countEmulators))
    allocate(used(size(modelParametersActive_)))
    used=.false.
    k        =  0
    emulator_ => self%emulators
    do while (associated(emulator_))
       k=k+1
       associate (state => self%states(k))
         ! Map each input of the emulator to the active parameter of the same name.
         call emulator_%emulator_%inputNames(names)
         allocate(state%indexParameter(size(names)))
         state%indexParameter=0
         do i=1,size(names)
            do j=1,size(modelParametersActive_)
               if (modelParametersActive_(j)%modelParameter_%name() == names(i)) state%indexParameter(i)=j
            end do
            if (state%indexParameter(i) == 0) call Error_Report(var_str("emulator input '")//names(i)//"' is not an active parameter"//{introspection:location})
            used(state%indexParameter(i))=.true.
         end do
         ! Check that the prior of each parameter is that under which the emulator was trained.
         call emulator_%emulator_%trainingPoints(quantiles,values)
         do i=1,size(names)
            differenceMaximum=0.0d0
            do j=1,size(quantiles,dim=2)
               differenceMaximum=max(                                                                                                  &
                    &                differenceMaximum                                                                               , &
                    &                abs(                                                                                              &
                    &                    +modelParametersActive_(state%indexParameter(i))%modelParameter_%priorCumulative(values(i,j)) &
                    &                    -                                                                                quantiles(i,j)  &
                    &                   )                                                                                              &
                    &               )
            end do
            if (differenceMaximum > emulatedToleranceQuantile) then
               write (label,'(e12.4)') differenceMaximum
               call Error_Report(                                                                                                &
                    &            var_str("the prior of parameter '")//names(i)//"' differs from that under which the emulator " // &
                    &            "was trained (the prior quantiles of training points differ by up to "//trim(adjustl(label))// &
                    &            ")"//{introspection:location}                                                                     &
                    &           )
            end if
         end do
         ! Find the bins of the target data to include - those with finite data and variance.
         call emulator_%emulator_%target(state%yTarget,state%covarianceTarget)
         n=size(state%yTarget)
         if (n /= emulator_%emulator_%countOutputs()) call Error_Report('target data has the wrong number of bins'//{introspection:location})
         allocate(state%included(n))
         do i=1,n
            state%included(i)=emulatedIsFinite(state%yTarget(i)) .and. emulatedIsFinite(state%covarianceTarget(i,i))
         end do
         if (self%likelihoodForms(k) == emulatedLikelihoodFormGaussianCovariance) then
            do i=1,n
               if (.not.state%included(i)) cycle
               do j=1,n
                  if (state%included(j) .and. .not.emulatedIsFinite(state%covarianceTarget(i,j))) &
                       & call Error_Report(var_str('target covariance of emulator ')//k//' is not finite between included bins'//{introspection:location})
               end do
            end do
         end if
       end associate
       emulator_ => emulator_%next
    end do
    ! Every active parameter must be used by some emulator.
    do j=1,size(modelParametersActive_)
       if (.not.used(j)) call Error_Report(var_str("active parameter '")//modelParametersActive_(j)%modelParameter_%name()//"' is not an input to any emulator"//{introspection:location})
    end do
    self%initialized=.true.
    return
  end subroutine emulatedInitialize

  double precision function emulatedEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !!{RST
    Return the log-likelihood for the emulated likelihood function.
    !!}
    use :: Linear_Algebra                , only : assignment(=)                  , matrix, matrixCholesky, vector
    use :: Numerical_Constants_Math      , only : Pi
    use :: Posterior_Sampling_Convergence, only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State      , only : posteriorSampleStateClass
    implicit none
    class           (posteriorSampleLikelihoodEmulated), intent(inout), target         :: self
    class           (posteriorSampleStateClass        ), intent(inout)                 :: simulationState
    type            (modelParameterList               ), intent(inout), dimension(:  ) :: modelParametersActive_, modelParametersInactive_
    class           (posteriorSampleConvergenceClass  ), intent(inout)                 :: simulationConvergence
    double precision                                   , intent(in   )                 :: temperature           , logLikelihoodCurrent    , &
         &                                                                                logPriorCurrent       , logPriorProposed
    real                                               , intent(inout)                 :: timeEvaluate
    double precision                                   , intent(  out), optional       :: logLikelihoodVariance
    logical                                            , intent(inout), optional       :: forceAcceptance
    type            (emulatedEmulatorList             ), pointer                       :: emulator_
    double precision                                   , allocatable  , dimension(:  ) :: stateArray            , quantiles               , &
         &                                                                                mean                  , variance                , &
         &                                                                                residual              , gradient
    double precision                                   , allocatable  , dimension(:,:) :: covariance
    integer                                            , allocatable  , dimension(:  ) :: indices
    double precision                                                                   :: varianceTotal         , varianceLogLikelihood
    integer                                                                            :: i                     , k
    !$GLC attributes unused :: timeEvaluate, temperature, simulationConvergence, logPriorProposed, logPriorCurrent, logLikelihoodCurrent, modelParametersInactive_, forceAcceptance

    if (.not.self%initialized) call self%initialize(modelParametersActive_)
    ! Find the physical values of the parameters.
    allocate(stateArray(simulationState%dimension()))
    stateArray=simulationState%get()
    do i=1,size(stateArray)
       stateArray(i)=modelParametersActive_(i)%modelParameter_%unmap(stateArray(i))
    end do
    emulatedEvaluate     =0.0d0
    varianceLogLikelihood=0.0d0
    k                    =0
    emulator_ => self%emulators
    do while (associated(emulator_))
       k=k+1
       associate (state => self%states(k))
         ! Predict the observable at the prior quantiles of the parameters.
         allocate(quantiles(size(state%indexParameter)))
         do i=1,size(state%indexParameter)
            quantiles(i)=modelParametersActive_(state%indexParameter(i))%modelParameter_%priorCumulative(stateArray(state%indexParameter(i)))
         end do
         allocate(mean    (size(state%yTarget)))
         allocate(variance(size(state%yTarget)))
         call emulator_%emulator_%predict(quantiles,mean,variance)
         if (self%likelihoodForms(k) == emulatedLikelihoodFormGaussianDiagonal) then
            do i=1,size(state%yTarget)
               if (.not.state%included(i)) cycle
               varianceTotal=state%covarianceTarget(i,i)
               if (self%includeEmulatorVariance) varianceTotal=varianceTotal+variance(i)
               if (varianceTotal <= 0.0d0) cycle
               emulatedEvaluate     =+emulatedEvaluate                                    &
                    &                -0.5d0                                               &
                    &                *(                                                   &
                    &                  +(state%yTarget(i)-mean(i))**2/varianceTotal       &
                    &                  +log(2.0d0*Pi*varianceTotal)                       &
                    &                 )
               varianceLogLikelihood=+varianceLogLikelihood                               &
                    &                +((state%yTarget(i)-mean(i))/varianceTotal)**2       &
                    &                *variance(i)
            end do
         else
            ! Select the bins to include - those with finite target data, and positive total variance.
            allocate(indices(0))
            do i=1,size(state%yTarget)
               if (.not.state%included(i)) cycle
               varianceTotal=state%covarianceTarget(i,i)
               if (self%includeEmulatorVariance) varianceTotal=varianceTotal+variance(i)
               if (varianceTotal <= 0.0d0) cycle
               indices=[indices,i]
            end do
            if (size(indices) > 0) then
               covariance=state%covarianceTarget(indices,indices)
               if (self%includeEmulatorVariance) then
                  do i=1,size(indices)
                     covariance(i,i)=covariance(i,i)+variance(indices(i))
                  end do
               end if
               residual=state%yTarget(indices)-mean(indices)
               allocate(gradient(size(indices)))
               block
                 type(matrixCholesky) :: decomposition
                 decomposition=matrixCholesky(matrix(covariance))
                 gradient     =decomposition%squareSystemSolve(vector(residual))
                 emulatedEvaluate=+emulatedEvaluate                                &
                      &           -0.5d0                                           &
                      &           *(                                               &
                      &             +dot_product(residual,gradient)                &
                      &             +decomposition%logarithmicDeterminant()        &
                      &             +dble(size(indices))*log(2.0d0*Pi)             &
                      &            )
               end block
               varianceLogLikelihood=+varianceLogLikelihood                        &
                    &                +sum(gradient**2*variance(indices))
               deallocate(covariance,residual,gradient)
            end if
            deallocate(indices)
         end if
         deallocate(quantiles,mean,variance)
       end associate
       emulator_ => emulator_%next
    end do
    if (present(logLikelihoodVariance)) logLikelihoodVariance=varianceLogLikelihood
    return
  end function emulatedEvaluate

  elemental logical function emulatedIsFinite(x)
    !!{RST
    Return true if ``x`` is finite. Galacticus is compiled with ``-ffinite-math-only``, which allows the compiler to assume that
    ``ieee_is_finite()`` is always true, so the exponent bits are examined directly instead: they are all set only for
    infinities and not-a-number values.
    !!}
    use, intrinsic :: ISO_Fortran_Env, only : int64
    implicit none
    double precision, intent(in   ) :: x

    emulatedIsFinite=ibits(transfer(x,0_int64),52,11) /= 2047_int64
    return
  end function emulatedIsFinite

  subroutine emulatedFunctionChanged(self)
    !!{RST
    Respond to possible changes in the likelihood function.
    !!}
    implicit none
    class(posteriorSampleLikelihoodEmulated), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine emulatedFunctionChanged
