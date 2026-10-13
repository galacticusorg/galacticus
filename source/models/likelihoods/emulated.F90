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
     its target data (and, for the ``poissonCounts`` form, target counts and the factors converting values to counts), whether
     it emulates the logarithm of the observable, and which bins of the target data are included in the likelihood.
     !!}
     integer         , allocatable, dimension(:  ) :: indexParameter
     double precision, allocatable, dimension(:  ) :: yTarget         , countTarget, &
          &                                           countPerValue
     double precision, allocatable, dimension(:,:) :: covarianceTarget
     logical         , allocatable, dimension(:  ) :: included
     logical                                       :: transformLog
     double precision                              :: floor           , rootVarianceFloored
     integer                                       :: unitSamples
     type            (varying_string)              :: fileNameSamples
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
   <entry label="poissonCounts"     />
   <entry label="logNormalFloored"  />
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

   ``poissonCounts``
     A likelihood of the target counts, :math:`n_i`, in each bin (as used by the
     :galacticus-class:`outputAnalysisProgenitorMassFunction` analysis with ``likelihoodInCounts=true``). The emulator must
     have been trained on an observable whose training set records target counts, and the factor, :math:`f_i`, converting
     the (untransformed) value of the observable to an expected count. The expected count is :math:`\lambda_i = f_i
     \bar{\phi}_i`, where :math:`\bar{\phi}_i` is the mean of the emulator's prediction of the untransformed observable. The
     likelihood is a negative binomial (an over-dispersed Poisson distribution),

     .. math::

        \log \mathcal{L}_i = n_i \log \lambda_i - \log n_i! + \log \Gamma(r_i+n_i) - \log \Gamma(r_i) - n_i \log
        (r_i+\lambda_i) - r_i \log (1+\lambda_i/r_i),

     with :math:`r_i = 1/c_i^2`, which reduces to a Poisson distribution if :math:`c_i=0`. Here :math:`c_i^2 = C^2 + v_i + C^2
     v_i`, where :math:`C` is ``[rootVarianceTargetFractional]`` (the fractional root variance of the expected count), and
     :math:`v_i` is the fractional variance of the emulator's prediction of the untransformed observable (:math:`v_i =
     \exp[(\sigma_i \log 10)^2]-1` if the logarithm of the observable is emulated, or :math:`v_i = \sigma_i^2/\mu_i^2`
     otherwise). This approximates the marginalization of the likelihood over the uncertainty of the emulation by matching
     the first two moments of the distribution of counts. If the expected count is not positive, a bin with no target counts
     is ignored, while a bin with target counts makes the model impossible. Bins whose target count is not finite are
     excluded.

     ``[rootVarianceTargetFractional]`` gives :math:`C` for each emulator (a single value applies to every emulator). If
     ``[rootVarianceTargetFractionalParameter]`` is given, :math:`C` is instead the value of the active parameter of that
     name, which is then a parameter of the likelihood rather than an input to any emulator.

   ``logNormalFloored``
     A Gaussian likelihood, neglecting covariance between bins, of an observable emulated in its logarithm with a floor,
     :math:`F`, at which values were set in training (with root variance :math:`\sigma_F`; both are recorded in the training
     set). Bins at or below the floor, in the target data or in the model, are treated as in training: the target datum is
     :math:`y_i = \max(\log_{10} D_i, F)` and the model prediction :math:`m_i = \max(\mu_i, F)`, so that

     .. math::

        \log \mathcal{L} = -{1 \over 2} \sum_i \left[ {(y_i - m_i)^2 \over s_i^2} + \log (2 \pi s_i^2) \right],

     with :math:`s_i^2 = \mathsf{C}_{ii} + \sigma_i^2`, plus :math:`\sigma_F^2` if the model prediction is at the floor (an
     "empty" prediction is only an upper limit). A target datum at or below the floor (including a bin in which nothing was
     observed, whose logarithm is :math:`-\infty`) is given variance :math:`\sigma_F^2` in place of :math:`\mathsf{C}_{ii}`.
     Bins with no target datum (not-a-number) are excluded. Empty bins therefore contribute softly, rather than being
     excluded (as in ``gaussianDiagonal``) or dominating the likelihood.

   If ``[includeEmulatorVariance]`` is false, the variance of the emulator's prediction, :math:`\sigma_i^2`, is omitted from
   :math:`s_i^2` and :math:`\mathsf{S}` (and :math:`v_i=0`). Bins whose target datum or variance is not finite (for example,
   bins of a mass function in which no galaxies were observed, after a logarithmic transformation), or whose total variance
   is not positive, are excluded.

   The variance of the log-likelihood due to the uncertainty of the emulation is estimated to first order as
   :math:`\sum_i g_i^2 \sigma_i^2`, where :math:`g_i` is the derivative of the log-likelihood with respect to the predicted
   mean in bin :math:`i` (:math:`\mathbf{g} = \mathsf{S}^{-1} \Delta` for the Gaussian forms).

   The likelihood is not evaluated for states outside the support of the prior (at which the prior quantiles, the inputs to
   the emulators, are undefined).

   If ``[pathSamples]`` is not ``none``, the predictions of each emulator at every evaluation are written to a file
   ``[pathSamples]/&lt;label&gt;_&lt;process&gt;.txt`` (one per emulator, labeled by its observable, per :term:`MPI` process), one row
   per evaluation: the simulation step and chain index, then the predicted mean of each bin, then its variance (in the space
   in which the observable is emulated). Rows are written for every evaluated proposal, not only for accepted states (the
   state proposed at a rejected step appears in no other output unless ``[logProposals]`` is set on the simulation), and
   must be attributed using the chain index rather than the file name.
   </description>
   <linkedList type="emulatedEmulatorList" variable="emulators" next="next" object="emulator_" objectType="emulatorClass"/>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodEmulated
     !!{RST
     Implementation of a posterior sampling likelihood class which evaluates the likelihood using emulators of model predictions.
     !!}
     private
     type   (emulatedEmulatorList                 ), pointer                   :: emulators                             => null()
     type   (enumerationEmulatedLikelihoodFormType), allocatable, dimension(:) :: likelihoodForms
     type   (emulatedEmulatorState                ), allocatable, dimension(:) :: states
     double precision                              , allocatable, dimension(:) :: rootVarianceTargetFractional
     type   (varying_string                       )                            :: rootVarianceTargetFractionalParameter
     integer                                                                   :: indexRootVarianceTargetFractional
     logical                                                                   :: includeEmulatorVariance                        , initialized          , &
          &                                                                       appendSamples
     type   (varying_string                       )                            :: pathSamples
   contains
     !![
     <methods docformat="rst">
       <method method="initialize" description="Map the inputs of each emulator to the active parameters, and check the priors of those parameters."/>
     </methods>
     !!]
     final     ::                    emulatedDestructor
     procedure :: evaluate        => emulatedEvaluate
     procedure :: willEvaluate    => emulatedWillEvaluate
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
    type            (posteriorSampleLikelihoodEmulated    )                              :: self
    type            (inputParameters                      ), intent(inout)               :: parameters
    type            (emulatedEmulatorList                 ), pointer                     :: emulators                            , emulator_
    type            (varying_string                       ), allocatable  , dimension(:) :: likelihoodForms
    type            (enumerationEmulatedLikelihoodFormType), allocatable  , dimension(:) :: likelihoodForms_
    double precision                                       , allocatable  , dimension(:) :: rootVarianceTargetFractional
    type            (varying_string                       )                              :: rootVarianceTargetFractionalParameter, pathSamples
    logical                                                                              :: includeEmulatorVariance              , appendSamples
    integer                                                                              :: countEmulators                       , countForms, &
         &                                                                                  countRootVariances                   , i

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
         The form of the likelihood (``gaussianDiagonal``, ``gaussianCovariance``, ``poissonCounts``, or ``logNormalFloored``)
         for each emulator, in order. A single value applies to every emulator. If not given, ``gaussianDiagonal`` is used for
         every emulator.
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
    countRootVariances=parameters%count('rootVarianceTargetFractional',zeroIfNotPresent=.true.)
    if (countRootVariances == 0) then
       allocate(rootVarianceTargetFractional(countEmulators))
       rootVarianceTargetFractional=0.0d0
    else
       if (countRootVariances /= 1 .and. countRootVariances /= countEmulators) call Error_Report('the number of `rootVarianceTargetFractional` must be 1 or equal to the number of emulators'//{introspection:location})
       allocate(rootVarianceTargetFractional(countRootVariances))
       !![
       <inputParameter docformat="rst">
         <name>rootVarianceTargetFractional</name>
         <description>
         The fractional root variance of the expected count, :math:`C`, for each emulator using the ``poissonCounts`` form, in
         order. A single value applies to every emulator. If not given, zero is used for every emulator.
         </description>
         <source>parameters</source>
       </inputParameter>
       !!]
       if (countRootVariances == 1) rootVarianceTargetFractional=[(rootVarianceTargetFractional(1),i=1,countEmulators)]
    end if
    !![
    <inputParameter docformat="rst">
      <name>rootVarianceTargetFractionalParameter</name>
      <defaultValue>var_str('')</defaultValue>
      <description>
      If not empty, the name of the active parameter whose value gives the fractional root variance of the expected count,
      :math:`C`, for every emulator using the ``poissonCounts`` form (in place of ``[rootVarianceTargetFractional]``).
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>pathSamples</name>
      <defaultValue>var_str('none')</defaultValue>
      <description>
      The path into which the predictions of the emulators at each evaluation should be written. If ``none``, they are not
      written.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>appendSamples</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, append to existing files of predictions (see ``[pathSamples]``) rather than replacing them.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=posteriorSampleLikelihoodEmulated(emulators,likelihoodForms_,includeEmulatorVariance,rootVarianceTargetFractional,rootVarianceTargetFractionalParameter,pathSamples,appendSamples)
    !![
    <inputParametersValidate source="parameters" multiParameters="emulator"/>
    !!]
    return
  end function emulatedConstructorParameters

  function emulatedConstructorInternal(emulators,likelihoodForms,includeEmulatorVariance,rootVarianceTargetFractional,rootVarianceTargetFractionalParameter,pathSamples,appendSamples) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`posteriorSampleLikelihoodEmulated` posterior sampling likelihood class.
    !!}
    use :: Error             , only : Error_Report
    use :: File_Utilities    , only : Directory_Make
    use :: ISO_Varying_String, only : operator(//)  , operator(/=), char, &
         &                            assignment(=)
    use :: MPI_Utilities     , only : mpiSelf
    implicit none
    type            (posteriorSampleLikelihoodEmulated    )                              :: self
    type            (emulatedEmulatorList                 ), intent(in   ), target       :: emulators
    type            (enumerationEmulatedLikelihoodFormType), intent(in   ), dimension(:) :: likelihoodForms
    logical                                                , intent(in   )               :: includeEmulatorVariance
    double precision                                       , intent(in   ), dimension(:) :: rootVarianceTargetFractional
    type            (varying_string                       ), intent(in   )               :: rootVarianceTargetFractionalParameter
    type            (varying_string                       ), intent(in   ), optional     :: pathSamples
    logical                                                , intent(in   ), optional     :: appendSamples
    type            (emulatedEmulatorList                 ), pointer                     :: emulator_
    double precision                                       , allocatable  , dimension(:) :: outputs
    type            (varying_string                       )                              :: rankLabel
    integer                                                                              :: countEmulators                       , k   , &
         &                                                                                  unitSamples
    !![
    <constructorAssign variables="likelihoodForms, includeEmulatorVariance, rootVarianceTargetFractional, rootVarianceTargetFractionalParameter"/>
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
    if (size(likelihoodForms             ) /= countEmulators ) call Error_Report('a likelihood form must be given for each emulator'                   //{introspection:location})
    if (size(rootVarianceTargetFractional) /= countEmulators ) call Error_Report('a fractional root variance of counts must be given for each emulator'//{introspection:location})
    if (any (rootVarianceTargetFractional  <  0.0d0         )) call Error_Report('fractional root variances of counts must be non-negative'            //{introspection:location})
    self%indexRootVarianceTargetFractional=0
    self%initialized=.false.
    self%pathSamples  ='none'
    self%appendSamples=.false.
    if (present(pathSamples  )) self%pathSamples  =pathSamples
    if (present(appendSamples)) self%appendSamples=appendSamples
    allocate(self%states(countEmulators))
    ! Initialize the files of predictions, if needed.
    if (self%pathSamples /= 'none') then
       call Directory_Make(self%pathSamples)
       rankLabel=mpiSelf%rankLabel()
       emulator_ => emulators
       k         =  0
       do while (associated(emulator_))
          k=k+1
          self%states(k)%fileNameSamples=self%pathSamples//'/'//emulator_%emulator_%labelObservable()//'_'//rankLabel//'.txt'
          if (.not.self%appendSamples) then
             outputs=emulator_%emulator_%outputs()
             open(newUnit=unitSamples,file=char(self%states(k)%fileNameSamples),form='formatted',status='unknown')
             write (unitSamples,'(a,a,a,a)') '# Emulated predictions of ',char(emulator_%emulator_%labelObservable()),' written by process ',char(rankLabel)
             write (unitSamples,'(a,i6,a)' ) '# One row per likelihood evaluation. Columns: simulation step, chain index, then the predicted mean of each of the ',size(outputs),' bins below,'
             write (unitSamples,'(a)'      ) '# then the variance of each, in the space in which the observable is emulated. Records must be attributed using the'
             write (unitSamples,'(a)'      ) '# chain index, not the file name, as any process may evaluate any chain.'
             write (unitSamples,*          ) '# Abscissae: ',outputs
             close(unitSamples)
          end if
          emulator_ => emulator_%next
       end do
    end if
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
    use :: ISO_Varying_String, only : operator(//), operator(==), var_str, len
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
    type            (varying_string                   ), allocatable  , dimension(:  ) :: fileNamesSamples

    ! Count emulators.
    countEmulators =  0
    emulator_      => self%emulators
    do while (associated(emulator_))
       countEmulators=countEmulators+1
       emulator_ => emulator_%next
    end do
    ! Reset the states of the emulators, retaining the names of their files of predictions (set on construction).
    allocate(fileNamesSamples(countEmulators))
    if (allocated(self%states)) then
       if (size(self%states) == countEmulators) fileNamesSamples=self%states%fileNameSamples
       deallocate(self%states)
    end if
    allocate(self%states(countEmulators))
    self%states%fileNameSamples=fileNamesSamples
    allocate(used(size(modelParametersActive_)))
    used=.false.
    k         =  0
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
               differenceMaximum=max(                                                                                                     &
                    &                differenceMaximum                                                                                  , &
                    &                abs(                                                                                                 &
                    &                    +modelParametersActive_(state%indexParameter(i))%modelParameter_%priorCumulative(values   (i,j)) &
                    &                    -                                                                                quantiles(i,j)  &
                    &                   )                                                                                                 &
                    &               )
            end do
            if (differenceMaximum > emulatedToleranceQuantile) then
               write (label,'(e12.4)') differenceMaximum
               call Error_Report(                                                                                                 &
                    &            var_str("the prior of parameter '")//names(i)//"' differs from that under which the emulator "// &
                    &            "was trained (the prior quantiles of training points differ by up to "//trim(adjustl(label))  // &
                    &            ")"//{introspection:location}                                                                    &
                    &           )
            end if
         end do
         ! Find the bins of the target data to include - those with finite data and variance.
         call emulator_%emulator_%target(state%yTarget,state%covarianceTarget)
         n=size(state%yTarget)
         if (n /= emulator_%emulator_%countOutputs()) call Error_Report('target data has the wrong number of bins'//{introspection:location})
         allocate(state%included(n))
         state%transformLog=emulator_%emulator_%transform() == 'log10'
         if (self%likelihoodForms(k) == emulatedLikelihoodFormPoissonCounts) then
            ! For the counts form, the included bins are those with finite target counts and conversion factors.
            call emulator_%emulator_%counts(state%countTarget,state%countPerValue)
            if (size(state%countTarget) /= n .or. size(state%countPerValue) /= n) call Error_Report('target counts have the wrong number of bins'//{introspection:location})
            do i=1,n
               state%included(i)=emulatedIsFinite(state%countTarget(i)) .and. emulatedIsFinite(state%countPerValue   (i  ))
            end do
         else if (self%likelihoodForms(k) == emulatedLikelihoodFormLogNormalFloored) then
            ! For the floored form, the observable must be emulated in its logarithm with a floor. Bins whose target is at or below
            ! the floor (including those in which nothing was observed, whose logarithm is -infinity) are included - only bins
            ! with no target datum (not-a-number) are excluded.
            if (.not.state%transformLog) call Error_Report(var_str('the `logNormalFloored` form requires an observable emulated in its logarithm, but emulator ')//k//' is not'//{introspection:location})
            call emulator_%emulator_%floor(state%floor,state%rootVarianceFloored)
            do i=1,n
               state%included(i)=.not.emulatedIsNotANumber(state%yTarget(i)) .and. emulatedIsFinite(state%covarianceTarget(i,i))
            end do
         else
            do i=1,n
               state%included(i)=emulatedIsFinite(state%yTarget    (i)) .and. emulatedIsFinite(state%covarianceTarget(i,i))
            end do
         end if
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
    ! Find the active parameter giving the fractional root variance of counts, if any. It is a parameter of the likelihood, so
    ! must not also be an input to an emulator.
    self%indexRootVarianceTargetFractional=0
    if (len(self%rootVarianceTargetFractionalParameter) > 0) then
       do j=1,size(modelParametersActive_)
          if (modelParametersActive_(j)%modelParameter_%name() == self%rootVarianceTargetFractionalParameter) self%indexRootVarianceTargetFractional=j
       end do
       if (self%indexRootVarianceTargetFractional == 0) call Error_Report(var_str("parameter '")//self%rootVarianceTargetFractionalParameter//"' (giving the fractional root variance of counts) is not an active parameter"      //{introspection:location})
       if (used(self%indexRootVarianceTargetFractional)) call Error_Report(var_str("parameter '")//self%rootVarianceTargetFractionalParameter//"' gives the fractional root variance of counts, so can not also be an input to an emulator"//{introspection:location})
       if (.not.any([(self%likelihoodForms(i) == emulatedLikelihoodFormPoissonCounts,i=1,size(self%likelihoodForms))])) call Error_Report('`rootVarianceTargetFractionalParameter` is given, but no emulator uses the `poissonCounts` likelihood form'//{introspection:location})
       used(self%indexRootVarianceTargetFractional)=.true.
    end if
    ! Every active parameter must be used by some emulator, or by the likelihood.
    do j=1,size(modelParametersActive_)
       if (.not.used(j)) call Error_Report(var_str("active parameter '")//modelParametersActive_(j)%modelParameter_%name()//"' is not an input to any emulator, or a parameter of the likelihood"//{introspection:location})
    end do
    self%initialized=.true.
    return
  end subroutine emulatedInitialize

  double precision function emulatedEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !!{RST
    Return the log-likelihood for the emulated likelihood function.
    !!}
    use :: Factorials                    , only : Logarithmic_Factorial
    use :: Gamma_Functions               , only : Gamma_Function_Logarithmic
    use :: ISO_Varying_String            , only : char                           , operator(/=)
    use :: Linear_Algebra                , only : assignment(=)                  , matrix, matrixCholesky, vector
    use :: Models_Likelihoods_Constants  , only : logImpossible                  , logImprobable
    use :: Numerical_Constants_Math      , only : Pi                             , ln10
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
    double precision                                                                   :: varianceTotal         , varianceLogLikelihood   , &
         &                                                                                rootVarianceFractional, countExpected           , &
         &                                                                                varianceFractional    , dispersion              , &
         &                                                                                valueMean             , derivativeCount         , &
         &                                                                                countTarget           , gradientCount           , &
         &                                                                                valueTarget           , valueModel              , &
         &                                                                                derivativeModel
    logical                                                                            :: targetFloored
    integer                                                                            :: i                     , k
    !$GLC attributes unused :: timeEvaluate, temperature, simulationConvergence, logPriorCurrent, logLikelihoodCurrent, modelParametersInactive_, forceAcceptance

    ! States outside the support of the prior are not evaluated - the prior quantiles, the inputs to the emulators, are undefined
    ! there.
    if (.not.self%willEvaluate(simulationState,modelParametersActive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)) then
       emulatedEvaluate=logImpossible
       if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
       return
    end if
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
         ! Write the predictions to file if needed.
         if (self%pathSamples /= 'none') then
            open(newUnit=state%unitSamples,file=char(state%fileNameSamples),form='formatted',status='unknown',position='append')
            write (state%unitSamples,*) simulationState%count(),simulationState%chainIndex(),mean,variance
            close(state%unitSamples)
         end if
         if (self%likelihoodForms(k) == emulatedLikelihoodFormGaussianDiagonal) then
            do i=1,size(state%yTarget)
               if (.not.state%included(i)) cycle
               varianceTotal=state%covarianceTarget(i,i)
               if (self%includeEmulatorVariance) varianceTotal=varianceTotal+variance(i)
               if (varianceTotal <= 0.0d0) cycle
               emulatedEvaluate     =+emulatedEvaluate                              &
                    &                -0.5d0                                         &
                    &                *(                                             &
                    &                  +(state%yTarget(i)-mean(i))**2/varianceTotal &
                    &                  +log(2.0d0*Pi*varianceTotal)                 &
                    &                 )
               varianceLogLikelihood=+varianceLogLikelihood                         &
                    &                +((state%yTarget(i)-mean(i))/varianceTotal)**2 &
                    &                *variance(i)
            end do
         else if (self%likelihoodForms(k) == emulatedLikelihoodFormLogNormalFloored) then
            do i=1,size(state%yTarget)
               if (.not.state%included(i)) cycle
               ! The target datum, floored. (The target may be -infinity, so is compared with the floor only once known to be
               ! finite - with -ffinite-math-only the comparison of an infinity can not be relied upon.)
               targetFloored=.true.
               if (emulatedIsFinite(state%yTarget(i))) targetFloored=state%yTarget(i) <= state%floor
               if (targetFloored) then
                  valueTarget   =state%floor
                  varianceTotal =state%rootVarianceFloored**2
               else
                  valueTarget   =state%yTarget(i)
                  varianceTotal =state%covarianceTarget(i,i)
               end if
               ! The model prediction, floored. A prediction at the floor is only an upper limit, so has the variance of floored
               ! values in training.
               if (mean(i) > state%floor) then
                  valueModel    =mean(i)
                  derivativeModel=1.0d0
               else
                  valueModel    =state%floor
                  derivativeModel=0.0d0
                  varianceTotal =varianceTotal+state%rootVarianceFloored**2
               end if
               if (self%includeEmulatorVariance) varianceTotal=varianceTotal+variance(i)
               if (varianceTotal <= 0.0d0) cycle
               emulatedEvaluate     =+emulatedEvaluate                                    &
                    &                -0.5d0                                               &
                    &                *(                                                   &
                    &                  +(valueTarget-valueModel)**2/varianceTotal          &
                    &                  +log(2.0d0*Pi*varianceTotal)                       &
                    &                 )
               varianceLogLikelihood=+varianceLogLikelihood                               &
                    &                +((valueTarget-valueModel)/varianceTotal*derivativeModel)**2 &
                    &                *variance(i)
            end do
         else if (self%likelihoodForms(k) == emulatedLikelihoodFormPoissonCounts) then
            ! Find the fractional root variance of the expected counts - from the active parameter if one was given.
            if (self%indexRootVarianceTargetFractional > 0) then
               rootVarianceFractional=stateArray                       (self%indexRootVarianceTargetFractional)
            else
               rootVarianceFractional=self%rootVarianceTargetFractional(                                     k)
            end if
            do i=1,size(state%yTarget)
               if (.not.state%included(i)) cycle
               ! Find the mean of the untransformed observable, the fractional variance of the emulation, and the derivative of the
               ! expected count with respect to the emulated mean.
               if (state%transformLog) then
                  if (self%includeEmulatorVariance) then
                     varianceFractional=exp((ln10**2)*variance(i))-1.0d0
                  else
                     varianceFractional=0.0d0
                  end if
                  valueMean         =10.0d0**mean(i)*sqrt(1.0d0+varianceFractional)
                  countExpected     =valueMean*state%countPerValue(i)
                  derivativeCount   =countExpected*ln10
               else
                  valueMean         =mean(i)
                  countExpected     =valueMean*state%countPerValue(i)
                  derivativeCount   =          state%countPerValue(i)
                  if (self%includeEmulatorVariance .and. valueMean > 0.0d0) then
                     varianceFractional=variance(i)/valueMean**2
                  else
                     varianceFractional=0.0d0
                  end if
               end if
               countTarget=state%countTarget(i)
               if (countExpected <= 0.0d0) then
                  ! No counts are expected in this bin. If there are also no target counts, the bin can not be populated, so does not
                  ! affect the likelihood. Otherwise, the model is impossible.
                  if (countTarget <= 0.0d0) cycle
                  emulatedEvaluate=logImprobable
                  if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
                  return
               end if
               ! Combine the fractional variances of the expected count, and of the emulation. (Written so that the result is
               ! exactly the square of the fractional root variance of the expected count if the emulation has no variance.)
               dispersion=+rootVarianceFractional**2                    &
                    &     +varianceFractional                           &
                    &     +rootVarianceFractional**2*varianceFractional
               if (dispersion > 0.0d0) then
                  ! Negative binomial distribution, with stopping-time parameter r=1/dispersion.
                  emulatedEvaluate=+emulatedEvaluate                                                                             &
                       &           +countTarget                     *log                       (                countExpected  ) &
                       &           -                                 Logarithmic_Factorial     (nint(countTarget)              ) &
                       &           +                                 Gamma_Function_Logarithmic(1.0d0/dispersion+countTarget   ) &
                       &           -                                 Gamma_Function_Logarithmic(1.0d0/dispersion               ) &
                       &           -countTarget                     *log                       (1.0d0/dispersion+countExpected ) &
                       &           -(1.0d0/dispersion)              *log                       (countExpected*dispersion+1.0d0 )
                  gradientCount   =+countTarget/countExpected-(countTarget+1.0d0/dispersion)/(countExpected+1.0d0/dispersion)
               else
                  ! Poisson distribution.
                  emulatedEvaluate=+emulatedEvaluate                                        &
                       &           +countTarget*log                  (     countExpected  ) &
                       &           -            Logarithmic_Factorial(nint(countTarget   )) &
                       &           -                                       countExpected
                  gradientCount   =+countTarget/countExpected-1.0d0
               end if
               varianceLogLikelihood=+varianceLogLikelihood              &
                    &                +(gradientCount*derivativeCount)**2 &
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
                 emulatedEvaluate=+emulatedEvaluate                         &
                      &           -0.5d0                                    &
                      &           *(                                        &
                      &             +dot_product(residual,gradient)         &
                      &             +decomposition%logarithmicDeterminant() &
                      &             +dble(size(indices))*log(2.0d0*Pi)      &
                      &            )
               end block
               varianceLogLikelihood=+varianceLogLikelihood                 &
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

  elemental logical function emulatedIsNotANumber(x)
    !!{RST
    Return true if ``x`` is not-a-number, examining its bits directly (see :func:`emulatedIsFinite`): its exponent bits are all
    set, and its mantissa is non-zero.
    !!}
    use, intrinsic :: ISO_Fortran_Env, only : int64
    implicit none
    double precision, intent(in   ) :: x

    emulatedIsNotANumber=ibits(transfer(x,0_int64),52,11) == 2047_int64 .and. ibits(transfer(x,0_int64),0,52) /= 0_int64
    return
  end function emulatedIsNotANumber

  logical function emulatedWillEvaluate(self,simulationState,modelParameters_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
    !!{RST
    Return true if the log-likelihood will be evaluated: only states within the support of the prior are evaluated, as the
    prior quantiles (the inputs to the emulators) are undefined outside it.
    !!}
    use :: Models_Likelihoods_Constants  , only : logImpossible
    use :: Posterior_Sampling_Convergence, only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State      , only : posteriorSampleStateClass
    implicit none
    class           (posteriorSampleLikelihoodEmulated), intent(inout)               :: self
    class           (posteriorSampleStateClass        ), intent(inout)               :: simulationState
    type            (modelParameterList               ), intent(in   ), dimension(:) :: modelParameters_
    class           (posteriorSampleConvergenceClass  ), intent(inout)               :: simulationConvergence
    double precision                                   , intent(in   )               :: temperature          , logLikelihoodCurrent, &
         &                                                                              logPriorCurrent      , logPriorProposed
    !$GLC attributes unused :: self, simulationState, modelParameters_, simulationConvergence, temperature, logLikelihoodCurrent, logPriorCurrent

    emulatedWillEvaluate=logPriorProposed > logImpossible
    return
  end function emulatedWillEvaluate

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
