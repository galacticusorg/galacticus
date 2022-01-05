!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  Implementation of a posterior sampling likelihood class which implements a likelihood for halo mass functions.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass, haloEnvironmentClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOClass

  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodHaloMassFunction">
   <description>A posterior sampling likelihood class which implements a likelihood for halo mass functions.</description>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodHaloMassFunction
     !!{
     Implementation of a posterior sampling likelihood class which implements a likelihood for halo mass functions.
     !!}
     private
     double precision                               , dimension(:  ), allocatable :: mass                               , massFunction                                    , &
          &                                                                          massMinimum                        , massMaximum
     double precision                               , dimension(:,:), allocatable :: covarianceMatrix
     class           (cosmologyFunctionsClass      ), pointer                     :: cosmologyFunctions_       => null()
     class           (cosmologyParametersClass     ), pointer                     :: cosmologyParameters_      => null()
     class           (cosmologicalMassVarianceClass), pointer                     :: cosmologicalMassVariance_ => null(), cosmologicalMassVarianceUnconditioned_ => null()
     class           (criticalOverdensityClass     ), pointer                     :: criticalOverdensity_      => null(), criticalOverdensityUnconditioned_      => null()
     class           (darkMatterHaloScaleClass     ), pointer                     :: darkMatterHaloScale_      => null()
     class           (darkMatterProfileDMOClass    ), pointer                     :: darkMatterProfileDMO_     => null()
     class           (haloEnvironmentClass         ), pointer                     :: haloEnvironment_          => null()
     double precision                                                             :: time                               , massParticle                                    , &
          &                                                                          massRangeMinimum                   , redshift
     type            (vector                       )                              :: means
     type            (matrix                       )                              :: covariance
     integer                                                                      :: errorModel                         , haloMassFunctionType
     type            (varying_string               )                              :: fileName
     logical                                                                      :: environmentAveraged
   contains
     final     ::                    haloMassFunctionDestructor
     procedure :: evaluate        => haloMassFunctionEvaluate
     procedure :: functionChanged => haloMassFunctionFunctionChanged
  end type posteriorSampleLikelihoodHaloMassFunction

  interface posteriorSampleLikelihoodHaloMassFunction
     !!{
     Constructors for the {\normalfont \ttfamily haloMassFunction} posterior sampling convergence class.
     !!}
     module procedure haloMassFunctionConstructorParameters
     module procedure haloMassFunctionConstructorInternal
  end interface posteriorSampleLikelihoodHaloMassFunction

  ! Mass function enumeration.
  !![
  <enumeration>
   <name>haloMassFunctionType</name>
   <description>Used to specify the halo mass function for likelihoods.</description>
   <visibility>public</visibility>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <entry label="shethTormen"     />
   <entry label="bhattacharya2011"/>
  </enumeration>
  !!]

  ! Mass function error model enumeration.
  !![
  <enumeration>
   <name>haloMassFunctionErrorModel</name>
   <description>Used to specify the error model to use for halo mass function likelihoods.</description>
   <visibility>public</visibility>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <entry label="none"                />
   <entry label="powerLaw"            />
   <entry label="sphericalOverdensity"/>
   <entry label="trenti2010"          />
  </enumeration>
  !!]

contains

  function haloMassFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily haloMassFunction} posterior sampling convergence class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleLikelihoodHaloMassFunction)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    type            (inputParameters                          )                :: parametersUnconditioned
    type            (varying_string                           )                :: fileName                 , errorModel                            , &
         &                                                                        haloMassFunctionType
    double precision                                                           :: redshift                 , massRangeMinimum                      , &
         &                                                                        massParticle
    integer                                                                    :: binCountMinimum
    logical                                                                    :: environmentAveraged
    class           (cosmologyFunctionsClass                  ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass                 ), pointer       :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass            ), pointer       :: cosmologicalMassVariance_, cosmologicalMassVarianceUnconditioned_
    class           (criticalOverdensityClass                 ), pointer       :: criticalOverdensity_     , criticalOverdensityUnconditioned_
    class           (darkMatterHaloScaleClass                 ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass                ), pointer       :: darkMatterProfileDMO_
    class           (haloEnvironmentClass                     ), pointer       :: haloEnvironment_

    !![
    <inputParameter>
      <name>fileName</name>
      <description>The name of the file containing the halo mass function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshift</name>
      <description>The redshift at which to evaluate the halo mass function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massRangeMinimum</name>
      <description>The minimum halo mass to include in the likelihood evaluation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>binCountMinimum</name>
      <description>The minimum number of halos per bin required to permit bin to be included in likelihood evaluation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>haloMassFunctionType</name>
      <description>The type of halo mass function to use in likelihood calculations.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>errorModel</name>
      <description>The error model to use for the halo mass function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massParticle</name>
      <description>The N-body particle mass.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>environmentAveraged</name>
      <description>If true, the mass function will ve averaged over all environments.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"      name="darkMatterHaloScale_"      source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"     name="darkMatterProfileDMO_"     source="parameters"/>
    <objectBuilder class="haloEnvironment"          name="haloEnvironment_"          source="parameters"/>
    !!]
    if (environmentAveraged) then
       parametersUnconditioned=parameters%subParameters("unconditioned",requireValue=.false.)
       !![
       <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVarianceUnconditioned_" source="parametersUnconditioned"/>
       <objectBuilder class="criticalOverdensity"      name="criticalOverdensityUnconditioned_"      source="parametersUnconditioned"/>
       !!]
    else
       cosmologicalMassVarianceUnconditioned_ => null()
       criticalOverdensityUnconditioned_      => null()
    end if
    self=posteriorSampleLikelihoodHaloMassFunction(char(fileName),redshift,massRangeMinimum,binCountMinimum,enumerationHaloMassFunctionTypeEncode(char(haloMassFunctionType),includesPrefix=.false.),enumerationHaloMassFunctionErrorModelEncode(char(errorModel),includesPrefix=.false.),massParticle,environmentAveraged,cosmologyFunctions_,cosmologyParameters_,cosmologicalMassVariance_,criticalOverdensity_,cosmologicalMassVarianceUnconditioned_,criticalOverdensityUnconditioned_,darkMatterHaloScale_,darkMatterProfileDMO_,haloEnvironment_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="darkMatterHaloScale_"     />
    <objectDestructor name="darkMatterProfileDMO_"    />
    <objectDestructor name="haloEnvironment_"         />
    !!]
    if (environmentAveraged) then
       !![
       <objectDestructor name="cosmologicalMassVarianceUnconditioned_"/>
       <objectDestructor name="criticalOverdensityUnconditioned_"     />
       !!]
    end if
    return
  end function haloMassFunctionConstructorParameters

  function haloMassFunctionConstructorInternal(fileName,redshift,massRangeMinimum,binCountMinimum,haloMassFunctionType,errorModel,massParticle,environmentAveraged,cosmologyFunctions_,cosmologyParameters_,cosmologicalMassVariance_,criticalOverdensity_,cosmologicalMassVarianceUnconditioned_,criticalOverdensityUnconditioned_,darkMatterHaloScale_,darkMatterProfileDMO_,haloEnvironment_) result(self)
    !!{
    Constructor for ``haloMassFunction'' posterior sampling likelihood class.
    !!}
    use :: Display          , only : displayMessage         , displayMagenta, displayReset
    use :: File_Utilities   , only : File_Name_Expand
    use :: Galacticus_Error , only : Galacticus_Error_Report
    use :: HDF5_Access      , only : hdf5Access
    use :: IO_HDF5          , only : hdf5Object
    use :: Linear_Algebra   , only : assignment(=)
    use :: Memory_Management, only : allocateArray
    implicit none
    type            (posteriorSampleLikelihoodHaloMassFunction)                                :: self
    character       (len=*                                    ), intent(in   )                 :: fileName
    double precision                                           , intent(in   )                 :: redshift                      , massRangeMinimum                      , &
         &                                                                                        massParticle
    integer                                                    , intent(in   )                 :: binCountMinimum               , errorModel                            , &
         &                                                                                        haloMassFunctionType
    logical                                                    , intent(in   )                 :: environmentAveraged
    class           (cosmologyFunctionsClass                  ), intent(in   ), target         :: cosmologyFunctions_
    class           (cosmologyParametersClass                 ), intent(in   ), target         :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass            ), intent(in   ), target         :: cosmologicalMassVariance_     , cosmologicalMassVarianceUnconditioned_
    class           (criticalOverdensityClass                 ), intent(in   ), target         :: criticalOverdensity_          , criticalOverdensityUnconditioned_
    class           (darkMatterHaloScaleClass                 ), intent(in   ), target         :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass                ), intent(in   ), target         :: darkMatterProfileDMO_
    class           (haloEnvironmentClass                     ), intent(in   ), target         :: haloEnvironment_
    double precision                                           , allocatable  , dimension(:  ) :: eigenValueArray               , massOriginal                          , &
         &                                                                                        massFunctionOriginal
    integer         (c_size_t                                 ), allocatable  , dimension(:  ) :: massFunctionCountOriginal
    double precision                                           , allocatable  , dimension(:,:) :: massFunctionCovarianceOriginal
    character       (len=12                                   )                                :: redshiftLabel
    type            (hdf5Object                               )                                :: massFunctionFile              , simulationGroup
    integer                                                                                    :: i                             , j                                     , &
         &                                                                                        ii                            , jj                                    , &
         &                                                                                        massCountReduced
    double precision                                                                           :: massIntervalLogarithmic
    type            (matrix                                   )                                :: eigenVectors
    type            (vector                                   )                                :: eigenValues
    !![
    <constructorAssign variables="fileName, redshift, massRangeMinimum, haloMassFunctionType, errorModel, massParticle, environmentAveraged, *cosmologyFunctions_, *cosmologyParameters_, *cosmologicalMassVariance_, *criticalOverdensity_, *cosmologicalMassVarianceUnconditioned_, *criticalOverdensityUnconditioned_, *darkMatterHaloScale_, *darkMatterProfileDMO_, *haloEnvironment_"/>
    !!]

    ! Convert redshift to time.
    self%time=self%cosmologyFunctions_ %cosmicTime                (          &
         &    self%cosmologyFunctions_%expansionFactorFromRedshift (         &
         &                                                          redshift &
         &                                                         )         &
         &                                                        )
    ! Read the halo mass function file.
    write (redshiftLabel,'(f6.3)') redshift
    !$ call hdf5Access%set()
    call massFunctionFile%openFile(char(File_Name_Expand(trim(fileName))),readOnly=.true.)
    simulationGroup=massFunctionFile%openGroup('simulation0001')
    call simulationGroup %readDataset("mass"        ,massOriginal             )
    call simulationGroup %readDataset("massFunction",massFunctionOriginal     )
    call simulationGroup %readDataset("count"       ,massFunctionCountOriginal)
    call simulationGroup %close      (                                        )
    call massFunctionFile%close      (                                        )
    !$ call hdf5Access%unset()
    ! Construct the covariance matrix.
    allocate(massFunctionCovarianceOriginal(size(massOriginal),size(massOriginal)))
    massFunctionCovarianceOriginal=0.0d0
    do i=1,size(massOriginal)
       if (massFunctionCountOriginal(i) > 0_c_size_t) &
            &  massFunctionCovarianceOriginal(i,i)=+     massFunctionOriginal     (i) **2 &
            &                                      /dble(massFunctionCountOriginal(i))
    end do
    ! Find a reduced mass function excluding any empty bins.
    massCountReduced=0
    do i=1,size(massOriginal)
       if (massFunctionOriginal     (i) <= 0.0d0           ) cycle
       if (massOriginal             (i) <= massRangeMinimum) cycle
       if (massFunctionCountOriginal(i) <= binCountMinimum ) cycle
       massCountReduced=massCountReduced+1
    end do
    if (massCountReduced == 0) call Galacticus_Error_Report('no usable bins in mass function'//{introspection:location})
    call allocateArray(self%mass            ,[massCountReduced                 ])
    call allocateArray(self%massFunction    ,[massCountReduced                 ])
    call allocateArray(self%covarianceMatrix,[massCountReduced,massCountReduced])
    ii=0
    do i=1,size(massOriginal)
       if (     massFunctionOriginal          (i  )  <= 0.0d0                                              ) cycle
       if (     massOriginal                  (i  )  <= massRangeMinimum                                   ) cycle
       if (sqrt(massFunctionCovarianceOriginal(i,i)) >= massFunctionOriginal(i)/sqrt(dble(binCountMinimum))) cycle
       ii=ii+1
       self%mass        (ii)=massOriginal        (i)
       self%massFunction(ii)=massFunctionOriginal(i)
       jj=0
       do j=1,size(massOriginal)
          if (     massFunctionOriginal          (j  )  <= 0.0d0                                              ) cycle
          if (     massOriginal                  (j  )  <= massRangeMinimum                                   ) cycle
          if (sqrt(massFunctionCovarianceOriginal(j,j)) >= massFunctionOriginal(j)/sqrt(dble(binCountMinimum))) cycle
          jj=jj+1
          self%covarianceMatrix(ii,jj)=massFunctionCovarianceOriginal(i,j)
       end do
    end do
    ! Compute mass ranges for bins.
    massIntervalLogarithmic=+log(                                  &
         &                       +massOriginal(size(massOriginal)) &
         &                       /massOriginal(                 1) &
         &                      )                                  &
         &                  /dble(                                 &
         &                        +size(massOriginal)              &
         &                        -1                               &
         &                       )
    call allocateArray(self%massMinimum,shape(self%mass))
    call allocateArray(self%massMaximum,shape(self%mass))
    do i=1,size(self%mass)
       self%massMinimum(i)=self%mass(i)*exp(-0.5d0*massIntervalLogarithmic)
       self%massMaximum(i)=self%mass(i)*exp(+0.5d0*massIntervalLogarithmic)
    end do
    ! Find the covariance matrix.
    self%covariance=self%covarianceMatrix
    ! Get eigenvalues and vectors of the covariance matrix.
    allocate(eigenValueArray(size(self%mass)))
    call self%covariance%eigenSystem(eigenVectors,eigenValues)
    eigenValueArray=eigenValues
    if (any(eigenValueArray < 0.0d0)) call displayMessage(displayMagenta()//'WARNING:'//displayReset()//' inverse covariance matrix is not semi-positive definite')
    deallocate(eigenValueArray               )
    deallocate(massOriginal                  )
    deallocate(massFunctionOriginal          )
    deallocate(massFunctionCovarianceOriginal)
    return
  end function haloMassFunctionConstructorInternal

  subroutine haloMassFunctionDestructor(self)
    !!{
    Destructor for ``haloMassFunction'' posterior sampling likelihood class.
    !!}
    implicit none
    type(posteriorSampleLikelihoodHaloMassFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%darkMatterHaloScale_"     />
    <objectDestructor name="self%darkMatterProfileDMO_"    />
    <objectDestructor name="self%haloEnvironment_"         />
    !!]
    if (self%environmentAveraged) then
       !![
       <objectDestructor name="self%cosmologicalMassVarianceUnconditioned_"/>
       <objectDestructor name="self%criticalOverdensityUnconditioned_"     />
       !!]
    end if
    return
  end subroutine haloMassFunctionDestructor

  double precision function haloMassFunctionEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !!{
    Return the log-likelihood for the halo mass function likelihood function.
    !!}
    use :: Galacticus_Error                 , only : Galacticus_Error_Report
    use :: Halo_Mass_Functions              , only : haloMassFunctionClass           , haloMassFunctionEnvironmentAveraged, haloMassFunctionErrorConvolved, haloMassFunctionShethTormen, &
         &                                           haloMassFunctionBhattacharya2011
    use :: Linear_Algebra                   , only : assignment(=)                   , operator(*)
    use :: Models_Likelihoods_Constants     , only : logImpossible
    use :: Posterior_Sampling_Convergence   , only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State         , only : posteriorSampleStateClass
    use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass         , nbodyHaloMassErrorPowerLaw         , nbodyHaloMassErrorSOHaloFinder, nbodyHaloMassErrorTrenti2010
    implicit none
    class           (posteriorSampleLikelihoodHaloMassFunction), intent(inout)               :: self
    class           (posteriorSampleStateClass                ), intent(inout)               :: simulationState
    type            (modelParameterList                       ), intent(in   ), dimension(:) :: modelParametersActive_          , modelParametersInactive_
    class           (posteriorSampleConvergenceClass          ), intent(inout)               :: simulationConvergence
    double precision                                           , intent(in   )               :: temperature                     , logLikelihoodCurrent          , &
         &                                                                                      logPriorCurrent                 , logPriorProposed
    real                                                       , intent(inout)               :: timeEvaluate
    double precision                                           , intent(  out), optional     :: logLikelihoodVariance
    logical                                                    , intent(inout), optional     :: forceAcceptance
    double precision                                           , allocatable  , dimension(:) :: stateVector                     , massFunction
    double precision                                           , parameter                   :: errorFractionalMaximum    =1.0d1
    class           (nbodyHaloMassErrorClass                  ), pointer                     :: nbodyHaloMassError_
    class           (haloMassFunctionClass                    ), pointer                     :: haloMassFunctionRaw_            , haloMassFunctionAveraged_     , &
         &                                                                                      haloMassFunctionConvolved_      , haloMassFunctionUnconditioned_
    type            (vector                                   )                              :: difference
    integer                                                                                  :: i                               , haloMassFunctionParameterCount
    !$GLC attributes unused :: simulationConvergence, temperature, timeEvaluate, logLikelihoodCurrent, logPriorCurrent, modelParametersInactive_, forceAcceptance

    ! There is no variance in our likelihood estimate.
    if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
    ! Do not evaluate if the proposed prior is impossible.
    if (logPriorProposed <= logImpossible) then
       haloMassFunctionEvaluate=0.0d0
       return
    end if
    ! Build the halo mass function object.
    stateVector=simulationState%get()
    do i=1,size(stateVector)
       stateVector(i)=modelParametersActive_(i)%modelParameter_%unmap(stateVector(i))
    end do
    ! Construct the raw mass function.
    select case (self%haloMassFunctionType)
    case (haloMassFunctionTypeShethTormen     )
       haloMassFunctionParameterCount=3
       if (size(stateVector) < haloMassFunctionParameterCount  )                                               &
            & call Galacticus_Error_Report(                                                                    &
            &                              'at least 3 parameters are required for this likelihood function'// &
            &                              {introspection:location}                                            &
            &                             )
       allocate(haloMassFunctionShethTormen      :: haloMassFunctionRaw_)
       select type (haloMassFunctionRaw_)
       type is (haloMassFunctionShethTormen     )
          haloMassFunctionRaw_=haloMassFunctionShethTormen     (                                   &
               &                                                self%cosmologyParameters_        , &
               &                                                self%cosmologicalMassVariance_   , &
               &                                                self%criticalOverdensity_        , &
               &                                                stateVector                   (1), &
               &                                                stateVector                   (2), &
               &                                                stateVector                   (3)  &
               &                                               )
       end select
    case (haloMassFunctionTypeBhattacharya2011)
       haloMassFunctionParameterCount=4
       if (size(stateVector) < haloMassFunctionParameterCount)                                                 &
            & call Galacticus_Error_Report(                                                                    &
            &                              'at least 4 parameters are required for this likelihood function'// &
            &                              {introspection:location}                                            &
            &                             )
       allocate(haloMassFunctionBhattacharya2011 :: haloMassFunctionRaw_)
       select type (haloMassFunctionRaw_)
       type is (haloMassFunctionBhattacharya2011)
          haloMassFunctionRaw_=haloMassFunctionBhattacharya2011(                                   &
               &                                                self%cosmologyParameters_        , &
               &                                                self%cosmologicalMassVariance_   , &
               &                                                self%criticalOverdensity_        , &
               &                                                stateVector                   (1), &
               &                                                stateVector                   (2), &
               &                                                stateVector                   (3), &
               &                                                stateVector                   (4)  &
               &                                               )
       end select
    case default
       haloMassFunctionParameterCount=-1
       call Galacticus_Error_Report('unknown halo mass function'//{introspection:location})
    end select
    ! If averaging over environment, build the averager.
    if (self%environmentAveraged) then
       ! First build an unconditioned (on environment) mass function.
       select case (self%haloMassFunctionType)
       case (haloMassFunctionTypeShethTormen     )
          allocate(haloMassFunctionShethTormen      :: haloMassFunctionUnconditioned_)
          select type (haloMassFunctionUnconditioned_)
          type is (haloMassFunctionShethTormen     )
             haloMassFunctionUnconditioned_=haloMassFunctionShethTormen     (                                                &
                  &                                                     self%cosmologyParameters_                     , &
                  &                                                          self%cosmologicalMassVarianceUnconditioned_   , &
                  &                                                          self%criticalOverdensityUnconditioned_        , &
                  &                                                          stateVector                                (1), &
                  &                                                          stateVector                                (2), &
                  &                                                          stateVector                                (3)  &
                  &                                                         )
          end select
       case (haloMassFunctionTypeBhattacharya2011     )
          allocate(haloMassFunctionBhattacharya2011 :: haloMassFunctionUnconditioned_)
          select type (haloMassFunctionUnconditioned_)
          type is (haloMassFunctionBhattacharya2011)
             haloMassFunctionUnconditioned_=haloMassFunctionBhattacharya2011(                                                &
                  &                                                          self%cosmologyParameters_                     , &
                  &                                                          self%cosmologicalMassVarianceUnconditioned_   , &
                  &                                                          self%criticalOverdensityUnconditioned_        , &
                  &                                                          stateVector                                (1), &
                  &                                                          stateVector                                (2), &
                  &                                                          stateVector                                (3), &
                  &                                                          stateVector                                (4)  &
                  &                                                         )
          end select
       case default
          call Galacticus_Error_Report('unknown halo mass function'//{introspection:location})
       end select
       ! Now build the environment averaged mass function.
       allocate(haloMassFunctionEnvironmentAveraged :: haloMassFunctionAveraged_)
       select type (haloMassFunctionAveraged_)
       type is (haloMassFunctionEnvironmentAveraged)
          haloMassFunctionAveraged_=haloMassFunctionEnvironmentAveraged(                                     &
               &                                                             haloMassFunctionRaw_          , &
               &                                                             haloMassFunctionUnconditioned_, &
               &                                                        self%haloEnvironment_              , &
               &                                                        self%cosmologyParameters_            &
               &                                                       )
       end select
    else
       haloMassFunctionAveraged_ => haloMassFunctionRaw_
    end if
    ! If convolving with an error distribution, build the error model.
    select case (self%errorModel)
    case (haloMassFunctionErrorModelNone                )
       if (size(stateVector) /= haloMassFunctionParameterCount  )                                          &
            & call Galacticus_Error_Report(                                                                &
            &                              'incorrect number of parameters for this likelihood function'// &
            &                              {introspection:location}                                        &
            &                             )
       haloMassFunctionConvolved_ => haloMassFunctionAveraged_
    case (haloMassFunctionErrorModelPowerLaw            )
       if (size(stateVector) /= haloMassFunctionParameterCount+3)                                          &
            & call Galacticus_Error_Report(                                                                &
            &                              'incorrect number of parameters for this likelihood function'// &
            &                              {introspection:location}                                        &
            &                             )
       allocate(nbodyHaloMassErrorPowerLaw :: nbodyHaloMassError_)
       select type (nbodyHaloMassError_)
       type is (nbodyHaloMassErrorPowerLaw)
          nbodyHaloMassError_=nbodyHaloMassErrorPowerLaw(                &
               &                                         stateVector(4), &
               &                                         stateVector(5), &
               &                                         stateVector(6)  &
               &                                        )
       end select
       allocate(haloMassFunctionErrorConvolved :: haloMassFunctionConvolved_)
       select type (haloMassFunctionConvolved_)
       type is (haloMassFunctionErrorConvolved)
          haloMassFunctionConvolved_   =haloMassFunctionErrorConvolved(                              &
               &                                                       haloMassFunctionAveraged_   , &
               &                                                       self%cosmologyParameters_   , &
               &                                                       nbodyHaloMassError_         , &
               &                                                       errorFractionalMaximum        &
               &                                                      )
       end select
       nullify(nbodyHaloMassError_)
    case (haloMassFunctionErrorModelSphericalOverdensity)
       if (size(stateVector) /= haloMassFunctionParameterCount  )                                          &
            & call Galacticus_Error_Report(                                                                &
            &                              'incorrect number of parameters for this likelihood function'// &
            &                              {introspection:location}                                        &
            &                             )
       ! Use a mass function convolved with an error model for spherical overdensity algorithm errors.
       allocate(nbodyHaloMassErrorSOHaloFinder :: nbodyHaloMassError_)
       select type (nbodyHaloMassError_)
       type is (nbodyHaloMassErrorSOHaloFinder)
          nbodyHaloMassError_=nbodyHaloMassErrorSOHaloFinder(                            &
               &                                             self%darkMatterHaloScale_ , &
               &                                             self%darkMatterProfileDMO_, &
               &                                             self%massParticle           &
               &                                            )
       end select
       allocate(haloMassFunctionErrorConvolved :: haloMassFunctionConvolved_)
       select type (haloMassFunctionConvolved_)
       type is (haloMassFunctionErrorConvolved)
          haloMassFunctionConvolved_=haloMassFunctionErrorConvolved(                              &
               &                                                    haloMassFunctionAveraged_   , &
               &                                                    self%cosmologyParameters_   , &
               &                                                    nbodyHaloMassError_         , &
               &                                                    errorFractionalMaximum        &
               &                                                   )
       end select
       nullify(nbodyHaloMassError_)
    case (haloMassFunctionErrorModelTrenti2010          )
       if (size(stateVector) /= haloMassFunctionParameterCount  )                                          &
            & call Galacticus_Error_Report(                                                                &
            &                              'incorrect number of parameters for this likelihood function'// &
            &                              {introspection:location}                                        &
            &                             )
       ! Use a mass function convolved with the error model from Trenti et al. (2010).
       allocate(nbodyHaloMassErrorTrenti2010 :: nbodyHaloMassError_)
       select type (nbodyHaloMassError_)
       type is (nbodyHaloMassErrorTrenti2010)
          nbodyHaloMassError_=nbodyHaloMassErrorTrenti2010(                  &
               &                                           self%massParticle &
               &                                          )
       end select
       allocate(haloMassFunctionErrorConvolved :: haloMassFunctionConvolved_)
       select type (haloMassFunctionConvolved_)
       type is (haloMassFunctionErrorConvolved)
          haloMassFunctionConvolved_=haloMassFunctionErrorConvolved(                              &
               &                                                    haloMassFunctionAveraged_   , &
               &                                                    self%cosmologyParameters_   , &
               &                                                    nbodyHaloMassError_         , &
               &                                                    errorFractionalMaximum        &
               &                                                   )
       end select
       nullify(nbodyHaloMassError_)
    end select
    ! Compute the mass function.
    allocate(massFunction(size(self%mass)))
    do i=1,size(self%mass)
       massFunction(i)=+haloMassFunctionConvolved_%integrated(                     &
            &                                                 self%time          , &
            &                                                 self%massMinimum(i), &
            &                                                 self%massMaximum(i)  &
            &                                                )                     &
            &          /log(                                                       &
            &                                                +self%massMaximum(i)  &
            &                                                /self%massMinimum(i)  &
            &              )
    end do
    ! Evaluate the log-likelihood.
    difference              =massFunction-self%massFunction
    haloMassFunctionEvaluate=-0.5d0*self%covariance%covarianceProduct(difference)
    ! Clean up.
    deallocate(stateVector )
    deallocate(massFunction)
    return
  end function haloMassFunctionEvaluate

  subroutine haloMassFunctionFunctionChanged(self)
    !!{
    Respond to possible changes in the likelihood function.
    !!}
    implicit none
    class(posteriorSampleLikelihoodHaloMassFunction), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine haloMassFunctionFunctionChanged
