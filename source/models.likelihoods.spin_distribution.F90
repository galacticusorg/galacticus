!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implementation of a posterior sampling likelihood class which implements a likelihood for halo spin distributions.
  
  use Cosmology_Functions
  use Halo_Mass_Functions
  use Statistics_NBody_Halo_Mass_Errors
  use Dark_Matter_Profiles
  use Dark_Matter_Halo_Scales
  use Dark_Matter_Profile_Scales       , only : darkMatterProfileScaleRadius, darkMatterProfileScaleRadiusClass

  !# <posteriorSampleLikelihood name="posteriorSampleLikelihoodSpinDistribution">
  !#  <description>A posterior sampling likelihood class which implements a likelihood for halo spin distributions.</description>
  !# </posteriorSampleLikelihood>
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodSpinDistribution
     !% Implementation of a posterior sampling likelihood class which implements a likelihood for SED fitting.
     private
     class           (cosmologyFunctionsClass          ), pointer                     :: cosmologyFunctions_
     class           (haloMassFunctionClass            ), pointer                     :: haloMassFunction_
     class           (nbodyHaloMassErrorClass          ), pointer                     :: nbodyHaloMassError_
     class           (darkMatterProfileClass           ), pointer                     :: darkMatterProfile_
     class           (darkMatterHaloScaleClass         ), pointer                     :: darkMatterHaloScale_
     class           (darkMatterProfileScaleRadiusClass), pointer                     :: darkMatterProfileScaleRadius_
     double precision                                   , dimension(:  ), allocatable :: spin                         , distribution                      , &
          &                                                                              spinMinimum                  , spinMaximum                       , &
          &                                                                              distributionError
     double precision                                                                 :: time                         , massParticle                      , &
          &                                                                              massHaloMinimum              , energyEstimateParticleCountMaximum, &
          &                                                                              redshift
     integer                                                                          :: particleCountMinimum         , distributionType
     type            (varying_string                   )                              :: fileName
   contains
     final     ::                    spinDistributionDestructor
     procedure :: evaluate        => spinDistributionEvaluate
     procedure :: functionChanged => spinDistributionFunctionChanged
  end type posteriorSampleLikelihoodSpinDistribution

  interface posteriorSampleLikelihoodSpinDistribution
     !% Constructors for the {\normalfont \ttfamily spinDistribution} posterior sampling convergence class.
     module procedure spinDistributionConstructorParameters
     module procedure spinDistributionConstructorInternal
  end interface posteriorSampleLikelihoodSpinDistribution

  !# <enumeration>
  !#  <name>spinDistributionType</name>
  !#  <description>Used to specify the type of intrinsic spin distribution.</description>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <entry label="logNormal" />
  !#  <entry label="bett2007"  />
  !# </enumeration>

contains

  function spinDistributionConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily spinDistribution} posterior sampling convergence class which builds the object from a
    !% parameter set.
    use Input_Parameters
    implicit none
    type            (posteriorSampleLikelihoodSpinDistribution)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                  ), pointer       :: cosmologyFunctions_
    class           (haloMassFunctionClass                    ), pointer       :: haloMassFunction_
    class           (nbodyHaloMassErrorClass                  ), pointer       :: nbodyHaloMassError_
    class           (darkMatterProfileClass                   ), pointer       :: darkMatterProfile_
    class           (darkMatterHaloScaleClass                 ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass        ), pointer       :: darkMatterProfileScaleRadius_
    type            (varying_string                           )                :: fileName                     , distributionType
    double precision                                                           :: redshift                     , massHaloMinimum                   , &
         &                                                                        massParticle                 , energyEstimateParticleCountMaximum
    integer                                                                    :: particleCountMinimum

    !# <inputParameter>
    !#   <name>fileName</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The name of the file containing the target spin distribution.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>distributionType</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The name of the spin distribution to use.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>redshift</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The redshift at which to compute the spin distribution.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massParticle</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The mass of a particle in the N-body simulation.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massHaloMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The minimum halo mass over which to integrate.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>particleCountMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The minimum particle count used in N-body halos.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>energyEstimateParticleCountMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The maximum number of N-body particles used in estimating halo energies.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    !# <objectBuilder class="haloMassFunction"             name="haloMassFunction_"             source="parameters"/>
    !# <objectBuilder class="nbodyHaloMassError"           name="nbodyHaloMassError_"           source="parameters"/>
    !# <objectBuilder class="darkMatterProfile"            name="darkMatterProfile_"            source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    !# <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    self=posteriorSampleLikelihoodSpinDistribution(char(fileName),enumerationSpinDistributionTypeEncode(char(distributionType),includesPrefix=.false.),redshift,massHaloMinimum,massParticle,particleCountMinimum,energyEstimateParticleCountMaximum,cosmologyFunctions_,haloMassFunction_,nbodyHaloMassError_,darkMatterProfile_,darkMatterHaloScale_,darkMatterProfileScaleRadius_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"          />
    !# <objectDestructor name="haloMassFunction_"            />
    !# <objectDestructor name="nbodyHaloMassError_"          />
    !# <objectDestructor name="darkMatterProfile_"           />
    !# <objectDestructor name="darkMatterHaloScale_"         />
    !# <objectDestructor name="darkMatterProfileScaleRadius_"/>
    return
  end function spinDistributionConstructorParameters

  function spinDistributionConstructorInternal(fileName,distributionType,redshift,massHaloMinimum,massParticle,particleCountMinimum,energyEstimateParticleCountMaximum,cosmologyFunctions_,haloMassFunction_,nbodyHaloMassError_,darkMatterProfile_,darkMatterHaloScale_,darkMatterProfileScaleRadius_) result(self)
    !% Constructor for ``spinDistribution'' posterior sampling likelihood class.
    use IO_HDF5
    use Memory_Management
    use Cosmology_Functions
    use Galacticus_Error
    implicit none
    type            (posteriorSampleLikelihoodSpinDistribution)                        :: self
    character       (len=*                                    ), intent(in   )         :: fileName
    double precision                                           , intent(in   )         :: redshift                     , massHaloMinimum                   , &
         &                                                                                massParticle                 , energyEstimateParticleCountMaximum
    integer                                                    , intent(in   )         :: particleCountMinimum         , distributionType
    class           (cosmologyFunctionsClass                  ), intent(in   ), target :: cosmologyFunctions_
    class           (haloMassFunctionClass                    ), intent(in   ), target :: haloMassFunction_
    class           (nbodyHaloMassErrorClass                  ), intent(in   ), target :: nbodyHaloMassError_
    class           (darkMatterProfileClass                   ), intent(in   ), target :: darkMatterProfile_
    class           (darkMatterHaloScaleClass                 ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass        ), intent(in   ), target :: darkMatterProfileScaleRadius_
    type            (hdf5Object                               )                        :: spinDistributionFile
    double precision                                                                   :: spinIntervalLogarithmic
    integer                                                                            :: i
    !# <constructorAssign variables="fileName, distributionType, redshift, massHaloMinimum, massParticle, particleCountMinimum, energyEstimateParticleCountMaximum, *cosmologyFunctions_, *haloMassFunction_, *nbodyHaloMassError_, *darkMatterProfile_, *darkMatterHaloScale_, *darkMatterProfileScaleRadius_"/>

    ! Convert redshift to time.
    self%time=self%cosmologyFunctions_ %cosmicTime                 (          &
         &     self%cosmologyFunctions_%expansionFactorFromRedshift (         &
         &                                                           redshift &
         &                                                          )         &
         &                                                         )
    ! Read the target spin distribution from file.
    !$ call hdf5Access%set()
    call spinDistributionFile%openFile   (trim(fileName),readOnly=.true.)
    call spinDistributionFile%readDataset("spinParameter"    ,self%spin             )
    call spinDistributionFile%readDataset("distribution"     ,self%distribution     )
    call spinDistributionFile%readDataset("distributionError",self%distributionError)
    call spinDistributionFile%close()
    !$ call hdf5Access%unset()
    ! Compute spin ranges for bins.
    spinIntervalLogarithmic=+log(                            &
         &                       +self%spin(size(self%spin)) &
         &                       /self%spin(              1) &
         &                      )                            &
         &                  /dble(                           &
         &                        +size(self%spin)           &
         &                        -1                         &
         &                       )
    call allocateArray(self%spinMinimum,shape(self%spin))
    call allocateArray(self%spinMaximum,shape(self%spin))
    do i=1,size(self%spin)
       self%spinMinimum(i)=self%spin(i)*exp(-0.5d0*spinIntervalLogarithmic)
       self%spinMaximum(i)=self%spin(i)*exp(+0.5d0*spinIntervalLogarithmic)
    end do
    return
  end function spinDistributionConstructorInternal

  subroutine spinDistributionDestructor(self)
    !% Destructor for ``spinDistribution'' posterior sampling likelihood class.
    implicit none
    type(posteriorSampleLikelihoodSpinDistribution), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"          />
    !# <objectDestructor name="self%haloMassFunction_"            />
    !# <objectDestructor name="self%nbodyHaloMassError_"          />
    !# <objectDestructor name="self%darkMatterProfile_"           />
    !# <objectDestructor name="self%darkMatterHaloScale_"         />
    !# <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    return
  end subroutine spinDistributionDestructor

  double precision function spinDistributionEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !% Return the log-likelihood for the halo spin distribution likelihood function.
    use, intrinsic :: ISO_C_Binding
    use               Numerical_Integration
    use               MPI_Utilities
    use               Posterior_Sampling_State
    use               Models_Likelihoods_Constants
    use               Posterior_Sampling_Convergence
    use               Galacticus_Error
    use               Halo_Spin_Distributions
    use               Galacticus_Nodes              , only : treeNode     , nodeComponentBasic        , nodeComponentSpin
    use               FGSL                          , only : fgsl_function, fgsl_integration_workspace
    implicit none
    class           (posteriorSampleLikelihoodSpinDistribution), intent(inout)               :: self
    class           (posteriorSampleStateClass                ), intent(inout)               :: simulationState
    type            (modelParameterList                       ), intent(in   ), dimension(:) :: modelParametersActive_, modelParametersInactive_
    class           (posteriorSampleConvergenceClass          ), intent(inout)               :: simulationConvergence
    double precision                                           , intent(in   )               :: temperature           , logLikelihoodCurrent    , &
         &                                                                                      logPriorCurrent       , logPriorProposed
    real                                                       , intent(inout)               :: timeEvaluate
    double precision                                           , intent(  out), optional     :: logLikelihoodVariance
    logical                                                    , intent(inout), optional     :: forceAcceptance
    double precision                                           , allocatable  , dimension(:) :: stateVector           , distribution
    type            (treeNode                                 ), pointer                     :: node
    class           (nodeComponentBasic                       ), pointer                     :: nodeBasic
    class           (nodeComponentSpin                        ), pointer                     :: nodeSpin
    class           (haloSpinDistributionLogNormal            ), pointer                     :: distributionLogNormal
    class           (haloSpinDistributionBett2007             ), pointer                     :: distributionBett2007
    type            (haloSpinDistributionNbodyErrors          )                              :: distributionNbody
    type            (fgsl_function                            )                              :: integrandFunction
    type            (fgsl_integration_workspace               )                              :: integrationWorkspace
    integer                                                                                  :: i
    !GCC$ attributes unused :: simulationConvergence, temperature, timeEvaluate, logLikelihoodCurrent, logPriorCurrent, modelParametersInactive_, forceAcceptance

    ! There is no variance in our likelihood estimate.
    if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
    ! Do not evaluate if the proposed prior is impossible.
    if (logPriorProposed <= logImpossible) then
       spinDistributionEvaluate=0.0d0
       return
    end if
    ! Extract parameters of the spin distribution.
    stateVector=simulationState%get()
    do i=1,size(stateVector)
       stateVector(i)=modelParametersActive_(i)%modelParameter_%unmap(stateVector(i))
    end do
    ! Build the halo spin distribution object.
    select case (self%distributionType)
    case (spinDistributionTypeLogNormal)
       if (size(stateVector) /= 2)                                                                             &
            & call Galacticus_Error_Report(                                                                    & 
            &                              '2 parameters are required for the "lognormal" spin distribution'// &
            &                              {introspection:location}                                            &
            &                             )
       allocate(haloSpinDistributionLogNormal :: distributionLogNormal)
       select type (distributionLogNormal)
       type is (haloSpinDistributionLogNormal)
          distributionLogNormal=haloSpinDistributionLogNormal  (                                             &
               &                                                stateVector(1)                             , &
               &                                                stateVector(2)                               &
               &                                               )
          distributionNbody    =haloSpinDistributionNbodyErrors(                                             &
               &                                                distributionLogNormal                      , &
               &                                                self%massParticle                          , &
               &                                                self%particleCountMinimum                  , &
               &                                                self%energyEstimateParticleCountMaximum    , &
               &                                                self%time                                  , &
               &                                                self%nbodyHaloMassError_                   , &
               &                                                self%haloMassFunction_                     , &
               &                                                self%darkMatterHaloScale_                  , &
               &                                                self%darkMatterProfile_                    , &
               &                                                self%darkMatterProfileScaleRadius_           &
               &                                               )
       end select
    case (spinDistributionTypeBett2007)
       if (size(stateVector) /= 2)                                                                            &
            & call Galacticus_Error_Report(                                                                   & 
            &                              '2 parameters are required for the "bett2007" spin distribution'// &
            &                              {introspection:location}                                           &
            &                             )
       allocate(haloSpinDistributionBett2007 :: distributionBett2007)
       select type (distributionBett2007)
       type is (haloSpinDistributionBett2007)
          distributionBett2007=haloSpinDistributionBett2007    (                                             &
               &                                                stateVector(1)                             , &
               &                                                stateVector(2)                               &
               &                                               )
          distributionNbody    =haloSpinDistributionNbodyErrors(                                             &
               &                                                distributionBett2007                       , &
               &                                                self%massParticle                          , &
               &                                                self%particleCountMinimum                  , &
               &                                                self%energyEstimateParticleCountMaximum    , &
               &                                                self%time                                  , &
               &                                                self%nbodyHaloMassError_                   , &
               &                                                self%haloMassFunction_                     , &
               &                                                self%darkMatterHaloScale_                  , &
               &                                                self%darkMatterProfile_                    , &
               &                                                self%darkMatterProfileScaleRadius_           &
               &                                               )
       end select
    end select
    ! Nullify objects to avoid them being destroyed.
    nullify(distributionLogNormal)
    ! Create a work node and set the appropriate cosmological time.
    node      => treeNode      (                 )
    nodeBasic => node    %basic(autoCreate=.true.)
    nodeSpin  => node    %spin (autoCreate=.true.)
    call nodeBasic%timeSet(self%time)
    ! Compute the spin distribution. The spin distribution is averaged over the width of each bin.
    allocate(distribution(size(self%spin)))
    do i=1,size(self%spin)
       distribution(i)=+Integrate(                                  &
            &                     self%spinMinimum(i)             , &
            &                     self%spinMaximum(i)             , &
            &                     spinDistributionIntegrate       , &
            &                     integrandFunction               , &
            &                     integrationWorkspace            , &
            &                     toleranceAbsolute        =0.0d+0, &
            &                     toleranceRelative        =1.0d-6  &
            &                    )                                  &
            &          /log10(                                      &
            &                 +self%spinMaximum(i)                  &
            &                 /self%spinMinimum(i)                  &
            &                )
       call Integrate_Done(integrandFunction,integrationWorkspace)
    end do
    ! Evaluate the log-likelihood.
    spinDistributionEvaluate=-0.5d0*sum(((distribution-self%distribution)/self%distributionError)**2,mask=self%distributionError > 0.0d0)
    ! Clean up.
    call node%destroy()
    deallocate(node        )
    deallocate(stateVector )
    deallocate(distribution)  
    return

  contains

    double precision function spinDistributionIntegrate(spinPrime)
      !% Integrand function used to find cumulative spin distribution over a bin.
      implicit none
      double precision, intent(in   ) :: spinPrime

      call nodeSpin%spinSet(spinPrime)
      spinDistributionIntegrate=distributionNbody%distributionAveraged(node,self%massHaloMinimum)
      return
    end function spinDistributionIntegrate

  end function spinDistributionEvaluate

  subroutine spinDistributionFunctionChanged(self)
    !% Respond to possible changes in the likelihood function.
    implicit none
    class(posteriorSampleLikelihoodSpinDistribution), intent(inout) :: self
    !GCC$ attributes unused :: self

    return
  end subroutine spinDistributionFunctionChanged
