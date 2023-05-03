!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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

  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Biases   , only : darkMatterHaloBiasClass
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Power_Spectra             , only : powerSpectrumClass
  use :: Numerical_Integration     , only : integrator
  
  ! Enumeration for likelihood model.
  !![
  <enumeration>
   <name>haloMassFunctionLikelihoodModel</name>
   <description>Specifies the likelihood model to use.</description>
   <visibility>public</visibility>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <entry label="uncorrelated"    />
   <entry label="simulationCube"  />
   <entry label="simulationSphere"/>
  </enumeration>
  !!]

  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodHaloMassFunction">
   <description>A posterior sampling likelihood class which implements a likelihood for halo mass functions.</description>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodBaseParameters) :: posteriorSampleLikelihoodHaloMassFunction
     !!{
     Implementation of a posterior sampling likelihood class which implements a likelihood for halo mass functions.
     !!}
     private
     class           (cosmologyParametersClass                      ), pointer                     :: cosmologyParameters_      => null()
     class           (cosmologyFunctionsClass                       ), pointer                     :: cosmologyFunctions_       => null()
     class           (darkMatterHaloBiasClass                       ), pointer                     :: darkMatterHaloBias_       => null()
     class           (powerSpectrumClass                            ), pointer                     :: powerSpectrum_            => null()
     class           (cosmologicalMassVarianceClass                 ), pointer                     :: cosmologicalMassVariance_ => null()
     double precision                                                , dimension(:  ), allocatable :: mass                               , massFunction                      , &
          &                                                                                           massMinimum                        , massMaximum
     double precision                                                , dimension(:,:), allocatable :: covarianceMatrix
     integer         (c_size_t                                      ), dimension(:  ), allocatable :: countHalos
     double precision                                                                              :: time                               , massParticle                      , &
          &                                                                                           massRangeMinimum                   , massRangeMaximum                  , &
          &                                                                                           countConversionFactor              , redshift                          , &
          &                                                                                           varianceSimulation                 , lengthSimulationCube              , &
          &                                                                                           massSphere                         , varianceFractionalModelDiscrepancy
     logical                                                                                       :: likelihoodPoisson                  , truncatePower                     , &
          &                                                                                           includeDiscrepancyChecked
     integer                                                                                       :: binCountMinimum                    , indexDiscrepancy
     type            (vector                                        )                              :: means
     type            (matrix                                        )                              :: covariance
     type            (varying_string                                )                              :: fileName
     type            (enumerationHaloMassFunctionLikelihoodModelType)                              :: likelihoodModel
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

  class           (powerSpectrumClass), pointer :: powerSpectrum__
  type            (integrator        )          :: integratorX            , integratorY      , &
       &                                           integratorZ            , integratorK
  double precision                              :: wavenumberX            , windowFunctionX  , &
       &                                           wavenumberY            , windowFunctionY  , &
       &                                           wavenumberZ            , windowFunctionZ  , &
       &                                           wavenumberMinimum      , wavenumberMaximum, &
       &                                           lengthSimulationCube_  , time_            , &
       &                                           radiusSimulationSphere_
  logical                                       :: truncatePower_
  !$omp threadprivate(powerSpectrum__,wavenumberX,windowFunctionX,wavenumberY,windowFunctionY,wavenumberZ,windowFunctionZ,lengthSimulationCube_,radiusSimulationSphere_,truncatePower_,integratorX,integratorY,integratorZ,integratorK,wavenumberMinimum,wavenumberMaximum,time_)
  
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
    type            (varying_string                           )                :: fileName                 , baseParametersFileName
    double precision                                                           :: redshift                 , massRangeMinimum                  , &
         &                                                                        massRangeMaximum         , lengthSimulationCube              , &
         &                                                                        massSphere               , varianceFractionalModelDiscrepancy
    integer                                                                    :: binCountMinimum
    logical                                                                    :: likelihoodPoisson        , truncatePower
    type            (inputParameters                          ), pointer       :: parametersModel
    class           (cosmologyParametersClass                 ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                  ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloBiasClass                  ), pointer       :: darkMatterHaloBias_
    class           (powerSpectrumClass                       ), pointer       :: powerSpectrum_
    class           (cosmologicalMassVarianceClass            ), pointer       :: cosmologicalMassVariance_
    type            (varying_string                           )                :: likelihoodModel

    !![
    <inputParameter>
      <name>baseParametersFileName</name>
      <description>The base set of parameters to use.</description>
      <source>parameters</source>
    </inputParameter>
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
      <name>massRangeMaximum</name>
      <description>The maximum halo mass to include in the likelihood evaluation.</description>
      <source>parameters</source>
      <defaultValue>huge(0.0d0)</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>binCountMinimum</name>
      <description>The minimum number of halos per bin required to permit bin to be included in likelihood evaluation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>likelihoodPoisson</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, likelihood is computed assuming a Poisson distribution for the number of halos in each bin (with no covariance between bins). Otherwise a multivariate normal is assumed when computing likelihood.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>likelihoodModel</name>
      <defaultValue>var_str('uncorrelated')</defaultValue>
      <description>The likelihood model to use for this halo mass function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>varianceFractionalModelDiscrepancy</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The fractional variance due to model discrepancy.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    if      (likelihoodModel == 'simulationCube'  ) then
       !![
       <inputParameter>
	 <name>lengthSimulationCube</name>
	 <description>The length of the simulation cube from which the target mass function was derived.</description>
	 <source>parameters</source>
       </inputParameter>
       !!]
    else if (likelihoodModel == 'simulationSphere') then
       !![
       <inputParameter>
	 <name>massSphere</name>
	 <description>The mass of the sphere from which the target mass function was derived.</description>
	 <source>parameters</source>
       </inputParameter>
       !!]
    end if
    if     (                                       &
         &   likelihoodModel == 'simulationCube'   &
         &  .or.                                   &
         &   likelihoodModel == 'simulationSphere' &
         & ) then
       !![
       <inputParameter>
	 <name>truncatePower</name>
	 <description>If true, truncate power on scales larger than the simulation volume when computing variance.</description>
	 <source>parameters</source>
       </inputParameter>
       !!]
    end if
    allocate(parametersModel)
    parametersModel=inputParameters(baseParametersFileName,noOutput=.true.)
    if (parametersModel%isPresent('powerSpectrumSimulation')) then
       !![
       <objectBuilder class="powerSpectrum" name="powerSpectrum_" source="parametersModel" parameterName="powerSpectrumSimulation"/>
       !!]
    else
       !![
       <objectBuilder class="powerSpectrum" name="powerSpectrum_" source="parametersModel"                                        />
       !!]
    end if
    !![
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parametersModel"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parametersModel"/>
    <objectBuilder class="darkMatterHaloBias"       name="darkMatterHaloBias_"       source="parametersModel"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parametersModel"/>
    <conditionalCall>
      <call>self=posteriorSampleLikelihoodHaloMassFunction(char(fileName),redshift,massRangeMinimum,massRangeMaximum,binCountMinimum,likelihoodPoisson,varianceFractionalModelDiscrepancy,enumerationHaloMassFunctionLikelihoodModelEncode(char(likelihoodModel),includesPrefix=.false.),parametersModel,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloBias_,powerSpectrum_,cosmologicalMassVariance_{conditions})</call>
      <argument name="lengthSimulationCube" value="lengthSimulationCube" condition="likelihoodModel == 'simulationCube'"                                           />
      <argument name="massSphere"           value="massSphere"           condition="likelihoodModel == 'simulationSphere'"                                         />
      <argument name="truncatePower"        value="truncatePower"        condition="likelihoodModel == 'simulationCube' .or. likelihoodModel == 'simulationSphere'"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="darkMatterHaloBias_"      />
    <objectDestructor name="powerSpectrum_"           />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    nullify(parametersModel)
    return
  end function haloMassFunctionConstructorParameters

  function haloMassFunctionConstructorInternal(fileName,redshift,massRangeMinimum,massRangeMaximum,binCountMinimum,likelihoodPoisson,varianceFractionalModelDiscrepancy,likelihoodModel,parametersModel,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloBias_,powerSpectrum_,cosmologicalMassVariance_,lengthSimulationCube,massSphere,truncatePower) result(self)
    !!{
    Constructor for ``haloMassFunction'' posterior sampling likelihood class.
    !!}
    use :: Display                 , only : displayMessage  , displayMagenta     , displayReset
    use :: File_Utilities          , only : File_Name_Expand, File_Exists        , File_Lock   , File_Unlock, &
         &                                  lockDescriptor
    use :: Error                   , only : Error_Report
    use :: HDF5_Access             , only : hdf5Access
    use :: IO_HDF5                 , only : hdf5Object
    use :: Input_Paths             , only : inputPath       , pathTypeDataDynamic
    use :: ISO_Varying_String      , only : char
    use :: Linear_Algebra          , only : assignment(=)
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (posteriorSampleLikelihoodHaloMassFunction     )                                :: self
    character       (len=*                                         ), intent(in   )                 :: fileName
    double precision                                                , intent(in   )                 :: redshift                            , massRangeMinimum                  , &
         &                                                                                             massRangeMaximum                    , varianceFractionalModelDiscrepancy
    integer                                                         , intent(in   )                 :: binCountMinimum
    logical                                                         , intent(in   )                 :: likelihoodPoisson
    type            (enumerationHaloMassFunctionLikelihoodModelType), intent(in   )                 :: likelihoodModel
    type            (inputParameters                               ), intent(inout), target         :: parametersModel
    class           (cosmologyParametersClass                      ), intent(inout), target         :: cosmologyParameters_
    class           (cosmologyFunctionsClass                       ), intent(inout), target         :: cosmologyFunctions_
    class           (darkMatterHaloBiasClass                       ), intent(inout), target         :: darkMatterHaloBias_
    class           (powerSpectrumClass                            ), intent(inout), target         :: powerSpectrum_
    class           (cosmologicalMassVarianceClass                 ), intent(inout), target         :: cosmologicalMassVariance_
    double precision                                                , intent(in   ), optional       :: lengthSimulationCube                , massSphere
    logical                                                         , intent(in   ), optional       :: truncatePower
    double precision                                                , allocatable  , dimension(:  ) :: eigenValueArray                     , massOriginal                      , &
         &                                                                                             massFunctionOriginal
    integer         (c_size_t                                      ), allocatable  , dimension(:  ) :: massFunctionCountOriginal
    double precision                                                , allocatable  , dimension(:,:) :: massFunctionCovarianceOriginal
    double precision                                                , parameter                     :: wavenumberMaximumFractional   =1.0d2
    character       (len=12                                        )                                :: redshiftLabel                       , lengthCubeLabel                   , &
         &                                                                                             timeLabel
    type            (hdf5Object                                    )                                :: massFunctionFile                    , simulationGroup                   , &
         &                                                                                             varianceFile
    integer                                                                                         :: i                                   , j                                 , &
         &                                                                                             ii                                  , jj                                , &
         &                                                                                             massCountReduced
    double precision                                                                                :: massIntervalLogarithmic
    type            (matrix                                        )                                :: eigenVectors
    type            (vector                                        )                                :: eigenValues
    type            (varying_string                                )                                :: fileNameVariance
    type            (lockDescriptor                                )                                :: fileLock
    !![
    <constructorAssign variables="fileName, redshift, binCountMinimum, massRangeMinimum, massRangeMaximum, likelihoodPoisson, varianceFractionalModelDiscrepancy, likelihoodModel, *parametersModel, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterHaloBias_, *powerSpectrum_, *cosmologicalMassVariance_, lengthSimulationCube, massSphere, truncatePower"/>
    !!]

    ! Convert redshift to time.
    self%time=cosmologyFunctions_%cosmicTime                 (          &
         &    cosmologyFunctions_%expansionFactorFromRedshift (         &
         &                                                     redshift &
         &                                                    )         &
         &                                                   )
    ! Record that we have not yet checked if model discrepancy terms should be included.
    self%includeDiscrepancyChecked=.false.
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
    ! Find the variance of the simulation.
    if (self%likelihoodModel == haloMassFunctionLikelihoodModelSimulationCube) then
       ! Find the variance of a cosmological box simulation.
       write (lengthCubeLabel,'(e12.6)') self%lengthSimulationCube
       write (      timeLabel,'(e12.6)') self%time
       fileNameVariance=inputPath(pathTypeDataDynamic)                                  // &
         &              'largeScaleStructure/'                                          // &
         &              self%objectType()                                               // &
         &              'CubeVariance_'                                                 // &
         &              'lengthSimulationCube_'//trim(lengthCubeLabel)//'_'             // &
         &              'time_'                //trim(      timeLabel)//'_'             // &
         &              self%powerSpectrum_%hashedDescriptor(includeSourceDigest=.true.)// &
         &              '.hdf5'
       if (File_Exists(fileNameVariance)) then
          call File_Lock(char(fileNameVariance),fileLock,lockIsShared=.true.)
          !$ call hdf5Access%set()
          call varianceFile%openFile     (char(fileNameVariance),readOnly=.true.                 )
          call varianceFile%readAttribute('variance'             ,        self%varianceSimulation)
          call varianceFile%close        (                                                       )
          !$ call hdf5Access%unset()
          call File_Unlock(fileLock)
       else
          powerSpectrum__         => powerSpectrum_
          truncatePower_          =  self%truncatePower
          time_                   =  self%time
          lengthSimulationCube_   =  self%lengthSimulationCube
          integratorX             =  integrator(integrandVarianceSimulationX,toleranceRelative=1.0d-6)
          integratorY             =  integrator(integrandVarianceSimulationY,toleranceRelative=1.0d-6)
          integratorZ             =  integrator(integrandVarianceSimulationZ,toleranceRelative=1.0d-6)
          wavenumberMinimum       =  +0.0d0
          wavenumberMaximum       =  +wavenumberMaximumFractional/self%lengthSimulationCube
          self%varianceSimulation =  integratorX%integrate(wavenumberMinimum,wavenumberMaximum)
          call File_Lock(char(fileNameVariance),fileLock,lockIsShared=.false.)
          !$ call hdf5Access%set()
          call varianceFile%openFile      (char(fileNameVariance) ,readOnly=.false.   ,overWrite=.true.)
          call varianceFile%writeAttribute(self%varianceSimulation,         'variance'                 )
          call varianceFile%close         (                                                            )
          !$ call hdf5Access%unset()
          call File_Unlock(fileLock)
       end if
    else if (self%likelihoodModel == haloMassFunctionLikelihoodModelSimulationSphere) then
       ! Find the variance in a spherical region.
       powerSpectrum__         => powerSpectrum_
       truncatePower_          =  self%truncatePower
       time_                   =  self%time
       radiusSimulationSphere_ =  (                                             &
            &                      +3.0d0                                       &
            &                      *self                     %massSphere        &
            &                      /4.0d0                                       &
            &                      /Pi                                          &
            &                      /self%cosmologyParameters_%OmegaMatter    () &
            &                      /self%cosmologyParameters_%densityCritical() &
            &                     )**(1.0d0/3.0d0)
       integratorK             =  integrator(integrandVarianceSimulationK,toleranceRelative=1.0d-6)
       wavenumberMinimum       =  +0.0d0
       wavenumberMaximum       =  +wavenumberMaximumFractional/radiusSimulationSphere_
       self%varianceSimulation =  integratorK%integrate(wavenumberMinimum,wavenumberMaximum)
    else
       self%varianceSimulation=0.0d0
    end if
    ! Compute quantities needed for likelihood calculations.
    if (self%likelihoodPoisson) then
       ! Find a reduced mass function excluding bins below the mass threshold.
       massCountReduced=0
       do i=1,size(massOriginal)       
          if     (                                    &
               &   massOriginal(i) < massRangeMinimum &
               &  .or.                                &
               &   massOriginal(i) > massRangeMaximum &
               & ) cycle
          massCountReduced=massCountReduced+1
       end do
       if (massCountReduced == 0) call Error_Report("no usable bins in mass function from file '"//trim(fileName)//"'"//{introspection:location})
       allocate(self%mass        (massCountReduced))
       allocate(self%massFunction(massCountReduced))
       allocate(self%countHalos  (massCountReduced))
       ii=0
       do i=1,size(massOriginal)
          if     (                                    &
               &   massOriginal(i) < massRangeMinimum &
               &  .or.                                &
               &   massOriginal(i) > massRangeMaximum &
               & ) cycle
          ii=ii+1
          self%mass        (ii)=massOriginal             (i)
          self%massFunction(ii)=massFunctionOriginal     (i)
          self%countHalos  (ii)=massFunctionCountOriginal(i)
       end do
       ! Compute the conversion factor between halo count per bin and the mass function.
       self%countConversionFactor=+     sum  (dble(self%countHalos)/self%massFunction,mask=self%massFunction > 0.0d0)  &
            &                     /dble(count(                                        mask=self%massFunction > 0.0d0))
    else
       ! Construct the covariance matrix.
       allocate(massFunctionCovarianceOriginal(size(massOriginal),size(massOriginal)))
       massFunctionCovarianceOriginal=0.0d0
       do i=1,size(massOriginal)
          do j=1,size(massOriginal)
             if   (                                           &
                &   massFunctionCountOriginal(i) > 0_c_size_t &
                &  .and.                                      &
                &   massFunctionCountOriginal(j) > 0_c_size_t &
                & ) then
                ! Compute the cosmic variance contribution.
                massFunctionCovarianceOriginal(i,j)=+                              massFunctionOriginal(i)           *                              massFunctionOriginal(j)            &
                     &                              *self%darkMatterHaloBias_%bias(massOriginal        (i),self%time)*self%darkMatterHaloBias_%bias(massOriginal        (j),self%time) &
                     &                              *self%varianceSimulation
                ! Compute the Poisson contribution.
                if (i == j) massFunctionCovarianceOriginal(i,j)=+     massFunctionCovarianceOriginal(i,j)     &
                     &                                          +     massFunctionOriginal          (i  ) **2 &
                     &                                          /dble(massFunctionCountOriginal     (i  ))
             end if
          end do
       end do
       ! Find a reduced mass function excluding any empty bins.
       massCountReduced=0
       do i=1,size(massOriginal)       
          if (massFunctionOriginal     (i) <= 0.0d0           ) cycle
          if (massOriginal             (i) <  massRangeMinimum) cycle
          if (massOriginal             (i) >  massRangeMaximum) cycle
          if (massFunctionCountOriginal(i) <  binCountMinimum ) cycle
          massCountReduced=massCountReduced+1
       end do
       if (massCountReduced == 0) call Error_Report("no usable bins in mass function from file '"//trim(fileName)//"'"//{introspection:location})
       allocate(self%mass            (massCountReduced                 ))
       allocate(self%massFunction    (massCountReduced                 ))
       allocate(self%covarianceMatrix(massCountReduced,massCountReduced))
       ii=0
       do i=1,size(massOriginal)
          if (massFunctionOriginal     (i) <= 0.0d0           ) cycle
          if (massOriginal             (i) <  massRangeMinimum) cycle
          if (massOriginal             (i) >  massRangeMaximum) cycle
          if (massFunctionCountOriginal(i) <  binCountMinimum ) cycle
          ii=ii+1
          self%mass        (ii)=massOriginal        (i)
          self%massFunction(ii)=massFunctionOriginal(i)
          jj=0
          do j=1,size(massOriginal)
             if (massFunctionOriginal     (j) <= 0.0d0           ) cycle
             if (massOriginal             (j) <  massRangeMinimum) cycle
             if (massOriginal             (j) >  massRangeMaximum) cycle
             if (massFunctionCountOriginal(j) <  binCountMinimum ) cycle
             jj=jj+1
             self%covarianceMatrix(ii,jj)=massFunctionCovarianceOriginal(i,j)
          end do
       end do
       ! Find the covariance matrix.
       self%covariance=self%covarianceMatrix
       ! Get eigenvalues and vectors of the covariance matrix.
       allocate(eigenValueArray(size(self%mass)))
       call self%covariance%eigenSystem(eigenVectors,eigenValues)
       eigenValueArray=eigenValues
       if (any(eigenValueArray < 0.0d0)) call displayMessage(displayMagenta()//'WARNING:'//displayReset()//' inverse covariance matrix is not semi-positive definite')
       deallocate(eigenValueArray               )
    end if
    ! Compute mass ranges for bins.
    massIntervalLogarithmic=+log(                                  &
         &                       +massOriginal(size(massOriginal)) &
         &                       /massOriginal(                 1) &
         &                      )                                  &
         &                  /dble(                                 &
         &                        +size(massOriginal)              &
         &                        -1                               &
         &                       )
    allocate(self%massMinimum,mold=self%mass)
    allocate(self%massMaximum,mold=self%mass)
    do i=1,size(self%mass)
       self%massMinimum(i)=self%mass(i)*exp(-0.5d0*massIntervalLogarithmic)
       self%massMaximum(i)=self%mass(i)*exp(+0.5d0*massIntervalLogarithmic)
    end do
    if (allocated(massOriginal                  )) deallocate(massOriginal                  )
    if (allocated(massFunctionOriginal          )) deallocate(massFunctionOriginal          )
    if (allocated(massFunctionCovarianceOriginal)) deallocate(massFunctionCovarianceOriginal)
    return
  end function haloMassFunctionConstructorInternal

  double precision function integrandVarianceSimulationX(wavenumber)
    !!{
    Integrand function for cosmic variance in a simulation cube: $x$-axis.
    !!}
    implicit none
    double precision, intent(in   ) :: wavenumber

    wavenumberX                 =wavenumber
    windowFunctionX             =windowFunctionCube(wavenumberX*lengthSimulationCube_)**2
    integrandVarianceSimulationX=integratorY%integrate(wavenumberMinimum,wavenumberMaximum)
    return
  end function integrandVarianceSimulationX
  
  double precision function integrandVarianceSimulationY(wavenumber)
    !!{
    Integrand function for cosmic variance in a simulation cube: $y$-axis.
    !!}
    implicit none
    double precision, intent(in   ) :: wavenumber

    wavenumberY                 =wavenumber
    windowFunctionY             =windowFunctionCube(wavenumberY*lengthSimulationCube_)**2
    integrandVarianceSimulationY=integratorZ%integrate(wavenumberMinimum,wavenumberMaximum)
    return
  end function integrandVarianceSimulationY
  
  double precision function integrandVarianceSimulationZ(wavenumber)
    !!{
    Integrand function for cosmic variance in a simulation cube: $z$-axis.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: wavenumber
    double precision                :: wavenumberMagnitude
    
    wavenumberZ                 =wavenumber
    ! Truncate the power spectrum on scales larger than the simulation volume if necessary.
    if     (                                             &
         &   truncatePower_                              &
         &  .and.                                        &
         &   abs(wavenumberX) < Pi/lengthSimulationCube_ &
         &  .and.                                        &
         &   abs(wavenumberX) < Pi/lengthSimulationCube_ &
         &  .and.                                        &
         &   abs(wavenumberX) < Pi/lengthSimulationCube_ &
         & ) then
       ! Truncate power missing from the cube.
       integrandVarianceSimulationZ=+0.0d0
    else
       windowFunctionZ             =+windowFunctionCube(wavenumberZ*lengthSimulationCube_)**2
       wavenumberMagnitude         =+sqrt(+wavenumberX**2+wavenumberY**2+wavenumberZ**2)    
       if (wavenumberMagnitude == 0.0d0) then
          integrandVarianceSimulationZ=+0.0d0
       else
          ! We integrate only over positive wavenumbers, so we include a factor 8 here to account for the other octants.
          integrandVarianceSimulationZ=+8.0d0                                            &
               &                       *powerSpectrum__%power(wavenumberMagnitude,time_) &
               &                       *windowFunctionX                                  &
               &                       *windowFunctionY                                  &
               &                       *windowFunctionZ
       end if
    end if
    return
  end function integrandVarianceSimulationZ

  double precision function integrandVarianceSimulationK(wavenumber)
    !!{
    Integrand function for cosmic variance in a simulation sphere.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: wavenumber
    double precision                :: windowFunction
    
    ! Truncate the power spectrum on scales larger than the simulation volume if necessary.
    if     (                                         &
         &   truncatePower_                          &
         &  .and.                                    &
         &   wavenumber < Pi/radiusSimulationSphere_ &        
         & ) then
       ! Truncate power missing from the cube.
       integrandVarianceSimulationK=+0.0d0
    else
       windowFunction              =+windowFunctionSphere(wavenumber*radiusSimulationSphere_)**2
       integrandVarianceSimulationK=+4.0d0                                      &
            &                       *Pi                                         &
            &                       *wavenumber                             **2 &
            &                       *powerSpectrum__%power(wavenumber,time_)    &
            &                       *windowFunction
    end if
    return
  end function integrandVarianceSimulationK

  double precision function windowFunctionCube(x)
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: x

    if (x == 0.0d0) then
       windowFunctionCube=+1.0d0               /sqrt(2.0d0*Pi)
    else
       windowFunctionCube=+2.0d0*sin(0.5d0*x)/x/sqrt(2.0d0*Pi)
    end if
    return
  end function windowFunctionCube
  
  double precision function windowFunctionSphere(x)
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: x

    if (x == 0.0d0) then
       windowFunctionSphere=+1.0d0                         /sqrt(2.0d0*Pi)**3
    else
       windowFunctionSphere=+3.0d0*(sin(x)-x*cos(x))/(x**3)/sqrt(2.0d0*Pi)**3
    end if
    return
  end function windowFunctionSphere
  
  subroutine haloMassFunctionDestructor(self)
    !!{
    Destructor for ``haloMassFunction'' posterior sampling likelihood class.
    !!}
    implicit none
    type(posteriorSampleLikelihoodHaloMassFunction), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%darkMatterHaloBias_"      />
    <objectDestructor name="self%powerSpectrum_"           />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    if (associated(self%parametersModel)) then
       call self%parametersModel%destroy()
       deallocate(self%parametersModel)
    end if
    return
  end subroutine haloMassFunctionDestructor

  double precision function haloMassFunctionEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !!{
    Return the log-likelihood for the halo mass function likelihood function. Several options are supported here:
    \begin{itemize}
    \item Poisson or Gaussian statistics:
      \begin{itemize}
      \item If {\normalfont \ttfamily [likelihoodPoisson]=false} then Gaussian statistics are assumed;
      \item otherwise, Poisson statistics are assumed.
      \end{itemize}
    \item Covariance model:
      \begin{itemize}
      \item If {\normalfont \ttfamily [likelihoodModel]=uncorrelated}, then assume no covariance between bins;
      \item otherwise, we assume that there is some fractional variance in the cosmological density field in the simulation volume
        of $\sigma_\mathrm{V}^2$. This is computed as
        \begin{equation}
          \sigma_\mathrm{V}^2 = \frac{1}{(2\pi)^3} \int \mathrm{d}^3\mathbf{k} |W(\mathbf{k})|^2 P(k).
        \end{equation}
        If
        \begin{itemize}
        \item {\normalfont \ttfamily [likelihoodModel]=simulationCube} then the simulation is assumed to be a cubic region with
          length {\normalfont \ttfamily [lengthSimulationCube]}, and the window function in the above equation is given by
          \begin{equation}
            W(\mathbf{k}) = \frac{1}{(2\pi)^{3/2}} \hbox{sinc}\left(\frac{L k_x}{2}\right) \hbox{sinc}\left(\frac{L k_y}{2}\right) \hbox{sinc}\left(\frac{L k_z}{2}\right);
          \end{equation}
        \item {\normalfont \ttfamily [likelihoodModel]=simulationSphere} then the simulation is assumed to be a spherical region with
          mean mass given by {\normalfont \ttfamily [massSphere]} and the window function is therefore:
          \begin{equation}
            W(k) = \frac{3}{(2\pi)^{3/2}} \frac{\sin(k R)- k R \cos (k R)}{(k R)^3},
          \end{equation}
          where $R = (3 M/4 \pi \bar{\rho})^{1/3}$ is the radius of the simulated region.
        \end{itemize}
        We model the distribution of this random variable as a log-normal, and define a variable $y$
        which is log-normally distributed such that the mean of $\log y$ is 0 and the variance of $\log y$ is 1 (so, a standard
        normal in the log).
        \begin{itemize}
        \item For Gaussian statistics, there is then a contribution to the covariance matrix:
          $C_{ij} = \phi_i \phi_j b_i b_j \sigma_\mathrm{V}^2$.
        \item For Poisson statistics the joint probability distribution is given by
          \begin{equation}
            p(N_1\ldots N_n) = \int_0^\infty \left( \prod_{i=1}^n \frac{\lambda_i^{N_i} \mathrm{e}^{-\lambda_i}}{N_i!} \right)\frac{\mathrm{e}^{-(\log y + 1)^2/2}}{\sqrt{2\pi}} \frac{\mathrm{d}y}{y}.
          \end{equation}
          The integral over the common log-normal distribution will then introduce covariance between the bins.

          Given the large values of the terms it is numerically convenient to write this integral as:
          \begin{equation}
            p(N_1\ldots N_n) = \frac{1}{\sqrt{2\pi}} \int_0^\infty \frac{\mathrm{d}y}{y} \exp\left[ f(y) \right],
          \end{equation}
          where
          \begin{equation}
            f(y) =  - \frac{(\log y  + 1)^2}{2} + \sum_{i=1}^n -\lambda_i + N_i \log \lambda_i - \log N_i!.
          \end{equation}
          This is still difficult to evaluate because $f(y)$ can be a large negative number, which results in floating point
          underflow when we compute $\exp[f(y)]$.

          To avoid this we first numerically solve for the value of $y$ that maximizes $f(y)$ - call this $y_\mathrm{max}$. Then,
          we compute
          \begin{equation}
            p(N_1\ldots N_n) = \frac{1}{\sqrt{2\pi}} \exp[f(y_\mathrm{max})] \int_0^\infty \frac{\mathrm{d}y}{y} \exp\left[ f(y) -f(y_\mathrm{max}) \right],
          \end{equation}
          which ensures that the argument of the exponential in the integral is always less than zero (so it cannot overflow), and
          will be $\approx 0$ close to the maximum of the function, so underflow will not be a problem. We never actually need to
          evaluate $\exp[f(y_\mathrm{max})]$ because we take the logarithm of this entire expression anyway to get
          $\log\mathcal{L}$.
        \end{itemize}
      \end{itemize}
    \end{itemize}
    !!}
    use :: Error                            , only : errorStatusSuccess
    use :: Halo_Mass_Functions              , only : haloMassFunctionClass
    use :: Linear_Algebra                   , only : assignment(=)                   , operator(*)
    use :: Models_Likelihoods_Constants     , only : logImpossible                   , logImprobable
    use :: Posterior_Sampling_Convergence   , only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State         , only : posteriorSampleStateClass
    use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass         , nbodyHaloMassErrorPowerLaw, nbodyHaloMassErrorSOHaloFinder, nbodyHaloMassErrorTrenti2010
    use :: Factorials                       , only : Logarithmic_Factorial
    use :: Gamma_Functions                  , only : Gamma_Function_Logarithmic
    use :: Root_Finder                      , only : rootFinder                      , rangeExpandMultiplicative , rangeExpandSignExpectPositive , rangeExpandSignExpectNegative
    implicit none
    class           (posteriorSampleLikelihoodHaloMassFunction), intent(inout)               :: self
    class           (posteriorSampleStateClass                ), intent(inout)               :: simulationState
    type            (modelParameterList                       ), intent(in   ), dimension(:) :: modelParametersActive_                       , modelParametersInactive_
    class           (posteriorSampleConvergenceClass          ), intent(inout)               :: simulationConvergence
    double precision                                           , intent(in   )               :: temperature                                  , logLikelihoodCurrent    , &
         &                                                                                      logPriorCurrent                              , logPriorProposed
    real                                                       , intent(inout)               :: timeEvaluate
    double precision                                           , intent(  out), optional     :: logLikelihoodVariance
    logical                                                    , intent(inout), optional     :: forceAcceptance
    double precision                                           , allocatable  , dimension(:) :: stateVector                                  , massFunction
    double precision                                           , parameter                   :: errorFractionalMaximum                =1.0d+1
    double precision                                           , parameter                   :: amplitudeFractionalPerturbationMinimum=1.0d-6, amplitudeFractionalPerturbationMaximum=1.0d+2
    class           (haloMassFunctionClass                    ), pointer                     :: haloMassFunction_
    type            (vector                                   )                              :: difference
    type            (integrator                               ), allocatable                 :: integrator_
    type            (rootFinder                               ), allocatable                 :: finder_
    logical                                                                                  :: evaluationFailed
    integer                                                                                  :: i                                            , status
    double precision                                                                         :: countHalosMean                               , likelihood                                   , &
         &                                                                                      argumentOffset                               , amplitudeFractionalPerturbationPeak          , &
         &                                                                                      varianceFractionalModelDiscrepancy           , stoppingTimeParameter
    !$GLC attributes unused :: simulationConvergence, temperature, timeEvaluate, logLikelihoodCurrent, logPriorCurrent, modelParametersInactive_, forceAcceptance

    ! There is no variance in our likelihood estimate.
    if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
    ! Do not evaluate if the proposed prior is impossible.
    if (logPriorProposed <= logImpossible) then
       haloMassFunctionEvaluate=0.0d0
       return
    end if
    ! Ensure pointers into the base parameters are initialized.
    call self%initialize(modelParametersActive_,modelParametersInactive_)
    ! Get states for all chains.
    allocate(stateVector(simulationState%dimension()))
    stateVector=simulationState%get()
    ! Update parameter values.
    call self%update(simulationState,modelParametersActive_,modelParametersInactive_,stateVector)
    if (.not.self%includeDiscrepancyChecked) then
       self%includeDiscrepancyChecked=.true.
       self%indexDiscrepancy         =-1
       do i=1,size(self%modelParametersActive_)
          if (self%modelParametersActive_(i)%definition == "varianceFractionalModelDiscrepancy") then
             if (self%indexDiscrepancy > 0) call Error_Report('multiple instances of parameter [varianceFractionalModelDiscrepancy] found'//{introspection:location})
             self%indexDiscrepancy=i
          end if
       end do
    end if
    ! Get the halo mass function object.
    !![
    <objectBuilder class="haloMassFunction" name="haloMassFunction_" source="self%parametersModel"/>
    !!]
    ! Determine the model discrepancy variance term.
    if (self%indexDiscrepancy > 0) then
       varianceFractionalModelDiscrepancy=stateVector(self%indexDiscrepancy)
    else
       varianceFractionalModelDiscrepancy=self%varianceFractionalModelDiscrepancy
    end if
    if (varianceFractionalModelDiscrepancy > 0.0d0) then
       stoppingTimeParameter=1.0d0/varianceFractionalModelDiscrepancy
    else
       stoppingTimeParameter=0.0d0
    end if
    ! Compute the mass function.
    allocate(massFunction(size(self%mass)))
    evaluationFailed=.false.
    do i=1,size(self%mass)
       massFunction(i)=+haloMassFunction_%integrated(                            &
            &                                               self%time          , &
            &                                               self%massMinimum(i), &
            &                                               self%massMaximum(i), &
            &                                        status=     status          &
            &                                       )                            &
            &          /log(                                                     &
            &                                       +self%massMaximum(i)         &
            &                                       /self%massMinimum(i)         &
            &              )
       if (status /= errorStatusSuccess) then
          haloMassFunctionEvaluate=logImprobable
          evaluationFailed        =.true.
          exit
       end if
    end do
    ! Evaluate the log-likelihood.
    if (.not.evaluationFailed) then
       if (self%likelihoodPoisson) then
          ! Assume Poisson statistics.
          if     (                                                                         &
               &   self%likelihoodModel == haloMassFunctionLikelihoodModelSimulationCube   &
               &  .or.                                                                     &
               &   self%likelihoodModel == haloMassFunctionLikelihoodModelSimulationSphere &
               & ) then
             ! Simulation variance is to be included. We model the likelihood as a mixture of Poisson with a log-normal distribution
             ! on the mean, with that log-normal correlated across all bins.
             if (any(massFunction <= 0.0d0 .and. self%countHalos > 0)) then
                ! If any bin is predicted to have a zero mass function, but the target has halos present, the model is improbable.
                haloMassFunctionEvaluate=logImprobable
             else
                ! First find the maximum of the argument appearing in the exponential term in the integrand. This allows us to avoid
                ! floating point over/underflow. We avoid amplitudes that are too small or too large as these are extremely unlikely -
                ! they can, however, be favored when the model mass function is extremely far from the target dataset. Such models
                ! will have extremely low likelihoods - it does not matter if their likelihood is poorly estimated.
                if      (argumentMaximumRoot(amplitudeFractionalPerturbationMinimum) < 0.0d0) then
                   amplitudeFractionalPerturbationPeak=amplitudeFractionalPerturbationMinimum
                else if (argumentMaximumRoot(amplitudeFractionalPerturbationMaximum) > 0.0d0) then
                   amplitudeFractionalPerturbationPeak=amplitudeFractionalPerturbationMaximum
                else
                   allocate(finder_)
                   finder_                            =rootFinder  (                                                             &
                        &                                           argumentMaximumRoot                                        , &
                        &                                           rangeExpandUpward            =2.0d0                        , &
                        &                                           rangeExpandDownward          =0.5d0                        , &
                        &                                           rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
                        &                                           rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
                        &                                           rangeExpandType              =rangeExpandMultiplicative    , &
                        &                                           toleranceRelative            =1.0d-3                         &
                        &                                          )
                   amplitudeFractionalPerturbationPeak=finder_%find(                                                             &
                        &                                           rootGuess                    =1.0d0                          &
                        &                                          )
                   deallocate(finder_)
                end if
                ! Compute the argument of the exponential at this maximum.
                argumentOffset=-0.5d0*(log(amplitudeFractionalPerturbationPeak)+1.0d0)**2
                do i=1,size(self%mass)
                   if (massFunction(i) <= 0.0d0) cycle
                   countHalosMean=+self%countConversionFactor                                                                                                      &
                        &         *     massFunction                                                                 (i)                                           &
                        &         *     amplitudeFractionalPerturbationPeak**(self%darkMatterHaloBias_%bias(self%mass(i),self%time)*sqrt(self%varianceSimulation))
                   if (varianceFractionalModelDiscrepancy <= 0.0d0) then
                      ! Evaluate the Poisson likelihood (zero model discrepancy term).
                      argumentOffset=+argumentOffset                                     &
                           &         +dble                 (    self%countHalos    (i))  &
                           &         *log                  (         countHalosMean   )  &
                           &         -                               countHalosMean      &
                           &         -Logarithmic_Factorial(int(self%countHalos    (i)))
                   else
                      ! Evaluate the negative binomial likelihood (non-zero model discrepancy term).
                      argumentOffset=+argumentOffset                                     &
                           &         +dble(self%countHalos           (i))*log                       (                                               countHalosMean    ) &
                           &         -                                    Logarithmic_Factorial     (                                    +int (self%countHalos    (i))) &
                           &         +                                    Gamma_Function_Logarithmic(               stoppingTimeParameter+dble(self%countHalos    (i))) &
                           &         -                                    Gamma_Function_Logarithmic(               stoppingTimeParameter                             ) &
                           &         -dble(self%countHalos           (i))*Gamma_Function_Logarithmic(               stoppingTimeParameter+          countHalosMean    ) &
                           &         -          stoppingTimeParameter    *log                       (countHalosMean/stoppingTimeParameter+1.0d0                       )
                   end if
                end do
                ! Integrate the Poisson distribution over the distribution of mean parameters.
                allocate(integrator_)             
                integrator_=integrator           (                                                          &
                     &                                              integrandLikelihoodPoisson            , &
                     &                            toleranceRelative=1.0d-3                                  &
                     &                           )
                likelihood =integrator_%integrate(                                                          &
                     &                                              amplitudeFractionalPerturbationMinimum, &
                     &                                              amplitudeFractionalPerturbationMaximum  &
                     &                           )
                deallocate(integrator_)
                ! Evaluate the log-likelihood, correcting for the maximum of the argument of the exponential which we factored out
                ! of the integrand.
                if (likelihood > 0.0d0) then
                   haloMassFunctionEvaluate=+log(likelihood    ) &
                        &                   +    argumentOffset
                else
                   haloMassFunctionEvaluate=logImprobable
                end if
             end if
          else
             ! No simulation variance is to be included. We treat each bin as independent with a pure Poisson/negative binomial distribution.
             haloMassFunctionEvaluate=0.0d0
             do i=1,size(self%mass)
                ! Find the mean number of halos expected in this bin based on our model mass function.
                countHalosMean          =+                          self%countConversionFactor         &
                     &                   *                               massFunction            (i)
                ! If the expected mean is zero, and the measured number is non-zero, this is impossible.
                if (countHalosMean <= 0.0d0) then
                   if (self%countHalos(i) > 0) then
                      haloMassFunctionEvaluate=logImprobable
                      exit
                   end if
                else
                   if (varianceFractionalModelDiscrepancy <= 0.0d0) then
                      ! Evaluate the Poisson likelihood (zero model discrepancy term).
                      haloMassFunctionEvaluate=+                               haloMassFunctionEvaluate      &
                           &                   +dble                 (    self%countHalos              (i))  &
                           &                   *log                  (         countHalosMean             )  &
                           &                   -                               countHalosMean                &
                           &                   -Logarithmic_Factorial(int(self%countHalos              (i)))
                   else
                      ! Evaluate the negative binomial likelihood (non-zero model discrepancy term).
                      haloMassFunctionEvaluate=+haloMassFunctionEvaluate                                                                                                          &
                           &                   +dble(self%countHalos           (i))*log                       (                                               countHalosMean    ) &
                           &                   -                                    Logarithmic_Factorial     (                                    +int (self%countHalos    (i))) &
                           &                   +                                    Gamma_Function_Logarithmic(               stoppingTimeParameter+dble(self%countHalos    (i))) &
                           &                   -                                    Gamma_Function_Logarithmic(               stoppingTimeParameter                             ) &
                           &                   -dble(self%countHalos           (i))*Gamma_Function_Logarithmic(               stoppingTimeParameter+          countHalosMean    ) &
                           &                   -          stoppingTimeParameter    *log                       (countHalosMean/stoppingTimeParameter+1.0d0                       )


write (0,*) "WTF ",i,haloMassFunctionEvaluate

                   end if
                end if
             end do
          end if
       else
          ! Assume Gaussian statistics.
          difference              =+massFunction                                  &
               &                   -self%massFunction
          haloMassFunctionEvaluate=-0.5d0                                         &
               &                   *self%covariance%covarianceProduct(difference)
       end if
    end if
    ! Clean up.
    !![
    <objectDestructor name="haloMassFunction_"/>
    !!]
    deallocate(stateVector )
    deallocate(massFunction)
    return

  contains

    double precision function integrandLikelihoodPoisson(amplitudeFractionalPerturbation)
      !!{
      Integrand function for Poission likelihood.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: amplitudeFractionalPerturbation
      double precision                :: argument
      integer                         :: i
      
      ! Term related to the log-normal distribution of the perturbation amplitude.
      argument=-0.5d0*(log(amplitudeFractionalPerturbation)+1.0d0)**2
      ! Sum terms corresponding to each mass bin.
      do i=1,size(self%mass)
         ! Skip empty bins.
         if (massFunction(i) <= 0.0d0) cycle
         ! Find the mean number of halos in this bin according to the model.
         countHalosMean=+self%countConversionFactor                                                                                                  &
              &         *     massFunction                                                             (i)                                           &
              &         *     amplitudeFractionalPerturbation**(self%darkMatterHaloBias_%bias(self%mass(i),self%time)*sqrt(self%varianceSimulation))
         if (varianceFractionalModelDiscrepancy <= 0.0d0) then
            ! Evaluate the Poisson likelihood (zero model discrepancy term).
            argument      =+argument                                           &
                 &         +dble                 (    self%countHalos    (i))  &
                 &         *log                  (         countHalosMean   )  &
                 &         -                               countHalosMean      &
                 &         -Logarithmic_Factorial(int(self%countHalos    (i)))
         else
            ! Evaluate the negative binomial likelihood (non-zero model discrepancy term).
            argument      =+argument                                                                                                                          &
                 &         +dble(self%countHalos           (i))*log                       (                                               countHalosMean    ) &
                 &         -                                    Logarithmic_Factorial     (                                    +int (self%countHalos    (i))) &
                 &         +                                    Gamma_Function_Logarithmic(               stoppingTimeParameter+dble(self%countHalos    (i))) &
                 &         -                                    Gamma_Function_Logarithmic(               stoppingTimeParameter                             ) &
                 &         -dble(self%countHalos           (i))*Gamma_Function_Logarithmic(               stoppingTimeParameter+          countHalosMean    ) &
                 &         -          stoppingTimeParameter    *log                       (countHalosMean/stoppingTimeParameter+1.0d0                       )
         end if
      end do
      integrandLikelihoodPoisson=+exp(                            &
           &                          +argument                   &
           &                          -argumentOffset             &
           &                         )                            &
           &                     /amplitudeFractionalPerturbation &
           &                     /sqrt(2.0d0*Pi)
      return
    end function integrandLikelihoodPoisson
    
    double precision function argumentMaximumRoot(amplitudeFractionalPerturbation)
      !!{
      Root function used to find the value of the perturbation amplitude which maximizes the integrand for the Poisson likelihood.
      !!}
      use :: Gamma_Functions, only : Digamma_Function
      implicit none
      double precision, intent(in   ) :: amplitudeFractionalPerturbation
      integer                         :: i
      double precision                :: amplitudeScaling

      ! Term related to the log-normal distribution of the perturbation amplitude.
      argumentMaximumRoot=-(+log(amplitudeFractionalPerturbation)+1.0d0) &
           &              /      amplitudeFractionalPerturbation
      ! Sum terms corresponding to each mass bin.
      do i=1,size(self%mass)
         ! Skip empty bins.
         if (massFunction(i) <= 0.0d0) cycle
         ! Find the mean number of halos in this bin according to the model.
         countHalosMean  =+     self                     %countConversionFactor                         &
              &           *                              massFunction          (          i           )
         ! Compute the contribution to the gradient of the argument.
         amplitudeScaling=+     self%darkMatterHaloBias_%bias                  (self%mass(i),self%time) &
              &           *sqrt(self                    %varianceSimulation                           )
         if (varianceFractionalModelDiscrepancy <= 0.0d0) then
            ! Evaluate the Poisson likelihood (zero model discrepancy term).
             argumentMaximumRoot=+argumentMaximumRoot                                                                                                                                                                &
              &                  +amplitudeScaling                  *dble(self%countHalos(i))*amplitudeFractionalPerturbation**(                -1.0d0)                                                              &              
              &                  -amplitudeScaling   *countHalosMean                         *amplitudeFractionalPerturbation**(amplitudeScaling-1.0d0)
         else
             argumentMaximumRoot=+argumentMaximumRoot                                                                                                                                                                &
                 &               +amplitudeScaling                  *dble(self%countHalos(i))*amplitudeFractionalPerturbation**(                -1.0d0)                                                              &
                 &               -amplitudeScaling   *countHalosMean                         *amplitudeFractionalPerturbation**(                -1.0d0)/                (1.0d0+countHalosMean/stoppingTimeParameter) &
                 &               -amplitudeScaling   *countHalosMean*dble(self%countHalos(i))*amplitudeFractionalPerturbation**(                -1.0d0)*Digamma_Function(      countHalosMean+stoppingTimeParameter)
         end if
      end do
      return
    end function argumentMaximumRoot
    
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
