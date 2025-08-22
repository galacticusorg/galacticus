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
  An implementation of the lightcone geometry class which assumes a cylindrical ``cone''.
  !!}

  use :: Cosmology_Functions            , only : cosmologyFunctionsClass
  use :: Correlation_Functions_Two_Point, only : correlationFunctionTwoPointClass
  use :: Dark_Matter_Halo_Biases        , only : darkMatterHaloBiasClass
  use :: Dark_Matter_Halo_Scales        , only : darkMatterHaloScaleClass
  use :: Kind_Numbers                   , only : kind_int8
  use :: Linear_Growth                  , only : linearGrowthClass
  use :: Numerical_Random_Numbers       , only : randomNumberGeneratorClass
  use :: Output_Times                   , only : outputTimesClass
  use :: Power_Spectra                  , only : powerSpectrumClass
  
  !![
  <geometryLightcone name="geometryLightconeCylindrical">
   <description>
    A lightcone geometry class which assumes a cylindrical ``cone'', i.e. defined such that a point $(x,y,z)$ is in the survey
    if $\sqrt{x^2+y^2} &lt; r$, where $r$ is the radius of the ``cone''.
   </description>
  </geometryLightcone>
  !!]
  type, extends(geometryLightconeClass) :: geometryLightconeCylindrical
     !!{
     A lightcone geometry class which assumes a cylindrical ``cone''.
     !!}
     private
     class           (cosmologyFunctionsClass         ), pointer                     :: cosmologyFunctions_          => null()
     class           (outputTimesClass                ), pointer                     :: outputTimes_                 => null()
     class           (randomNumberGeneratorClass      ), pointer                     :: randomNumberGenerator_       => null()
     class           (powerSpectrumClass              ), pointer                     :: powerSpectrum_               => null()
     class           (linearGrowthClass               ), pointer                     :: linearGrowth_                => null()
     class           (darkMatterHaloBiasClass         ), pointer                     :: darkMatterHaloBias_          => null()
     class           (darkMatterHaloScaleClass        ), pointer                     :: darkMatterHaloScale_         => null()
     class           (correlationFunctionTwoPointClass), pointer                     :: correlationFunctionTwoPoint_ => null()
     double precision                                  , dimension(:  ), allocatable :: distanceMinimum                       , distanceMaximum       , &
          &                                                                             volume                                , densityContrast
     double precision                                                                :: radiusCylinderComoving                , radiusBufferComoving  , &
          &                                                                             massHaloLens                          , redshiftLens          , &
          &                                                                             distanceComovingLens                  , timeLens              , &
          &                                                                             radiusVirialComovingLens              , correlationLensMaximum
     integer         (kind_int8                       )                              :: activeCenterUniqueID                  , activeUniqueID
     integer                                                                         :: activeCenterCount                     , activeCount
     logical                                           , dimension(:  ), allocatable :: activeInCylinder
     integer                                           , dimension(:  ), allocatable :: activeInstances
     double precision                                  , dimension(:,:), allocatable :: activeCenterPosition                  , activePosition
     integer         (c_size_t                        )                              :: activeOutput
   contains
     !![
     <methods>
      <method method="sampleNode">
       <description>Determine if, and how many times, the given node appears in the ``lightcone'', and choose positions for each instance.</description>
      </method>
     </methods>
     !!]
     final     ::                              cylindricalDestructor
     procedure :: timeMinimum               => cylindricalTimeMinimum
     procedure :: timeMaximum               => cylindricalTimeMaximum
     procedure :: isInLightcone             => cylindricalIsInLightcone
     procedure :: replicationCount          => cylindricalReplicationCount
     procedure :: solidAngle                => cylindricalSolidAngle
     procedure :: position                  => cylindricalPosition
     procedure :: velocity                  => cylindricalVelocity
     procedure :: timeLightconeCrossing     => cylindricalTimeLightconeCrossing
     procedure :: positionLightconeCrossing => cylindricalPositionLightconeCrossing
     procedure :: sampleNode                => cylindricalSampleNode
  end type geometryLightconeCylindrical

  interface geometryLightconeCylindrical
     !!{
     Constructors for the \refClass{geometryLightconeCylindrical} dark matter halo spin distribution class.
     !!}
     module procedure cylindricalConstructorParameters
     module procedure cylindricalConstructorInternal
  end interface geometryLightconeCylindrical

  ! Tolerance for matching to output times.
  double precision, parameter :: timeTolerance =1.0d-3

contains

  function cylindricalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{geometryLightconeCylindrical} lightcone geometry distribution class which takes a parameter list as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (geometryLightconeCylindrical    )                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (outputTimesClass                ), pointer       :: outputTimes_
    class           (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    class           (randomNumberGeneratorClass      ), pointer       :: randomNumberGenerator_
    class           (powerSpectrumClass              ), pointer       :: powerSpectrum_
    class           (linearGrowthClass               ), pointer       :: linearGrowth_
    class           (darkMatterHaloBiasClass         ), pointer       :: darkMatterHaloBias_
    class           (darkMatterHaloScaleClass        ), pointer       :: darkMatterHaloScale_
    class           (correlationFunctionTwoPointClass), pointer       :: correlationFunctionTwoPoint_
    double precision                                                  :: radiusCylinderComoving      , radiusBufferComoving, &
         &                                                               massHaloLens                , redshiftLens
    
    !![
    <inputParameter>
      <name>radiusCylinderComoving</name>
      <source>parameters</source>
      <description>The comoving radius of the cylinder to populate.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusBufferComoving</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The comoving buffer radius to add around the cylinder. This is used to ensure that the sample within the cylinder is complete.</description>
    </inputParameter>
    <inputParameter>
      <name>massHaloLens</name>
      <source>parameters</source>
      <defaultValue>-1.0d0</defaultValue>
      <description>The mass of the primary lens halo (or a negative value for no lens).</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftLens</name>
      <source>parameters</source>
      <defaultValue>-1.0d0</defaultValue>
      <description>The redshift of the primary lens halo (or a negative value for no lens).</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"          source="parameters"/>
    <objectBuilder class="outputTimes"                  name="outputTimes_"                 source="parameters"/>
    <objectBuilder class="randomNumberGenerator"        name="randomNumberGenerator_"       source="parameters"/>
    <objectBuilder class="powerSpectrum"                name="powerSpectrum_"               source="parameters"/>
    <objectBuilder class="linearGrowth"                 name="linearGrowth_"                source="parameters"/>
    <objectBuilder class="darkMatterHaloBias"           name="darkMatterHaloBias_"          source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"         source="parameters"/>
    <objectBuilder class="correlationFunctionTwoPoint"  name="correlationFunctionTwoPoint_" source="parameters"/>
    !!]
    self=geometryLightconeCylindrical(radiusCylinderComoving,radiusBufferComoving,massHaloLens,redshiftLens,cosmologyFunctions_,powerSpectrum_,linearGrowth_,outputTimes_,darkMatterHaloBias_,darkMatterHaloScale_,correlationFunctionTwoPoint_,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"                />
    <objectDestructor name="cosmologyFunctions_"         />
    <objectDestructor name="randomNumberGenerator_"      />
    <objectDestructor name="powerSpectrum_"              />
    <objectDestructor name="linearGrowth_"               />
    <objectDestructor name="darkMatterHaloBias_"         />
    <objectDestructor name="darkMatterHaloScale_"        />
    <objectDestructor name="correlationFunctionTwoPoint_"/>
    !!]
    return
  end function cylindricalConstructorParameters

  function cylindricalConstructorInternal(radiusCylinderComoving,radiusBufferComoving,massHaloLens,redshiftLens,cosmologyFunctions_,powerSpectrum_,linearGrowth_,outputTimes_,darkMatterHaloBias_,darkMatterHaloScale_,correlationFunctionTwoPoint_,randomNumberGenerator_) result(self)
    !!{
    Internal constructor for the \refClass{geometryLightconeCylindrical} lightcone geometry distribution class.
    !!}
    use :: File_Utilities          , only : File_Exists, File_Lock          , File_Unlock  , lockDescriptor
    use :: Galacticus_Nodes        , only : treeNode   , nodeComponentBasic
    use :: Input_Paths             , only : inputPath  , pathTypeDataDynamic
    use :: HDF5_Access             , only : hdf5Access
    use :: IO_HDF5                 , only : hdf5Object
    use :: Linear_Algebra          , only : matrix     , vector             , assignment(=), operator(*)
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator , GSL_Integ_Gauss15
    use :: Hashes_Cryptographic    , only : Hash_MD5
    use :: MPI_Utilities           , only : mpiSelf
    implicit none
    type            (geometryLightconeCylindrical    )                              :: self
    class           (cosmologyFunctionsClass         ), target     , intent(in   )  :: cosmologyFunctions_
    class           (outputTimesClass                ), target     , intent(in   )  :: outputTimes_
    class           (randomNumberGeneratorClass      ), target     , intent(in   )  :: randomNumberGenerator_
    class           (powerSpectrumClass              ), target     , intent(in   )  :: powerSpectrum_
    class           (linearGrowthClass               ), target     , intent(in   )  :: linearGrowth_
    class           (darkMatterHaloBiasClass         ), target     , intent(in   )  :: darkMatterHaloBias_
    class           (darkMatterHaloScaleClass        ), target     , intent(in   )  :: darkMatterHaloScale_
    class           (correlationFunctionTwoPointClass), target     , intent(in   )  :: correlationFunctionTwoPoint_
    double precision                                               , intent(in   )  :: radiusCylinderComoving       , radiusBufferComoving       , &
         &                                                                             massHaloLens                 , redshiftLens
    double precision                                  , allocatable, dimension(:  ) :: deviates                     , densityContrastLogarithmic , &
         &                                                                             logNormalSigma               , logNormalMu
    double precision                                  , allocatable, dimension(:,:) :: covariance
     type            (treeNode                       ), pointer                     :: node
     class           (nodeComponentBasic             ), pointer                     :: basic
    double precision                                  , parameter                   :: wavenumberMaximumFactor=1.0d2
    integer          (c_size_t                       )                              :: output                       , output1                    , &
         &                                                                             output2
    double precision                                                                :: timeMinimum                  , timeMaximum                , &
         &                                                                             heightRegion1Lower           , heightRegion1Upper         , &
         &                                                                             heightRegion2Lower           , heightRegion2Upper         , &         
         &                                                                             wavenumberVertical           , wavenumberMaximumRadial    , &
         &                                                                             wavenumberMinimumVertical    , wavenumberMaximumVertical  , & 
         &                                                                             time
    type            (integrator                      )                              :: integratorVertical           , integratorRadial
    type            (matrix                          )                              :: covarianceMatrix
    type            (vector                          )                              :: deviateVector
    type            (hdf5Object                      )                              :: file
    type            (lockDescriptor                  )                              :: fileLock
    type            (varying_string                  )                              :: fileName
    character       (len=18                          )                              :: label
    !![
    <constructorAssign variables="radiusCylinderComoving, radiusBufferComoving, massHaloLens, redshiftLens, *cosmologyFunctions_, *powerSpectrum_, *linearGrowth_, *outputTimes_, *darkMatterHaloBias_, *darkMatterHaloScale_, *correlationFunctionTwoPoint_, *randomNumberGenerator_"/>
    !!]
    
    ! Find the minimum and maximum distance associated with each output time.
    allocate(self%distanceMinimum(self%outputTimes_%count()))
    allocate(self%distanceMaximum(self%outputTimes_%count()))
    allocate(self%volume         (self%outputTimes_%count()))
    allocate(self%densityContrast(self%outputTimes_%count()))
    do output=1,self%outputTimes_%count()
       if (output == 1                        ) then
          timeMinimum=                                      self%outputTimes_%time(output)
       else
          timeMinimum=sqrt(self%outputTimes_%time(output-1)*self%outputTimes_%time(output))
       end if
       if (output == self%outputTimes_%count()) then
          timeMaximum=                                      self%outputTimes_%time(output)
       else
          timeMaximum=sqrt(self%outputTimes_%time(output+1)*self%outputTimes_%time(output))
       end if
       self%distanceMinimum(output)=self%cosmologyFunctions_%distanceComoving(timeMaximum)
       self%distanceMaximum(output)=self%cosmologyFunctions_%distanceComoving(timeMinimum)
    end do
    ! Compute the comoving volume associated with each output time, including the buffer region around our cylinder.
    self%volume=+Pi                            &
         &      *(                             &
         &        +self%radiusCylinderComoving &
         &        +self%radiusBufferComoving   &
         &       )**2                          &
         &      *(                             &
         &        +self%distanceMaximum        &
         &        -self%distanceMinimum        &
         &      )
    ! Compute cosmic variance within each cylinder slice. Note that we construct a descriptor for the file name directly - we do
    ! not use the automatically generated one as it will include the random number generator hash which does not affect the
    ! covariance.
    write (label,'(e17.10)') self%radiusCylinderComoving
    fileName=         inputPath                                (pathTypeDataDynamic                                           )// &
         &            'largeScaleStructure/'                                                                                   // &
         &            self%objectType                          (                                                              )// &
         &            'CosmicVariance_'                                                                                        // &
         &   Hash_MD5(                                                                                                            &
         &            self%cosmologyFunctions_%hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)// &
         &            self%outputTimes_       %hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)// &
         &            self%powerSpectrum_     %hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)// &
         &            self%linearGrowth_      %hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)// &
         &            trim(adjustl(label))                                                                                        &
         &           )                                                                                                         // &
         &            '.hdf5'
    allocate(covariance(self%outputTimes_%count(),self%outputTimes_%count()))
    !! Read the covariance matrix from file if possible.
    if (File_Exists(fileName)) then
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
       !$ call hdf5Access%set()
       file=hdf5Object(char(fileName))
       call file%readDataset('covariance',covariance)
       !$ call hdf5Access%unset()
       call File_Unlock(fileLock)
    else
       !! Compute the covariance matrix.
       integratorVertical     = integrator(cosmicVarianceIntegrandVertical,integrationRule=GSL_Integ_Gauss15,toleranceRelative=1.0d-2,toleranceAbsolute=1.0d-4)
       integratorRadial       = integrator(cosmicVarianceIntegrandRadial  ,integrationRule=GSL_Integ_Gauss15,toleranceRelative=1.0d-2                         )
       wavenumberMaximumRadial=+wavenumberMaximumFactor &
            &                  /radiusCylinderComoving
       ! Set the time to the present epoch - we will apply our own linear growth to the power spectrum below.
       time                   = self%cosmologyFunctions_%cosmicTime(1.0d0)
       do output1=1,self%outputTimes_%count()
          heightRegion1Lower                    =self%distanceMinimum(output1)
          heightRegion1Upper                    =self%distanceMaximum(output1)
          do output2=output1,self%outputTimes_%count()
             heightRegion2Lower                        =self%distanceMinimum(output2)
             heightRegion2Upper                        =self%distanceMaximum(output2)
             wavenumberMinimumVertical                 =1.0d0/wavenumberMaximumFactor/max((heightRegion1Upper-heightRegion1Lower),(heightRegion2Upper-heightRegion2Lower))
             wavenumberMaximumVertical                 =1.0d0*wavenumberMaximumFactor/min((heightRegion1Upper-heightRegion1Lower),(heightRegion2Upper-heightRegion2Lower))
             covariance               (output1,output2)=+2.0d0                                                                                       &
                  &                                     *integratorVertical%integrate(log(wavenumberMinimumVertical),log(wavenumberMaximumVertical)) &
                  &                                     *self%linearGrowth_%value(self%outputTimes_%time(output1))                                   &
                  &                                     *self%linearGrowth_%value(self%outputTimes_%time(output2))
             covariance               (output2,output1)=+covariance(output1,output2          )
          end do
       end do
       !! Store the covariance matrix to file.
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
       !$ call hdf5Access%set()
       file=hdf5Object(char(fileName))
       call file%writeDataset(covariance,'covariance')
       !$ call hdf5Access%unset()
       call File_Unlock(fileLock)
    end if
    ! Generate a realization of the cosmic density field.
    !! Realization is generated only on the master process so that we guarantee that it is the same across all processes.
    if (mpiSelf%isMaster()) then
       ! We assume that the density field follows a log-normal distribution. Note
       ! that we don't explicitly transform the covariance (computed in terms of Δ=ρ/ρ̅) to that for logΔ because (since the
       ! unperturbed Δ̅=1) the transformed covariance is just equal to the original covariance.
       allocate(deviates      (self%outputTimes_%count()))    
       allocate(logNormalSigma(self%outputTimes_%count()))    
       allocate(logNormalMu   (self%outputTimes_%count()))    
       do output=1_c_size_t,self%outputTimes_%count()
          !! Generate a deviate.
          deviates      (output)=self%randomNumberGenerator_%standardNormalSample()
          !! Compute the parameters of the log-normal distribution.
          logNormalSigma(output)=+sqrt(log(1.0d0+covariance(output,output)))
          logNormalMu   (output)=-logNormalSigma(output)**2 &
               &                 /2.0d0
       end do
       deviateVector=vector(deviates)
       !! Cholesky decompose the covariance matrix.
       covarianceMatrix=matrix(covariance)
       call covarianceMatrix%choleskyDecomposition()
       !! Generate the realization.
       allocate(densityContrastLogarithmic(self%outputTimes_%count()))
       densityContrastLogarithmic=covarianceMatrix*deviateVector
       self%densityContrast      =exp(densityContrastLogarithmic)
    else
       self%densityContrast      =0.0d0
    end if
#ifdef USEMPI
    self%densityContrast=mpiSelf%sum(self%densityContrast)
#endif
    ! Initialize the unique ID of the active node to an impossible value.
    self%activeCenterUniqueID=-1_kind_int8
    self%activeUniqueID      =-1_kind_int8
    self%activeOutput        =-1_c_size_t
    ! Determine properties of the lens halo.
    if (self%massHaloLens > 0.0d0) then
       allocate(node)
       node                           =>  treeNode                                     (                                                                       )
       basic                          =>  node    %basic                               (autoCreate=.true.                                                      )
       self %timeLens                 =   self    %cosmologyFunctions_%cosmicTime      (self%cosmologyFunctions_%expansionFactorFromRedshift(     redshiftLens))
       self %distanceComovingLens     =   self    %cosmologyFunctions_%distanceComoving(                                                     self%timeLens     )
       call basic%massSet            (self%massHaloLens)
       call basic%timeSet            (self%timeLens    )
       call basic%timeLastIsolatedSet(self%timeLens    )
       self %radiusVirialComovingLens =  +self    %darkMatterHaloScale_%radiusVirial   (                                                     node              ) &
            &                            /self    %cosmologyFunctions_ %expansionFactor(                                                     self%timeLens     )
       self %correlationLensMaximum   =  +self%correlationFunctionTwoPoint_%correlation(self%radiusVirialComovingLens,self%timeLens) &
            &                            *self%darkMatterHaloBias_         %bias       (node                                       )
       call node%destroy()
       deallocate(node)
    end if
    return

  contains

    double precision function cosmicVarianceIntegrandVertical(wavenumberVerticalLogarithmic)
      !!{
      Vertical integrand used in evaluating the cosmic variance.
      !!}
      use :: Bessel_Functions, only : Bessel_Function_J1_Zero
      implicit none
      double precision, intent(in   ) :: wavenumberVerticalLogarithmic
      double precision                :: wavenumberLower              , wavenumberUpper
      integer                         :: besselZero

      wavenumberVertical             =exp(wavenumberVerticalLogarithmic)
      besselZero                     =0
      wavenumberUpper                =0.0d0
      cosmicVarianceIntegrandVertical=0.0d0
      do while (wavenumberUpper < wavenumberMaximumRadial)
         besselZero                     =+besselZero                                                                              &
              &                          +1
         wavenumberLower                =                                                                wavenumberUpper
         wavenumberUpper                = min(Bessel_Function_J1_Zero(besselZero)/radiusCylinderComoving,wavenumberMaximumRadial)
         cosmicVarianceIntegrandVertical=+cosmicVarianceIntegrandVertical                                                         &
              &                          +integratorRadial%integrate(wavenumberLower,wavenumberUpper)
      end do
      cosmicVarianceIntegrandVertical=+cosmicVarianceIntegrandVertical                                                               &
           &                          *wavenumberVertical                                                                            &
           &                          *real(                                                                                         &
           &                                +      windowFunctionVertical(wavenumberVertical,heightRegion1Lower,heightRegion1Upper)  &
           &                                *conjg(windowFunctionVertical(wavenumberVertical,heightRegion2Lower,heightRegion2Upper)) &
           &                               )
      return
    end function cosmicVarianceIntegrandVertical

    double precision function cosmicVarianceIntegrandRadial(wavenumberRadial)
      !!{
      Radial integrand used in evaluating the cosmic variance.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: wavenumberRadial
      double precision                :: wavenumber

      wavenumber                   =+sqrt(                       &
           &                              +wavenumberVertical**2 &
           &                              +wavenumberRadial  **2 &
           &                             )
      cosmicVarianceIntegrandRadial=+self%powerSpectrum_%power(wavenumber,time)                       &
           &                        *windowFunctionRadial(wavenumberRadial,radiusCylinderComoving)**2 &
           &                        * 2.0d0*Pi                                                        &
           &                        *wavenumberRadial                                                 &
           &                        /(2.0d0*Pi)**3
      return
    end function cosmicVarianceIntegrandRadial

    double complex function windowFunctionVertical(wavenumberVertical,heightLower,heightUpper)
      !!{
      Window function used in evaluating the cosmic variance.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: wavenumberVertical, heightLower, &
           &                             heightUpper

      windowFunctionVertical=+2.0d0                                                                  &
           &                 *      dcmplx(0.0d0,1.0d0)                                              &
           &                 *(                                                                      &
           &                   +exp(dcmplx(0.0d0,1.0d0)*wavenumberVertical*             heightLower) &
           &                   -exp(dcmplx(0.0d0,1.0d0)*wavenumberVertical* heightUpper            ) &
           &                 )                                                                       &
           &                 /                          wavenumberVertical/(heightUpper-heightLower) &
           &                 /sqrt(2.0d0*Pi)
      return
    end function windowFunctionVertical
    
    double precision function windowFunctionRadial(wavenumberRadial,radius)
      !!{
      Window function used in evaluating the cosmic variance.
      !!}
      use :: Bessel_Functions, only : Bessel_Function_J1
      implicit none
      double precision, intent(in   ) :: wavenumberRadial, radius

      if (wavenumberRadial == 0.0d0) then
         windowFunctionRadial=+0.5d0
      else
         windowFunctionRadial=+Bessel_Function_J1(wavenumberRadial*radius) &
              &               /                   wavenumberRadial/radius
      end if
      return
    end function windowFunctionRadial
    
  end function cylindricalConstructorInternal

  subroutine cylindricalDestructor(self)
    !!{
    Destructor for the \refClass{geometryLightconeCylindrical} lightcone geometry distribution class.
    !!}
    implicit none
    type(geometryLightconeCylindrical), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"         />
    <objectDestructor name="self%outputTimes_"                />
    <objectDestructor name="self%randomNumberGenerator_"      />
    <objectDestructor name="self%powerSpectrum_"              />
    <objectDestructor name="self%linearGrowth_"               />
    <objectDestructor name="self%darkMatterHaloBias_"         />
    <objectDestructor name="self%darkMatterHaloScale_"        />
    <objectDestructor name="self%correlationFunctionTwoPoint_"/>
    !!]
    return
  end subroutine cylindricalDestructor

  function cylindricalReplicationCount(self,node)
    !!{
    Determine the number of times {\normalfont \ttfamily node} appears in the lightcone.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    integer(c_size_t                    )                :: cylindricalReplicationCount
    class  (geometryLightconeCylindrical), intent(inout) :: self
    type   (treeNode                    ), intent(inout) :: node

    call self%sampleNode(node)
    cylindricalReplicationCount=self%activeCount
    return
  end function cylindricalReplicationCount
  
  double precision function cylindricalTimeMinimum(self)
    !!{
    Return the minimum time in the lightcone.
    !!}
    implicit none
    class(geometryLightconeCylindrical), intent(inout) :: self

    cylindricalTimeMinimum=self%outputTimes_%time(1_c_size_t)
    return
  end function cylindricalTimeMinimum

  double precision function cylindricalTimeMaximum(self)
    !!{
    Return the minimum time in the lightcone.
    !!}
    implicit none
    class(geometryLightconeCylindrical), intent(inout) :: self

    cylindricalTimeMaximum=self%outputTimes_%time(self%outputTimes_%count())
    return
  end function cylindricalTimeMaximum

  logical function cylindricalIsInLightcone(self,node,atPresentEpoch,radiusBuffer)
    !!{
    Determine if the given {\normalfont \ttfamily node} lies within the lightcone.
    !!}
    use :: Galacticus_Nodes    , only : nodeComponentBasic 
    use :: Numerical_Comparison, only : Values_Agree
    use :: Error               , only : Error_Report
    use :: String_Handling     , only : operator(//)
    implicit none
    class           (geometryLightconeCylindrical), intent(inout)            :: self
    type            (treeNode                    ), intent(inout)            :: node
    logical                                       , intent(in   ) , optional :: atPresentEpoch
    double precision                              , intent(in   ) , optional :: radiusBuffer
    class           (nodeComponentBasic          )                , pointer  :: basic 
    character       (len=10                      )                           :: label
    type            (varying_string              )                           :: message
    integer         (c_size_t                    )                           :: output
    !![
    <optionalArgument name="atPresentEpoch" defaultsTo=".true." />
    !!]

    ! Get the basic component.
    basic => node%basic()
    ! Assume not in the lightcone by default.
    cylindricalIsInLightcone=.false.
    ! Check if this node exists prior to any lightcone time. If it does it will not be output.
    if (basic%time() < self%outputTimes_%time(1_c_size_t)*(1.0d0-timeTolerance)) return
    ! Determine to which output this galaxy corresponds.
    output=self%outputTimes_%index(basic%time(),findClosest=.true.)
    if (atPresentEpoch_) then
       ! We want to check only the current time for this node. Check that the node exists precisely at a lightcone snapshot time,
       ! and then set the maximum output to check to equal to minimum, such that we test only the current time.
       if (.not.Values_Agree(self%outputTimes_%time(output),basic%time(),relTol=timeTolerance)) then
          message=         'failed to find matching time in lightcone'                      //char(10)
          write (label,'(f6.3)') basic%time()
          message=message//'               node time = '//trim(label)//' Gyr'               //char(10)
          write (label,'(f6.3)') self%outputTimes_%time(output)
          message=message//'  closest lightcone time = '//trim(label)//' Gyr ['//output//']'
          call Error_Report(message//{introspection:location})
       end if
       call self%sampleNode(node)
       cylindricalIsInLightcone=self%activeCount > 0_c_size_t
    else
       cylindricalIsInLightcone=.false.
       call Error_Report('not well-defined'//{introspection:location})  
    end if
    return
  end function cylindricalIsInLightcone

  double precision function cylindricalSolidAngle(self)
    !!{
    Return the solid angle (in steradians) of a cylindrical lightcone.
    !!}
    implicit none
    class(geometryLightconeCylindrical), intent(inout) :: self
    !$GLC attributes unused :: self

    ! Solid angle is not well-defined for this "lightcone" class. Simply return zero.
    cylindricalSolidAngle=0.0d0
    return
  end function cylindricalSolidAngle

  function cylindricalPosition(self,node,instance)
    !!{
    Return the position of the node in lightcone coordinates.
    !!}
    use            :: Error               , only : Error_Report
    use            :: Galacticus_Nodes    , only : nodeComponentBasic
    use, intrinsic :: ISO_C_Binding       , only : c_size_t
    use            :: ISO_Varying_String  , only : varying_string
    use            :: Numerical_Comparison, only : Values_Agree
    use            :: String_Handling     , only : operator(//)
    implicit none
    double precision                              , dimension(3)           :: cylindricalPosition
    class           (geometryLightconeCylindrical), intent(inout)          :: self
    type            (treeNode                    ), intent(inout), target  :: node
    integer         (c_size_t                    ), intent(in   )          :: instance
    class           (nodeComponentBasic          )               , pointer :: basic
    integer         (c_size_t                    )                         :: output
    character       (len=10                      )                         :: label
    type            (varying_string              )                         :: message

    ! Get the basic component.
    basic => node%basic()
    ! Determine to which output this node corresponds.
    output=self%outputTimes_%index(basic%time(),findClosest=.true.)
    if (.not.Values_Agree(self%outputTimes_%time(output),basic%time(),relTol=timeTolerance)) then
       message=         'failed to find matching time in lightcone'                      //char(10)
       write (label,'(f6.3)') basic%time()
       message=message//'               node time = '//trim(label)//' Gyr'               //char(10)
       write (label,'(f6.3)') self%outputTimes_%time(output)
       message=message//'  closest lightcone time = '//trim(label)//' Gyr ['//output//']'
       call Error_Report(message//{introspection:location})
    end if
    call self%sampleNode(node)
    cylindricalPosition=self%activePosition(:,self%activeInstances(instance))
    return
  end function cylindricalPosition

  function cylindricalVelocity(self,node,instance)
    !!{
    Return the velocity of the node in lightcone coordinates.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentPosition, treeNode
    implicit none
    double precision                              , dimension(3)  :: cylindricalvelocity
    class           (geometryLightconeCylindrical), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    integer         (c_size_t                    ), intent(in   ) :: instance
    !$GLC attributes unused :: self, node, instance

    ! Currently this class provides no model for velocities.
    cylindricalVelocity=0.0d0
    return
  end function cylindricalVelocity

  double precision function cylindricalTimeLightconeCrossing(self,node,timeStart,timeEnd,timesCrossing)
    !!{
    Return the time of the next lightcone crossing for this node.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (geometryLightconeCylindrical), intent(inout), target                              :: self
    type            (treeNode                    ), intent(inout), target                              :: node
    double precision                              , intent(in   )                                      :: timeStart    , timeEnd
    double precision                              , intent(inout), dimension(:), allocatable, optional :: timesCrossing
   !$GLC attributes unused :: self, node, timeStart, timeEnd, timesCrossing

    cylindricalTimeLightconeCrossing=0.0d0
    call Error_Report('not implemented'//{introspection:location})
    return
  end function cylindricalTimeLightconeCrossing
  
  function cylindricalPositionLightconeCrossing(self,node)
    !!{
    Return the position at the next lightcone crossing for this node.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision                              , dimension(3)  :: cylindricalPositionLightconeCrossing
    class           (geometryLightconeCylindrical), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    cylindricalPositionLightconeCrossing=0.0d0
    call Error_Report('not implemented'//{introspection:location})
    return
  end function cylindricalPositionLightconeCrossing
  
  subroutine cylindricalSampleNode(self,node)
    !!{
    Determine how many times the given node appears in the ``lightcone''.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Comparison    , only : Values_Agree
    use :: Error                   , only : Error_Report
    use :: Galacticus_Nodes        , only : nodeComponentBasic, nodeComponentSatellite
    use :: String_Handling         , only : operator(//)
    implicit none
    class           (geometryLightconeCylindrical), intent(inout)          :: self
    type            (treeNode                    ), intent(inout), target  :: node
    class           (nodeComponentBasic          )               , pointer :: basicHost
    class           (nodeComponentSatellite      )               , pointer :: satellite
    type            (treeNode                    )               , pointer :: nodeHost
    double precision                              , dimension(3)           :: positionOffset
    double precision                                                       :: distanceComoving  , radiusComoving   , &
         &                                                                    theta             , numberMean       , &
         &                                                                    separationLens    , probabilityAccept, &
         &                                                                    linearGrowthFactor, bias
    logical                                                                :: accept
    integer                                                                :: i                 , j
    character       (len=10                      )                         :: label
    type            (varying_string              )                         :: message
    integer         (c_size_t                    )                         :: output
    
    ! Find the isolated node.
    nodeHost => node
    do while (nodeHost%isSatellite())
       nodeHost => nodeHost%parent
    end do
    ! Get the basic component.
    basicHost => nodeHost%basic()
    ! Determine to which output this isolated node corresponds.
    output=self%outputTimes_%index(basicHost%time(),findClosest=.true.)
    if (.not.Values_Agree(self%outputTimes_%time(output),basicHost%time(),relTol=timeTolerance)) then
       message=         'failed to find matching time in lightcone'                      //char(10)
       write (label,'(f6.3)') basicHost%time()
       message=message//'               node time = '//trim(label)//' Gyr'               //char(10)
       write (label,'(f6.3)') self%outputTimes_%time(output)
       message=message//'  closest lightcone time = '//trim(label)//' Gyr ['//output//']'
       call Error_Report(message//{introspection:location})
    end if
    ! Determine if the output has changed.
    if (output /= self%activeOutput) then
       self%activeOutput        =output
       self%activeCenterUniqueID=-1_c_size_t
       self%activeUniqueID      =-1_c_size_t
    end if
    ! Check if we already have results for this isolated node stored.
    if (nodeHost%uniqueID() /= self%activeCenterUniqueID) then
       self%activeCenterUniqueID=nodeHost%uniqueID()
       self%activeUniqueID      =-1_c_size_t
       ! Compute the number of instances of this halo which appear in the cylinder.
       numberMean            =+self                           %volume         (output    ) &
            &                 *nodeHost%hostTree              %volumeWeight                &
            &                 *self                           %densityContrast(output    )
       self%activeCenterCount=+self    %randomNumberGenerator_%poissonSample  (numberMean)
       ! Generate random positions for each instance.
       if (allocated(self%activeCenterPosition)) deallocate(self%activeCenterPosition)
       if (allocated(self%activePosition      )) deallocate(self%activePosition      )
       if (allocated(self%activeInCylinder    )) deallocate(self%activeInCylinder    )
       if (allocated(self%activeInstances     )) deallocate(self%activeInstances     )
       allocate(self%activeCenterPosition(3,self%activeCenterCount))
       allocate(self%activePosition      (3,self%activeCenterCount))
       allocate(self%activeInCylinder    (  self%activeCenterCount))
       allocate(self%activeInstances     (  self%activeCenterCount))
       do i=1,self%activeCenterCount
          accept=.false.
          do while (.not.accept)
             radiusComoving  =sqrt(                                                                                  &
                  &                self%randomNumberGenerator_%uniformSample()*(                                     &
                  &                                                             +self%radiusCylinderComoving         &
                  &                                                             +self%radiusBufferComoving           &
                  &                                                            )**2                                  &
                  &               )
             theta           =     self%randomNumberGenerator_%uniformSample()*  2.0d0                               &
                  &                                                           *  Pi
             distanceComoving=     self%randomNumberGenerator_%uniformSample()*(                                     &
                  &                                                             +self%distanceMaximum       (output) &
                  &                                                             -self%distanceMinimum       (output) &
                  &                                                            )                                     &
                  &                                                             +self%distanceMinimum       (output)
             self%activeCenterPosition(1,i)=radiusComoving  *cos(theta)
             self%activeCenterPosition(2,i)=radiusComoving  *sin(theta)
             self%activeCenterPosition(3,i)=distanceComoving
             ! Decide whether to accept this based on the correlation with the lens.
             if (self%massHaloLens > 0.0d0) then
                separationLens    =+sqrt(radiusComoving**2+(distanceComoving-self%distanceComovingLens)**2)
                if (separationLens < self%radiusVirialComovingLens) then
                   accept=.false.
                else
                   linearGrowthFactor= self%linearGrowth_      %value(self%cosmologyFunctions_%timeAtDistanceComoving(sqrt(sum(self%activeCenterPosition(:,i)**2))))
                   bias              = self%darkMatterHaloBias_%bias (nodeHost                                                                                     )
                   probabilityAccept =+(1.0d0+self%correlationFunctionTwoPoint_%correlation(separationLens,self%timeLens)*linearGrowthFactor*bias) &
                        &             /(1.0d0+self%correlationLensMaximum                                                *linearGrowthFactor*bias)
                   accept            = self%randomNumberGenerator_%uniformSample() <= probabilityAccept
                end if
             else
                accept               =.true.
             end if
          end do
       end do
    end if
    ! Check if we already have results for this specific node stored.
    if (node%uniqueID() /= self%activeUniqueID) then
       self%activeUniqueID=node%uniqueID()
       ! Calculate offset from the isolated node center.
       positionOffset=0.0d0
       nodeHost => node
       do while (nodeHost%isSatellite())
          satellite      =>  nodeHost      %satellite                          (                )
          positionOffset =  +positionOffset                                                       &
               &            +satellite                         %position       (                ) &
               &            /self          %cosmologyFunctions_%expansionFactor(basicHost%time())
          nodeHost       =>  nodeHost      %parent
       end do
       do i=1,self%activeCenterCount          
          self%activePosition(:,i)=self%activeCenterPosition(:,i)+positionOffset
       end do
       self%activeInCylinder=sum  (self%activePosition  (1:2,:)**2,dim=1) <= self%radiusCylinderComoving**2
       self%activeCount     =count(self%activeInCylinder                )
       j                    =0
       do i=1,self%activeCenterCount
          if (self%activeInCylinder(i)) then
             j                      =j+1
             self%activeInstances(j)=i
          end if
       end do
    end if
    return
  end subroutine cylindricalSampleNode
