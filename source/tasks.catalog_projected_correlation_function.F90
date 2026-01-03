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

  use :: Cosmology_Functions     , only : cosmologyFunctions        , cosmologyFunctionsClass
  use :: Cosmology_Parameters    , only : cosmologyParameters       , cosmologyParametersClass
  use :: Geometry_Surveys        , only : surveyGeometry            , surveyGeometryClass
  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !![
  <enumeration>
   <name>randomSampleCountType</name>
   <description>Used to specify the type of random sample count.</description>
   <visibility>private</visibility>
   <entry label="fixed"        />
   <entry label="multiplicative"/>
  </enumeration>
  !!]
  
  !![
  <task name="taskCatalogProjectedCorrelationFunction">
   <description>A task which computes projected correlation functions based on a simple halo model approach.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskCatalogProjectedCorrelationFunction
     !!{
     Implementation of a task which computes projected correlation functions based on a simple halo model approach.
     !!}
     private
     class           (cosmologyParametersClass            ), pointer      :: cosmologyParameters_      => null()
     class           (cosmologyFunctionsClass             ), pointer      :: cosmologyFunctions_       => null()
     class           (surveyGeometryClass                 ), pointer      :: surveyGeometry_           => null()
     class           (randomNumberGeneratorClass          ), pointer      :: randomNumberGenerator_    => null()
     type            (varying_string                      )               :: galaxyCatalogFileName
     integer                                                              :: separationCount                     , randomSampleCountNumber
     type            (enumerationRandomSampleCountTypeType)               :: randomSampleCountType
     type            (varying_string                      )               :: randomSampleCount
     double precision                                                     :: massMinimum                         , massMaximum            , &
          &                                                                  separationMinimum                   , separationMaximum      , &
          &                                                                  separationRadialMaximum,              widthBuffer            , &
          &                                                                  angleRotation
     double precision                                      , dimension(3) :: origin
     double precision                                      , dimension(2) :: vectorRotation
     logical                                                              :: nodeComponentsInitialized =  .false., halfIntegral
     ! Pointer to the parameters for this task.
     type            (inputParameters                     )               :: parameters
   contains
     final     ::            catalogProjectedCorrelationFunctionDestructor
     procedure :: perform => catalogProjectedCorrelationFunctionPerform
  end type taskCatalogProjectedCorrelationFunction

  interface taskCatalogProjectedCorrelationFunction
     !!{
     Constructors for the \refClass{taskCatalogProjectedCorrelationFunction} task.
     !!}
     module procedure catalogProjectedCorrelationFunctionConstructorParameters
     module procedure catalogProjectedCorrelationFunctionConstructorInternal
  end interface taskCatalogProjectedCorrelationFunction

contains

  function catalogProjectedCorrelationFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskCatalogProjectedCorrelationFunction} task class which takes a parameter set as input.
    !!}
    use :: Galacticus_Nodes        , only : nodeClassHierarchyInitialize
    use :: Input_Parameters        , only : inputParameter              , inputParameters
    use :: Node_Components         , only : Node_Components_Initialize
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (taskCatalogProjectedCorrelationFunction)                        :: self
    type            (inputParameters                        ), intent(inout), target :: parameters
    class           (cosmologyFunctionsClass                ), pointer               :: cosmologyFunctions_
    class           (cosmologyParametersClass               ), pointer               :: cosmologyParameters_
    class           (surveyGeometryClass                    ), pointer               :: surveyGeometry_
    class           (randomNumberGeneratorClass             ), pointer               :: randomNumberGenerator_
    type            (inputParameters                        ), pointer               :: parametersRoot
    integer                                                                          :: separationCount        , randomSampleCountNumber
    type            (enumerationRandomSampleCountTypeType   )                        :: randomSampleCountType
    double precision                                                                 :: massMinimum            , massMaximum            , &
         &                                                                              separationMinimum      , separationMaximum      , &
         &                                                                              separationRadialMaximum, widthBuffer            , &
         &                                                                              angleRotation
    double precision                                         , dimension(3)          :: origin
    double precision                                         , dimension(2)          :: vectorRotation
    logical                                                                          :: halfIntegral
    type            (varying_string                         )                        :: galaxyCatalogFileName  , randomSampleCount
    character       (len=128                                )                        :: label

    ! Ensure the nodes objects are initialized.
    if (associated(parameters%parent)) then
       parametersRoot => parameters%parent
       do while (associated(parametersRoot%parent))
          parametersRoot => parametersRoot%parent
       end do
       call nodeClassHierarchyInitialize(parametersRoot)
       call Node_Components_Initialize  (parametersRoot)
    else
       parametersRoot => parameters
       call nodeClassHierarchyInitialize(parameters    )
       call Node_Components_Initialize  (parameters    )
    end if
    self%nodeComponentsInitialized=.true.
    !![
    <inputParameter>
      <name>galaxyCatalogFileName</name>
      <description>The file name from which the galaxy catalog should be read.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMinimum</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum mass galaxy to include in a mock catalog correlation function calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum mass galaxy to include in a mock catalog correlation function calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>separationMinimum</name>
      <defaultValue>0.1d0</defaultValue>
      <description>The minimum separation to compute in a mock catalog correlation function calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>separationMaximum</name>
      <defaultValue>30.0d0</defaultValue>
      <description>The maximum separation to compute in a mock catalog correlation function calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>separationCount</name>
      <defaultValue>15</defaultValue>
      <description>The number of bins in separation to compute in a mock catalog correlation function calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>randomSampleCount</name>
      <defaultValue>var_str('*10')</defaultValue>
      <description>The number of random points to use when constructing random catalogs. Can be either a fixed number or, if prefixed with ``{\normalfont \ttfamily *}'', a multiplicative factor.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    if (extract(randomSampleCount,1,1) == "*") then
       randomSampleCountType=randomSampleCountTypeMultiplicative
       label=char(extract(randomSampleCount,2))
       read (label,*) randomSampleCountNumber
    else
       randomSampleCountType=randomSampleCountTypeFixed
       label=char(randomSampleCount)
       read (label,*) randomSampleCountNumber
    end if
    !![
    <inputParameter>
      <name>separationRadialMaximum</name>
      <defaultValue>40.0d0</defaultValue>
      <description>The maximum radial separation of galaxies to consider when computing projected correlation functions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>halfIntegral</name>
      <defaultValue>.false.</defaultValue>
      <description>Set to {\normalfont \ttfamily true} if the projected correlation function is computed as $w_\mathrm{p}(r_\mathrm{p})=\int_0^{+\pi_\mathrm{max}} \xi(r_\mathrm{p},\pi) \mathrm{d} \pi$, instead of the usual $w_\mathrm{p}(r_\mathrm{p})=\int_{-\pi_\mathrm{max}}^{+\pi_\mathrm{max}} \xi(r_\mathrm{p},\pi) \mathrm{d} \pi$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>widthBuffer</name>
      <defaultValue>30.0d0</defaultValue>
      <description>The width of the buffer region around survey geometry to ensure galaxies are not lost when moving to redshift space.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="surveyGeometry"        name="surveyGeometry_"        source="parameters"/>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    <inputParameter>
      <name>origin</name>
      <defaultValue>[randomNumberGenerator_%uniformSample(),randomNumberGenerator_%uniformSample(),randomNumberGenerator_%uniformSample()]</defaultValue>
      <defaultSource>Uniformly random distribution within the box.</defaultSource>
      <description>The vector (in units of the box length) giving the origin of the coordinate system to use in mock catalog construction.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>vectorRotation</name>
      <defaultValue>[acos(2.0d0*randomNumberGenerator_%uniformSample()-1.0d0),2.0d0*Pi*randomNumberGenerator_%uniformSample()]</defaultValue>
      <defaultSource>Isotropically random on the unit sphere.</defaultSource>
      <description>The vector, in spherical coordinates $(\theta,\phi)$, about which the mock catalog should be rotated.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>angleRotation</name>
      <defaultValue>2.0d0*Pi*randomNumberGenerator_%uniformSample()</defaultValue>
      <defaultSource>Uniformly random distribution between $0$ and $2\pi$.</defaultSource>
      <description>The angle through which the mock catalog should be rotated.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=taskCatalogProjectedCorrelationFunction(galaxyCatalogFileName,massMinimum,massMaximum,separationCount,separationMinimum, separationMaximum, separationRadialMaximum,widthBuffer,origin,vectorRotation,angleRotation,randomSampleCountNumber,randomSampleCountType,halfIntegral,cosmologyParameters_,cosmologyFunctions_,surveyGeometry_,randomNumberGenerator_,parametersRoot)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="surveyGeometry_"       />
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function catalogProjectedCorrelationFunctionConstructorParameters

  function catalogProjectedCorrelationFunctionConstructorInternal(galaxyCatalogFileName,massMinimum,massMaximum,separationCount,separationMinimum, separationMaximum, separationRadialMaximum,widthBuffer,origin,vectorRotation,angleRotation,randomSampleCountNumber,randomSampleCountType,halfIntegral,cosmologyParameters_,cosmologyFunctions_,surveyGeometry_,randomNumberGenerator_,parameters) result(self)
    !!{
    Constructor for the \refClass{taskCatalogProjectedCorrelationFunction} task class which takes a parameter set as input.
    !!}
    use :: String_Handling, only : operator(//)
    implicit none
    type            (taskCatalogProjectedCorrelationFunction)                              :: self
    type            (varying_string                         ), intent(in   )               :: galaxyCatalogFileName
    integer                                                  , intent(in   )               :: separationCount        , randomSampleCountNumber
    type            (enumerationRandomSampleCountTypeType   ), intent(in   )               :: randomSampleCountType
    double precision                                         , intent(in   )               :: massMinimum            , massMaximum            , &
         &                                                                                    separationMinimum      , separationMaximum      , &
         &                                                                                    separationRadialMaximum, widthBuffer            , &
         &                                                                                    angleRotation
    double precision                                         , intent(in   ), dimension(3) :: origin
    double precision                                         , intent(in   ), dimension(2) :: vectorRotation
    logical                                                  , intent(in   )               :: halfIntegral
    class           (cosmologyParametersClass               ), intent(in   ), target       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                ), intent(in   ), target       :: cosmologyFunctions_
    class           (surveyGeometryClass                    ), intent(in   ), target       :: surveyGeometry_
    class           (randomNumberGeneratorClass             ), intent(in   ), target       :: randomNumberGenerator_
    type            (inputParameters                        ), intent(in   ), target :: parameters
    !![
    <constructorAssign variables="galaxyCatalogFileName,massMinimum,massMaximum,separationCount,separationMinimum, separationMaximum, separationRadialMaximum,widthBuffer,origin,vectorRotation,angleRotation,randomSampleCountNumber, randomSampleCountType, halfIntegral, *cosmologyParameters_, *cosmologyFunctions_, *surveyGeometry_, *randomNumberGenerator_"/>
    !!]

    if      (randomSampleCountType == randomSampleCountTypeMultiplicative) then
       self%randomSampleCount=var_str("*")//randomSampleCountNumber
    else if (randomSampleCountType == randomSampleCountTypeFixed         ) then
       self%randomSampleCount=var_str("" )//randomSampleCountNumber
    else
       call Error_Report('unknown randomSampleCountType'//{introspection:location})
    end if
    self%parameters=inputParameters(parameters)
    call self%parameters%parametersGroupCopy(parameters)
    return
  end function catalogProjectedCorrelationFunctionConstructorInternal

  subroutine catalogProjectedCorrelationFunctionDestructor(self)
    !!{
    Destructor for the \refClass{taskCatalogProjectedCorrelationFunction} task class.
    !!}
    use :: Node_Components, only : Node_Components_Uninitialize
    implicit none
    type(taskCatalogProjectedCorrelationFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%surveyGeometry_"       />
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    if (self%nodeComponentsInitialized) call Node_Components_Uninitialize()
    return
  end subroutine catalogProjectedCorrelationFunctionDestructor

  subroutine catalogProjectedCorrelationFunctionPerform(self,status)
    !!{
    Compute the projected correlation function from a galaxy catalog.
    !!}
    use :: Display                         , only : displayIndent                    , displayMessage                     , displayUnindent
    use :: Error                           , only : Error_Report                     , errorStatusSuccess
    use :: Output_HDF5                     , only : outputFile
    use :: IO_HDF5                         , only : hdf5Object
    use :: IO_IRATE                        , only : irate
    use :: ISO_Varying_String              , only : varying_string
    use :: Node_Components                 , only : Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Points                          , only : Points_Prune                     , Points_Replicate                   , Points_Rotate  , Points_Survey_Geometry, &
          &                                         Points_Translate
    use :: Statistics_Points_Correlations  , only : Statistics_Points_Correlation
    use :: String_Handling                 , only : operator(//)
    implicit none
    class           (taskCatalogProjectedCorrelationFunction), intent(inout), target         :: self
    integer                                                  , intent(  out), optional       :: status
    double precision                                         , pointer      , dimension(:,:) :: galaxyPosition_      , galaxyVelocity_
    double precision                                         , allocatable  , dimension(:,:) :: galaxyPosition       , galaxyVelocity
    double precision                                         , allocatable  , dimension(:,:) :: randomPosition
    double precision                                         , pointer      , dimension(:  ) :: galaxyMass
    double precision                                         , allocatable  , dimension(:  ) :: correlation          , separation              , &
         &                                                                                      correlationSurvey
    double precision                                                        , dimension(3  ) :: rotationAxis
    type            (varying_string                         )                                :: message
    type            (hdf5Object                             )                                :: dataset              , correlationFunctionGroup
    type            (irate                                  )                                :: galaxyFile
    double precision                                                                         :: simulationBoxSize    , time                    , &
         &                                                                                      redshift
    integer                                                                                  :: randomPointCount     , replications            , &
         &                                                                                      i                    , j                       , &
         &                                                                                      replicatedGalaxyCount

    call displayIndent('Begin task: catalog projected correlation function')
    ! Call routines to perform initialization which must occur for all threads if run in parallel.
    call Node_Components_Thread_Initialize(self%parameters)
    ! Read the galaxy catalog.
    call displayIndent("Reading galaxy catalog")
    galaxyFile=irate(char(self%galaxyCatalogFileName),self%cosmologyParameters_,self%cosmologyFunctions_)
    call galaxyFile%readHalos     (                            &
         &                         snapshot=1                , &
         &                         redshift=redshift         , &
         &                         center  =galaxyPosition_  , &
         &                         velocity=galaxyVelocity_  , &
         &                         mass    =galaxyMass         &
         &                         )
    call galaxyFile%readSimulation(                            &
         &                         boxSize =simulationBoxSize  &
         &                         )
    message="Read "
    message=message//size(galaxyPosition,dim=2)//" galaxies"
    call displayMessage(message)
    call displayUnindent("done")
    ! Copy data.
    allocate(galaxyPosition,source=galaxyPosition_)
    allocate(galaxyVelocity,source=galaxyVelocity_)
    deallocate(galaxyPosition_)
    deallocate(galaxyVelocity_)
    ! Get cosmic time.
    time=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshift))
    ! Prune points below our mass threshold.
    call Points_Prune(galaxyPosition,galaxyMass >= self%massMinimum .and. galaxyMass < self%massMaximum)
    call Points_Prune(galaxyVelocity,galaxyMass >= self%massMinimum .and. galaxyMass < self%massMaximum)
    message="Pruned on mass leaving "
    message=message//size(galaxyPosition,dim=2)//" galaxies"
    call displayMessage(message)
    ! Generate random points.
    select case (self%randomSampleCountType%ID)
    case (randomSampleCountTypeFixed         %ID)
       randomPointCount=int(self%randomSampleCountNumber)
    case (randomSampleCountTypeMultiplicative%ID)
       randomPointCount=int(self%randomSampleCountNumber*dble(size(galaxyPosition,dim=2)))
    case default
       randomPointCount=0
       call Error_Report('unknown random sample count type'//{introspection:location})
    end select
    allocate(randomPosition(3,randomPointCount))
    do i=1,3
       do j=1,randomPointCount
          randomPosition(i,j)=self%randomNumberGenerator_%uniformSample()*simulationBoxSize
       end do
    end do
    ! Shift origin to a random point in the box.
    call Points_Translate(galaxyPosition,-      simulationBoxSize*self%origin,periodicLength=simulationBoxSize)
    ! Shift points so origin is at the center of the box.
    call Points_Translate(galaxyPosition,-0.5d0*simulationBoxSize*[1.0d0,1.0d0,1.0d0])
    call Points_Translate(randomPosition,-0.5d0*simulationBoxSize*[1.0d0,1.0d0,1.0d0])
    ! Compute correlation function.
    call Statistics_Points_Correlation(                                                                       &
         &                             galaxyPosition                                                       , &
         &                             randomPosition                                                       , &
         &                             self%separationMinimum                                               , &
         &                             self%separationMaximum                                               , &
         &                             self%separationCount                                                 , &
         &                             separation                                                           , &
         &                             correlation                                                          , &
         &                             projected                               =.true.                      , &
         &                             radialSeparationMaximum                 =self%separationRadialMaximum  &
         &                            )
    ! Replicate points to encompass survey geometry.
    replications=int((self%surveyGeometry_%distanceMaximum(self%massMaximum)+self%widthBuffer)/simulationBoxSize+0.5d0)
    call Points_Replicate(galaxyPosition,simulationBoxSize,-replications*[1,1,1],+replications*[1,1,1])
    message="Replicated to cover survey volume giving "
    message=message//size(galaxyPosition,dim=2)//" galaxies"
    call displayMessage(message)
    replicatedGalaxyCount=size(galaxyPosition,dim=2)
    ! Shift points into redshift space.
    call pointsToRedshiftSpace()
    ! Rotate the points.
    rotationAxis =[                                                         &
         &         sin(self%vectorRotation(1))*cos(self%vectorRotation(2)), &
         &         sin(self%vectorRotation(1))*sin(self%vectorRotation(2)), &
         &         cos(self%vectorRotation(1))                              &
         &        ]
    call Points_Rotate(galaxyPosition,rotationAxis,self%angleRotation)
    ! Limit points to survey geometry.
    call Points_Survey_Geometry(galaxyPosition,self%surveyGeometry_,self%massMaximum)
    message="Pruned on survey geometry leaving "
    message=message//size(galaxyPosition,dim=2)//" galaxies"
    call displayMessage(message)
    ! Generate random points.
    select case (self%randomSampleCountType%ID)
    case (randomSampleCountTypeFixed         %ID)
       randomPointCount=int(self%randomSampleCountNumber*dble(replicatedGalaxyCount)/dble(size(galaxyPosition,dim=2)))
    case (randomSampleCountTypeMultiplicative%ID)
       randomPointCount=int(self%randomSampleCountNumber*dble(replicatedGalaxyCount))
    case default
       randomPointCount=0
       call Error_Report('unknown random sample count type'//{introspection:location})
    end select
    deallocate(randomPosition                     )
    allocate(randomPosition(3,randomPointCount))
    do i=1,3
       do j=1,randomPointCount
          randomPosition(i,j)=(self%randomNumberGenerator_%uniformSample()-0.5d0)*dble(2*replications+1)*simulationBoxSize
       end do
    end do
    message="Generated "
    message=message//size(randomPosition,dim=2)//" random points"
    call displayMessage(message)
    call Points_Survey_Geometry(randomPosition,self%surveyGeometry_,self%massMaximum)
    message="Pruned on survey geometry leaving "
    message=message//size(randomPosition,dim=2)//" random points"
    call displayMessage(message)
    ! Compute the correlation function.
    call Statistics_Points_Correlation(                                                                       &
         &                             galaxyPosition                                                       , &
         &                             randomPosition                                                       , &
         &                             self%separationMinimum                                               , &
         &                             self%separationMaximum                                               , &
         &                             self%separationCount                                                 , &
         &                             separation                                                           , &
         &                             correlationSurvey                                                    , &
         &                             projected                               =.true.                      , &
         &                             radialSeparationMaximum                 =self%separationRadialMaximum, &
         &                             halfIntegral                            =self%halfIntegral             &
         &                            )
    ! Write correlations to file.
    correlationFunctionGroup=outputFile%openGroup('projectedCorrelationFunction')
    call correlationFunctionGroup%writeDataset(separation       ,'separation'                ,comment='Galaxy separation.'                                ,datasetReturned=dataset)
    call dataset%writeAttribute(megaParsec,'unitsInSI')
    call dataset%close         (                      )
    call correlationFunctionGroup%writeDataset(correlation      ,'projectedCorrelation'      ,comment='Projected correlation function from the full mock.',datasetReturned=dataset)
    call dataset%writeAttribute(megaParsec,'unitsInSI')
    call dataset%close         (                      )
    call correlationFunctionGroup%writeDataset(correlationSurvey,'projectedCorrelationSurvey',comment='Projected correlation function from survey region.',datasetReturned=dataset)
    call dataset%writeAttribute(megaParsec,'unitsInSI')
    call dataset%close         (                      )
    call correlationFunctionGroup%close()
    ! Clean up.
    deallocate(galaxyPosition)
    deallocate(galaxyVelocity)
    deallocate(galaxyMass    )
    call Node_Components_Thread_Uninitialize()
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: catalog projected correlation function' )

  contains

    subroutine pointsToRedshiftSpace()
      !!{
      Shift the set of points into redshift space.
      !!}
      use :: Vectors, only : Vector_Magnitude
      implicit none
      integer          :: galaxyCount       , replicantCount, i, j, k, l
      double precision :: velocityToDistance

      velocityToDistance=(1.0d0+redshift)/self%cosmologyFunctions_%hubbleParameterEpochal(time)
      galaxyCount       =size(galaxyVelocity,dim=2)
      replicantCount    =0
      do i=-replications,+replications
         do j=-replications,+replications
            do k=-replications,+replications
               forall(l=1:galaxyCount)
                  galaxyPosition                         (:,replicantCount*galaxyCount+l)=   &
                       & +galaxyPosition                 (:,replicantCount*galaxyCount+l)    &
                       & *(                                                                  &
                       &   +1.0d0                                                            &
                       &   +Dot_Product     (                                                &
                       &                     galaxyVelocity(:,                           l), &
                       &                     galaxyPosition(:,replicantCount*galaxyCount+l)  &
                       &                    )                                                &
                       &   /Vector_Magnitude(                                                &
                       &                     galaxyPosition(:,replicantCount*galaxyCount+l)  &
                       &                    )**2                                             &
                       &   *velocityToDistance                                               &
                       &  )
               end forall
               replicantCount=replicantCount+1
            end do
         end do
      end do
      return
    end subroutine pointsToRedshiftSpace

  end subroutine catalogProjectedCorrelationFunctionPerform
