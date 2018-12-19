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
  
  use Cosmology_Functions , only : cosmologyFunctions , cosmologyFunctionsClass
  use Cosmology_Parameters, only : cosmologyParameters, cosmologyParametersClass
  use Geometry_Surveys    , only : surveyGeometry     , surveyGeometryClass

  !# <task name="taskCatalogProjectedCorrelationFunction" defaultThreadPrivate="yes">
  !#  <description>A task which computes projected correlation functions based on a simple halo model approach.</description>
  !# </task>
  type, extends(taskClass) :: taskCatalogProjectedCorrelationFunction
     !% Implementation of a task which computes projected correlation functions based on a simple halo model approach.
     private 
     class           (cosmologyParametersClass         ), pointer      :: cosmologyParameters_
     class           (cosmologyFunctionsClass          ), pointer      :: cosmologyFunctions_
     class           (surveyGeometryClass              ), pointer      :: surveyGeometry_
     type            (varying_string                   )               :: galaxyCatalogFileName
     integer                                                           :: separationCount        , randomSampleCount, &
          &                                                               randomSampleCountType
     double precision                                                  :: massMinimum            , massMaximum      , &
          &                                                               separationMinimum      , separationMaximum, &
          &                                                               separationRadialMaximum, widthBuffer      , &
          &                                                               angleRotation
     double precision                                   , dimension(3) :: origin
     double precision                                   , dimension(2) :: vectorRotation
     logical                                                           :: halfIntegral
   contains
     final     ::            catalogProjectedCorrelationFunctionDestructor
     procedure :: perform => catalogProjectedCorrelationFunctionPerform
  end type taskCatalogProjectedCorrelationFunction

  interface taskCatalogProjectedCorrelationFunction
     !% Constructors for the {\normalfont \ttfamily catalogProjectedCorrelationFunction} task.
     module procedure catalogProjectedCorrelationFunctionConstructorParameters
     module procedure catalogProjectedCorrelationFunctionConstructorInternal
  end interface taskCatalogProjectedCorrelationFunction

  !# <enumeration>
  !#  <name>randomSampleCountType</name>
  !#  <description>Used to specify the type of random sample count.</description>
  !#  <visibility>private</visibility>
  !#  <entry label="fixed"        />
  !#  <entry label="multiplicative"/>
  !# </enumeration>
  
contains
  
  function catalogProjectedCorrelationFunctionConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily catalogProjectedCorrelationFunction} task class which takes a parameter set as input.
    use Input_Parameters        , only : inputParameters             , inputParameter
    use Galacticus_Nodes        , only : nodeClassHierarchyInitialize
    use Node_Components         , only : Node_Components_Initialize  , Node_Components_Thread_Initialize
    use Pseudo_Random           , only : pseudoRandom
    use Numerical_Constants_Math, only : Pi
    implicit none
    type            (taskCatalogProjectedCorrelationFunction)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass               ), pointer       :: cosmologyParameters_
    class           (surveyGeometryClass                    ), pointer       :: surveyGeometry_
    integer                                                                  :: separationCount        , randomSampleCount    , &
         &                                                                      randomSampleCountType
    double precision                                                         :: massMinimum            , massMaximum          , &
         &                                                                      separationMinimum      , separationMaximum    , &
         &                                                                      separationRadialMaximum, widthBuffer          , &
         &                                                                      angleRotation
    double precision                                         , dimension(3)  :: origin
    double precision                                         , dimension(2)  :: vectorRotation
    logical                                                                  :: halfIntegral
    type            (varying_string                         )                :: galaxyCatalogFileName  , randomSampleCountText
    type            (pseudoRandom                           )                :: randomSequence
    character       (len=128                                )                :: label

    call nodeClassHierarchyInitialize     (parameters)
    call Node_Components_Initialize       (parameters)
    call Node_Components_Thread_Initialize(parameters)
    !# <inputParameter>
    !#   <name>galaxyCatalogFileName</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The file name from which the galaxy catalog should be read.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The minimum mass galaxy to include in a mock catalog correlation function calculation.</description>
    !#   <source>parameters</source>
    !#   <type>float</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d16</defaultValue>
    !#   <description>The maximum mass galaxy to include in a mock catalog correlation function calculation.</description>
    !#   <source>parameters</source>
    !#   <type>float</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>separationMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.1d0</defaultValue>
    !#   <description>The minimum separation to compute in a mock catalog correlation function calculation.</description>
    !#   <source>parameters</source>
    !#   <type>float</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>separationMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>30.0d0</defaultValue>
    !#   <description>The maximum separation to compute in a mock catalog correlation function calculation.</description>
    !#   <source>parameters</source>
    !#   <type>float</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>separationCount</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>15</defaultValue>
    !#   <description>The number of bins in separation to compute in a mock catalog correlation function calculation.</description>
    !#   <source>parameters</source>
    !#   <type>float</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomSampleCount</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>var_str('*10')</defaultValue>
    !# <variable>randomSampleCountText</variable>
    !#   <description>The number of random points to use when constructing random catalogs. Can be either a fixed number or, if prefixed with ``{\normalfont \ttfamily *}'', a multiplicative factor.</description>
    !#   <source>parameters</source>
    !#   <type>float</type>
    !# </inputParameter>
    if (extract(randomSampleCountText,1,1) == "*") then
       randomSampleCountType=randomSampleCountTypeMultiplicative
       label=char(extract(randomSampleCountText,2))
       read (label,*) randomSampleCount
    else
       randomSampleCountType=randomSampleCountTypeFixed
       label=char(randomSampleCountText)
       read (label,*) randomSampleCount
    end if
    !# <inputParameter>
    !#   <name>separationRadialMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>40.0d0</defaultValue>
    !#   <description>The maximum radial separation of galaxies to consider when computing projected correlation functions.</description>
    !#   <source>parameters</source>
    !#   <type>float</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>halfIntegral</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Set to {\normalfont \ttfamily true} if the projected correlation function is computed as $w_\mathrm{p}(r_\mathrm{p})=\int_0^{+\pi_\mathrm{max}} \xi(r_\mathrm{p},\pi) \mathrm{d} \pi$, instead of the usual $w_\mathrm{p}(r_\mathrm{p})=\int_{-\pi_\mathrm{max}}^{+\pi_\mathrm{max}} \xi(r_\mathrm{p},\pi) \mathrm{d} \pi$.</description>
    !#   <source>parameters</source>
    !#   <type>float</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>widthBuffer</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>30.0d0</defaultValue>
    !#   <description>The width of the buffer region around survey geometry to ensure galaxies are not lost when moving to redshift space.</description>
    !#   <source>parameters</source>
    !#   <type>float</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>origin</name>
    !#   <defaultValue>[randomSequence%uniformSample(),randomSequence%uniformSample(),randomSequence%uniformSample()]</defaultValue>
    !#   <defaultSource>Uniformly random distribution within the box.</defaultSource>
    !#   <description>The vector (in units of the box length) giving the origin of the coordinate system to use in mock catalog construction.</description>
    !#   <type>float</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>vectorRotation</name>
    !#   <defaultValue>[acos(2.0d0*randomSequence%uniformSample()-1.0d0),2.0d0*Pi*randomSequence%uniformSample()]</defaultValue>
    !#   <defaultSource>Isotropically random on the unit sphere.</defaultSource>
    !#   <description>The vector, in spherical coordinates $(\theta,\phi)$, about which the mock catalog should be rotated.</description>
    !#   <type>float</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>angleRotation</name>
    !#   <defaultValue>2.0d0*Pi*randomSequence%uniformSample()</defaultValue>
    !#   <defaultSource>Uniformly random distribution between $0$ and $2\pi$.</defaultSource>
    !#   <description>The angle through which the mock catalog should be rotated.</description>
    !#   <type>float</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !# <objectBuilder class="surveyGeometry"      name="surveyGeometry_"      source="parameters"/>
    self=taskCatalogProjectedCorrelationFunction(galaxyCatalogFileName,massMinimum,massMaximum,separationCount,separationMinimum, separationMaximum, separationRadialMaximum,widthBuffer,origin,vectorRotation,angleRotation,randomSampleCount,randomSampleCountType,halfIntegral,cosmologyParameters_,cosmologyFunctions_,surveyGeometry_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function catalogProjectedCorrelationFunctionConstructorParameters
  
  function catalogProjectedCorrelationFunctionConstructorInternal(galaxyCatalogFileName,massMinimum,massMaximum,separationCount,separationMinimum, separationMaximum, separationRadialMaximum,widthBuffer,origin,vectorRotation,angleRotation,randomSampleCount,randomSampleCountType,halfIntegral,cosmologyParameters_,cosmologyFunctions_,surveyGeometry_) result(self)
    !% Constructor for the {\normalfont \ttfamily catalogProjectedCorrelationFunction} task class which takes a parameter set as input.
    implicit none
    type            (taskCatalogProjectedCorrelationFunction)                              :: self
    type            (varying_string                         ), intent(in   )               :: galaxyCatalogFileName
    integer                                                  , intent(in   )               :: separationCount        , randomSampleCount, &
         &                                                                                    randomSampleCountType
    double precision                                         , intent(in   )               :: massMinimum            , massMaximum      , &
         &                                                                                    separationMinimum      , separationMaximum, &
         &                                                                                    separationRadialMaximum, widthBuffer      , &
         &                                                                                    angleRotation
    double precision                                         , intent(in   ), dimension(3) :: origin
    double precision                                         , intent(in   ), dimension(2) :: vectorRotation
    logical                                                  , intent(in   )               :: halfIntegral
    class           (cosmologyParametersClass               ), intent(in   ), target       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                ), intent(in   ), target       :: cosmologyFunctions_
    class           (surveyGeometryClass                    ), intent(in   ), target       :: surveyGeometry_
    !# <constructorAssign variables="galaxyCatalogFileName,massMinimum,massMaximum,separationCount,separationMinimum, separationMaximum, separationRadialMaximum,widthBuffer,origin,vectorRotation,angleRotation,randomSampleCount, randomSampleCountType, halfIntegral, *cosmologyParameters_, *cosmologyFunctions_, *surveyGeometry_"/>

    return
  end function catalogProjectedCorrelationFunctionConstructorInternal

  subroutine catalogProjectedCorrelationFunctionDestructor(self)
    !% Destructor for the {\normalfont \ttfamily catalogProjectedCorrelationFunction} task class.
    implicit none
    type(taskCatalogProjectedCorrelationFunction), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%cosmologyFunctions_" />
    !# <objectDestructor name="self%surveyGeometry_"     />
    return
  end subroutine catalogProjectedCorrelationFunctionDestructor

  subroutine catalogProjectedCorrelationFunctionPerform(self)
    !% Compute the projected correlation function from a galaxy catalog.
    use Galacticus_Error
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Display
    use Memory_Management
    use IO_HDF5
    use IO_IRATE
    use String_Handling
    use Statistics_Points_Correlations
    use Pseudo_Random
    use Points
    use Numerical_Constants_Astronomical
    use Galacticus_HDF5
    implicit none
    class           (taskCatalogProjectedCorrelationFunction), intent(inout)                 :: self
    double precision                                         , allocatable  , dimension(:,:) :: galaxyPosition       , galaxyVelocity          , &
         &                                                                                      randomPosition
    double precision                                         , allocatable  , dimension(:  ) :: correlation          , separation              , &
         &                                                                                      galaxyMass           , correlationSurvey
    double precision                                                        , dimension(3  ) :: rotationAxis
    type            (varying_string                         )                                :: message
    type            (hdf5Object                             )                                :: thisDataset          , correlationFunctionGroup
    type            (irate                                  )                                :: galaxyFile
    double precision                                                                         :: simulationBoxSize    , time                    , &
         &                                                                                      redshift
    integer                                                                                  :: randomPointCount     , replications            , &
         &                                                                                      i                    , j                       , &
         &                                                                                      replicatedGalaxyCount
    type            (pseudoRandom                           )                                :: randomSequence
    
    call Galacticus_Display_Indent('Begin task: catalog projected correlation function')
    ! Read the galaxy catalog.
    call Galacticus_Display_Indent("Reading galaxy catalog")
    galaxyFile=irate(char(self%galaxyCatalogFileName),self%cosmologyParameters_,self%cosmologyFunctions_)
    call galaxyFile%readHalos     (                            &
         &                         snapshot=1                , &
         &                         redshift=redshift         , &
         &                         center  =galaxyPosition   , &
         &                         velocity=galaxyVelocity   , &
         &                         mass    =galaxyMass         &
         &                         )
    call galaxyFile%readSimulation(                            &
         &                         boxSize =simulationBoxSize  &
         &                         )
    message="Read "
    message=message//size(galaxyPosition,dim=2)//" galaxies"
    call Galacticus_Display_Message(message)
    call Galacticus_Display_Unindent("done")
    ! Get cosmic time.
    time=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshift))
    ! Prune points below our mass threshold.
    call Points_Prune(galaxyPosition,galaxyMass >= self%massMinimum .and. galaxyMass < self%massMaximum)
    call Points_Prune(galaxyVelocity,galaxyMass >= self%massMinimum .and. galaxyMass < self%massMaximum)
    message="Pruned on mass leaving "
    message=message//size(galaxyPosition,dim=2)//" galaxies"
    call Galacticus_Display_Message(message)
    ! Generate random points.
    select case (self%randomSampleCountType)
    case (randomSampleCountTypeFixed         )
       randomPointCount=int(self%randomSampleCount)
    case (randomSampleCountTypeMultiplicative)
       randomPointCount=int(self%randomSampleCount*dble(size(galaxyPosition,dim=2)))
    case default
       randomPointCount=0
       call Galacticus_Error_Report('unknown random sample count type'//{introspection:location})
    end select
    call allocateArray(randomPosition,[3,randomPointCount])
    do i=1,3
       do j=1,randomPointCount
          randomPosition(i,j)=randomSequence%uniformSample()*simulationBoxSize
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
    call Galacticus_Display_Message(message)
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
    call Galacticus_Display_Message(message)
    ! Generate random points.
    select case (self%randomSampleCountType)
    case (randomSampleCountTypeFixed         )
       randomPointCount=int(self%randomSampleCount*dble(replicatedGalaxyCount)/dble(size(galaxyPosition,dim=2)))
    case (randomSampleCountTypeMultiplicative)
       randomPointCount=int(self%randomSampleCount*dble(replicatedGalaxyCount))
    case default
       randomPointCount=0
       call Galacticus_Error_Report('unknown random sample count type'//{introspection:location})
    end select
    call deallocateArray(randomPosition                     )
    call allocateArray  (randomPosition,[3,randomPointCount])
    do i=1,3
       do j=1,randomPointCount
          randomPosition(i,j)=(randomSequence%uniformSample()-0.5d0)*dble(2*replications+1)*simulationBoxSize
       end do
    end do
    message="Generated "
    message=message//size(randomPosition,dim=2)//" random points"
    call Galacticus_Display_Message(message)
    call Points_Survey_Geometry(randomPosition,self%surveyGeometry_,self%massMaximum)
    message="Pruned on survey geometry leaving "
    message=message//size(randomPosition,dim=2)//" random points"
    call Galacticus_Display_Message(message)
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
    correlationFunctionGroup=galacticusOutputFile%openGroup('projectedCorrelationFunction')
    call correlationFunctionGroup%writeDataset(separation       ,'separation'                ,commentText='Galaxy separation.'                                ,datasetReturned=thisDataset)
    call thisDataset%writeAttribute(megaParsec,'unitsInSI')
    call thisDataset%close         (                      )
    call correlationFunctionGroup%writeDataset(correlation      ,'projectedCorrelation'      ,commentText='Projected correlation function from the full mock.',datasetReturned=thisDataset)
    call thisDataset%writeAttribute(megaParsec,'unitsInSI')
    call thisDataset%close         (                      )
    call correlationFunctionGroup%writeDataset(correlationSurvey,'projectedCorrelationSurvey',commentText='Projected correlation function from survey region.',datasetReturned=thisDataset)
    call thisDataset%writeAttribute(megaParsec,'unitsInSI')
    call thisDataset%close         (                      )
    call correlationFunctionGroup%close()
    call Galacticus_Display_Unindent('Done task: catalog projected correlation function' )

  contains

    subroutine pointsToRedshiftSpace()
      !% Shift the set of points into redshift space.
      use Vectors
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
