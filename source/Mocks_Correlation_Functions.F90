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

!% Generate a mock realization of a halo model.

program Mocks_Correlation_Functions
  !% Generates a mock realization from an input halo catalog and a halo model prescription.
  use Galacticus_Error
  use Command_Arguments
  use ISO_Varying_String
  use Input_Parameters
  use Galacticus_Display
  use Geometry_Surveys
  use Memory_Management
  use IO_HDF5
  use IO_IRATE
  use String_Handling
  use Statistics_Points_Correlations
  use Pseudo_Random
  use Points
  use Cosmology_Functions
  use Numerical_Constants_Math
  use Numerical_Constants_Astronomical
  implicit none
  double precision                         , allocatable, dimension(:,:) :: galaxyPosition                          , galaxyVelocity                                , &
       &                                                                    randomPosition
  double precision                         , allocatable, dimension(:  ) :: correlation                             , separation                                    , &
       &                                                                    galaxyMass                              , correlationSurvey
  double precision                                      , dimension(3  ) :: rotationAxis                            , mockCorrelationFunctionOrigin
  double precision                                      , dimension(2  ) :: mockCorrelationFunctionRotationVector
  class           (surveyGeometryClass    ), pointer                     :: surveyGeometry_
  class           (cosmologyFunctionsClass), pointer                     :: cosmologyFunctions_
  integer                                  , parameter                   :: mockCorrelationFunctionRandomSampleCountMultiplicative=0
  integer                                  , parameter                   :: mockCorrelationFunctionRandomSampleCountFixed         =1
  type            (varying_string         )                              :: parameterFileName                       , galaxyCatalogFileName                         , &
       &                                                                    correlationFunctionFileName             , message                                       , &
       &                                                                    mockCorrelationFunctionRandomSampleCount
  type            (hdf5Object             )                              :: thisDataset                             , correlationFunctionFile
  type            (irate                  )                              :: galaxyFile
  double precision                                                       :: simulationBoxSize                       , mockCorrelationFunctionRadialSeparationMaximum, &
       &                                                                    mockCorrelationFunctionMassMinimum      , mockCorrelationFunctionMassMaximum            , &
       &                                                                    mockCorrelationFunctionSeparationMinimum, mockCorrelationFunctionSeparationMaximum      , &
       &                                                                    redshift                                , time                                          , &
       &                                                                    mockCorrelationFunctionBufferWidth      , mockCorrelationFunctionRotationAngle          , &
       &                                                                    mockCorrelationFunctionRandomSampleCountValue
  integer                                                                :: mockCorrelationFunctionSeparationCount  , randomPointCount                              , &
       &                                                                    i                                       , j                                             , &
       &                                                                    mockCorrelationFunctionRandomSampleCountType, replicatedGalaxyCount                     , &
       &                                                                    replications
  type            (pseudoRandom           )                              :: randomSequence
  character       (len=128                )                              :: label
  type            (inputParameters        )                              :: parameters
  logical                                                                :: mockCorrelationFunctionHalfIntegral

  ! Read in basic code memory usage.
  call Code_Memory_Usage('Mocks_Correlation_Functions.size')
  ! Check that correct number of arguments have been supplied.
  if (Command_Argument_Count() /= 3) call Galacticus_Error_Report(message="Usage: Mocks_Correlation_Functions.exe <parameterFileName> <galaxyCatalog> <correlationFunctionFileName>")
  ! Get the arguments.
  call Get_Argument(1,parameterFileName          )
  call Get_Argument(2,galaxyCatalogFileName      )
  call Get_Argument(3,correlationFunctionFileName)
  ! Open the parameter file.
  parameters=inputParameters(parameterFileName)
  call parameters%markGlobal()
  !# <inputParameter>
  !#   <name>mockCorrelationFunctionMassMinimum</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>0.0d0</defaultValue>
  !#   <description>The minimum mass galaxy to include in a mock catalog correlation function calculation.</description>
  !#   <source>globalParameters</source>
  !#   <type>float</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>mockCorrelationFunctionMassMaximum</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>1.0d16</defaultValue>
  !#   <description>The maximum mass galaxy to include in a mock catalog correlation function calculation.</description>
  !#   <source>globalParameters</source>
  !#   <type>float</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>mockCorrelationFunctionSeparationMinimum</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>0.1d0</defaultValue>
  !#   <description>The minimum separation to compute in a mock catalog correlation function calculation.</description>
  !#   <source>globalParameters</source>
  !#   <type>float</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>mockCorrelationFunctionSeparationMaximum</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>30.0d0</defaultValue>
  !#   <description>The maximum separation to compute in a mock catalog correlation function calculation.</description>
  !#   <source>globalParameters</source>
  !#   <type>float</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>mockCorrelationFunctionSeparationCount</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>15</defaultValue>
  !#   <description>The number of bins in separation to compute in a mock catalog correlation function calculation.</description>
  !#   <source>globalParameters</source>
  !#   <type>float</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>mockCorrelationFunctionRandomSampleCount</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>var_str('*10')</defaultValue>
  !#   <description>The number of random points to use when constructing random catalogs. Can be either a fixed number or, if prefixed with ``{\normalfont \ttfamily *}'', a multiplicative factor.</description>
  !#   <source>globalParameters</source>
  !#   <type>float</type>
  !# </inputParameter>
  if (extract(mockCorrelationFunctionRandomSampleCount,1,1) == "*") then
     mockCorrelationFunctionRandomSampleCountType=mockCorrelationFunctionRandomSampleCountMultiplicative
     label=char(extract(mockCorrelationFunctionRandomSampleCount,2))
     read (label,*) mockCorrelationFunctionRandomSampleCountValue
  else
     mockCorrelationFunctionRandomSampleCountType=mockCorrelationFunctionRandomSampleCountFixed
     label=char(mockCorrelationFunctionRandomSampleCount)
     read (label,*) mockCorrelationFunctionRandomSampleCountValue
  end if
  !# <inputParameter>
  !#   <name>mockCorrelationFunctionRadialSeparationMaximum</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>40.0d0</defaultValue>
  !#   <description>The maximum radial separation of galaxies to consider when computing projected correlation functions.</description>
  !#   <source>globalParameters</source>
  !#   <type>float</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>mockCorrelationFunctionHalfIntegral</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>.false.</defaultValue>
  !#   <description>Set to {\tt true} if the projected correlation function is computed as $w_{\mathrm p}(r_{\mathrm p})=\int_0^{+\pi_{\mathrm max}} \xi(r_{\mathrm p},\pi) {\mathrm d} \pi$, instead of the usual $w_{\mathrm p}(r_{\mathrm p})=\int_{-\pi_{\mathrm max}}^{+\pi_{\mathrm max}} \xi(r_{\mathrm p},\pi) {\mathrm d} \pi$.</description>
  !#   <source>globalParameters</source>
  !#   <type>float</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>mockCorrelationFunctionBufferWidth</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>30.0d0</defaultValue>
  !#   <description>The width of the buffer region around survey geometry to ensure galaxies are not lost when moving to redshift space.</description>
  !#   <source>globalParameters</source>
  !#   <type>float</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>mockCorrelationFunctionOrigin</name>
  !#   <defaultValue>[randomSequence%sample(),randomSequence%sample(),randomSequence%sample()]</defaultValue>
  !#   <defaultSource>Uniformly random distribution within the box.</defaultSource>
  !#   <description>The vector (in units of the box length) giving the origin of the coordinate system to use in mock catalog construction.</description>
  !#   <type>float</type>
  !#   <cardinality>1</cardinality>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>mockCorrelationFunctionRotationVector</name>
  !#   <defaultValue>[acos(2.0d0*randomSequence%sample()-1.0d0),2.0d0*Pi*randomSequence%sample()]</defaultValue>
  !#   <defaultSource>Isotropically random on the unit sphere.</defaultSource>
  !#   <description>The vector, in spherical coordinates $(\theta,\phi)$, about which the mock catalog should be rotated.</description>
  !#   <type>float</type>
  !#   <cardinality>1</cardinality>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>mockCorrelationFunctionRotationAngle</name>
  !#   <defaultValue>2.0d0*Pi*randomSequence%sample()</defaultValue>
  !#   <defaultSource>Uniformly random distribution between $0$ and $2\pi$.</defaultSource>
  !#   <description>The angle through which the mock catalog should be rotated.</description>
  !#   <type>float</type>
  !#   <cardinality>1</cardinality>
  !# </inputParameter>
  ! Get required objects.
  cosmologyFunctions_ => cosmologyFunctions()
  surveyGeometry_     => surveyGeometry    ()
  ! Read the galaxy catalog.
  call Galacticus_Display_Indent("Reading galaxy catalog")
  galaxyFile=irate(char(galaxyCatalogFileName))
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
  time=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift))
  ! Prune points below our mass threshold.
  call Points_Prune(galaxyPosition,galaxyMass >= mockCorrelationFunctionMassMinimum .and. galaxyMass < mockCorrelationFunctionMassMaximum)
  call Points_Prune(galaxyVelocity,galaxyMass >= mockCorrelationFunctionMassMinimum .and. galaxyMass < mockCorrelationFunctionMassMaximum)
  message="Pruned on mass leaving "
  message=message//size(galaxyPosition,dim=2)//" galaxies"
  call Galacticus_Display_Message(message)
  ! Generate random points.
  select case (mockCorrelationFunctionRandomSampleCountType)
  case (mockCorrelationFunctionRandomSampleCountFixed         )
     randomPointCount=int(mockCorrelationFunctionRandomSampleCountValue)
  case (mockCorrelationFunctionRandomSampleCountMultiplicative)
     randomPointCount=int(mockCorrelationFunctionRandomSampleCountValue*dble(size(galaxyPosition,dim=2)))
  end select
  call allocateArray(randomPosition,[3,randomPointCount])
  do i=1,3
     do j=1,randomPointCount
        randomPosition(i,j)=randomSequence%sample()*simulationBoxSize
     end do
  end do
  ! Shift origin to a random point in the box.
  call Points_Translate(galaxyPosition,-      simulationBoxSize*mockCorrelationFunctionOrigin,periodicLength=simulationBoxSize)
  ! Shift points so origin is at the center of the box.
  call Points_Translate(galaxyPosition,-0.5d0*simulationBoxSize*[1.0d0,1.0d0,1.0d0]                                           )
  call Points_Translate(randomPosition,-0.5d0*simulationBoxSize*[1.0d0,1.0d0,1.0d0]                                           )
  ! Compute correlation function.
  call Statistics_Points_Correlation(                                                                                         &
       &                             galaxyPosition                                                                         , &
       &                             randomPosition                                                                         , &
       &                             mockCorrelationFunctionSeparationMinimum                                               , &
       &                             mockCorrelationFunctionSeparationMaximum                                               , &
       &                             mockCorrelationFunctionSeparationCount                                                 , &
       &                             separation                                                                             , &
       &                             correlation                                                                            , &
       &                             projected                               =.true.                                        , &
       &                             radialSeparationMaximum                 =mockCorrelationFunctionRadialSeparationMaximum  &
       &                            )
  ! Replicate points to encompass survey geometry.
  replications=int((surveyGeometry_%distanceMaximum(mockCorrelationFunctionMassMaximum)+mockCorrelationFunctionBufferWidth)/simulationBoxSize+0.5d0)
  call Points_Replicate(galaxyPosition,simulationBoxSize,-replications*[1,1,1],+replications*[1,1,1])
  message="Replicated to cover survey volume giving "
  message=message//size(galaxyPosition,dim=2)//" galaxies"
  call Galacticus_Display_Message(message)
  replicatedGalaxyCount=size(galaxyPosition,dim=2)
  ! Shift points into redshift space.
  call Points_To_Redshift_Space()
  ! Rotate the points.
  rotationAxis =[                     &
       &         sin(mockCorrelationFunctionRotationVector(1))*cos(mockCorrelationFunctionRotationVector(2)), &
       &         sin(mockCorrelationFunctionRotationVector(1))*sin(mockCorrelationFunctionRotationVector(2)), &
       &         cos(mockCorrelationFunctionRotationVector(1))                                                &
       &        ]
  call Points_Rotate(galaxyPosition,rotationAxis,mockCorrelationFunctionRotationAngle)
  ! Limit points to survey geometry.
  call Points_Survey_Geometry(galaxyPosition,surveyGeometry_,mockCorrelationFunctionMassMaximum)
  message="Pruned on survey geometry leaving "
  message=message//size(galaxyPosition,dim=2)//" galaxies"
  call Galacticus_Display_Message(message)
  ! Generate random points.
  select case (mockCorrelationFunctionRandomSampleCountType)
  case (mockCorrelationFunctionRandomSampleCountFixed         )
     randomPointCount=int(mockCorrelationFunctionRandomSampleCountValue*dble(replicatedGalaxyCount)/dble(size(galaxyPosition,dim=2)))
  case (mockCorrelationFunctionRandomSampleCountMultiplicative)
     randomPointCount=int(mockCorrelationFunctionRandomSampleCountValue*dble(replicatedGalaxyCount))
  end select
  call deallocateArray(randomPosition                     )
  call allocateArray  (randomPosition,[3,randomPointCount])
  do i=1,3
     do j=1,randomPointCount
        randomPosition(i,j)=(randomSequence%sample()-0.5d0)*dble(2*replications+1)*simulationBoxSize
     end do
  end do
  message="Generated "
  message=message//size(randomPosition,dim=2)//" random points"
  call Galacticus_Display_Message(message)
  call Points_Survey_Geometry(randomPosition,surveyGeometry_,mockCorrelationFunctionMassMaximum)
  message="Pruned on survey geometry leaving "
  message=message//size(randomPosition,dim=2)//" random points"
  call Galacticus_Display_Message(message)
  ! Compute the correlation function.
  call Statistics_Points_Correlation(                                                                                         &
       &                             galaxyPosition                                                                         , &
       &                             randomPosition                                                                         , &
       &                             mockCorrelationFunctionSeparationMinimum                                               , &
       &                             mockCorrelationFunctionSeparationMaximum                                               , &
       &                             mockCorrelationFunctionSeparationCount                                                 , &
       &                             separation                                                                             , &
       &                             correlationSurvey                                                                      , &
       &                             projected                               =.true.                                        , &
       &                             radialSeparationMaximum                 =mockCorrelationFunctionRadialSeparationMaximum, &
       &                             halfIntegral                            =mockCorrelationFunctionHalfIntegral             &
       &                            )
  ! Write correlations to file.
  call correlationFunctionFile%openFile(char(correlationFunctionFileName))
  call correlationFunctionFile%writeDataset(separation       ,'separation'                ,commentText='Galaxy separation.'                                ,datasetReturned=thisDataset)
  call thisDataset%writeAttribute(megaParsec,'unitsInSI')
  call thisDataset%close         (                      )
  call correlationFunctionFile%writeDataset(correlation      ,'projectedCorrelation'      ,commentText='Projected correlation function from the full mock.',datasetReturned=thisDataset)
  call thisDataset%writeAttribute(megaParsec,'unitsInSI')
  call thisDataset%close         (                      )
  call correlationFunctionFile%writeDataset(correlationSurvey,'projectedCorrelationSurvey',commentText='Projected correlation function from survey region.',datasetReturned=thisDataset)
  call thisDataset%writeAttribute(megaParsec,'unitsInSI')
  call thisDataset%close         (                      )
  call correlationFunctionFile%close()

contains

  subroutine Points_To_Redshift_Space()
    !% Shift the set of points into redshift space.
    use Vectors
    implicit none
    integer          :: galaxyCount       , replicantCount, i, j, k, l
    double precision :: velocityToDistance

    velocityToDistance=(1.0d0+redshift)/cosmologyFunctions_%hubbleParameterEpochal(time)
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
  end subroutine Points_To_Redshift_Space

end program Mocks_Correlation_Functions
