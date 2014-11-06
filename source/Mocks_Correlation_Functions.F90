!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  type            (hdf5Object             )                              :: galaxyCatalogFile                       , snapshotGroup                                 , &
       &                                                                    galaxyCatalogGroup                      , simulationGroup                               , &
       &                                                                    correlationFunctionFile                 , thisDataset
  type            (irate                  )                              :: galaxyFile
  double precision                                                       :: simulationBoxSize                       , separationMinimum                             , &
       &                                                                    separationMaximum                       , mockCorrelationFunctionRadialSeparationMaximum, &
       &                                                                    mockCorrelationFunctionMassMinimum      , mockCorrelationFunctionMassMaximum            , &
       &                                                                    mockCorrelationFunctionSeparationMinimum, mockCorrelationFunctionSeparationMaximum      , &
       &                                                                    redshift                                , time                                          , &
       &                                                                    mockCorrelationFunctionBufferWidth      , mockCorrelationFunctionRotationAngle          , &
       &                                                                    mockCorrelationFunctionRandomSampleCountValue
  integer                                                                :: separationCount                         , randomPointCount                              , &
       &                                                                    i                                       , j                                             , &
       &                                                                    replications                            , mockCorrelationFunctionSeparationCount        , &
       &                                                                    mockCorrelationFunctionRandomSampleCountType, replicatedGalaxyCount
  type            (pseudoRandom  )                                       :: randomSequence
  character       (len=128       )                                       :: label

  ! Read in basic code memory usage.
  call Code_Memory_Usage('Mocks_Correlation_Functions.size')
  ! Check that correct number of arguments have been supplied.
  if (Command_Argument_Count() /= 3) call Galacticus_Error_Report(message="Usage: Mocks_Correlation_Functions.exe <parameterFileName> <galaxyCatalog> <correlationFunctionFileName>")
  ! Get the arguments.
  call Get_Argument              (1,parameterFileName          )
  call Get_Argument              (2,galaxyCatalogFileName      )
  call Get_Argument              (3,correlationFunctionFileName)
  ! Open the parameter file.
  call Input_Parameters_File_Open(  parameterFileName          )
  !@ <inputParameter>
  !@   <name>mockCorrelationFunctionMassMinimum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$0$</defaultValue>
  !@   <description>
  !@     The minimum mass galaxy to include in a mock catalog correlation function calculation.
  !@   </description>
  !@   <type>float</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('mockCorrelationFunctionMassMinimum',mockCorrelationFunctionMassMinimum,defaultValue=0.0d0)
  !@ <inputParameter>
  !@   <name>mockCorrelationFunctionMassMaximum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$10^{16}M_\odot$</defaultValue>
  !@   <description>
  !@     The maximum mass galaxy to include in a mock catalog correlation function calculation.
  !@   </description>
  !@   <type>float</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('mockCorrelationFunctionMassMaximum',mockCorrelationFunctionMassMaximum,defaultValue=1.0d16)
  !@ <inputParameter>
  !@   <name>mockCorrelationFunctionSeparationMinimum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$0.1$ Mpc</defaultValue>
  !@   <description>
  !@     The minimum separation to compute in a mock catalog correlation function calculation.
  !@   </description>
  !@   <type>float</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('mockCorrelationFunctionSeparationMinimum',mockCorrelationFunctionSeparationMinimum,defaultValue=0.1d0)
  !@ <inputParameter>
  !@   <name>mockCorrelationFunctionSeparationMaximum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$30$ Mpc</defaultValue>
  !@   <description>
  !@     The maximum separation to compute in a mock catalog correlation function calculation.
  !@   </description>
  !@   <type>float</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('mockCorrelationFunctionSeparationMaximum',mockCorrelationFunctionSeparationMaximum,defaultValue=30.0d0)
  !@ <inputParameter>
  !@   <name>mockCorrelationFunctionSeparationCount</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>15</defaultValue>
  !@   <description>
  !@     The number of bins in separation to compute in a mock catalog correlation function calculation.
  !@   </description>
  !@   <type>float</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('mockCorrelationFunctionSeparationCount',mockCorrelationFunctionSeparationCount,defaultValue=15)
  !@ <inputParameter>
  !@   <name>mockCorrelationFunctionRandomSampleCount</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>*10</defaultValue>
  !@   <description>
  !@     The number of random points to use when constructing random catalogs. Can be either a fixed number or, if prefixed with ``{\tt *}'', a multiplicative factor.
  !@   </description>
  !@   <type>float</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('mockCorrelationFunctionRandomSampleCount',mockCorrelationFunctionRandomSampleCount,defaultValue='*10')
  if (extract(mockCorrelationFunctionRandomSampleCount,1,1) == "*") then
     mockCorrelationFunctionRandomSampleCountType=mockCorrelationFunctionRandomSampleCountMultiplicative
     label=char(extract(mockCorrelationFunctionRandomSampleCount,2))
     read (label,*) mockCorrelationFunctionRandomSampleCountValue
  else
     mockCorrelationFunctionRandomSampleCountType=mockCorrelationFunctionRandomSampleCountFixed
     label=char(mockCorrelationFunctionRandomSampleCount)
     read (label,*) mockCorrelationFunctionRandomSampleCountValue
  end if
  !@ <inputParameter>
  !@   <name>mockCorrelationFunctionRadialSeparationMaximum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>40 Mpc</defaultValue>
  !@   <description>
  !@     The maximum radial separation of galaxies to consider when computing projected correlation functions.
  !@   </description>
  !@   <type>float</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('mockCorrelationFunctionRadialSeparationMaximum',mockCorrelationFunctionRadialSeparationMaximum,defaultValue=40.0d0)
  !@ <inputParameter>
  !@   <name>mockCorrelationFunctionBufferWidth</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>30 Mpc</defaultValue>
  !@   <description>
  !@     The width of the buffer region around survey geometry to ensure galaxies are not lost when moving to redshift space.
  !@   </description>
  !@   <type>float</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('mockCorrelationFunctionBufferWidth',mockCorrelationFunctionBufferWidth,defaultValue=30.0d0)
  !@ <inputParameter>
  !@   <name>mockCorrelationFunctionOrigin</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>random</defaultValue>
  !@   <description>
  !@     The vector (in units of the box length) giving the origin of the coordinate system to use in mock catalog construction.
  !@   </description>
  !@   <type>float</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('mockCorrelationFunctionOrigin',mockCorrelationFunctionOrigin,defaultValue=[randomSequence%sample(),randomSequence%sample(),randomSequence%sample()])
  !@ <inputParameter>
  !@   <name>mockCorrelationFunctionRotationVector</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>random</defaultValue>
  !@   <description>
  !@     The vector, in spherical coordinates $(\theta,\phi)$, about which the mock catalog should be rotated.
  !@   </description>
  !@   <type>float</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('mockCorrelationFunctionRotationVector',mockCorrelationFunctionRotationVector,defaultValue=[acos(2.0d0   *randomSequence%sample()-1.0d0),2.0d0*Pi*randomSequence%sample()])
  !@ <inputParameter>
  !@   <name>mockCorrelationFunctionRotationAngle</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>random</defaultValue>
  !@   <description>
  !@     The angle through which the mock catalog should be rotated.
  !@   </description>
  !@   <type>float</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('mockCorrelationFunctionRotationAngle',mockCorrelationFunctionRotationAngle,defaultValue=2.0d0*Pi*randomSequence%sample())
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
  call Alloc_Array(randomPosition,[3,randomPointCount])
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
  call Dealloc_Array(randomPosition                     )
  call Alloc_Array  (randomPosition,[3,randomPointCount])
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
       &                             radialSeparationMaximum                 =mockCorrelationFunctionRadialSeparationMaximum  &
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
  ! Close the parameter file.
  call Input_Parameters_File_Close()

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
