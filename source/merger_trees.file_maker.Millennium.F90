!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which handles reading of data from CSV files extracted from the Millennium Simulation database.

module Merger_Trees_Millennium
  !% Handles reading of data from CSV files extracted from the Millennium Simulation database.
  implicit none
  private
  public :: Merger_Trees_Millennium_Process

contains

  subroutine Merger_Trees_Millennium_Process(nodesFile,particlesFile,mergerTrees,generation)
    !% Read and process a CSV file of merger trees extracted from the Millennium Simulation database.
    use ISO_Varying_String
    use Merger_Tree_Data_Structure
    use Dates_and_Times
    use Numerical_Constants_Astronomical
    use File_Utilities
    use Galacticus_Error
    implicit none
    character       (len=*         ), intent(in   ) :: nodesFile     , particlesFile
    type            (mergerTreeData), intent(inout) :: mergerTrees
    integer                         , intent(in   ) :: generation
    integer                                         :: fileUnit      , lineCountData, lineCountTotal, lineNumberStart, &
         &                                             lineNumberStop
    character       (len=1024      )                :: sqlQuery
    logical                                         :: traceParticles
    double precision                                :: particleMass

    ! Process the nodes file.
    ! Determine if particles are being traced.
    traceParticles=(trim(particlesFile) /= "none")

    ! Retrieve the SQL query used to generate this file.
    open(newunit=fileUnit,file=nodesFile,status='old',form='formatted')
    read (fileUnit,'(a)') sqlQuery
    read (fileUnit,'(a)') sqlQuery
    close(fileUnit)

    ! Find number of lines in file, with and without comments.
    lineCountTotal=Count_Lines_in_File(nodesFile    )
    lineCountData =Count_Lines_in_File(nodesFile,"#")-1

    ! Find lines number ranges with data.
    lineNumberStart=lineCountTotal-lineCountData
    lineNumberStop =lineCountTotal-1

    ! Set columns to read.
    call                     mergerTrees%setPropertyColumn(propertyTypeTreeIndex             , 1)
    call                     mergerTrees%setPropertyColumn(propertyTypeNodeIndex             , 2)
    call                     mergerTrees%setPropertyColumn(propertyTypeDescendentIndex       , 3)
    call                     mergerTrees%setPropertyColumn(propertyTypeHostIndex             , 4)
    call                     mergerTrees%setPropertyColumn(propertyTypeSnapshot              , 5)
    call                     mergerTrees%setPropertyColumn(propertyTypeRedshift              , 6)
    call                     mergerTrees%setPropertyColumn(propertyTypeNodeMass              , 7)
    call                     mergerTrees%setPropertyColumn(propertyTypeParticleCount         , 8)
    call                     mergerTrees%setPropertyColumn(propertyTypePositionX             , 9)
    call                     mergerTrees%setPropertyColumn(propertyTypePositionY             ,10)
    call                     mergerTrees%setPropertyColumn(propertyTypePositionZ             ,11)
    call                     mergerTrees%setPropertyColumn(propertyTypeVelocityX             ,12)
    call                     mergerTrees%setPropertyColumn(propertyTypeVelocityY             ,13)
    call                     mergerTrees%setPropertyColumn(propertyTypeVelocityZ             ,14)
    call                     mergerTrees%setPropertyColumn(propertyTypeAngularMomentumX      ,15)
    call                     mergerTrees%setPropertyColumn(propertyTypeAngularMomentumY      ,16)
    call                     mergerTrees%setPropertyColumn(propertyTypeAngularMomentumZ      ,17)
    call                     mergerTrees%setPropertyColumn(propertyTypeHalfMassRadius        ,18)
    if (traceParticles) call mergerTrees%setPropertyColumn(propertyTypeMostBoundParticleIndex,19)

    ! Read in the data.
    call mergerTrees%readASCII(nodesFile,lineNumberStart=lineNumberStart,lineNumberStop=lineNumberStop,separator=",")


    ! Process the particles file if one is specified.
    if (traceParticles) then

       ! Find number of lines in file, with and without comments.
       lineCountTotal=Count_Lines_in_File(particlesFile    )
       lineCountData =Count_Lines_in_File(particlesFile,"#")-1

       ! Find lines number ranges with data.
       lineNumberStart=lineCountTotal-lineCountData+1
       lineNumberStop =lineCountTotal

       ! Set columns to read.
       call mergerTrees%setParticlePropertyColumn(propertyTypeParticleIndex,1)
       call mergerTrees%setParticlePropertyColumn(propertyTypeRedshift     ,2)
       call mergerTrees%setParticlePropertyColumn(propertyTypeSnapshot     ,3)
       call mergerTrees%setParticlePropertyColumn(propertyTypePositionX    ,4)
       call mergerTrees%setParticlePropertyColumn(propertyTypePositionY    ,5)
       call mergerTrees%setParticlePropertyColumn(propertyTypePositionZ    ,6)
       call mergerTrees%setParticlePropertyColumn(propertyTypeVelocityX    ,7)
       call mergerTrees%setParticlePropertyColumn(propertyTypeVelocityY    ,8)
       call mergerTrees%setParticlePropertyColumn(propertyTypeVelocityZ    ,9)

       ! Read in the data.
       call mergerTrees%readParticlesASCII(particlesFile,lineNumberStart=lineNumberStart,lineNumberStop=lineNumberStop,separator=",")

    end if

    ! Specify that we do not want to create individual merger tree reference datasets.
    call mergerTrees%makeReferences          (.false.)

    ! Specify particle mass (in whatever mass units the trees use - in this case 1e10 Msun/h).
    select case (generation)
    case (1)
       particleMass=8.60d-2
    case (2)
       particleMass=6.89d-4
    case default
       call Galacticus_Error_Report('Merger_Trees_Millennium_Process','unknown Millennium Simulation generation')
    end select
    call mergerTrees%setParticleMass         (particleMass)

    ! Specify that trees are self-contained (i.e. nodes never move from one tree to another).
    call mergerTrees%setSelfContained        (.true. )

    ! Specify that velocities do not include the Hubble flow.
    call mergerTrees%setIncludesHubbleFlow   (.false.)

    ! Specify that positions are periodic.
    call mergerTrees%setPositionsArePeriodic (.true. )

    ! Specify that halo masses do include subhalo contributions. (Note that the expectation is that masses of isolated halos have
    ! been derived from something like the "m_tophat" column of the Millennium Database.)
    call mergerTrees%setIncludesSubhaloMasses(.true. )

    ! Specify units system used.
    call mergerTrees%setUnits(unitsMass    ,unitsInSI=1.0d10*massSolar,hubbleExponent=-1,scaleFactorExponent=0,name="1e10 Msolar/h" )
    call mergerTrees%setUnits(unitsLength  ,unitsInSI=megaParsec      ,hubbleExponent=-1,scaleFactorExponent=1,name="comoving Mpc/h")
    call mergerTrees%setUnits(unitsVelocity,unitsInSI=kilo            ,hubbleExponent= 0,scaleFactorExponent=0,name="km/s"          )

    ! Set cosmology metadata.
    call mergerTrees%addMetadata(metaDataCosmology  ,'OmegaMatter'               ,0.250d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'OmegaBaryon'               ,0.045d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'OmegaLambda'               ,0.750d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'HubbleParam'               ,0.730d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'sigma_8'                   ,0.900d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'powerSpectrumIndex'        ,1.000d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'transferFunction'          ,'CMBFast'                          )

    ! Set simulation metadata.
    call mergerTrees%addMetadata(metaDataSimulation ,'code'                      ,'GADGET-2'                         )
    call mergerTrees%addMetadata(metaDataSimulation ,'boxSize'                   , 5.000d2                           )
    call mergerTrees%addMetadata(metaDataSimulation ,'startRedshift'             , 1.270d2                           )
    call mergerTrees%addMetadata(metaDataSimulation ,'initialConditions'         ,'glass'                            )
    call mergerTrees%addMetadata(metaDataSimulation ,'softeningKernel'           ,'spline'                           )
    call mergerTrees%addMetadata(metaDataSimulation ,'softeningPlummerEquivalent', 5.0d-3                            )
    call mergerTrees%addMetadata(metaDataSimulation ,'TypeOfTimestepCriterion'   , 0                                 )
    call mergerTrees%addMetadata(metaDataSimulation ,'ErrTolIntAccuracy'         , 0.02d0                            )

    ! Set group finder metadata.
    call mergerTrees%addMetadata(metaDataGroupFinder,'code'                      ,'SUBFIND'                          )
    call mergerTrees%addMetadata(metaDataGroupFinder,'minimumParticleNumber'     ,20                                 )
    call mergerTrees%addMetadata(metaDataGroupFinder,'linkingLength'             , 0.2d0                             )

    ! Set provenance metadata.
    call mergerTrees%addMetadata(metaDataProvenance ,'fileBuiltBy'               ,'Galacticus'                       )
    call mergerTrees%addMetadata(metaDataProvenance ,'fileTimestamp'             ,char(Formatted_Date_and_Time())    )
    call mergerTrees%addMetadata(metaDataProvenance ,'source'                    ,'http://www.g-vo.org/MyMillennium3')
    call mergerTrees%addMetadata(metaDataProvenance ,'sqlQuery'                  ,sqlQuery(7:len_trim(sqlQuery))     )

    return
  end subroutine Merger_Trees_Millennium_Process

end module Merger_Trees_Millennium
