!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which handles reading of data from CSV files extracted from the Millennium Simulation database.

module Merger_Trees_Millennium
  !% Handles reading of data from CSV files extracted from the Millennium Simulation database.
  private
  public :: Merger_Trees_Millennium_Process

contains

  subroutine Merger_Trees_Millennium_Process(nodesFile,particlesFile,mergerTrees)
    !% Read and process a CSV file of merger trees extracted from the Millennium Simulation database.
    use ISO_Varying_String
    use Merger_Tree_Data_Structure
    use File_Utilities
    use Dates_and_Times
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    implicit none
    character(len=*),     intent(in)    :: nodesFile,particlesFile
    type(mergerTreeData), intent(inout) :: mergerTrees
    integer                             :: lineCountTotal,lineCountData,lineNumberStart,lineNumberStop,fileUnit
    character(len=1024)                 :: sqlQuery

    ! Process the nodes file.

    ! Retrieve the SQL query used to generate this file.
    fileUnit=File_Units_Get()
    open(fileUnit,file=nodesFile,status='old',form='formatted')
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
    call mergerTrees%setProperty(propertyTypeTreeIndex             , 1)
    call mergerTrees%setProperty(propertyTypeNodeIndex             , 2)
    call mergerTrees%setProperty(propertyTypeDescendentIndex       , 3)
    call mergerTrees%setProperty(propertyTypeHostIndex             , 4)
    call mergerTrees%setProperty(propertyTypeRedshift              , 6)
    call mergerTrees%setProperty(propertyTypeNodeMass              , 7)
    call mergerTrees%setProperty(propertyTypeParticleCount         , 8)
    call mergerTrees%setProperty(propertyTypePositionX             , 9)
    call mergerTrees%setProperty(propertyTypePositionY             ,10)
    call mergerTrees%setProperty(propertyTypePositionZ             ,11)
    call mergerTrees%setProperty(propertyTypeVelocityX             ,12)
    call mergerTrees%setProperty(propertyTypeVelocityY             ,13)
    call mergerTrees%setProperty(propertyTypeVelocityZ             ,14)
    call mergerTrees%setProperty(propertyTypeAngularMomentumX      ,15)
    call mergerTrees%setProperty(propertyTypeAngularMomentumY      ,16)
    call mergerTrees%setProperty(propertyTypeAngularMomentumZ      ,17)
    call mergerTrees%setProperty(propertyTypeHalfMassRadius        ,18)
    call mergerTrees%setProperty(propertyTypeMostBoundParticleIndex,19)

    ! Read in the data.
    call mergerTrees%readASCII(nodesFile,lineNumberStart=lineNumberStart,lineNumberStop=lineNumberStop,separator=",")


    ! Process the particles file.

    ! Find number of lines in file, with and without comments.
    lineCountTotal=Count_Lines_in_File(particlesFile    )
    lineCountData =Count_Lines_in_File(particlesFile,"#")-1

    ! Find lines number ranges with data.
    lineNumberStart=lineCountTotal-lineCountData
    lineNumberStop =lineCountTotal-1

    ! Set columns to read.
    call mergerTrees%setParticleProperty(propertyTypeParticleIndex,1)
    call mergerTrees%setParticleProperty(propertyTypeRedshift     ,2)
    call mergerTrees%setParticleProperty(propertyTypePositionX    ,3)
    call mergerTrees%setParticleProperty(propertyTypePositionY    ,4)
    call mergerTrees%setParticleProperty(propertyTypePositionZ    ,5)
    call mergerTrees%setParticleProperty(propertyTypeVelocityX    ,6)
    call mergerTrees%setParticleProperty(propertyTypeVelocityY    ,7)
    call mergerTrees%setParticleProperty(propertyTypeVelocityZ    ,8)

    ! Read in the data.
    call mergerTrees%readParticlesASCII(particlesFile,lineNumberStart=lineNumberStart,lineNumberStop=lineNumberStop,separator=",")

    ! Specify particle mass (in whatever mass units the trees use - in this case 1e10 Msun/h).
    call mergerTrees%setParticleMass         (8.6d-2 )

    ! Specify that trees are self-contained (i.e. nodes never move from one tree to another).
    call mergerTrees%setSelfContained        (.true. )

    ! Specify that velocities do not include the Hubble flow.
    call mergerTrees%setIncludesHubbleFlow   (.false.)

    ! Specify that halo masses do not include subhalo contributions.
    call mergerTrees%setIncludesSubhaloMasses(.false.)

    ! Specify units system used.
    call mergerTrees%setUnits(unitsMass    ,unitsInSI=1.0d10*massSolar,hubbleExponent=-1,scaleFactorExponent=0)
    call mergerTrees%setUnits(unitsLength  ,unitsInSI=megaParsec      ,hubbleExponent=-1,scaleFactorExponent=1)
    call mergerTrees%setUnits(unitsVelocity,unitsInSI=kilo            ,hubbleExponent= 0,scaleFactorExponent=0)

    ! Set cosmology metadata.
    call mergerTrees%addMetadata(metaDataCosmology  ,'Omega0'                    ,0.250d0                            )
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
