!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!+ Contributions to this file made by: Stephanie Dörschner.

!% Contains a module which handles reading the Rockstar merger trees from the \href{http://hipacc.ucsc.edu/Bolshoi/MergerTrees.html}{Bolshoi} simulation.

module Merger_Trees_Bolshoi
  !% Handles reading the Rockstar merger trees from the \href{http://hipacc.ucsc.edu/Bolshoi/MergerTrees.html}{Bolshoi} simulation.
  implicit none
  private
  public :: Merger_Trees_Bolshoi_Process

contains

  subroutine Merger_Trees_Bolshoi_Process(nodesFile,mergerTrees)
    !% Reads the Rockstar merger trees from the \href{http://hipacc.ucsc.edu/Bolshoi/MergerTrees.html}{Bolshoi} simulation.
    use ISO_Varying_String
    use Merger_Tree_Data_Structure
    use File_Utilities
    use Input_Parameters
    use Dates_and_Times
    use Numerical_Constants_Astronomical   
    implicit none
    character(len=*         ), intent(in   ) :: nodesFile
    type     (mergerTreeData), intent(inout) :: mergerTrees
    integer                                  :: lineCountData, lineCountTotal, lineNumberStart, lineNumberStop


    ! Process the file.

    ! Find number of lines in file, with and without comments.
    lineCountTotal=Count_Lines_in_File(nodesFile    )
    lineCountData =Count_Lines_in_File(nodesFile,"#")

    ! Find lines number ranges with data.
    lineNumberStart=lineCountTotal-lineCountData + 1
    lineNumberStop =lineCountTotal

    ! Set particle mass.
    call mergerTrees%setParticleMass         (1.35d8)

    ! Specify that trees are self-contained.
    call mergerTrees%setSelfContained        (.true.)

    ! Specify Hubble flow inclusion
    call mergerTrees%setIncludesHubbleFlow   (.false.)
   
    ! Specify subhalo mass inclusion
    call mergerTrees%setIncludesSubhaloMasses(.true.)
    
    ! Specify dummyHostId for self-hosting halos. As default Galacticus assumes hostIndex == nodeIndex in case of self-hosting halos.
    call mergerTrees%setDummyHostId          (-1)

    ! Specify properties that need conversion due to inconsistent units.
    call mergerTrees%setConversionFactor     (propertyTypeScaleRadius, 1.0d-3)

    ! Specify that we do not want to create individual merger tree reference datasets.
    call mergerTrees%makeReferences          (.false.)

    ! Specify units system used.
    call mergerTrees%setUnits(unitsMass    ,unitsInSI=massSolar ,hubbleExponent=-1,scaleFactorExponent=0, name = "Msolar/h"      )
    call mergerTrees%setUnits(unitsLength  ,unitsInSI=megaParsec,hubbleExponent=-1,scaleFactorExponent=1, name = "Mpc/h comoving")
    call mergerTrees%setUnits(unitsVelocity,unitsInSI=kilo      ,hubbleExponent= 0,scaleFactorExponent=0, name = "km/s"          )

    ! Set cosmology metadata.
    call mergerTrees%addMetadata(metaDataCosmology  ,'OmegaMatter'               , 0.2700d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'OmegaBaryon'               , 0.0469d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'OmegaLambda'               , 0.7300d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'HubbleParam'               , 0.7000d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'sigma_8'                   , 0.8200d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'powerSpectrumIndex'        , 0.9500d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'transferFunction'          ,'CAMB'                               )

    ! Set simulation metadata.
    call mergerTrees%addMetadata(metaDataSimulation ,'code'                      ,'ART (Kravtsov et al. 1997, Kravtsov 1999)')
    call mergerTrees%addMetadata(metaDataSimulation ,'boxSize'                   , 2.5000d2                           )
    call mergerTrees%addMetadata(metaDataSimulation ,'startRedshift'             , 8.0000d1                           )
    call mergerTrees%addMetadata(metaDataSimulation ,'initialConditions'         ,"Zel\'dovich approximation"         )
    call mergerTrees%addMetadata(metaDataSimulation ,'softeningKernel'           ,'n/a'                               )
    call mergerTrees%addMetadata(metaDataSimulation ,'softeningPlummerEquivalent', 1.0d-3                             )
    call mergerTrees%addMetadata(metaDataSimulation ,'TypeOfTimestepCriterion'   ,'a'                                 )

    ! Set group finder metadata.
    call mergerTrees%addMetadata(metaDataGroupFinder,'code'                      ,'ROCKSTAR'                          )
    call mergerTrees%addMetadata(metaDataGroupFinder,'minimumParticleNumber'     , 20                                 )
    call mergerTrees%addMetadata(metaDataGroupFinder,'linkingLength'             , 0.2800d0                           )

    ! Set provenance metadata.
    call mergerTrees%addMetadata(metaDataProvenance ,'fileBuiltBy'               ,'Galacticus'                        )
    call mergerTrees%addMetadata(metaDataProvenance ,'fileTimestamp'             ,char(Formatted_Date_and_Time())     )
    call mergerTrees%addMetadata(metaDataProvenance ,'source'                    ,'http://hipacc.ucsc.edu/Bolshoi/MergerTrees.html')

    ! Set columns to read. 
    call mergerTrees%setPropertyColumn(propertyTypeTreeIndex         , 1)
    call mergerTrees%setPropertyColumn(propertyTypeScaleFactor       , 2)
    call mergerTrees%setPropertyColumn(propertyTypeNodeIndex         , 3)
    call mergerTrees%setPropertyColumn(propertyTypeDescendentIndex   , 5)
    call mergerTrees%setPropertyColumn(propertyTypeHostIndex         , 8)
    call mergerTrees%setPropertyColumn(propertyTypeNodeMass          ,11)
    call mergerTrees%setPropertyColumn(propertyTypeScaleRadius       ,14)
    call mergerTrees%setPropertyColumn(propertyTypeVelocityDispersion,15)
    call mergerTrees%setPropertyColumn(propertyTypeVelocityMaximum   ,18)
    call mergerTrees%setPropertyColumn(propertyTypePositionX         ,19)
    call mergerTrees%setPropertyColumn(propertyTypePositionY         ,20)
    call mergerTrees%setPropertyColumn(propertyTypePositionZ         ,21)
    call mergerTrees%setPropertyColumn(propertyTypeVelocityX         ,22)
    call mergerTrees%setPropertyColumn(propertyTypeVelocityY         ,23)
    call mergerTrees%setPropertyColumn(propertyTypeVelocityZ         ,24)
    call mergerTrees%setPropertyColumn(propertyTypeAngularMomentumX  ,25)
    call mergerTrees%setPropertyColumn(propertyTypeAngularMomentumY  ,26)
    call mergerTrees%setPropertyColumn(propertyTypeAngularMomentumZ  ,27)
    call mergerTrees%setPropertyColumn(propertyTypeSpin              ,28)
    call mergerTrees%setPropertyColumn(propertyTypeSnapshot          ,33)

    ! Read in the data.
    call mergerTrees%readASCII(nodesFile,lineNumberStart=lineNumberStart,lineNumberStop=lineNumberStop,separator=" ",maximumRedshift=80.0d0)   
    return
  end subroutine Merger_Trees_Bolshoi_Process

end module Merger_Trees_Bolshoi
