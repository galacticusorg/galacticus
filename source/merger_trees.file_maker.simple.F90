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


!% Contains a module which handles reading of data from CSV files of simple merger trees.

module Merger_Trees_Simple
  !% Handles reading of data from CSV files of simple merger trees.
  implicit none
  private
  public :: Merger_Trees_Simple_Process

contains

  subroutine Merger_Trees_Simple_Process(nodesFile,mergerTrees)
    !% Read and process a CSV file of simple merger trees.
    use ISO_Varying_String
    use Merger_Tree_Data_Structure
    use File_Utilities
    use Dates_and_Times
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    implicit none
    character(len=*),     intent(in)    :: nodesFile
    type(mergerTreeData), intent(inout) :: mergerTrees
    integer                             :: lineCountTotal,lineCountData,lineNumberStart,lineNumberStop

    ! Process the file.

    ! Find number of lines in file, with and without comments.
    lineCountTotal=Count_Lines_in_File(nodesFile    )
    lineCountData =Count_Lines_in_File(nodesFile,"#")-1

    ! Find lines number ranges with data.
    lineNumberStart=lineCountTotal-lineCountData
    lineNumberStop =lineCountTotal-1

    ! Set columns to read.
    call mergerTrees%setProperty(propertyTypeTreeIndex      ,1)
    call mergerTrees%setProperty(propertyTypeNodeIndex      ,2)
    call mergerTrees%setProperty(propertyTypeHostIndex      ,3)
    call mergerTrees%setProperty(propertyTypeDescendentIndex,4)
    call mergerTrees%setProperty(propertyTypeNodeMass       ,5)
    call mergerTrees%setProperty(propertyTypeRedshift       ,6)

    ! Read in the data.
    call mergerTrees%readASCII(nodesFile,lineNumberStart=lineNumberStart,lineNumberStop=lineNumberStop,separator=",")
    
    ! Specify that we do not want to create individual merger tree reference datasets.
    call mergerTrees%makeReferences          (.false.)

    ! Specify that trees are self-contained (i.e. nodes never move from one tree to another).
    call mergerTrees%setSelfContained        (.true. )

    ! Specify units system used.
    call mergerTrees%setUnits(unitsMass,unitsInSI=massSolar,hubbleExponent=0,scaleFactorExponent=0)

    ! Set cosmology metadata.
    call mergerTrees%addMetadata(metaDataCosmology  ,'Omega0'                    ,0.300d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'OmegaBaryon'               ,0.040d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'OmegaLambda'               ,0.700d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'HubbleParam'               ,0.700d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'sigma_8'                   ,0.930d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'powerSpectrumIndex'        ,1.000d0                            )
    call mergerTrees%addMetadata(metaDataCosmology  ,'transferFunction'          ,'BBKS'                             )

    ! Set provenance metadata.
    call mergerTrees%addMetadata(metaDataProvenance ,'fileBuiltBy'               ,'Galacticus'                       )
    call mergerTrees%addMetadata(metaDataProvenance ,'fileTimestamp'             ,char(Formatted_Date_and_Time())    )
    call mergerTrees%addMetadata(metaDataProvenance ,'source'                    ,'???'                              )

    return
  end subroutine Merger_Trees_Simple_Process
  
end module Merger_Trees_Simple
