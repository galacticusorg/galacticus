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
    use Input_Parameters
    use Dates_and_Times
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Cosmological_Parameters
    implicit none
    character       (len=*         ), intent(in   ) :: nodesFile
    type            (mergerTreeData), intent(inout) :: mergerTrees
    integer                                         :: lineCountData            , lineCountTotal    , &
         &                                             lineNumberStart          , lineNumberStop
    double precision                                :: boxSize                  , powerSpectrumIndex, &
         &                                             sigma_8
    type            (varying_string)                :: source                   , transferFunction
    logical                                         :: haloMassesIncludeSubhalos

    ! Process the file.
    ! Find number of lines in file, with and without comments.
    lineCountTotal=Count_Lines_in_File(nodesFile    )
    lineCountData =Count_Lines_in_File(nodesFile,"#")

    ! Find lines number ranges with data.
    lineNumberStart=lineCountTotal-lineCountData+1
    lineNumberStop =lineCountTotal

    ! Set columns to read.
    call mergerTrees%setPropertyColumn(propertyTypeTreeIndex      ,1)
    call mergerTrees%setPropertyColumn(propertyTypeNodeIndex      ,2)
    call mergerTrees%setPropertyColumn(propertyTypeHostIndex      ,3)
    call mergerTrees%setPropertyColumn(propertyTypeDescendentIndex,4)
    call mergerTrees%setPropertyColumn(propertyTypeNodeMass       ,5)
    call mergerTrees%setPropertyColumn(propertyTypeRedshift       ,6)

    ! Read in the data.
    call mergerTrees%readASCII(nodesFile,lineNumberStart=lineNumberStart,lineNumberStop=lineNumberStop,separator=",")

    ! Specify that we do not want to create individual merger tree reference datasets.
    call mergerTrees%makeReferences          (.false.)

    ! Specify that trees are self-contained (i.e. nodes never move from one tree to another).
    call mergerTrees%setSelfContained        (.true. )

    ! Specify units system used.
    call mergerTrees%setUnits(unitsMass  ,unitsInSI=massSolar ,hubbleExponent=0,scaleFactorExponent=0)
    call mergerTrees%setUnits(unitsLength,unitsInSI=megaParsec,hubbleExponent=0,scaleFactorExponent=0)

    ! Set cosmology metadata.
    call    mergerTrees%addMetadata(metaDataCosmology  ,'OmegaMatter'               ,Omega_Matter()                     )
    call    mergerTrees%addMetadata(metaDataCosmology  ,'OmegaBaryon'               ,Omega_b     ()                     )
    call    mergerTrees%addMetadata(metaDataCosmology  ,'OmegaLambda'               ,Omega_DE    ()                     )
    call    mergerTrees%addMetadata(metaDataCosmology  ,'HubbleParam'               ,Little_H_0  ()                     )
    if (Input_Parameter_Is_Present("sigma_8"           )) then
       !@ <inputParameter>
       !@   <name>sigma_8</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The fractional mass fluctuation in the linear density field at the present day in spheres of radius 8~Mpc/h.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>cosmology</group>
       !@ </inputParameter>
       call Get_Input_Parameter    (                    'sigma_8'                   ,sigma_8                            )
       call mergerTrees%addMetadata(metaDataCosmology  ,'sigma_8'                   ,sigma_8                            )
    end if
    if (Input_Parameter_Is_Present("powerSpectrumIndex")) then
       !@ <inputParameter>
       !@   <name>powerSpectrumIndex</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The index of the power-law primordial power spectrum.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>cosmology</group>
       !@ </inputParameter>
       call Get_Input_Parameter    (                    'powerSpectrumIndex'        ,powerSpectrumIndex                 )
       call mergerTrees%addMetadata(metaDataCosmology  ,'powerSpectrumIndex'        ,powerSpectrumIndex                 )
    end if
    if (Input_Parameter_Is_Present("transferFunction"  )) then
       !@ <inputParameter>
       !@   <name>transferFunction</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The type of transfer function used.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>cosmology</group>
       !@ </inputParameter>
       call Get_Input_Parameter    (                    'transferFunction'          ,     transferFunction              )
       call mergerTrees%addMetadata(metaDataCosmology  ,'transferFunction'          ,char(transferFunction)             )
    end if

    ! Set provenance metadata.
    call    mergerTrees%addMetadata(metaDataProvenance ,'fileBuiltBy'               ,'Galacticus'                       )
    call    mergerTrees%addMetadata(metaDataProvenance ,'fileTimestamp'             ,char(Formatted_Date_and_Time())    )
    if (Input_Parameter_Is_Present("source"            )) then
       !@ <inputParameter>
       !@   <name>source</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The source of the merger trees.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>cosmology</group>
       !@ </inputParameter>
       call Get_Input_Parameter    (                    'source'                    ,     source                        )
       call mergerTrees%addMetadata(metaDataProvenance ,'source'                    ,char(source)                       )
    end if

    ! Set simulation metadata.
    if (Input_Parameter_Is_Present("boxSize"           )) then
       !@ <inputParameter>
       !@   <name>boxSize</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The box size of the simulation from which merger trees were extracted.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>cosmology</group>
       !@ </inputParameter>
       call Get_Input_Parameter    (                    'boxSize'                   ,boxSize                            )
       call mergerTrees%addMetadata(metaDataSimulation ,'boxSize'                   ,boxSize                            )
    end if

    ! Set halo properties.
    !@ <inputParameter>
    !@   <name>haloMassesIncludeSubhalos</name>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     Specifies whether or not halo masses include the masses of their subhalos.
    !@   </description>
    !@   <type>real</type>
    !@   <cardinality>1</cardinality>
    !@   <group>cosmology</group>
    !@ </inputParameter>
    call Get_Input_Parameter('haloMassesIncludeSubhalos',haloMassesIncludeSubhalos,defaultValue=.true.)
    call mergerTrees%setIncludesSubhaloMasses(haloMassesIncludeSubhalos)
    return
  end subroutine Merger_Trees_Simple_Process

end module Merger_Trees_Simple
