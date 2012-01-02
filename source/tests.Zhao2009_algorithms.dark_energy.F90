!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a program which tests the \cite{zhao_accurate_2009} halo mass formation history and halo concentration algorithms in a
!% dark energy Universe.

program Test_Zhao2009_Dark_Energy
  !% Tests the \cite{zhao_accurate_2009} halo mass formation history and halo concentration algorithms in a dark energy
  !% Universe. Comparisons are made to the \href{http://202.127.29.4/dhzhao/mandc_calculator.htm}{``{\tt mandc}''} Note that
  !% comparison tolerances are relatively large since we have not attempted to match details (such as critical density
  !% calculation) with ``{\tt mandc}''.
  use ISO_Varying_String
  use Memory_Management
  use Input_Parameters
  use Dark_Matter_Profiles_Concentrations
  use Dark_Matter_Halo_Mass_Accretion_Histories
  use Cosmology_Functions
  use Merger_Trees
  use Tree_Nodes
  use Unit_Tests
  use String_Handling
  use Galacticus_Input_Paths
  use File_Utilities
  implicit none
  type(mergerTree),    pointer                 :: thisTree
  type(treeNode),      pointer                 :: thisNode
  integer,             parameter, dimension(2) :: logarithmicHaloMasses           =[12    ,15    ]
  double precision,    parameter, dimension(2) :: timeDifferenceTolerance         =[2.5d-2,2.0d-2], &
       &                                          concentrationDifferenceTolerance=[3.1d-2,5.5d-4]
  type(varying_string)                         :: parameterFile,fileName,message
  integer                                      :: fUnit,totalLinesInFile,dataLinesInFile,iLine,iMass
  double precision                             :: redshift,haloMass,theirConcentration,theirTime,ourConcentration,ourTime&
       &,timeDifferenceMaximum,concentrationDifferenceMaximum

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.Zhao2009_algorithms.dark_energy.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Zhao et al. 2009 algorithms: dark energy cosmology")

  ! Test Zhao et al. 2009 algorithms in a dark energy universe.
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/Zhao2009Algorithms/darkEnergy.xml'
  call Input_Parameters_File_Open(parameterFile)
  
  ! Create a node.
  call thisTree%createNode(thisNode)
  
  ! Loop over halo masses to test.
  do iMass=1,size(logarithmicHaloMasses)

     ! Count lines in the "mandc" comparison file.
     fileName=char(Galacticus_Input_Path())//'testSuite/data/zhao2009MassAccretionHistories/mandcoutputDarkEnergylgM'
     fileName=fileName//logarithmicHaloMasses(iMass)//'.data'
     totalLinesInFile=Count_Lines_in_File(fileName    )
     dataLinesInFile =Count_Lines_in_File(fileName,'#')-1

     ! Discard file header.
     open(newunit=fUnit,file=char(fileName),status='old',form='formatted')
     do iLine=1,(totalLinesInFile-dataLinesInFile)
        read (fUnit,*)
     end do

     ! Initialize maximum differences to zero.
     timeDifferenceMaximum         =0.0d0
     concentrationDifferenceMaximum=0.0d0
     
     ! Read all data lines from the comparison file.
     do iLine=1,dataLinesInFile
        read (fUnit,*) redshift,haloMass,theirConcentration
        
        ! Compute the corresponding cosmological time.
        theirTime=Cosmology_Age(Expansion_Factor_From_Redshift(redshift))
        
        ! Set the mass and time of the original node.
        call Tree_Node_Mass_Set(thisNode,10.0d0**logarithmicHaloMasses(iMass))
        call Tree_Node_Time_Set(thisNode,Cosmology_Age                (1.0d0))

        ! Get the time corresponding to the current halo mass.
        ourTime=Dark_Matter_Halo_Mass_Accretion_Time(thisNode,haloMass)
        
        ! Set the node mass and time to the current values.
        call Tree_Node_Mass_Set(thisNode,              haloMass )
        call Tree_Node_Time_Set(thisNode,              theirTime)
        
        ! Get the corresponding halo concentration.
        ourConcentration=Dark_Matter_Profile_Concentration(thisNode)
        
        ! Compute the difference between our values and the comparison values.
        timeDifferenceMaximum         =max(                                                             &
             &                              timeDifferenceMaximum                                       &
             &                             ,dabs(ourTime         -theirTime         )/theirTime         &
             &                             )
        concentrationDifferenceMaximum=max(                                                             &
             &                              concentrationDifferenceMaximum                              &
             &                             ,dabs(ourConcentration-theirConcentration)/theirConcentration&
             &                            )

     end do
     close(fUnit)
     
     ! Perform the tests.
     message='10^'
     message=message//logarithmicHaloMasses(iMass)//' M⊙ halo mass accretion history'
     call Assert(char(message),timeDifferenceMaximum         ,0.0d0,absTol=timeDifferenceTolerance         (iMass))
     message='10^'
     message=message//logarithmicHaloMasses(iMass)//' M⊙ halo concentration history'
     call Assert(char(message),concentrationDifferenceMaximum,0.0d0,absTol=concentrationDifferenceTolerance(iMass))

  end do

  ! Close the input parameter file.
  call Input_Parameters_File_Close
 
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()
 
end program Test_Zhao2009_Dark_Energy
