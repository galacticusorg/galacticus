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

!% Contains a program which tests the \cite{zhao_accurate_2009} halo mass formation history and halo concentration algorithms in an
!% open Universe.

program Test_Zhao2009_Open
  !% Tests the \cite{zhao_accurate_2009} halo mass formation history and halo concentration algorithms in an open
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
  use Galacticus_Nodes
  use Unit_Tests
  use String_Handling
  use Galacticus_Input_Paths
  use File_Utilities
  implicit none
  type (mergerTree        ), pointer                 :: thisTree
  type (treeNode          ), pointer                 :: thisNode
  class(nodeComponentBasic), pointer                 :: thisBasicComponent
  integer                  , parameter, dimension(1) :: logarithmicHaloMasses           =[12    ]
  double precision         , parameter, dimension(1) :: timeDifferenceTolerance         =[2.3d-2],&
       & concentrationDifferenceTolerance=[3.6d-2]
  type (varying_string    )                          :: parameterFile,fileName,message
  integer                                            :: fUnit,totalLinesInFile,dataLinesInFile,iLine,iMass
  double precision                                   :: redshift,haloMass,theirConcentration,theirTime,ourConcentration,ourTime &
       &,timeDifferenceMaximum,concentrationDifferenceMaximum

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.Zhao2009_algorithms.open.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Zhao et al. 2009 algorithms: open cosmology")

  ! Test Zhao et al. 2009 algorithms in an open universe.
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/Zhao2009Algorithms/open.xml'
  call Input_Parameters_File_Open(parameterFile)
  
  ! Create a node.
  call thisTree%createNode(thisNode)

  ! Get the basic component.
  thisBasicComponent => thisNode%basic(autoCreate=.true.)
  
  ! Loop over halo masses to test.
  do iMass=1,size(logarithmicHaloMasses)

     ! Count lines in the "mandc" comparison file.
     fileName=char(Galacticus_Input_Path())//'testSuite/data/zhao2009MassAccretionHistories/mandcoutputOpenlgM'
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
        call thisBasicComponent%massSet(10.0d0**logarithmicHaloMasses(iMass))
        call thisBasicComponent%timeSet(Cosmology_Age                (1.0d0))

        ! Get the time corresponding to the current halo mass.
        ourTime=Dark_Matter_Halo_Mass_Accretion_Time(thisNode,haloMass)
        
        ! Set the node mass and time to the current values.
        call thisBasicComponent%massSet(haloMass )
        call thisBasicComponent%timeSet(theirTime)
        
        ! Get the corresponding halo concentration.
        ourConcentration=Dark_Matter_Profile_Concentration(thisNode)
        
        ! Compute the difference between our values and the comparison values.
        timeDifferenceMaximum         =max(                                                             &
             &                              timeDifferenceMaximum                                       &
             &                             ,abs(ourTime         -theirTime         )/theirTime         &
             &                             )
        concentrationDifferenceMaximum=max(                                                             &
             &                              concentrationDifferenceMaximum                              &
             &                             ,abs(ourConcentration-theirConcentration)/theirConcentration&
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
 
end program Test_Zhao2009_Open
