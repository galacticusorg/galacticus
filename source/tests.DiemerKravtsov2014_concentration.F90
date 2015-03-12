!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a program which tests the \cite{diemer_universal_2014} halo concentration algorithm.

program Test_DiemerKravtsov2014_Concentration
  !% Tests the \cite{diemer_universal_2014} halo concentration algorithm. Values of concentration are taken from their \href{http://www.benediktdiemer.com/wp-content/uploads/2014/07/Concentration_WMAP7_median.txt}{website}.
  use ISO_Varying_String
  use Memory_Management
  use Input_Parameters
  use Dark_Matter_Profiles_Concentration
  use Cosmology_Functions
  use Cosmology_Parameters
  use Galacticus_Nodes
  use Unit_Tests
  use System_Command
  use Galacticus_Error
  use File_Utilities
  use Galacticus_Input_Paths
  implicit none
  type            (treeNode                                        ), pointer :: node
  class           (nodeComponentBasic                              ), pointer :: basic
  class           (cosmologyFunctionsClass                         ), pointer :: cosmologyFunctions_
  class           (cosmologyParametersClass                        ), pointer :: cosmologyParameters_
  type            (darkMatterProfileConcentrationDiemerKravtsov2014)          :: darkMatterProfileConcentration_
  type            (varying_string                                  )          :: message                        , parameterFile
  double precision                                                            :: ourConcentration               , differenceFractional, &
       &                                                                         concentration                  , mass                , &
       &                                                                         redshift                       , nu                  , &
       &                                                                         differenceFractionalMaximum
  integer                                                                     :: iMass                          , ioStatus            , &
       &                                                                         referenceUnit

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.DiemerKravtsov2014_concentration.size')
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("DiemerKravtsov2014 halo concentration algorithm")
  ! Test DiemerKravtsov2014 halo concentration algorithm.
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/DiemerKravtsov2014HaloConcentration/testParameters.xml'
  call Input_Parameters_File_Open(parameterFile)

  ! Get the data file if we don't have it.
  if (.not.File_Exists(Galacticus_Input_Path()//"testSuite/data/diemerKravtsov2014Concentration.txt")) then
     call System_Command_Do(                                                                                                     &
          &                 "wget http://www.benediktdiemer.com/wp-content/uploads/2014/07/Concentration_WMAP7_median.txt -O "// &
          &                 Galacticus_Input_Path()                                                                           // &
          &                 "testSuite/data/diemerKravtsov2014Concentration.txt"                                                 &
          &                )
     if (.not.File_Exists(Galacticus_Input_Path()//"testSuite/data/diemerKravtsov2014Concentration.txt")) &
          & call Galacticus_Error_Report('Test_DiemerKravtsov2014_Concentration','unable to retrieve reference dataset')
  end if
  ! Create a node.
  node                            => treeNode                                        (                 )
  ! Get the basic component.
  basic                           => node%basic                                      (autoCreate=.true.)
  ! Get required objects.
  cosmologyFunctions_             => cosmologyFunctions                              (                 )
  cosmologyParameters_            => cosmologyParameters                             (                 )  
  darkMatterProfileConcentration_ =  darkMatterProfileConcentrationDiemerKravtsov2014(                 )
  ! Read the reference file.
  differenceFractionalMaximum=0.0d0
  open(newUnit=referenceUnit,file=char(Galacticus_Input_Path()//"testSuite/data/diemerKravtsov2014Concentration.txt"),status='old',form='formatted',iostat=ioStatus)
  read (referenceUnit,*,ioStat=ioStatus) ! Skip header.
  read (referenceUnit,*,ioStat=ioStatus)
  read (referenceUnit,*,ioStat=ioStatus)
  do while (ioStatus == 0)
     read (referenceUnit,*,ioStat=ioStatus) redshift,nu,mass,concentration
     if (ioStatus /= 0) exit
     ! Set the time for the node.
     call basic%timeSet            (cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)))
     call basic%timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)))
     ! Set the mass of the original node (Diemer & Kravtsov masses are in units of M☉/h, so
     ! we convert from that system).
     call basic%massSet(mass/cosmologyParameters_%HubbleConstant(units=unitsLittleH))
     ! Compute and compare concentration at z=0.
     ourConcentration           =darkMatterProfileConcentration_%concentration(node)
     differenceFractional       =abs(ourConcentration-concentration)/concentration
     differenceFractionalMaximum=max(differenceFractionalMaximum,differenceFractional)
  end do
  ! Assert that the maximum fractional difference is not too large. The ~2% differences are
  ! presumably because we don't use precisely the same transfer function as do Diemer &
  ! Kravtsov.
  call Assert("Halo concentration in WMAP7 reference model",differenceFractionalMaximum,0.0d0,absTol=2.1d-2)
  ! Close the input parameter file.
  call Input_Parameters_File_Close
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()
end program Test_DiemerKravtsov2014_Concentration
