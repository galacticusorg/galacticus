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

!% Contains a program which tests the \cite{prada_halo_2011} halo concentration algorithm.

program Test_Prada2011_Concentration
  !% Tests the \cite{prada_halo_2011} halo concentration algorithm. Values of concentration were read from their Figure~12.
  use ISO_Varying_String
  use Memory_Management
  use Input_Parameters
  use Dark_Matter_Profiles_Concentration
  use Cosmology_Functions
  use Cosmology_Parameters
  use Galacticus_Nodes
  use Unit_Tests
  implicit none
  type            (treeNode                           )                                 , pointer :: thisNode
  class           (nodeComponentBasic                 )                                 , pointer :: thisBasicComponent
  class           (cosmologyParametersClass           )                                 , pointer :: thisCosmologyParameters
  class           (darkMatterProfileConcentrationClass)                                 , pointer :: darkMatterProfileConcentration_
  type            (varying_string                     )                                           :: message                                                        , parameterFile
  integer                                                                    , parameter          :: massCount                =4
  double precision                                     , dimension(massCount), parameter          :: logMass                  =[11.000d0,12.000d0,13.000d0,14.000d0]
  double precision                                     , dimension(massCount), parameter          :: pradaLogConcentration    =[0.966d0,0.887d0,0.804d0,0.728d0]
  double precision                                     , dimension(massCount)                     :: ourLogConcentration
  class           (cosmologyFunctionsClass )                                            , pointer :: cosmologyFunctionsDefault
  integer                                                                                         :: iMass

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.Prada2011_concentration.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Prada2011 halo concentration algorithm")

  ! Test Prada2011 halo concentration algorithm.
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/Prada2011HaloConcentration/testParameters.xml'
  call Input_Parameters_File_Open(parameterFile)

  ! Create a node.
  thisNode => treeNode()

  ! Get the basic component.
  thisBasicComponent              => thisNode%basic(autoCreate=.true.)
  ! Get the default cosmology functions object.
  cosmologyFunctionsDefault       => cosmologyFunctions            ()
  ! Get the default concentrations object.
  darkMatterProfileConcentration_ => darkMatterProfileConcentration()

  ! Set the time for the node.
  call thisBasicComponent%timeSet(cosmologyFunctionsDefault%cosmicTime(1.00d0))

  ! Get the default cosmology.
  thisCosmologyParameters => cosmologyParameters()

  ! Loop over halo masses
  do iMass=1,massCount

     ! Set the mass of the original node.
     call thisBasicComponent%massSet(10.0d0**logMass(iMass)/thisCosmologyParameters%HubbleConstant(unitsLittleH))

     ! Compute and compare concentration at z=0.
     ourLogConcentration(iMass)=log10(darkMatterProfileConcentration_%concentration(thisNode))

  end do

  ! Check that results are as expected.
  message="Halo concentration at z=0"
  call Assert(char(message),ourLogConcentration,pradaLogConcentration,absTol=0.01d0)

  ! Close the input parameter file.
  call Input_Parameters_File_Close

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Prada2011_Concentration
