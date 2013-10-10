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

!% Contains a program which tests the \cite{navarro_structure_1996} halo concentration algorithm in a dark energy Universe.
!% Comparisons are made to the \href{http://www.astro.uvic.ca/~jfn/charden/charden.tar.gz}{``{\tt charden}''} code written by
!% Julio Navarro.

program Test_NFW96_Concentration_Dark_Energy
  !% Tests the \cite{navarro_structure_1996} halo concentration algorithm in a dark energy Universe. Comparisons are made to the
  !% \href{http://www.astro.uvic.ca/~jfn/charden/charden.tar.gz}{``{\tt charden}''} code written by Julio Navarro.
  use ISO_Varying_String
  use Memory_Management
  use Input_Parameters
  use Dark_Matter_Profiles_Concentrations
  use Cosmology_Functions
  use Cosmology_Parameters
  use Galacticus_Nodes
  use String_Handling
  use Unit_Tests
  implicit none
  type            (treeNode                )                         , pointer :: thisNode
  class           (nodeComponentBasic      )                         , pointer :: thisBasicComponent
  class           (cosmologyParametersClass)                         , pointer :: thisCosmologyParameters
  class           (cosmologyFunctionsClass )                         , pointer :: cosmologyFunctionsDefault
  integer                                   , dimension(6), parameter          :: chardenLogHaloMass       =[10,11,12,13,14,15]
  double precision                          , dimension(6), parameter          :: chardenConcentrationZ0   =[10.2700200d00,9.0204391d00,7.8041310d00,6.6154380d00,5.4956946d00,4.4538398d00], chardenConcentrationZ3=[5.8715897d00,5.4417138d00,5.0239682d00,4.6186433d00,4.2366042d00,3.8884208d00]
  type            (varying_string          )                                   :: message                                                                                                   , parameterFile
  integer                                                                      :: iMass
  double precision                                                             :: ourConcentration

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.NFW96_concentration.dark_energy.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("NFW96 halo concentration algorithm: dark energy cosmology")

  ! Test NFW96 halo concentration algorithm in a dark energy universe.
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/NFW96HaloConcentration/darkEnergy.xml'
  call Input_Parameters_File_Open(parameterFile)

  ! Create a node.
  thisNode => treeNode()

  ! Get the basic component.
  thisBasicComponent => thisNode%basic(autoCreate=.true.)

  ! Get the default cosmology.
  thisCosmologyParameters => cosmologyParameters()
  ! Get the default cosmology functions object.
  cosmologyFunctionsDefault => cosmologyFunctions()

  ! Loop over halo masses.
  do iMass=1,size(chardenLogHaloMass)

     ! Set the mass of the original node.
     call thisBasicComponent%massSet((10.0d0**chardenLogHaloMass(iMass))/thisCosmologyParameters%HubbleConstant(unitsLittleH))

     ! Compute and compare concentration at z=0.
     call thisBasicComponent%timeSet(cosmologyFunctionsDefault%cosmicTime(1.00d0))
     ourConcentration=Dark_Matter_Profile_Concentration(thisNode)
     message="10^"
     message=message//chardenLogHaloMass(iMass)//" M⊙/h halo concentration at z=0"
     call Assert(char(message),ourConcentration,chardenConcentrationZ0(iMass),relTol=0.02d0)

     ! Compute and compare concentration at z=3.
     call thisBasicComponent%timeSet(cosmologyFunctionsDefault%cosmicTime(0.25d0))
     ourConcentration=Dark_Matter_Profile_Concentration(thisNode)
     message="10^"
     message=message//chardenLogHaloMass(iMass)//" M⊙/h halo concentration at z=3"
     call Assert(char(message),ourConcentration,chardenConcentrationZ3(iMass),relTol=0.01d0)

  end do

  ! Close the input parameter file.
  call Input_Parameters_File_Close

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_NFW96_Concentration_Dark_Energy
