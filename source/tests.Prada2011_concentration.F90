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


!% Contains a program which tests the \cite{prada_halo_2011} halo concentration algorithm.

program Test_Prada2011_Concentration
  !% Tests the \cite{prada_halo_2011} halo concentration algorithm. Values of concentration were read from their Figure~12.
  use ISO_Varying_String
  use Memory_Management
  use Input_Parameters
  use Dark_Matter_Profiles_Concentrations
  use Cosmology_Functions
  use Cosmological_Parameters
  use Merger_Trees
  use Tree_Nodes
  use String_Handling
  use Unit_Tests
  implicit none
  type(mergerTree),     pointer                         :: thisTree
  type(treeNode),       pointer                         :: thisNode
  type(varying_string)                                  :: parameterFile,message
  integer             ,                       parameter :: massCount=4
  double precision    , dimension(massCount), parameter :: logMass              =[11.000d0,12.000d0,13.000d0,14.000d0]
  double precision    , dimension(massCount), parameter :: pradaLogConcentration=[ 0.966d0, 0.887d0, 0.804d0, 0.728d0]
  double precision    , dimension(massCount)            :: ourLogConcentration
  integer                                               :: iMass

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.Prada2011_concentration.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Prada2011 halo concentration algorithm")

  ! Test Prada2011 halo concentration algorithm.
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/Prada2011HaloConcentration/testParameters.xml'
  call Input_Parameters_File_Open(parameterFile)
  
  ! Create a node.
  call thisTree%createNode(thisNode)

  ! Set the time for the node.
  call Tree_Node_Time_Set(thisNode,Cosmology_Age(1.00d0))
  
  ! Loop over halo masses
  do iMass=1,massCount

     ! Set the mass of the original node.
     call Tree_Node_Mass_Set(thisNode,10.0d0**logMass(iMass)/Little_H_0())
     
     ! Compute and compare concentration at z=0.
     ourLogConcentration(iMass)=log10(Dark_Matter_Profile_Concentration(thisNode))

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
