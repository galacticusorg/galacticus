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

program Tests_Comoving_Distance_Dark_Energy
  !% Tests comoving distance calculations for a dark energy Universe. Distances calculated using Python implementation of Ned Wright's
  !% cosmology calculator available from: http://www.astro.ucla.edu/~wright/CC.python
  use Unit_Tests
  use Input_Parameters
  use ISO_Varying_String
  use Cosmology_Functions
  use Memory_Management
  implicit none
  double precision, parameter, dimension(8) :: redshift          =[0.1d0,1.0d0,3.0d0,9.0d0,30.0d0,100.0d0,300.0d0,1000.0d0]
  double precision, parameter, dimension(8) :: distanceDarkEnergy=[2.9291813d0,23.1267935d0,44.4897529d0,64.4722404d0,79.4221370d0,88.1893579d0,92.7669523d0,95.5884181d0]
  type(varying_string)                      :: parameterFile
  character(len=1024)                       :: message
  integer                                   :: iExpansion
  double precision                          :: time,timeLookup,distance

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.comoving_distance.dark_energy.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Comoving distance: dark energy cosmology")

  ! Test comoving distance in a dark energy universe.
  parameterFile='testSuite/parameters/comovingDistance/darkEnergy.xml'
  call Input_Parameters_File_Open(parameterFile)
  do iExpansion=1,size(redshift)
     time=Cosmology_Age(Expansion_Factor_From_Redshift(redshift(iExpansion)))
     distance=Comoving_Distance(time)
     write (message,'(a,f6.1,a)') "comoving distance [z=",redshift(iExpansion),"]"
     call Assert(trim(message),distance,distanceDarkEnergy(iExpansion),relTol=1.0d-3)
     timeLookup=Time_From_Comoving_Distance(distance)
     write (message,'(a,f6.1,a)') "cosmic time [distance=",distance," Mpc]"
     call Assert(trim(message),timeLookup,time,relTol=1.0d-3)
  end do
  call Input_Parameters_File_Close

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Comoving_Distance_Dark_Energy
