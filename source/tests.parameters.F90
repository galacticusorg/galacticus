!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a program which tests parameter input.

!: $(BUILDPATH)/tests.parameters.C.o

program Test_Parameters
  !% Test reading of input parameters.
  use, intrinsic :: ISO_C_Binding
  use Unit_Tests
  use IO_HDF5
  use ISO_Varying_String
  use Memory_Management
  use Input_Parameters2
  use Cosmology_Parameters
  implicit none
  type   (hdf5Object              )               :: outputFile
  type   (varying_string          )               :: parameterFile
  class  (cosmologyParametersClass), pointer      :: cosmologyParameters_
  type   (inputParameters         ), target       :: testParameters
  logical(kind=c_bool             ), dimension(3) :: results  
  ! Define an interface to the C-based portion of the testing code.
  interface
     subroutine testParametersC(results) bind(c,name='testParametersC')
       import c_bool
       logical(c_bool), intent(  out) :: results(3)
     end subroutine testParametersC
  end interface

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.parameters.size')
  ! Open an output file.
  call outputFile%openFile("testSuite/outputs/testParameters.hdf5",overWrite=.true.)
  parameterFile  ='testSuite/parameters/testsParameters.xml'
  testParameters=inputParameters(parameterFile,allowedParametersFile='tests.parameters.parameters.xml',outputParametersGroup=outputFile)
  call testParameters%markGlobal()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Parameter input")
  ! Test retrieval of cosmology parameters (simple).
  cosmologyParameters_ => cosmologyParameters()
  call Unit_Tests_Begin_Group("Retrieve cosmological parameters (simple)")
  call Assert('Ωₘ  ',cosmologyParameters_%OmegaMatter    (), 0.2725d0,relTol=1.0d-6)
  call Assert('Ωb  ',cosmologyParameters_%OmegaBaryon    (), 0.0455d0,relTol=1.0d-6)
  call Assert('ΩΛ  ',cosmologyParameters_%OmegaDarkEnergy(), 0.7275d0,relTol=1.0d-6)
  call Assert('H₀  ',cosmologyParameters_%HubbleConstant (),70.2000d0,relTol=1.0d-6)
  call Assert('TCMB',cosmologyParameters_%temperatureCMB (),2.72548d0,relTol=1.0d-6)
  call Unit_Tests_End_Group()
  ! Perform tests of the C interface to input parameter functions.
  call Unit_Tests_Begin_Group("C interface")
  call testParametersC(results)
  call Assert('Retrieve double'                   ,logical(results(1)),.true.)
  call Assert('Retrieve long'                     ,logical(results(2)),.true.)
  call Assert('Retrieve double from subparameters',logical(results(3)),.true.)
  call Unit_Tests_End_Group()
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  ! Close down.
  call testParameters%destroy()
  call outputFile     %close ()
end program Test_Parameters
