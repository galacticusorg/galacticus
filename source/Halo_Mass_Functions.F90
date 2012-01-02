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


!% Contains a code which computes dark matter halo mass functions and associated data.

program Halo_Mass_Functions
  !% Computes dark matter halo mass functions and associated data.
  use Galacticus_Error
  use ISO_Varying_String
  use Input_Parameters
  use Memory_Management
  use Halo_Mass_Function_Tasks
  implicit none
  integer,                             parameter :: fileNameLengthMaximum=1024
  character(len=fileNameLengthMaximum)           :: fileCharacter
  type(varying_string)                           :: parameterFile,haloMassFunctionOutputFileName

  ! Read in basic code memory usage.
  call Code_Memory_Usage('Halo_Mass_Functions.size')

  ! Get the name of the parameter file from the first command line argument.
  call Get_Command_Argument(1,fileCharacter)
  if (len_trim(fileCharacter) == 0) call Galacticus_Error_Report(message="Usage: Halo_Mass_Function.exe <parameterFile> <outputFile>")
  parameterFile=fileCharacter
  call Get_Command_Argument(2,fileCharacter)
  if (len_trim(fileCharacter) == 0) call Galacticus_Error_Report(message="Usage: Halo_Mass_Function.exe <parameterFile> <outputFile>")
  haloMassFunctionOutputFileName=fileCharacter

  ! Open the output file.
  call Halo_Mass_Function_Open_File(haloMassFunctionOutputFileName)
  
  ! Open the parameter file.
  call Input_Parameters_File_Open(parameterFile,haloMassFunctionOutputFile)

  ! Compute the mass function.
  call Halo_Mass_Function_Compute

  ! Output the mass function.
  call Halo_Mass_Function_Output

  ! Close the parameter file.
  call Input_Parameters_File_Close
  
  ! Close the output file.
  call Halo_Mass_Function_Close_File

  ! All done, finish.
end program Halo_Mass_Functions
