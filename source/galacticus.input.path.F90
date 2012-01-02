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


!% Contains a module which provides the path for \glc\ inputs and scripts.

module Galacticus_Input_Paths
  !% Provides the path for \glc\ inputs and scripts.
  use ISO_Varying_String
  implicit none
  private
  public :: Galacticus_Input_Path

  ! The input path.
  type(varying_string) :: galacticusInputPath

  ! Flag indicating if we've retrieved the path yet.
  logical              :: pathRetrieved=.false.

contains

  function Galacticus_Input_Path()
    !% Returns the path that \glc\ should use as the root for all input data reads.
    implicit none
    type(varying_string) :: Galacticus_Input_Path
    integer              :: pathLength,pathStatus

    ! Check if the path has been retrieved.
    !$omp critical (Galacticus_Path_Initialize)
    if (.not.pathRetrieved) then
       ! Get the status and path length.
       call Get_Environment_Variable("GALACTICUS_ROOT_V091",length=pathLength,status=pathStatus) 
       if (pathStatus == 0) then
          ! Path is defined, retrieve it.
          call Get_Path(pathLength)
       else
          ! No path is defined, default to current working directory.
          galacticusInputPath="./"
       end if
       ! Flag that the path has been retrieved.
       pathRetrieved=.true.
    end if
    !$omp end critical (Galacticus_Path_Initialize)

    ! Return the input path.
    Galacticus_Input_Path=char(galacticusInputPath)
    return
  end function Galacticus_Input_Path
  
  subroutine Get_Path(pathLength)
    !% Retrieve the \glc\ input data path from the environment.
    implicit none
    integer,                  intent(in) :: pathLength
    character(len=pathLength+1)            :: pathName
    
    ! Get the path.
    call Get_Environment_Variable("GALACTICUS_ROOT_V091",value=pathName)
    ! Append a final "/" if necessary.
    if (pathName(pathLength:pathLength) /= "/") pathName=trim(pathName)//"/"
    ! Store the path.
    galacticusInputPath=pathName
    return
  end subroutine Get_Path

end module Galacticus_Input_Paths
