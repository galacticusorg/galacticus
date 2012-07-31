!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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

!% Contains a module which provides estimates of the time taken to evolve a merger tree in \glc.

module Galacticus_Meta_Compute_Times_File
  !% Provides estimates of the time taken to evolve a merger tree in \glc.
  private
  public :: Galacticus_Time_Per_Tree_File_Initialize

  ! Coefficients of the fitting function.
  double precision, dimension(0:2) :: fitCoefficient

contains

  !# <timePerTreeMethod>
  !#  <unitName>Galacticus_Time_Per_Tree_File_Initialize</unitName>
  !# </timePerTreeMethod>
  subroutine Galacticus_Time_Per_Tree_File_Initialize(timePerTreeMethod,Galacticus_Time_Per_Tree_Get)
    !% Initializes the ``file'' time per tree module.
    use ISO_Varying_String
    use Input_Parameters
    use FoX_DOM
    use Galacticus_Error
    implicit none
    type(varying_string),                 intent(in)    :: timePerTreeMethod
    procedure(double precision), pointer, intent(inout) :: Galacticus_Time_Per_Tree_Get
    type(Node),                  pointer                :: doc,thisFit,thisCoefficient
    type(NodeList),              pointer                :: fitList,coefficientList
    integer                                             :: iTerm,ioErr
    double precision                                    :: coefficient(1)
    type(varying_string)                                :: timePerTreeFitFileName

    if (timePerTreeMethod == 'file') then
       Galacticus_Time_Per_Tree_Get => Galacticus_Time_Per_Tree_File

       ! Get the name of the file containing the fit coefficients.
       !@ <inputParameter>
       !@   <name>timePerTreeFitFileName</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The name of the file which contains fit coefficients for the time per tree fitting function.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('timePerTreeFitFileName',timePerTreeFitFileName)

       ! Parse the fit file.
       !$omp critical (FoX_DOM_Access)
       doc => parseFile(char(timePerTreeFitFileName),iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Galacticus_Time_Per_Tree_File_Initialize','Unable to find or parse tree timing file')
       fitList         => getElementsByTagname(doc,"fit")
       thisFit         => item(fitList,0)
       coefficientList => getElementsByTagname(thisFit,"coefficient")
       do iTerm=0,2
          thisCoefficient => item(coefficientList,iTerm)
          call extractDataContent(thisCoefficient,coefficient)
          fitCoefficient(iTerm)=coefficient(1)
       end do
       !$omp end critical (FoX_DOM_Access)
    end if
    return
  end subroutine Galacticus_Time_Per_Tree_File_Initialize

  double precision function Galacticus_Time_Per_Tree_File(treeRootMass)
    !% Provides estimates of the time taken to evolve a merger tree in \glc.
    implicit none
    double precision,     intent(in) :: treeRootMass
    integer                          :: iTerm
    double precision                 :: logTreeRootMass

    ! Find the logarithm of the tree mass.
    logTreeRootMass=dlog10(treeRootMass)
    
    Galacticus_Time_Per_Tree_File=0.0
    do iTerm=0,2
       Galacticus_Time_Per_Tree_File=Galacticus_Time_Per_Tree_File+fitCoefficient(iTerm)*logTreeRootMass**iTerm
    end do
    Galacticus_Time_Per_Tree_File=10.0d0**Galacticus_Time_Per_Tree_File

    return
  end function Galacticus_Time_Per_Tree_File
  
end module Galacticus_Meta_Compute_Times_File
