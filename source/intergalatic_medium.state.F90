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


!% Contains a module that implements calculations of the intergalactic medium thermal and ionization state.

module Intergalactic_Medium_State
  !% Implements calculations of the intergalactic medium thermal and ionization state.
  private
  public :: Intergalactic_Medium_Electron_Fraction, Intergalactic_Medium_Temperature

  ! Flag to indicate if this module has been initialized.  
  logical :: igmStateInitialized=.false.

  ! Pointer to the function that actually does the calculation.
  procedure(Intergalactic_Medium_State_Get_Template), pointer :: Intergalactic_Medium_Electron_Fraction_Get => null()
  procedure(Intergalactic_Medium_State_Get_Template), pointer :: Intergalactic_Medium_Temperature_Get       => null()
  abstract interface
     double precision function Intergalactic_Medium_State_Get_Template(time)
       double precision, intent(in) :: time
     end function Intergalactic_Medium_State_Get_Template
  end interface
  
contains

  subroutine Intergalactic_Medium_State_Initialize
    !% Initialize the intergalactic medium state module.
    use ISO_Varying_String 
    use Galacticus_Error
    use Input_Parameters
    use Memory_Management
    !# <include directive="intergalaticMediumStateMethod" type="moduleUse">
    include 'intergalactic_medium.state.modules.inc'
    !# </include>
    implicit none
    type(varying_string) :: intergalaticMediumStateMethod
    
    !$omp critical(Intergalactic_Medium_State_Initialization) 
    ! Initialize if necessary.
    if (.not.igmStateInitialized) then
       ! Get the cooling function method parameter.
       !@ <inputParameter>
       !@   <name>intergalaticMediumStateMethod</name>
       !@   <defaultValue>RecFast</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for computing the state of the intergalactic medium.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('intergalaticMediumStateMethod',intergalaticMediumStateMethod,defaultValue='RecFast')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="intergalaticMediumStateMethod" type="code" action="subroutine">
       !#  <subroutineArgs>intergalaticMediumStateMethod,Intergalactic_Medium_Electron_Fraction_Get,Intergalactic_Medium_Temperature_Get</subroutineArgs>
       include 'intergalactic_medium.state.inc'
       !# </include>
       if (.not.(associated(Intergalactic_Medium_Electron_Fraction_Get).and.associated(Intergalactic_Medium_Temperature_Get))) call&
            & Galacticus_Error_Report('Intergalactic_Medium_State_Initialize','method ' //char(intergalaticMediumStateMethod)//' is unrecognized')

       igmStateInitialized=.true.
    end if
    !$omp end critical(Intergalactic_Medium_State_Initialization) 
    return
  end subroutine Intergalactic_Medium_State_Initialize

  double precision function Intergalactic_Medium_Electron_Fraction(time)
    !% Return the electron fraction in the intergalactic medium at the specified {\tt time}.
    implicit none
    double precision, intent(in) :: time

    ! Initialize the module.
    call Intergalactic_Medium_State_Initialize
  
    Intergalactic_Medium_Electron_Fraction=Intergalactic_Medium_Electron_Fraction_Get(time)
    return
  end function Intergalactic_Medium_Electron_Fraction
  
  double precision function Intergalactic_Medium_Temperature(time)
    !% Return the temperature of the intergalactic medium at the specified {\tt time}.
    implicit none
    double precision, intent(in) :: time

    ! Initialize the module.
    call Intergalactic_Medium_State_Initialize
  
    Intergalactic_Medium_Temperature=Intergalactic_Medium_Temperature_Get(time)
    return
  end function Intergalactic_Medium_Temperature
  
end module Intergalactic_Medium_State
