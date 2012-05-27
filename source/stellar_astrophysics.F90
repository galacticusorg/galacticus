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


!% Contains a module which implements calculation of stellar astrophysics.

module Stellar_Astrophysics
  !% Implements calculation of stellar astrophysics.
  use ISO_Varying_String
  implicit none
  private
  public :: Star_Ejected_Mass, Star_Initial_Mass, Star_Metal_Yield_Mass, Star_Lifetime
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: stellarAstrophysicsInitialized=.false.

  ! Name of cosmology functions method used.
  type(varying_string) :: stellarAstrophysicsMethod

  ! Pointer to the functions that actually do the calculations.
  procedure(Stellar_Astrophysics_Template),       pointer :: Star_Ejected_Mass_Get     => null()
  procedure(Stellar_Astrophysics_Template),       pointer :: Star_Initial_Mass_Get     => null()
  procedure(Stellar_Astrophysics_Yield_Template), pointer :: Star_Metal_Yield_Mass_Get => null()
  procedure(Stellar_Astrophysics_Template),       pointer :: Star_Lifetime_Get         => null()

  abstract interface
     double precision function Stellar_Astrophysics_Template(inputParameter1,inputParameter2)
       double precision, intent(in) :: inputParameter1,inputParameter2
     end function Stellar_Astrophysics_Template
  end interface

  abstract interface
     double precision function Stellar_Astrophysics_Yield_Template(inputParameter1,inputParameter2,atomIndex)
       double precision, intent(in)           :: inputParameter1,inputParameter2
       integer,          intent(in), optional :: atomIndex
     end function Stellar_Astrophysics_Yield_Template
  end interface

contains

  subroutine Stellar_Astrophysics_Initialize
    !% Initialize the stellar astrophysics module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="stellarAstrophysicsMethod" type="moduleUse">
    include 'stellar_astrophysics.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.stellarAstrophysicsInitialized) then
       !$omp critical(Stellar_Astrophysics_Initialization) 
       if (.not.stellarAstrophysicsInitialized) then
          ! Get the stellar tracks method parameter.
          !@ <inputParameter>
          !@   <name>stellarAstrophysicsMethod</name>
          !@   <defaultValue>file</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for stellar astrophysics calculations.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('stellarAstrophysicsMethod',stellarAstrophysicsMethod,defaultValue='file')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="stellarAstrophysicsMethod" type="code" action="subroutine">
          !#  <subroutineArgs>stellarAstrophysicsMethod,Star_Ejected_Mass_Get,Star_Initial_Mass_Get,Star_Metal_Yield_Mass_Get,Star_Lifetime_Get</subroutineArgs>
          include 'stellar_astrophysics.inc'
          !# </include>
          if (.not.(associated(Star_Ejected_Mass_Get).and.associated(Star_Initial_Mass_Get).and.associated(Star_Metal_Yield_Mass_Get).and.associated(Star_Lifetime_Get))) &
               & call Galacticus_Error_Report('Stellar_Astrophysics','method '//char(stellarAstrophysicsMethod)//' is unrecognized')
          stellarAstrophysicsInitialized=.true.
       end if
       !$omp end critical(Stellar_Astrophysics_Initialization) 
    end if
    return
  end subroutine Stellar_Astrophysics_Initialize

  double precision function Star_Initial_Mass(lifetime,metallicity)
    !% Returns the initial mass of a star of given {\tt lifetime} and {\tt metallicity}.
    implicit none
    double precision, intent(in) :: lifetime,metallicity
    
    ! Initialize the module.
    call Stellar_Astrophysics_Initialize
    
    ! Get the answer using the selected method.
    Star_Initial_Mass=Star_Initial_Mass_Get(lifetime,metallicity)
    
    return
  end function Star_Initial_Mass
  
  double precision function Star_Ejected_Mass(initialMass,metallicity)
    !% Returns the mass ejected by a star of given {\tt initialMass} and {\tt metallicity}.
    implicit none
    double precision, intent(in) :: initialMass,metallicity
    
    ! Initialize the module.
    call Stellar_Astrophysics_Initialize
    
    ! Get the answer using the selected method.
    Star_Ejected_Mass=Star_Ejected_Mass_Get(initialMass,metallicity)
    
    return
  end function Star_Ejected_Mass
  
  double precision function Star_Metal_Yield_Mass(initialMass,metallicity,atomIndex)
    !% Returns the metal mass yielded by a star of given {\tt initialMass} and {\tt metallicity}.
    implicit none
    double precision, intent(in)           :: initialMass,metallicity
    integer,          intent(in), optional :: atomIndex
    
    ! Initialize the module.
    call Stellar_Astrophysics_Initialize
    
    ! Get the answer using the selected method.
    Star_Metal_Yield_Mass=Star_Metal_Yield_Mass_Get(initialMass,metallicity,atomIndex)
    
    return
  end function Star_Metal_Yield_Mass
  
  double precision function Star_Lifetime(initialMass,metallicity)
    !% Returns the lifetime of a star of given {\tt initialMass} and {\tt metallicity}.
    implicit none
    double precision, intent(in) :: initialMass,metallicity
    
    ! Initialize the module.
    call Stellar_Astrophysics_Initialize
    
    ! Get the answer using the selected method.
    Star_Lifetime=Star_Lifetime_Get(initialMass,metallicity)
    
    return
  end function Star_Lifetime
  
end module Stellar_Astrophysics
