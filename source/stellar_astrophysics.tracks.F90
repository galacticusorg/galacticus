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


!% Contains a module which implements calculation of stellar tracks.

module Stellar_Astrophysics_Tracks
  !% Implements stellar tracks.
  use ISO_Varying_String
  implicit none
  private
  public :: Stellar_Luminosity, Stellar_Effective_Temperature
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: stellarTracksInitialized=.false.

  ! Name of cosmology functions method used.
  type(varying_string) :: stellarTracksMethod

  ! Pointer to the functions that actually do the calculations.
  procedure(Stellar_Tracks_Template), pointer :: Stellar_Luminosity_Get            => null()
  procedure(Stellar_Tracks_Template), pointer :: Stellar_Effective_Temperature_Get => null()

  abstract interface
     double precision function Stellar_Tracks_Template(initialMass,metallicity,age)
       double precision, intent(in) :: initialMass,metallicity,age
     end function Stellar_Tracks_Template
  end interface

contains

  subroutine Stellar_Tracks_Initialize
    !% Initialize the cosmology functions module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="stellarTracksMethod" type="moduleUse">
    include 'stellar_astrophysics.tracks.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Stellar_Tracks_Initialization) 
    ! Initialize if necessary.
    if (.not.stellarTracksInitialized) then
       ! Get the stellar tracks method parameter.
       !@ <inputParameter>
       !@   <name>stellarTracksMethod</name>
       !@   <defaultValue>file</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for stellar tracks calculations.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('stellarTracksMethod',stellarTracksMethod,defaultValue='file')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="stellarTracksMethod" type="code" action="subroutine">
       !#  <subroutineArgs>stellarTracksMethod,Stellar_Luminosity_Get,Stellar_Effective_Temperature_Get</subroutineArgs>
       include 'stellar_astrophysics.tracks.inc'
       !# </include>
       if (.not.(associated(Stellar_Luminosity_Get).and.associated(Stellar_Effective_Temperature_Get))) &
            & call Galacticus_Error_Report('Stellar_Tracks','method '//char(stellarTracksMethod)//' is unrecognized')
       stellarTracksInitialized=.true.
    end if
    !$omp end critical(Stellar_Tracks_Initialization) 

    return
  end subroutine Stellar_Tracks_Initialize

  double precision function Stellar_Luminosity(initialMass,metallicity,age)
    !% Returns the bolometric luminosity of a star of given {\tt initialMass}, {\tt metallicity} and {\tt age}.
    implicit none
    double precision, intent(in) :: initialMass,metallicity,age
    
    ! Initialize the module.
    call Stellar_Tracks_Initialize
    
    ! Get the answer using the selected method.
    Stellar_Luminosity=Stellar_Luminosity_Get(initialMass,metallicity,age)
    
    return
  end function Stellar_Luminosity
  
  double precision function Stellar_Effective_Temperature(initialMass,metallicity,age)
    !% Returns the effective temperature of a star of given {\tt initialMass}, {\tt metallicity} and {\tt age}.
    implicit none
    double precision, intent(in) :: initialMass,metallicity,age
    
    ! Initialize the module.
    call Stellar_Tracks_Initialize
    
    ! Get the answer using the selected method.
    Stellar_Effective_Temperature=Stellar_Effective_Temperature_Get(initialMass,metallicity,age)
    
    return
  end function Stellar_Effective_Temperature
  
end module Stellar_Astrophysics_Tracks
