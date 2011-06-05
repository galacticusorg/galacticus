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


!% Contains a module which defines the radiation structure data type, used to describe radiation fields. (Currently only includes
!% CMB radiation temperature.)

module Radiation_Structure
  !% Defines the radiation structure data type, used to describe radiation fields. (Currently only includes
  !% CMB radiation temperature.)
  private
  public :: radiationStructure

  type radiationStructure
     !% The radiation structure data type, used to describe radiation fields. (Currently only includes
     !% CMB radiation temperature.)
     private
     double precision :: temperatureCosmicMicrowaveBackground
   contains
     procedure        :: set            => Radiation_Set
     procedure        :: setCMB         => Radiation_Set_CMB_Temperature
     procedure        :: temperatureCMB => Radiation_Get_CMB_Temperature
  end type radiationStructure

  ! Option labels.
  integer, public, parameter :: noRadiation=0

contains

  double precision function Radiation_Get_CMB_Temperature(radiation)
    !% Get the CMB temperature in a radiation structure.
    use Cosmology_Functions
#ifdef GCC45
    class(radiationStructure), intent(in) :: radiation
#else
    type(radiationStructure),  intent(in) :: radiation
#endif

#ifdef GCC45
    select type (radiation)
    type is (radiationStructure)
#endif
    Radiation_Get_CMB_Temperature=radiation%temperatureCosmicMicrowaveBackground
#ifdef GCC45
    end select
#endif
    return
  end function Radiation_Get_CMB_Temperature

  subroutine Radiation_Set_CMB_Temperature(radiation,cosmicTime)
    !% Set the CMB temperature in a radiation structure.
    use Cosmology_Functions
#ifdef GCC45
    class(radiationStructure), intent(inout) :: radiation
#else
    type(radiationStructure),  intent(inout) :: radiation
#endif
    double precision,          intent(in)    :: cosmicTime
    
#ifdef GCC45
    select type (radiation)
    type is (radiationStructure)
#endif
    radiation%temperatureCosmicMicrowaveBackground=CMB_Temperature(cosmicTime)
#ifdef GCC45
    end select
#endif
    return
  end subroutine Radiation_Set_CMB_Temperature

  subroutine Radiation_Set(radiation,setOption)
    !% Set the {\tt radiation} field as specified.
    implicit none
#ifdef GCC45
    class(radiationStructure), intent(inout)          :: radiation
#else
    type(radiationStructure),  intent(inout)          :: radiation
#endif
    integer,                   intent(in),   optional :: setOption

    ! AJB:: Currently does nothing as we don't support radiation structures yet.

    return
  end subroutine Radiation_Set

end module Radiation_Structure
