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


!% Contains a module which implements a barrier for excursion set calculations of dark matter halo formation which equals the
!% critical overdensity for collapse.

module Excursion_Sets_Barriers_Critical_Overdensity
  use FGSL
  !% Implements a barrier for excursion set calculations of dark matter halo formation which equals the critical overdensity for
  !% collapse.
  private
  public :: Excursion_Sets_Barriers_Critical_Overdensity_Initialize

contains

  !# <excursionSetBarrierMethod>
  !#  <unitName>Excursion_Sets_Barriers_Critical_Overdensity_Initialize</unitName>
  !# </excursionSetBarrierMethod>
  subroutine Excursion_Sets_Barriers_Critical_Overdensity_Initialize(excursionSetBarrierMethod,Excursion_Sets_Barrier_Get,Excursion_Sets_Barrier_Gradient_Get,barrierName)
    !% Initialize the critical overdensity excursion set barrier module.
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in)    :: excursionSetBarrierMethod
    procedure(double precision), pointer, intent(inout) :: Excursion_Sets_Barrier_Get,Excursion_Sets_Barrier_Gradient_Get
    type(varying_string),                 intent(inout) :: barrierName

    if (excursionSetBarrierMethod == 'criticalOverdensity') then
       Excursion_Sets_Barrier_Get          => Excursion_Sets_Barrier_Critical_Overdensity
       Excursion_Sets_Barrier_Gradient_Get => Excursion_Sets_Barrier_Gradient_Critical_Overdensity
       ! Construct a name for this barrier.
       barrierName=barrierName//":barrierCriticalOverdensity"
    end if
    return
  end subroutine Excursion_Sets_Barriers_Critical_Overdensity_Initialize

  double precision function Excursion_Sets_Barrier_Critical_Overdensity(variance,time)
    !% Return a critical overdensity barrier for excursion set calculations at the given {\tt variance}.
    use CDM_Power_Spectrum
    use Critical_Overdensity
    implicit none
    double precision, intent(in) :: variance,time
    double precision             :: mass

    if (variance <= 0.0d0) then
       ! Return the critical overdensity at this time for infinite mass.
       Excursion_Sets_Barrier_Critical_Overdensity=Critical_Overdensity_for_Collapse(time=time          )
    else
       ! Get the halo mass corresponding to this variance.
       mass=Mass_from_Sigma(dsqrt(variance))
       ! Return the critical overdensity at this time at the computed mass scale.
       Excursion_Sets_Barrier_Critical_Overdensity=Critical_Overdensity_for_Collapse(time=time,mass=mass)
    end if
    return
  end function Excursion_Sets_Barrier_Critical_Overdensity

  double precision function Excursion_Sets_Barrier_Gradient_Critical_Overdensity(variance,time)
    !% Return the gradient of a critical overdensity barrier for excursion set calculations at the given {\tt variance}.
    use CDM_Power_Spectrum
    use Critical_Overdensity
    implicit none
    double precision, intent(in) :: variance,time
    double precision             :: mass,alpha

    if (variance <= 0.0d0) then
       ! Return zero critical overdensity gradient at this time for infinite mass.
       Excursion_Sets_Barrier_Gradient_Critical_Overdensity=0.0d0
    else
       ! Get the halo mass corresponding to this variance.
       mass=Mass_from_Sigma(dsqrt(variance))
       ! Get the logarithmic slope of sigma(M).
       alpha=sigma_CDM_Logarithmic_Derivative(mass)
       ! Return the critical overdensity at this time at the computed mass scale.
       Excursion_Sets_Barrier_Gradient_Critical_Overdensity=(0.5d0*mass/variance/alpha)*Critical_Overdensity_for_Collapse(time=time,mass=mass)&
            &*Critical_Overdensity_Mass_Scaling_Gradient(mass)/Critical_Overdensity_Mass_Scaling(mass)
    end if
    return
  end function Excursion_Sets_Barrier_Gradient_Critical_Overdensity
  
end module Excursion_Sets_Barriers_Critical_Overdensity
