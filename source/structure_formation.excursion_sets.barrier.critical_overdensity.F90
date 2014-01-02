!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a barrier for excursion set calculations of dark matter halo formation which equals the
!% critical overdensity for collapse.

module Excursion_Sets_Barriers_Critical_Overdensity
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
    type     (varying_string  ), intent(in   )          :: excursionSetBarrierMethod
    procedure(Excursion_Sets_Barrier_Critical_Overdensity), intent(inout), pointer :: Excursion_Sets_Barrier_Get
    procedure(Excursion_Sets_Barrier_Gradient_Critical_Overdensity), intent(inout), pointer :: Excursion_Sets_Barrier_Gradient_Get
    type     (varying_string  ), intent(inout)          :: barrierName

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
    use Power_Spectra
    use Critical_Overdensity
    implicit none
    double precision, intent(in   ) :: time, variance
    double precision                :: mass

    if (variance <= 0.0d0) then
       ! Return the critical overdensity at this time for infinite mass.
       Excursion_Sets_Barrier_Critical_Overdensity=Critical_Overdensity_for_Collapse(time=time          )
    else
       ! Get the halo mass corresponding to this variance.
       mass=Mass_from_Cosmolgical_Root_Variance(sqrt(variance))
       ! Return the critical overdensity at this time at the computed mass scale.
       Excursion_Sets_Barrier_Critical_Overdensity=Critical_Overdensity_for_Collapse(time=time,mass=mass)
    end if
    return
  end function Excursion_Sets_Barrier_Critical_Overdensity

  double precision function Excursion_Sets_Barrier_Gradient_Critical_Overdensity(variance,time)
    !% Return the gradient of a critical overdensity barrier for excursion set calculations at the given {\tt variance}.
    use Power_Spectra
    use Critical_Overdensity
    implicit none
    double precision, intent(in   ) :: time , variance
    double precision                :: alpha, mass

    if (variance <= 0.0d0) then
       ! Return zero critical overdensity gradient at this time for infinite mass.
       Excursion_Sets_Barrier_Gradient_Critical_Overdensity=0.0d0
    else
       ! Get the halo mass corresponding to this variance.
       mass=Mass_from_Cosmolgical_Root_Variance(sqrt(variance))
       ! Get the logarithmic slope of sigma(M).
       alpha=Cosmological_Mass_Root_Variance_Logarithmic_Derivative(mass)
       ! Return the critical overdensity at this time at the computed mass scale.
       Excursion_Sets_Barrier_Gradient_Critical_Overdensity=(0.5d0*mass/variance/alpha)*Critical_Overdensity_for_Collapse(time=time,mass=mass)&
            &*Critical_Overdensity_Mass_Scaling_Gradient(mass)/Critical_Overdensity_Mass_Scaling(mass)
    end if
    return
  end function Excursion_Sets_Barrier_Gradient_Critical_Overdensity

end module Excursion_Sets_Barriers_Critical_Overdensity
