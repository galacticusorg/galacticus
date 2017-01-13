!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
!!    Andrew Benson <abenson@carnegiescience.edu>
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
    !% Return a critical overdensity barrier for excursion set calculations at the given {\normalfont \ttfamily variance}.
    use Cosmological_Mass_Variance
    use Critical_Overdensities
    implicit none
    double precision                               , intent(in   ) :: time                , variance
    class           (criticalOverdensityClass     ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    double precision                                               :: mass

    ! Get default objects.
    criticalOverdensity_      => criticalOverdensity     ()
    cosmologicalMassVariance_ => cosmologicalMassVariance()
    if (variance <= 0.0d0) then
       ! Return the critical overdensity at this time for infinite mass.
       Excursion_Sets_Barrier_Critical_Overdensity=criticalOverdensity_%value(time=time          )
    else
       ! Get the halo mass corresponding to this variance.
       mass=cosmologicalMassVariance_%mass(sqrt(variance))
       ! Return the critical overdensity at this time at the computed mass scale.
       Excursion_Sets_Barrier_Critical_Overdensity=criticalOverdensity_%value(time=time,mass=mass)
    end if
    return
  end function Excursion_Sets_Barrier_Critical_Overdensity

  double precision function Excursion_Sets_Barrier_Gradient_Critical_Overdensity(variance,time)
    !% Return the gradient of a critical overdensity barrier for excursion set calculations at the given {\normalfont \ttfamily variance}.
    use Cosmological_Mass_Variance
    use Critical_Overdensities
    implicit none
    double precision                               , intent(in   ) :: time                , variance
    class           (criticalOverdensityClass     ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    double precision                                               :: alpha               , mass

    if (variance <= 0.0d0) then
       ! Return zero critical overdensity gradient at this time for infinite mass.
       Excursion_Sets_Barrier_Gradient_Critical_Overdensity=0.0d0
    else
       ! Get default objects.
       criticalOverdensity_      => criticalOverdensity     ()
       cosmologicalMassVariance_ => cosmologicalMassVariance()
       ! Get the halo mass corresponding to this variance.
       mass=cosmologicalMassVariance_%mass(sqrt(variance))
       ! Get the logarithmic slope of sigma(M).
       alpha=cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass)
       ! Return the critical overdensity at this time at the computed mass scale.
       Excursion_Sets_Barrier_Gradient_Critical_Overdensity=+0.5d0                                                  &
            &                                               *mass                                                   &
            &                                               /variance                                               &
            &                                               /alpha                                                  &
            &                                               *criticalOverdensity_%gradientMass(time=time,mass=mass)
    end if
    return
  end function Excursion_Sets_Barrier_Gradient_Critical_Overdensity

end module Excursion_Sets_Barriers_Critical_Overdensity
