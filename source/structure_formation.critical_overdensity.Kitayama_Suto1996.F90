!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of critical overdensity using the fitting function of
!% \cite{kitayama_semianalytic_1996}.

module Critical_Overdensities_Kitayama_Suto1996
  implicit none
  private
  public :: Critical_Overdensity_Kitayama_Suto1996_Initialize

  ! Variables to hold the tabulated critical overdensity data.
  double precision            :: deltaTableTimeMaximum     =20.0d0, deltaTableTimeMinimum=1.0d0
  integer         , parameter :: deltaTableNPointsPerDecade=100

contains

  !# <criticalOverdensityMethod>
  !#  <unitName>Critical_Overdensity_Kitayama_Suto1996_Initialize</unitName>
  !# </criticalOverdensityMethod>
  subroutine Critical_Overdensity_Kitayama_Suto1996_Initialize(criticalOverdensityMethod,Critical_Overdensity_Contrast_Tabulate)
    !% Initializes the $\delta_{\rm c}$ calculation for the \cite{kitayama_semianalytic_1996} fitting function module.
    use ISO_Varying_String
    use Numerical_Comparison
    use Galacticus_Error
    use Cosmological_Parameters
   implicit none
    type     (varying_string                        ), intent(in   )          :: criticalOverdensityMethod
    procedure(Critical_Overdensity_Kitayama_Suto1996), intent(inout), pointer :: Critical_Overdensity_Contrast_Tabulate

    if (criticalOverdensityMethod == 'Kitayama-Suto1996') then
       Critical_Overdensity_Contrast_Tabulate => Critical_Overdensity_Kitayama_Suto1996
       ! Check that fitting formula is applicable to this cosmology.
       if (Values_Differ(Omega_Matter()+Omega_DE(),1.0d0,absTol=1.0d-6)) call Galacticus_Error_Report('Critical_Overdensity_Kitayama_Suto1996_Initialize','no fitting formula available for this cosmology')
    end if
    return
  end subroutine Critical_Overdensity_Kitayama_Suto1996_Initialize

  subroutine Critical_Overdensity_Kitayama_Suto1996(time,deltaTable)
    !% Tabulate the virial density contrast for the \cite{kitayama_semianalytic_1996} fitting function module.
    use Cosmology_Functions
    use Numerical_Constants_Math
    use Linear_Growth
    use Tables
    implicit none
    double precision                      , intent(in   ) :: time
    class           (table1D), allocatable, intent(inout) :: deltaTable
    integer                                               :: deltaTableNumberPoints, iTime

    ! Find minimum and maximum times to tabulate.
    deltaTableTimeMinimum=min(deltaTableTimeMinimum,time/2.0d0)
    deltaTableTimeMaximum=max(deltaTableTimeMaximum,time*2.0d0)

    ! Determine number of points to tabulate.
    deltaTableNumberPoints=int(log10(deltaTableTimeMaximum/deltaTableTimeMinimum)&
         &*dble(deltaTableNPointsPerDecade))

    ! Create the table.
    if (allocated(deltaTable)) then
       call deltaTable%destroy()
       deallocate(deltaTable)
    end if
    allocate(table1DLogarithmicLinear :: deltaTable)
    select type (deltaTable)
    type is (table1DLogarithmicLinear)
       ! Create the table.
       call deltaTable%create(deltaTableTimeMinimum,deltaTableTimeMaximum,deltaTableNumberPoints)
       ! Evaluate the fitting formula of Kitayama & Suto at each time to get the critical overdensity.
       do iTime=1,deltaTableNumberPoints
          call deltaTable%populate(                                                                  &
               &                   (3.0d0*(12.0d0*Pi)**(2.0d0/3.0d0)/20.0d0)                         &
               &                   *(1.0d0+0.0123d0*log10(Omega_Matter_Total(deltaTable%x(iTime))))  &
               &                   /Linear_Growth_Factor                    (deltaTable%x(iTime))  , &
               &                   iTime                                                             &
               &                  )
       end do
    end select
    return
  end subroutine Critical_Overdensity_Kitayama_Suto1996

end module Critical_Overdensities_Kitayama_Suto1996
