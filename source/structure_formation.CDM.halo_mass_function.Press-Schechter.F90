!! Copyright 2009, Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which generates a tabulated Press-Schechter halo mass function.

module Halo_Mass_Function_Press_Schechter
  !% Implements generation of a tabulated power-law primordial power spectrum.
  private
  public :: Halo_Mass_Function_Press_Schechter_Initialize

  ! Parameters controlling the gridding of the power spectrum and default wavenumber range.
  integer,          parameter :: nPointsPerDecade=1000
  double precision            :: logMassMinimum=dlog(1.0d9), logMassMaximum=dlog(1.0d15)

contains
  
  !# <haloMassFunctionMethod>
  !#  <unitName>Halo_Mass_Function_Press_Schechter_Initialize</unitName>
  !# </haloMassFunctionMethod>
  subroutine Halo_Mass_Function_Press_Schechter_Initialize(haloMassFunctionMethod,Halo_Mass_Function_Tabulate)
    !% Initializes the ``Press-Schechter mass functon'' module.
    use ISO_Varying_String
    implicit none
    type(varying_string),          intent(in)    :: haloMassFunctionMethod
    procedure(),          pointer, intent(inout) :: Halo_Mass_Function_Tabulate
    
    if (haloMassFunctionMethod.eq.'Press-Schechter') Halo_Mass_Function_Tabulate => Halo_Mass_Function_Press_Schechter_Tabulate
    return
  end subroutine Halo_Mass_Function_Press_Schechter_Initialize

  subroutine Halo_Mass_Function_Press_Schechter_Tabulate(time,logMass,haloMassFunctionNumberPoints,haloMassFunctionLogMass &
       &,haloMassFunctionLogAbundance)
    !% Tabulate a Press-Schechter halo mass function.
    use Memory_Management
    use Numerical_Ranges
    use Numerical_Constants_Math
    use CDM_Power_Spectrum
    use Critical_Overdensity
    use Cosmological_Parameters
    implicit none
    double precision,                            intent(in)    :: time,logMass
    double precision, allocatable, dimension(:), intent(inout) :: haloMassFunctionLogMass,haloMassFunctionLogAbundance
    integer,                                     intent(out)   :: haloMassFunctionNumberPoints
    integer                                                    :: iMass
    double precision                                           :: mass,nu,alpha

    ! Determine range of masss required.
    logMassMinimum=min(logMassMinimum,logMass-ln10)
    logMassMaximum=max(logMassMaximum,logMass+ln10)
    
    ! Determine number of points to tabulate.
    haloMassFunctionNumberPoints=int((logMassMaximum-logMassMinimum)*dble(nPointsPerDecade)/ln10)

    ! Deallocate arrays if currently allocated.
    if (allocated(haloMassFunctionLogMass))      call Dealloc_Array(haloMassFunctionLogMass     )
    if (allocated(haloMassFunctionLogAbundance)) call Dealloc_Array(haloMassFunctionLogAbundance)
    ! Allocate the arrays to current required size.
    call Alloc_Array(haloMassFunctionLogMass     ,haloMassFunctionNumberPoints,'haloMassFunctionLogMass'     )
    call Alloc_Array(haloMassFunctionLogAbundance,haloMassFunctionNumberPoints,'haloMassFunctionLogAbundance')

    ! Tabulate the function.
    haloMassFunctionLogMass=Make_Range(logMassMinimum,logMassMaximum,haloMassFunctionNumberPoints,rangeTypeLinear)
    do iMass=1,haloMassFunctionNumberPoints
       mass=dexp(haloMassFunctionLogMass(iMass))
       nu=Critical_Overdensity_for_Collapse(time)/sigma_CDM(mass)
       alpha=dabs(sigma_CDM_Logarithmic_Derivative(mass))
       haloMassFunctionLogAbundance(iMass)=(Omega_0()*Critical_Density()/mass**2)*alpha*dsqrt(2.0d0/Pi)*nu&
            &*dexp(-0.5d0*nu**2)
    end do
    haloMassFunctionLogAbundance=dlog(haloMassFunctionLogAbundance)
    
    return
  end subroutine Halo_Mass_Function_Press_Schechter_Tabulate
  
end module Halo_Mass_Function_Press_Schechter
