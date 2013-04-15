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

!% Contains a module which implements the nonlinear power spectrum using the algorithm of \cite{peacock_non-linear_1996}.

module Power_Spectra_Nonlinear_PeacockDodds1996
  !% Implements the nonlinear power spectrum using the algorithm of \cite{peacock_non-linear_1996}.
  implicit none
  private
  public :: Power_Spectrum_Nonlinear_PeacockDodds1996_Initialize

contains

  !# <powerSpectrumNonlinearMethod>
  !#  <unitName>Power_Spectrum_Nonlinear_PeacockDodds1996_Initialize</unitName>
  !# </powerSpectrumNonlinearMethod>
  subroutine Power_Spectrum_Nonlinear_PeacockDodds1996_Initialize(powerSpectrumNonlinearMethod,Power_Spectrum_Nonlinear_Get)
    !% Initializes the ``Peacock-Dodds1996'' nonlinear power spectrum module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in   ) :: powerSpectrumNonlinearMethod
    procedure(double precision), pointer, intent(inout) :: Power_Spectrum_Nonlinear_Get
    
    if (powerSpectrumNonlinearMethod == 'Peacock-Dodds1996') Power_Spectrum_Nonlinear_Get => Power_Spectrum_Nonlinear_PeacockDodds1996
    return
  end subroutine Power_Spectrum_Nonlinear_PeacockDodds1996_Initialize

  double precision function Power_Spectrum_Nonlinear_PeacockDodds1996(waveNumber,time)
    !% Return a nonlinear power spectrum equal using the algorithm of \cite{peacock_non-linear_1996}.
    use Numerical_Constants_Math
    use Power_Spectra
    use Cosmology_Functions
    use Galacticus_Error
    use Linear_Growth
    implicit none
    double precision, intent(in) :: waveNumber,time
    integer         , parameter  :: iterationCountMaximum=1000
    double precision, parameter  :: tolerance            =1.0d-3
    double precision, parameter  :: updateFraction       =0.1d0
    logical                      :: converged
    integer                      :: iterationCount
    double precision             :: fNL,fNLPrevious,waveNumberLinear,x,n,A,B,alpha,beta,V,g

    ! Make an initial guess that the nonlinear power spectrum equals the linear power spectrum.
    fNL        =Power_Spectrum_Dimensionless(waveNumber)
    fNLPrevious=fNL

    ! Iterate until a converged solution is found.
    converged     =.false.
    iterationCount=0
    do while (.not.converged .and. iterationCount < iterationCountMaximum)
       ! Find the corresponding linear wavenumber.
       waveNumberLinear=waveNumber/(1.0d0+fNL)**(1.0d0/3.0d0)
       ! Get the dimensionless linear power spectrum and its logarithmic slope.
       x=Power_Spectrum_Dimensionless         (      waveNumberLinear)*Linear_Growth_Factor(time)**2
       n=Power_Spectrum_Logarithmic_Derivative(0.5d0*waveNumberLinear)
       ! Compute parameters of the Peacock & Dodds fitting function.
       A    = 0.482d0/(1.0d0+n/3.0d0)**0.947d0
       B    = 0.226d0/(1.0d0+n/3.0d0)**1.778d0
       alpha= 3.310d0/(1.0d0+n/3.0d0)**0.244d0
       beta = 0.862d0/(1.0d0+n/3.0d0)**0.287d0
       V    =11.550d0/(1.0d0+n/3.0d0)**0.423d0
       ! Compute growth factor using same fitting function as Peacock & Dodds (from Carroll, Press & Turner 1992).
       g=2.5d0*Omega_Matter_Total(time)/(Omega_Matter_Total(time)**(4.0d0/7.0d0)-Omega_Dark_Energy(time)+(1.0d0+0.5d0&
            &*Omega_Matter_Total(time))*(1.0d0+Omega_Dark_Energy(time)/70.0d0))
       ! Compute new estimate of non-linear power-spectrum.
       fNL=    x                                             &
            & *(                                             &
            &    (1.0d0+B*beta*x+(A*x)**(alpha*beta))        &
            &   /(1.0d0+((A*x)**alpha*g**3/V/sqrt(x))**beta) &
            &  )**(1.0d0/beta)
       ! Update our estimate using a mixture of new and old results (this avoids oscillating solutions by preventing large initial
       ! jumps in the estimate).
       fNL=updateFraction*fNL+(1.0d0-updateFraction)*fNLPrevious
       ! Test for convergence.
       converged=(abs(fNL-fNLPrevious) < tolerance*0.5d0*(abs(fNL)+abs(fNLPrevious)))
       ! Move to next iteration.
       fNLPrevious=fNL
       iterationCount=iterationCount+1
    end do
    if (.not.converged) call Galacticus_Error_Report('Power_Spectrum_Nonlinear_PeacockDodds1996','nonlinear power spectrum calculation failed to converge')

    ! Convert to a dimensionful power spectrum.
    Power_Spectrum_Nonlinear_PeacockDodds1996=(2.0d0*Pi)**3*fNL/4.0d0/Pi/waveNumber**3
    return
  end function Power_Spectrum_Nonlinear_PeacockDodds1996
  
end module Power_Spectra_Nonlinear_PeacockDodds1996
