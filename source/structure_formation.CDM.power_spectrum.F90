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


!% Contains a module which implements the CDM power spectrum.

module CDM_Power_Spectrum
  !% Implements the CDM power spectrum.
  use FGSL
  use, intrinsic :: ISO_C_Binding                             
  private
  public :: Power_Spectrum_CDM, sigma_CDM, sigma_8, sigma_CDM_Logarithmic_Derivative, sigma_CDM_Plus_Logarithmic_Derivative,&
       & CDM_Power_Spectrum_State_Store, CDM_Power_Spectrum_State_Retrieve

  ! Flag to indicate if the power spectrum has been normalized.  
  logical                     :: sigmaInitialized   =.false.
  logical                     :: sigmaNormalized    =.false.
  logical                     :: normalizingSigma   =.false. ! Will be true during the normalization calculation.
  double precision, parameter :: radiusNormalization=8.0d0   ! Radius for sigma(M) normalization in Mpc/h.
  double precision            :: massNormalization           ! Mass for sigma(M) normalization in M_Solar.
  double precision            :: sigmaNormalization =1.0d0   ! Normalization for sigma(M).
  double precision            :: sigma_8_Value               ! Power spectrum normalization parameter.

  ! Variables to hold the tabulated sigma(M) data.
  integer                                        :: sigmaTableNPoints=-1
  double precision,    allocatable, dimension(:) :: sigmaTableLogMass, sigmaTable
  double precision                               :: sigmaTableLogMassMinimum=dlog(1.0d6), sigmaTableLogMassMaximum=dlog(1.0d15)
  integer,             parameter                 :: sigmaTableNPointsPerDecade=10
  type(fgsl_interp)                              :: interpolationObject
  type(fgsl_interp_accel)                        :: interpolationAccelerator
  logical                                        :: resetInterpolation=.true.

contains

  double precision function Power_Spectrum_CDM(wavenumber)
    !% Return the CDM power spectrum for $k=${\tt wavenumber} [Mpc$^{-1}$].
    use CDM_Transfer_Function
    use CDM_Primordial_Power_Spectrum
    use Numerical_Constants_Math
    use Cosmological_Parameters
    implicit none
    double precision, intent(in) :: wavenumber
    double precision             :: mass,logMass

    ! If this function is called not via the sigma(M) normalization routines, then ensure that sigma has been initialized so that
    ! we have the correct normalization.
    if (.not.normalizingSigma) then
       mass=(4.0d0*PI/3.0d0)*Omega_0()*Critical_Density()/waveNumber**3
       logMass=dlog(mass)
       call Initialize_Sigma(logMass)
    end if
    
    ! Compute the power spectrum.
    Power_Spectrum_CDM=(Transfer_Function_CDM(wavenumber)**2)*Primordial_Power_Spectrum_CDM(wavenumber)

    ! If this is not part of the normalization calculation, then scale by the normalization factor.
    if (.not.normalizingSigma) Power_Spectrum_CDM=Power_Spectrum_CDM*sigmaNormalization**2

    return
  end function Power_Spectrum_CDM

  double precision function sigma_CDM(mass)
    !% Computes the fractional mass fluctuation in real-space spherical top hats enclosing mass {\tt mass}.
    use Numerical_Interpolation
    implicit none
    double precision, intent(in) :: mass
    double precision             :: logMass

    ! Get logarithm of mass.
    logMass=dlog(mass)

    ! Check if we need to initialize this function.
    call Initialize_Sigma(logMass)
    
    ! Interpolate in tabulated function and return result.
    !$omp critical(Sigma_CDM_Interpolate)
    sigma_CDM=Interpolate(sigmaTableNPoints,sigmaTableLogMass,sigmaTable,interpolationObject,interpolationAccelerator,logMass&
         &,interpolationType=FGSL_Interp_CSpline,reset=resetInterpolation)
    !$omp end critical(Sigma_CDM_Interpolate)
    return
  end function sigma_CDM
  
  double precision function sigma_CDM_Logarithmic_Derivative(mass)
    !% Computes the logarithmic derivative in the fractional mass fluctuation in real-space spherical top hats enclosing mass {\tt
    !% mass}.
    use Numerical_Interpolation
    implicit none
    double precision, intent(in) :: mass
    double precision             :: sigma

    call sigma_CDM_Plus_Logarithmic_Derivative(mass,sigma,sigma_CDM_Logarithmic_Derivative)
    return
  end function sigma_CDM_Logarithmic_Derivative
  
  subroutine sigma_CDM_Plus_Logarithmic_Derivative(mass,sigma,sigmaLogarithmicDerivative)
    !% Returns both the fractional mass fluctuation in real-space spherical top hats enclosing mass {\tt mass} and its logarithmic derivative.
    use Numerical_Interpolation
    implicit none
    double precision, intent(in)  :: mass
    double precision, intent(out) :: sigma,sigmaLogarithmicDerivative
    double precision              :: logMass

    ! Get logarithm of mass.
    logMass=dlog(mass)

    ! Check if we need to initialize this function.
    call Initialize_Sigma(logMass)
    
    ! Interpolate in tabulated function and return result.
    !$omp critical(Sigma_CDM_Interpolate)
    sigma=Interpolate(sigmaTableNPoints,sigmaTableLogMass,sigmaTable,interpolationObject,interpolationAccelerator,logMass&
         &,interpolationType=FGSL_Interp_CSpline,reset=resetInterpolation)
    sigmaLogarithmicDerivative=Interpolate_Derivative(sigmaTableNPoints,sigmaTableLogMass,sigmaTable,interpolationObject&
         &,interpolationAccelerator,logMass,interpolationType=FGSL_Interp_CSpline,reset=resetInterpolation)/sigma
    !$omp end critical(Sigma_CDM_Interpolate)

    return
  end subroutine sigma_CDM_Plus_Logarithmic_Derivative

  double precision function sigma_8()
    !% Return the value of $\sigma_8$.
    implicit none

    ! Ensure the module has been initialized.
    call Initialize_Sigma(sigmaTableLogMassMinimum)

    sigma_8=sigma_8_Value
    return
  end function sigma_8

  subroutine Initialize_Sigma(logMass)
    !% Ensure that $\sigma(M)$ is tabulated over a range that includes {\tt logMass}. The default normalization, $\sigma_9=0.807$,
    !% is taken from \cite{komatsu_seven-year_2010}.
    use Input_Parameters
    use Memory_Management
    use Cosmological_Parameters
    use Numerical_Ranges
    use Numerical_Constants_Math
    use Numerical_Interpolation
    implicit none
    double precision, intent(in) :: logMass
    integer                      :: iMass
    double precision             :: mass
    
    !$omp critical (Sigma_CDM_Interpolate)
    if (.not.sigmaInitialized.or.logMass<sigmaTableLogMassMinimum.or.logMass>sigmaTableLogMassMaximum) then
       ! Compute the normalization if required.
       if (.not.sigmaNormalized) then
          !@ <inputParameter>
          !@   <name>sigma_8</name>
          !@   <defaultValue>0.807 \citep{komatsu_seven-year_2010}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The fractional mass fluctuation in the linear density field at the present day in spheres of radius 8~Mpc/h.
          !@   </description>
          !@ </inputParameter>
          call Get_Input_Parameter('sigma_8',sigma_8_Value,defaultValue=0.807d0)
          massNormalization=(4.0d0*PI/3.0d0)*Omega_0()*Critical_Density()*(radiusNormalization/Little_H_0())**3
          sigmaNormalization=sigma_8_Value/sigma_CDM_Integral(massNormalization)
          sigmaNormalized=.true.
       end if
       ! Find suitable range of masses to tabulate.
       sigmaTableLogMassMinimum=min(sigmaTableLogMassMinimum,logMass-ln10)
       sigmaTableLogMassMaximum=max(sigmaTableLogMassMaximum,logMass+ln10)
       sigmaTableNPoints=int((sigmaTableLogMassMaximum-sigmaTableLogMassMinimum)*dble(sigmaTableNPointsPerDecade)/ln10)
       ! Allocate arrays.
       if (allocated(sigmaTableLogMass)) call Dealloc_Array(sigmaTableLogMass)
       if (allocated(sigmaTable    ))    call Dealloc_Array(sigmaTable       )
       call Alloc_Array(sigmaTableLogMass,[sigmaTableNPoints])
       call Alloc_Array(sigmaTable       ,[sigmaTableNPoints])
       ! Generate a range of mass values to tabulate.
       sigmaTableLogMass=Make_Range(sigmaTableLogMassMinimum,sigmaTableLogMassMaximum,sigmaTableNPoints,rangeTypeLinear)
       ! Compute sigma(M) at each tabulated point.
       do iMass=1,sigmaTableNPoints
          mass=dexp(sigmaTableLogMass(iMass))
          sigmaTable(iMass)=sigma_CDM_Integral(mass)*sigmaNormalization
       end do
       ! Reset the interpolator.
       call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)
       resetInterpolation=.true.
       ! Flag that this module is now initialized.
       sigmaInitialized=.true.
    end if
    !$omp end critical (Sigma_CDM_Interpolate)
    return
  end subroutine Initialize_Sigma

  double precision function sigma_CDM_Integral(mass)
    !% Compute the root-variance of mass in (real space) top-hat spheres enclosing the given {\tt mass} from the power spectrum.
    use Numerical_Constants_Math
    use Numerical_Integration
    use Cosmological_Parameters
    implicit none
    double precision,                 intent(in) :: mass
    double precision,                 target     :: topHatRadius
    double precision                             :: wavenumberMinimum,wavenumberMaximum
    type(c_ptr)                                  :: parameterPointer
    type(fgsl_function)                          :: integrandFunction
    type(fgsl_integration_workspace)             :: integrationWorkspace

    parameterPointer=c_loc(topHatRadius)
    topHatRadius=((3.0d0/4.0d0/Pi)*mass/Omega_0()/Critical_Density())**(1.0d0/3.0d0)
    wavenumberMinimum=0.0d0/topHatRadius
    wavenumberMaximum=1.0d3/topHatRadius
    normalizingSigma=.true.
    sigma_CDM_Integral=Integrate(wavenumberMinimum,wavenumberMaximum,sigma_CDM_Integrand,parameterPointer,integrandFunction&
         &,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6)/2.0d0/Pi**2
    call Integrate_Done(integrandFunction,integrationWorkspace)
    normalizingSigma=.false.
    sigma_CDM_Integral=dsqrt(sigma_CDM_Integral)
    return
  end function sigma_CDM_Integral

  function sigma_CDM_Integrand(wavenumber,parameterPointer) bind(c)
    !% Integrand function used in compute the variance in (real space) top-hat spheres from the power spectrum.
    implicit none
    real(c_double)          :: sigma_CDM_Integrand
    real(c_double), value   :: wavenumber
    type(c_ptr),    value   :: parameterPointer
    real(c_double), pointer :: topHatRadius

    ! Extract integrand parameters.
    call c_f_pointer(parameterPointer,topHatRadius)
    ! Return power spectrum multiplied by window function and volume element in k-space. We don't include factors of 4 Pi here
    ! since this is unnormalized anyway.
    sigma_CDM_Integrand=Power_Spectrum_CDM(wavenumber)*(Window_Function_RealSpaceTopHat(wavenumber,topHatRadius)*wavenumber)**2
    return
  end function sigma_CDM_Integrand

  double precision function Window_Function_RealSpaceTopHat(wavenumber,topHatRadius)
    !% Top hat in real space window function Fourier transformed into $k$-space.
    implicit none
    double precision, intent(in) :: wavenumber,topHatRadius
    double precision, parameter  :: xSeriesMaximum=1.0d-3
    double precision             :: x,xSquared

    x=wavenumber*topHatRadius
    if      (x <= 0.0d0) then
       Window_Function_RealSpaceTopHat=0.0d0
    else if (x <= xSeriesMaximum) then 
       ! Use a series expansion of the window function for small x.
       xSquared=x**2
       Window_Function_RealSpaceTopHat=1.0d0 &
            & +xSquared*(-1.0d0/   10.0d0    &
            & +xSquared*(+1.0d0/  280.0d0    &
            & +xSquared*(-1.0d0/15120d0      &
            &           )))
    else
       ! For larger x, use the full expression.
       Window_Function_RealSpaceTopHat=3.0d0*(dsin(x)-x*dcos(x))/(x**3)
    end if
    return
  end function Window_Function_RealSpaceTopHat

  !# <galacticusStateStoreTask>
  !#  <unitName>CDM_Power_Spectrum_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine CDM_Power_Spectrum_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) sigmaTableLogMassMinimum,sigmaTableLogMassMaximum
    return
  end subroutine CDM_Power_Spectrum_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>CDM_Power_Spectrum_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine CDM_Power_Spectrum_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) sigmaTableLogMassMinimum,sigmaTableLogMassMaximum
    ! Flag that sigma tabulation needs to be reinitialized.
    sigmaInitialized=.false.
    return
  end subroutine CDM_Power_Spectrum_State_Retrieve
    
end module CDM_Power_Spectrum
