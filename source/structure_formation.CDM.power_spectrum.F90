!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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








!% Contains a module which implements the CDM power spectrum.

module CDM_Power_Spectrum
  !% Implements the CDM power spectrum.
  use FGSL
  use, intrinsic :: ISO_C_Binding                             
  private
  public :: Power_Spectrum_CDM, sigma_CDM, sigma_CDM_Logarithmic_Derivative, sigma_CDM_Plus_Logarithmic_Derivative,&
       & CDM_Power_Spectrum_State_Store, CDM_Power_Spectrum_State_Retrieve

  ! Flag to indicate if the power spectrum has been normalized.  
  logical                     :: sigmaInitialized   =.false.
  logical                     :: sigmaNormalized    =.false.
  double precision, parameter :: radiusNormalization=8.0d0 ! Radius for sigma(M) normalization in Mpc/h.
  double precision            :: massNormalization         ! Mass for sigma(M) normalization in M_Solar.
  double precision            :: sigmaNormalization =1.0d0 ! Normalization for sigma(M).
  double precision            :: sigma_8                   ! Power spectrum normalization parameter.

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
    implicit none
    double precision, intent(in) :: wavenumber

    ! Compute the power spectrum.
    Power_Spectrum_CDM=(Transfer_Function_CDM(wavenumber)**2)*Primordial_Power_Spectrum_CDM(wavenumber)

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
          call Get_Input_Parameter('sigma_8',sigma_8,defaultValue=0.807d0)
          massNormalization=(4.0d0*PI/3.0d0)*Omega_0()*Critical_Density()*(radiusNormalization/Little_H_0())**3
          sigmaNormalization=sigma_8/sigma_CDM_Integral(massNormalization)
          sigmaNormalized=.true.
       end if
       ! Find suitable range of masses to tabulate.
       sigmaTableLogMassMinimum=min(sigmaTableLogMassMinimum,logMass-ln10)
       sigmaTableLogMassMaximum=max(sigmaTableLogMassMaximum,logMass+ln10)
       sigmaTableNPoints=int((sigmaTableLogMassMaximum-sigmaTableLogMassMinimum)*dble(sigmaTableNPointsPerDecade)/ln10)
       ! Allocate arrays.
       if (allocated(sigmaTableLogMass)) call Dealloc_Array(sigmaTableLogMass)
       if (allocated(sigmaTable    ))    call Dealloc_Array(sigmaTable       )
       call Alloc_Array(sigmaTableLogMass,sigmaTableNPoints,'sigmaTableLogMass')
       call Alloc_Array(sigmaTable       ,sigmaTableNPoints,'sigmaTable'       )
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
    sigma_CDM_Integral=Integrate(wavenumberMinimum,wavenumberMaximum,sigma_CDM_Integrand,parameterPointer,integrandFunction&
         &,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    sigma_CDM_Integral=dsqrt(sigma_CDM_Integral)
    return
  end function sigma_CDM_Integral

  function sigma_CDM_Integrand(wavenumber,parameterPointer) bind(c)
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
    double precision             :: x

    x=wavenumber*topHatRadius
    if (x<=0.0d0) then
       Window_Function_RealSpaceTopHat=0.0d0
    else
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
