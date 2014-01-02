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

!% Contains a module which implements calculation of $\sigma(M)$ via filtering of the power spectrum.

module Cosmological_Mass_Variance_Filtered_Power_Spectrum
  !% Implements calculation of $\sigma(M)$ via filtering of the power spectrum.
  implicit none
  private
  public :: Cosmological_Mass_Variance_Filtered_Power_Spectrum_Initialize, &
       &    Cosmological_Mass_Variance_FPS_State_Store                   , &
       &    Cosmological_Mass_Variance_FPS_State_Retrieve

  ! Initial ranges and tabulation scale for the mass variance.
  double precision            :: sigmaTableMassMaximum     =1.0d15, sigmaTableMassMinimum=1.0d6
  integer         , parameter :: sigmaTableNPointsPerDecade=10

  ! Smoothing mass scale used in computing variance.
  double precision            :: smoothingMass

contains

  !# <cosmologicalMassVarianceMethod>
  !#  <unitName>Cosmological_Mass_Variance_Filtered_Power_Spectrum_Initialize</unitName>
  !# </cosmologicalMassVarianceMethod>
  subroutine Cosmological_Mass_Variance_Filtered_Power_Spectrum_Initialize(cosmologicalMassVarianceMethod&
       &,Cosmological_Mass_Variance_Tabulate)
    !% Initializes the $\sigma(M)$ calculation for the ``filtered power spectrum'' method.
    use ISO_Varying_String
    implicit none
    type     (varying_string                                             ), intent(in   )          :: cosmologicalMassVarianceMethod
    procedure(Cosmological_Mass_Variance_Filtered_Power_Spectrum_Tabulate), intent(inout), pointer :: Cosmological_Mass_Variance_Tabulate

    if (cosmologicalMassVarianceMethod == 'filteredPowerSpectrum') Cosmological_Mass_Variance_Tabulate =>&
         & Cosmological_Mass_Variance_Filtered_Power_Spectrum_Tabulate
    return
  end subroutine Cosmological_Mass_Variance_Filtered_Power_Spectrum_Initialize

  subroutine Cosmological_Mass_Variance_Filtered_Power_Spectrum_Tabulate(mass,massNormalization,sigmaNormalization,sigmaTable)
    !% Tabulate the virial density contrast for the \cite{kitayama_semianalytic_1996} fitting function module.
    use Tables
    double precision                      , intent(in   ) :: mass              , massNormalization
    double precision                      , intent(inout) :: sigmaNormalization
    class           (table1D), allocatable, intent(inout) :: sigmaTable
    integer                                               :: iMass             , sigmaTableNPoints
    double precision                                      :: sigma

    ! Create the table object.
    if (allocated(sigmaTable)) then
       call sigmaTable%destroy()
       deallocate(sigmaTable)
    end if
    allocate(table1DLogarithmicCSpline :: sigmaTable)
    select type (sigmaTable)
    type is (table1DLogarithmicCSpline)
       ! Determine the normalization of the power spectrum.
       sigmaNormalization=sigmaNormalization/Variance_Integral(massNormalization,useTopHat=.true.)
       ! Find suitable range of masses to tabulate.
       sigmaTableMassMinimum=min(sigmaTableMassMinimum,mass/10.0d0)
       sigmaTableMassMaximum=max(sigmaTableMassMaximum,mass*10.0d0)
       sigmaTableNPoints=int(log10(sigmaTableMassMaximum/sigmaTableMassMinimum)*dble(sigmaTableNPointsPerDecade))
       ! Allocate table grid.
       call sigmaTable%destroy()
       call sigmaTable%create(sigmaTableMassMinimum,sigmaTableMassMaximum,sigmaTableNPoints)
       ! Compute sigma(M) at each tabulated point.
       do iMass=1,sigmaTableNPoints
          sigma=Variance_Integral(sigmaTable%x(iMass),useTopHat=.false.)*sigmaNormalization
          ! Enforce monotonicity.
          if (iMass > 1) then
             if (sigma > sigmaTable%y(iMass-1)) sigma=sigmaTable%y(iMass-1)
          end if
          call sigmaTable%populate(sigma,iMass,computeSpline=(iMass == sigmaTableNPoints))
       end do
    end select
    return
  end subroutine Cosmological_Mass_Variance_Filtered_Power_Spectrum_Tabulate

  double precision function Variance_Integral(mass,useTopHat)
    !% Compute the root-variance of mass in spheres enclosing the given {\tt mass} from the power spectrum.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Math
    use Numerical_Integration
    use Cosmology_Parameters
    use Power_Spectrum_Window_Functions
    implicit none
    double precision                            , intent(in   ) :: mass
    logical                                     , intent(in   ) :: useTopHat
    class           (cosmologyParametersClass  ), pointer       :: thisCosmologyParameters
    double precision                                            :: topHatRadius           , wavenumberMaximum, &
         &                                                         wavenumberMinimum
    type            (c_ptr                     )                :: parameterPointer
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace

    ! Get the default cosmology.
    thisCosmologyParameters => cosmologyParameters()
    smoothingMass=mass
    topHatRadius=((3.0d0/4.0d0/Pi)*mass/thisCosmologyParameters%OmegaMatter()/thisCosmologyParameters%densityCritical())**(1.0d0/3.0d0)
    wavenumberMinimum=    0.0d0/topHatRadius
    wavenumberMaximum=min(1.0d3/topHatRadius,Power_Spectrum_Window_Function_Wavenumber_Maximum(smoothingMass))
    if (useTopHat) then
       Variance_Integral=Integrate(wavenumberMinimum,wavenumberMaximum,Variance_Integrand_TopHat,parameterPointer&
            &,integrandFunction ,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6,integrationRule=FGSL_Integ_Gauss15)/2.0d0/Pi**2
    else
       Variance_Integral=Integrate(wavenumberMinimum,wavenumberMaximum,Variance_Integrand,parameterPointer&
            &,integrandFunction ,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=4.0d-6,integrationRule=FGSL_Integ_Gauss15)/2.0d0/Pi**2
    end if
    call Integrate_Done(integrandFunction,integrationWorkspace)
    Variance_Integral=sqrt(Variance_Integral)
    return
  end function Variance_Integral

  function Variance_Integrand(wavenumber,parameterPointer) bind(c)
    !% Integrand function used in compute the variance in (real space) top-hat spheres from the power spectrum.
    use, intrinsic :: ISO_C_Binding
    use Power_Spectrum_Window_Functions
    use Primordial_Power_Spectra_Transferred
    implicit none
    real(kind=c_double)        :: Variance_Integrand
    real(kind=c_double), value :: wavenumber
    type(c_ptr        ), value :: parameterPointer

    ! Return power spectrum multiplied by window function and volume element in k-space. Factors of 2 and Pi are included
    ! elsewhere.
    Variance_Integrand=Primordial_Power_Spectrum_Transferred(wavenumber)*(Power_Spectrum_Window_Function(wavenumber,smoothingMass)*wavenumber)**2
    return
  end function Variance_Integrand

  function Variance_Integrand_TopHat(wavenumber,parameterPointer) bind(c)
    !% Integrand function used in compute the variance in (real space) top-hat spheres from the power spectrum.
    use, intrinsic :: ISO_C_Binding
    use Power_Spectrum_Window_Functions_Top_Hat
    use Primordial_Power_Spectra_Transferred
    implicit none
    real(kind=c_double)        :: Variance_Integrand_TopHat
    real(kind=c_double), value :: wavenumber
    type(c_ptr        ), value :: parameterPointer

    ! Return power spectrum multiplied by window function and volume element in k-space. Factors of 2 and Pi are included
    ! elsewhere.
    Variance_Integrand_TopHat=Primordial_Power_Spectrum_Transferred(wavenumber)*(Power_Spectrum_Window_Function_Top_Hat(wavenumber,smoothingMass)*wavenumber)**2
    return
  end function Variance_Integrand_TopHat

  !# <galacticusStateStoreTask>
  !#  <unitName>Cosmological_Mass_Variance_FPS_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Cosmological_Mass_Variance_FPS_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    write (stateFile) sigmaTableMassMinimum,sigmaTableMassMaximum
    return
  end subroutine Cosmological_Mass_Variance_FPS_State_Store

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Cosmological_Mass_Variance_FPS_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Cosmological_Mass_Variance_FPS_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) sigmaTableMassMinimum,sigmaTableMassMaximum
    return
  end subroutine Cosmological_Mass_Variance_FPS_State_Retrieve

end module Cosmological_Mass_Variance_Filtered_Power_Spectrum
