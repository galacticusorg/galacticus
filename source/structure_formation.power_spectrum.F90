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

!% Contains a module which implements the cosmological power spectrum.

module Power_Spectra
  !% Implements the cosmological power spectrum.
  use Tables
  use ISO_Varying_String
  implicit none
  private
  public :: Power_Spectrum, Power_Spectrum_Logarithmic_Derivative, Power_Spectrum_Dimensionless, Cosmological_Mass_Root_Variance, Mass_from_Cosmolgical_Root_Variance,&
       & sigma_8, Cosmological_Mass_Root_Variance_Logarithmic_Derivative, Cosmological_Mass_Root_Variance_Plus_Logarithmic_Derivative

  ! Flag to indicate if the power spectrum has been normalized.  
  logical                                       :: sigmaInitialized   =.false.
  logical                                       :: sigmaNormalized    =.false.
  double precision                , parameter   :: radiusNormalization=8.0d0   ! Radius for sigma(M) normalization in Mpc/h.
  double precision                              :: massNormalization           ! Mass for sigma(M) normalization in M_Solar.
  double precision                              :: sigmaNormalization =1.0d0   ! Normalization for sigma(M).
  double precision                              :: sigma_8_Value               ! Power spectrum normalization parameter.

  ! Variables to hold the tabulated sigma(M) data.
  class           (table1D       ), allocatable :: sigmaTable

  ! Name of mass variance method used.
  type            (varying_string)              :: cosmologicalMassVarianceMethod

  ! Pointer to the subroutine that tabulates the virial overdensity and template interface for that subroutine.
  procedure(Cosmological_Mass_Variance_Tabulate_Template), pointer :: Cosmological_Mass_Variance_Tabulate => null()
  abstract interface
     subroutine Cosmological_Mass_Variance_Tabulate_Template(mass,massNormalization,sigmaNormalization,sigmaTable)
       import table1D
       double precision         , intent(in   )              :: mass,massNormalization
       double precision         , intent(  out)              :: sigmaNormalization
       class           (table1D), intent(inout), allocatable :: sigmaTable
     end subroutine Cosmological_Mass_Variance_Tabulate_Template
  end interface
  
contains

  double precision function Power_Spectrum_Dimensionless(wavenumber)
    !% Return the dimensionless power spectrum, $\Delta^2(k)$, for $k=${\tt wavenumber} [Mpc$^{-1}$].
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in) :: wavenumber

    Power_Spectrum_Dimensionless=4.0d0*Pi*wavenumber**3*Power_Spectrum(wavenumber)/(2.0d0*Pi)**3
    return
  end function Power_Spectrum_Dimensionless

  double precision function Power_Spectrum_Logarithmic_Derivative(wavenumber)
    !% Return the logarithmic derivative of the power spectrum, ${\rm d}\ln P(k)/{\rm d}\ln k$, for $k=${\tt wavenumber} [Mpc$^{-1}$].
    use Transfer_Functions
    use Primordial_Power_Spectra
    implicit none
    double precision, intent(in) :: wavenumber

    Power_Spectrum_Logarithmic_Derivative=                                                                    &
         &                                       Primordial_Power_Spectrum_Logarithmic_Derivative(wavenumber) &
         &                                +2.0d0*        Transfer_Function_Logarithmic_Derivative(wavenumber)
    return
  end function Power_Spectrum_Logarithmic_Derivative

  double precision function Power_Spectrum(wavenumber)
    !% Return the cosmological power spectrum for $k=${\tt wavenumber} [Mpc$^{-1}$].
    use Primordial_Power_Spectra_Transferred
    use Numerical_Constants_Math
    use Cosmological_Parameters
    use Galacticus_Error
    use ISO_Varying_String
    implicit none
    double precision                , intent(in   ) :: wavenumber
    double precision                                :: mass
    character       (len=15        )                :: label
    type            (varying_string)                :: message

    ! Ensure that the normalization of the power spectrum has been computed.
    mass=(4.0d0*PI/3.0d0)*Omega_Matter()*Critical_Density()/waveNumber**3
    if (mass <= 0.0d0) then
       message="zero mass when trying to initialize power spectrum"//char(10)
       write (label,'(e12.6)') waveNumber
       message=message//"        waveNumber  : "//trim(label)//char(10)
       write (label,'(e12.6)') Omega_Matter()
       message=message//"      Omega_Matter(): "//trim(label)//char(10)
       write (label,'(e12.6)') Critical_Density()
       message=message//"  Critical_Density(): "//trim(label)
       call Galacticus_Error_Report("Power_Spectrum",message)
    end if
    call Initialize_Cosmological_Mass_Variance(mass)
    ! Compute the power spectrum.
    Power_Spectrum=Primordial_Power_Spectrum_Transferred(wavenumber)
    ! Scale by the normalization factor.
    Power_Spectrum=Power_Spectrum*sigmaNormalization**2
    return
  end function Power_Spectrum

  double precision function Mass_from_Cosmolgical_Root_Variance(sigma)
    !% Computes the mass corresponding to the given fractional mass fluctuation in real-space spherical top hats.
    implicit none
    double precision, intent(in) :: sigma
    integer                      :: iMass
    double precision             :: h,logMass
    
    ! Ensure that the sigma(M) tabulation exists.
    call Initialize_Cosmological_Mass_Variance()

    ! If the requested sigma is below the lowest value tabulated, attempt to tabulate to higher mass (lower sigma).
    do while (sigma < sigmaTable%y(-1))
       call Initialize_Cosmological_Mass_Variance(log(sigmaTable%x(-1))+1.0d0)
    end do

    ! If sigma exceeds the highest value tabulated, simply return the lowest tabulated mass.
    if (sigma > sigmaTable%y(1)) then
       Mass_from_Cosmolgical_Root_Variance=exp(sigmaTable%x(1))
       return
    end if

    ! Find the largest mass corresponding to this sigma.
    iMass=sigmaTable%size()
    do while (iMass > 1 .and. sigmaTable%y(iMass-1) < sigma)
       iMass=iMass-1
    end do

    h=(sigma-sigmaTable%y(iMass))/(sigmaTable%y(iMass-1)-sigmaTable%y(iMass))
    logMass=log(sigmaTable%x(iMass))*(1.0d0-h)+log(sigmaTable%x(iMass-1))*h
    Mass_from_Cosmolgical_Root_Variance=exp(logMass)
    return
  end function Mass_from_Cosmolgical_Root_Variance

  double precision function Cosmological_Mass_Root_Variance(mass)
    !% Computes the fractional mass fluctuation in real-space spherical top hats enclosing mass {\tt mass}.
    use Numerical_Interpolation
    implicit none
    double precision, intent(in) :: mass

    ! Check if we need to initialize this function.
    call Initialize_Cosmological_Mass_Variance(mass)
    
    ! Interpolate in tabulated function and return result.
    !$omp critical(Cosmological_Mass_Variance_Interpolate)
    Cosmological_Mass_Root_Variance=sigmaTable%interpolate(mass)
    !$omp end critical(Cosmological_Mass_Variance_Interpolate)
    return
  end function Cosmological_Mass_Root_Variance
  
  double precision function Cosmological_Mass_Root_Variance_Logarithmic_Derivative(mass)
    !% Computes the logarithmic derivative in the fractional mass fluctuation in real-space spherical top hats enclosing mass {\tt
    !% mass}.
    implicit none
    double precision, intent(in) :: mass
    double precision             :: sigma

    call Cosmological_Mass_Root_Variance_Plus_Logarithmic_Derivative(mass,sigma,Cosmological_Mass_Root_Variance_Logarithmic_Derivative)
    return
  end function Cosmological_Mass_Root_Variance_Logarithmic_Derivative
  
  subroutine Cosmological_Mass_Root_Variance_Plus_Logarithmic_Derivative(mass,sigma,sigmaLogarithmicDerivative)
    !% Returns both the fractional mass fluctuation in real-space spherical top hats enclosing mass {\tt mass} and its logarithmic derivative.
    use Numerical_Interpolation
    implicit none
    double precision, intent(in)  :: mass
    double precision, intent(out) :: sigma,sigmaLogarithmicDerivative

    ! Check if we need to initialize this function.
    call Initialize_Cosmological_Mass_Variance(mass)
    
    ! Interpolate in tabulated function and return result.
    !$omp critical(Cosmological_Mass_Variance_Interpolate)
    sigma                     =sigmaTable%interpolate        (mass)
    sigmaLogarithmicDerivative=sigmaTable%interpolateGradient(mass)*mass/sigma
    !$omp end critical(Cosmological_Mass_Variance_Interpolate)
    return
  end subroutine Cosmological_Mass_Root_Variance_Plus_Logarithmic_Derivative

  double precision function sigma_8()
    !% Return the value of $\sigma_8$.
    implicit none

    ! Ensure the module has been initialized.
    call Initialize_Cosmological_Mass_Variance()
    ! Return the value of sigma_8.
    sigma_8=sigma_8_Value
    return
  end function sigma_8

  subroutine Initialize_Cosmological_Mass_Variance(mass)
    !% Ensure that $\sigma(M)$ is tabulated over a range that includes {\tt logMass}. The default normalization, $\sigma_9=0.807$,
    !% is taken from \cite{komatsu_seven-year_2010}.
    use Input_Parameters
    use Memory_Management
    use Cosmological_Parameters
    use Numerical_Ranges
    use Numerical_Constants_Math
    use Numerical_Interpolation
    use Galacticus_Error
    !# <include directive="cosmologicalMassVarianceMethod" type="moduleUse">
    include 'structure_formation.cosmological_mass_variance.modules.inc'
    !# </include>
    implicit none
    double precision, intent(in   ), optional :: mass
    logical                                   :: remakeTable
    double precision                          :: massActual
    
    !$omp critical (Cosmological_Mass_Variance_Interpolate)
    ! Compute the normalization if required.
    remakeTable=.false.
    if (.not.sigmaInitialized) then
       !@ <inputParameter>
       !@   <name>sigma_8</name>
       !@   <defaultValue>0.817 (\citealt{hinshaw_nine-year_2012}; CMB$+H_0+$BAO)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The fractional mass fluctuation in the linear density field at the present day in spheres of radius 8~Mpc/h.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('sigma_8',sigma_8_Value,defaultValue=0.817d0)
       !@ <inputParameter>
       !@   <name>cosmologicalMassVarianceMethod</name>
       !@   <defaultValue>filteredPowerSpectrum</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Selects the method to be used for computing the cosmological mass variance.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('cosmologicalMassVarianceMethod',cosmologicalMassVarianceMethod,defaultValue='filteredPowerSpectrum')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="cosmologicalMassVarianceMethod" type="functionCall" functionType="void">
       !#  <functionArgs>cosmologicalMassVarianceMethod,Cosmological_Mass_Variance_Tabulate</functionArgs>
       include 'structure_formation.cosmological_mass_variance.inc'
       !# </include>
       if (.not.associated(Cosmological_Mass_Variance_Tabulate)) call Galacticus_Error_Report('Initialize_Cosmological_Mass_Variance','method ' &
            &//char(cosmologicalMassVarianceMethod)//' is unrecognized')
       ! Compute the mass at which the mass variance is normalized.
       massNormalization=(4.0d0*Pi/3.0d0)*Omega_Matter()*Critical_Density()*(radiusNormalization/Little_H_0())**3
       ! Flag that this module is now initialized.
       sigmaInitialized=.true.
       ! Table must be rebuilt.
       remakeTable=.true.
    end if
    ! Tabulate the mass variance.
    sigmaNormalization=sigma_8_Value
    if (present(mass)) then
       massActual=mass
    else
       massActual=massNormalization
    end if
    if (.not.remakeTable) remakeTable=(massActual < sigmaTable%x(1) .or. massActual > sigmaTable%x(-1))
    if (remakeTable) call Cosmological_Mass_Variance_Tabulate(massActual,massNormalization,sigmaNormalization,sigmaTable)
    !$omp end critical (Cosmological_Mass_Variance_Interpolate)
    return
  end subroutine Initialize_Cosmological_Mass_Variance
 
end module Power_Spectra
