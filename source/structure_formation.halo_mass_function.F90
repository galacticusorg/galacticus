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

module Halo_Mass_Function
  use ISO_Varying_String
  use FGSL
  use, intrinsic :: ISO_C_Binding                             
  implicit none
  private
  public :: Halo_Mass_Function_Differential, Halo_Mass_Function_Integrated, Halo_Mass_Fraction_Integrated

  ! Flag to indicate if this module has been initialized.  
  logical                                        :: haloMassFunctionInitialized=.false.

  ! Name of power spectrum method used.
  type(varying_string)                           :: haloMassFunctionMethod

  ! Pointer to the function that computes the mass function.
  procedure(Halo_Mass_Function_Differential), pointer :: Halo_Mass_Function_Differential_Get => null()
  
contains

  subroutine Halo_Mass_Function_Initialize()
    !% Initializes the halo mass function module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="haloMassFunctionMethod" type="moduleUse">
    include 'structure_formation.CDM.halo_mass_function.modules.inc'
    !# </include>
    implicit none

    if (.not.haloMassFunctionInitialized) then
       !$omp critical(Halo_Mass_Function_Initialize)
       if (.not.haloMassFunctionInitialized) then
          ! Get the transfer function method parameter.
          !@ <inputParameter>
          !@   <name>haloMassFunctionMethod</name>
          !@   <defaultValue>Tinker2008</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing the dark matter halo mass function.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('haloMassFunctionMethod',haloMassFunctionMethod,defaultValue='Tinker2008')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="haloMassFunctionMethod" type="functionCall" functionType="void">
          !#  <functionArgs>haloMassFunctionMethod,Halo_Mass_Function_Differential_Get</functionArgs>
          include 'structure_formation.CDM.halo_mass_function.inc'
          !# </include>
          if (.not.associated(Halo_Mass_Function_Differential_Get)) call Galacticus_Error_Report('Halo_Mass_Function_Initialize','method '&
               &//char(haloMassFunctionMethod)//' is unrecognized')
          ! Flag that the module is now initialized.
          haloMassFunctionInitialized=.true.
       end if
       !$omp end critical(Halo_Mass_Function_Initialize)
    end if
    return
  end subroutine Halo_Mass_Function_Initialize

  double precision function Halo_Mass_Function_Differential(time,mass)
    !% Return the differential halo mass function for {\tt mass} [$M_\odot$] at {\tt time}.
    use Numerical_Interpolation
    implicit none
    double precision, intent(in) :: time,mass

    ! Ensure that the module is initialized.
    call Halo_Mass_Function_Initialize()
    
    ! Call the function that does the work.
    Halo_Mass_Function_Differential=Halo_Mass_Function_Differential_Get(time,mass)
    return
  end function Halo_Mass_Function_Differential

  double precision function Halo_Mass_Function_Integrated(time,massLow,massHigh)
    !% Return tha halo mass function integrated between {\tt massLow} and {\tt massHigh}.
    use Numerical_Integration
    implicit none
    double precision,                intent(in) :: time,massLow,massHigh
    double precision,                target     :: timeTarget
    double precision                            :: logMassLow,logMassHigh
    type(c_ptr)                                 :: parameterPointer
    type(fgsl_function)                         :: integrandFunction
    type(fgsl_integration_workspace)            :: integrationWorkspace

    parameterPointer=c_loc(timeTarget)
    timeTarget=time
    logMassLow =dlog(massLow )
    logMassHigh=dlog(massHigh)
    Halo_Mass_Function_Integrated=Integrate(logMassLow,logMassHigh,Halo_Mass_Function_Integrand,parameterPointer &
         &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-5,integrationRule&
         &=FGSL_Integ_Gauss15)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return
  end function Halo_Mass_Function_Integrated

  function Halo_Mass_Function_Integrand(logMass,parameterPointer) bind(c)
    implicit none
    real(c_double)          :: Halo_Mass_Function_Integrand
    real(c_double), value   :: logMass
    type(c_ptr),    value   :: parameterPointer
    real(c_double), pointer :: time
    real(c_double)          :: mass

    ! Extract integrand parameters.
    call c_f_pointer(parameterPointer,time)
    mass=dexp(logMass)

    ! Return the differential mass function multiplied by mass since we are integrating over log(mass).
    Halo_Mass_Function_Integrand=Halo_Mass_Function_Differential(time,mass)*mass
    return
  end function Halo_Mass_Function_Integrand

  double precision function Halo_Mass_Fraction_Integrated(time,massLow,massHigh)
    !% Return the halo mass fraction integrated between {\tt massLow} and {\tt massHigh}.
    use Numerical_Integration
    use Cosmological_Parameters
    implicit none
    double precision,                intent(in) :: time,massLow,massHigh
    double precision,                target     :: timeTarget
    double precision                            :: logMassLow,logMassHigh
    type(c_ptr)                                 :: parameterPointer
    type(fgsl_function)                         :: integrandFunction
    type(fgsl_integration_workspace)            :: integrationWorkspace

    parameterPointer=c_loc(timeTarget)
    timeTarget=time
    logMassLow =dlog(massLow )
    logMassHigh=dlog(massHigh)
    Halo_Mass_Fraction_Integrated=Integrate(logMassLow,logMassHigh,Halo_Mass_Fraction_Integrand,parameterPointer &
         &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-5,integrationRule&
         &=FGSL_Integ_Gauss15)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    ! Convert to a mass fraction.
    Halo_Mass_Fraction_Integrated=Halo_Mass_Fraction_Integrated/Critical_Density()/Omega_Matter()
    return
  end function Halo_Mass_Fraction_Integrated

  function Halo_Mass_Fraction_Integrand(logMass,parameterPointer) bind(c)
    !% Integrand function used in computing the halo mass fraction.
    implicit none
    real(c_double)          :: Halo_Mass_Fraction_Integrand
    real(c_double), value   :: logMass
    type(c_ptr),    value   :: parameterPointer
    real(c_double), pointer :: time
    real(c_double)          :: mass

    ! Extract integrand parameters.
    call c_f_pointer(parameterPointer,time)
    mass=dexp(logMass)

    ! Return the differential mass function multiplied by mass since we are integrating over log(mass) and by mass again to get
    ! the mass fraction.
    Halo_Mass_Fraction_Integrand=Halo_Mass_Function_Differential(time,mass)*(mass**2)
    return
  end function Halo_Mass_Fraction_Integrand

end module Halo_Mass_Function
