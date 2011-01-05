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


module Halo_Mass_Function
  use ISO_Varying_String
  use FGSL
  use, intrinsic :: ISO_C_Binding                             
  private
  public :: Halo_Mass_Function_Differential, Halo_Mass_Function_Integrated, Halo_Mass_Fraction_Integrated

  ! Flag to indicate if this module has been initialized.  
  logical                                        :: haloMassFunctionInitialized=.false.

  ! Variables to hold the tabulated halo mass function data.
  integer                                        :: haloMassFunctionNumberPoints=-1
  double precision,    allocatable, dimension(:) :: haloMassFunctionLogMass,haloMassFunctionLogAbundance
  double precision                               :: haloMassFunctionTime
  type(fgsl_interp)                              :: interpolationObject
  type(fgsl_interp_accel)                        :: interpolationAccelerator
  logical                                        :: resetInterpolation=.true.

  ! Name of power spectrum method used.
  type(varying_string)                           :: haloMassFunctionMethod

  ! Pointer to the subroutine that tabulates the transfer function and template interface for that subroutine.
  procedure(Halo_Mass_Function_Tabulate_Template), pointer :: Halo_Mass_Function_Tabulate => null()
  abstract interface
     subroutine Halo_Mass_Function_Tabulate_Template(time,logMass,haloMassFunctionNumberPoints,haloMassFunctionLogMass &
          &,haloMassFunctionLogAbundance)
    double precision,                            intent(in)    :: time,logMass
    double precision, allocatable, dimension(:), intent(inout) :: haloMassFunctionLogMass,haloMassFunctionLogAbundance
    integer,                                     intent(out)   :: haloMassFunctionNumberPoints
  end subroutine Halo_Mass_Function_Tabulate_Template
 end interface
  
contains

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
    Halo_Mass_Function_Integrated=Integrate(logMassLow,logMassHigh,Halo_Mass_Function_Integrand,parameterPointer&
         &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-8)
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
         &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-8)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    ! Convert to a mass fraction.
    Halo_Mass_Fraction_Integrated=Halo_Mass_Fraction_Integrated/Critical_Density()/Omega_0()
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

  double precision function Halo_Mass_Function_Differential(time,mass)
    !% Return the differential halo mass function for {\tt mass} [$M_\odot$] at {\tt time}.
    use Numerical_Interpolation
    implicit none
    double precision, intent(in) :: time,mass
    double precision             :: logMass

    ! Get logarithm of mass.
    logMass=dlog(mass)

    !$omp critical(Halo_Mass_Function_Initialization) 
    ! Initialize if necessary.
    if (.not.haloMassFunctionInitialized) then
       call Halo_Mass_Function_Initialize(time,logMass)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)
       haloMassFunctionTime=time
       resetInterpolation=.true.
    end if

    ! If mass is out of range, attempt to remake the table.
    if (logMass<haloMassFunctionLogMass(1) .or. logMass>haloMassFunctionLogMass(haloMassFunctionNumberPoints) .or. time /=&
         & haloMassFunctionTime ) then
       call Halo_Mass_Function_Tabulate(time,logMass,haloMassFunctionNumberPoints,haloMassFunctionLogMass&
            &,haloMassFunctionLogAbundance)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)  
       haloMassFunctionTime=time
       resetInterpolation=.true.
    end if
    !$omp end critical(Halo_Mass_Function_Initialization)

    ! Interpolate in the tabulated function and return a value.
    Halo_Mass_Function_Differential=dexp(Interpolate(haloMassFunctionNumberPoints,haloMassFunctionLogMass &
         &,haloMassFunctionLogAbundance,interpolationObject,interpolationAccelerator,logMass,reset=resetInterpolation&
         &,interpolationType=fgsl_interp_cspline))
    return
  end function Halo_Mass_Function_Differential

  subroutine Halo_Mass_Function_Initialize(time,logMass)
    !% Initializes the transfer function module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="haloMassFunctionMethod" type="moduleUse">
    include 'structure_formation.CDM.halo_mass_function.modules.inc'
    !# </include>
    implicit none
    double precision, intent(in) :: time,logMass

    ! Get the transfer function method parameter.
    !@ <inputParameter>
    !@   <name>haloMassFunctionMethod</name>
    !@   <defaultValue>Tinker2008</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The name of the method to be used for computing the dark matter halo mass function.
    !@   </description>
    !@ </inputParameter>
    call Get_Input_Parameter('haloMassFunctionMethod',haloMassFunctionMethod,defaultValue='Tinker2008')
    ! Include file that makes calls to all available method initialization routines.
    !# <include directive="haloMassFunctionMethod" type="code" action="subroutine">
    !#  <subroutineArgs>haloMassFunctionMethod,Halo_Mass_Function_Tabulate</subroutineArgs>
    include 'structure_formation.CDM.halo_mass_function.inc'
    !# </include>
    if (.not.associated(Halo_Mass_Function_Tabulate)) call Galacticus_Error_Report('Halo_Mass_Function_Initialize','method '&
         &//char(haloMassFunctionMethod)//' is unrecognized')
    ! Call routine to initialize the transfer function.
    call Halo_Mass_Function_Tabulate(time,logMass,haloMassFunctionNumberPoints,haloMassFunctionLogMass&
         &,haloMassFunctionLogAbundance)
    ! Flag that the module is now initialized.
    haloMassFunctionInitialized=.true.
    return
  end subroutine Halo_Mass_Function_Initialize
  
end module Halo_Mass_Function
