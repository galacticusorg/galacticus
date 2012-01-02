!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which generates a tabulated transfer function using the Eisenstein \& Hu fitting formula.

module Transfer_Function_Eisenstein_Hu
  !% Implements generation of a tabulated transfer function using the Eisenstein \& Hu fitting formula.
  use ISO_Varying_String
  implicit none
  private
  public :: Transfer_Function_Eisenstein_Hu_Initialize, Transfer_Function_Eisenstein_Hu_State_Store,&
       & Transfer_Function_Eisenstein_Hu_State_Retrieve
  
  ! Flag to indicate if this module has been initialized.
  logical                     :: transferFunctionInitialized=.false.

  ! Wavenumber range and fineness of gridding.
  double precision            :: logWavenumberMaximum=dlog(10.0d0)
  double precision            :: logWavenumberMinimum=dlog(1.0d-5)
  integer,          parameter :: numberPointsPerDecade=1000

  ! Neutrino properties.
  double precision            :: effectiveNumberNeutrinos,summedNeutrinoMasses

contains
  
  !# <transferFunctionMethod>
  !#  <unitName>Transfer_Function_Eisenstein_Hu_Initialize</unitName>
  !# </transferFunctionMethod>
  subroutine Transfer_Function_Eisenstein_Hu_Initialize(transferFunctionMethod,Transfer_Function_Tabulate)
    !% Initializes the ``transfer function from Eisenstein \& Hu'' module.
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: transferFunctionMethod
    procedure(),          pointer, intent(inout) :: Transfer_Function_Tabulate
    
    if (transferFunctionMethod == 'Eisenstein-Hu1999') then
       Transfer_Function_Tabulate => Transfer_Function_Eisenstein_Hu_Make
       !@ <inputParameter>
       !@   <name>effectiveNumberNeutrinos</name>
       !@   <defaultValue>3.04 \citep{komatsu_seven-year_2010}</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The effective number of neutrino species as used in the \cite{eisenstein_power_1999} transfer function.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('effectiveNumberNeutrinos',effectiveNumberNeutrinos,defaultValue=3.04d0)
       !@ <inputParameter>
       !@   <name>summedNeutrinoMasses</name>
       !@   <defaultValue>0</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The summed mass (in electron volts) of all neutrino species.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('summedNeutrinoMasses'    ,summedNeutrinoMasses    ,defaultValue=0.00d0)
    end if
    return
  end subroutine Transfer_Function_Eisenstein_Hu_Initialize

  subroutine Transfer_Function_Eisenstein_Hu_Make(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber&
       &,transferFunctionLogT)
    !% Build a transfer function using the \cite{eisenstein_power_1999} fitting formula.
    use Memory_Management
    use Cosmological_Parameters
    use Numerical_Ranges
    use Numerical_Constants_Math
    implicit none
    double precision,                            intent(in)    :: logWavenumber
    double precision, allocatable, dimension(:), intent(inout) :: transferFunctionLogWavenumber,transferFunctionLogT
    integer,                                     intent(out)   :: transferFunctionNumberPoints
    integer                                                    :: iWavenumber
    double precision                                           :: qEH,wavenumber,Theta27,zeq,b1,b2,zd,yd,s,fv,Nv,fb ,fc,fcb,fvb&
         &,pc,pcb,alphav,Gammaeff,qeff,betac,L,C,Tsup,qv,Bk
 
    ! Set wavenumber range and number of points in table.
    logWavenumberMinimum=min(logWavenumberMinimum,logWavenumber-ln10)
    logWavenumberMaximum=max(logWavenumberMaximum,logWavenumber+ln10)
    transferFunctionNumberPoints=int((logWavenumberMaximum-logWavenumberMinimum)*dble(numberPointsPerDecade)/ln10)
    ! Deallocate arrays if currently allocated.
    if (allocated(transferFunctionLogWavenumber)) call Dealloc_Array(transferFunctionLogWavenumber)
    if (allocated(transferFunctionLogT))          call Dealloc_Array(transferFunctionLogT         )
    ! Allocate the arrays to current required size.
    call Alloc_Array(transferFunctionLogWavenumber,[transferFunctionNumberPoints])
    call Alloc_Array(transferFunctionLogT         ,[transferFunctionNumberPoints])
    ! Create range of wavenumbers.
    transferFunctionLogWavenumber=Make_Range(logWavenumberMinimum,logWavenumberMaximum,transferFunctionNumberPoints&
         &,rangeTypeLinear)
    ! Create transfer function.
    ! Present day CMB temperature [in units of 2.7K].
    Theta27=T_CMB()/2.7d0
    ! Redshift of matter-radiation equality.
    zeq=2.50d4*Omega_Matter()*(Little_H_0()**2)/(Theta27**4)
    ! Compute redshift at which baryons are released from Compton drag of photons (eqn. 2)
    b1=0.313d0*((Omega_Matter()*(Little_H_0()**2))**(-0.419d0))*(1.0d0+0.607d0*((Omega_Matter()*(Little_H_0()**2))**0.674d0))
    b2=0.238d0*((Omega_Matter()*(Little_H_0()**2))**0.223d0)
    zd=1291.0d0*((Omega_Matter()*(Little_H_0()**2))**0.251d0)*(1.0d0+b1*((Omega_b()*(Little_H_0()**2))**b2))/(1.0d0+0.659d0*((Omega_Matter()*(Little_H_0()**2))**0.828d0))
    ! Relative expansion factor between previous two computed redshifts.
    yd=(1.0d0+zeq)/(1.0d0+zd)
    ! Compute the comoving distance that a sound wave can propagate prior to zd (i.e. sound horizon; eq. 4)
    s=44.5d0*dlog(9.83d0/Omega_Matter()/(Little_H_0()**2))/dsqrt(1.0d0+10.0d0*((Omega_b()*(Little_H_0()**2))**0.75d0))
    ! Specify properties of neutrinos. Mass fraction formula is from Komatsu et al. (2007; http://adsabs.harvard.edu/abs/2010arXiv1001.4538K).
    fv=summedNeutrinoMasses/94.0d0/(Little_H_0()**2)/Omega_Matter()
    Nv=effectiveNumberNeutrinos
    ! Compute baryonic and cold dark matter fractions.
    fb=Omega_b()/Omega_Matter()
    fc=(Omega_Matter()-Omega_b())/Omega_Matter()
    ! Total matter fraction.
    fcb=fb+fc
    ! Baryonic + neutrino fraction.
    fvb=fv+fb
    ! Compute small scale suppression factor (eqn. 15).
    pc =0.25d0*(5.0d0-dsqrt(1.0d0+24.0d0*fc ))
    pcb=0.25d0*(5.0d0-dsqrt(1.0d0+24.0d0*fcb))
    alphav=(fc/fcb)*((5.0d0-2.0d0*(pc+pcb))/(5.0d0-4.0d0*pcb))*((1.0d0-0.533d0*fvb+0.126d0*(fvb**3))*((1.0d0+yd)**(pcb-pc))/(1.0d0-0.193d0&
         &*dsqrt(fv*Nv)+0.169d0*fv*(Nv**0.2d0)))*(1.0d0+0.5d0*(pc-pcb)*(1.0d0+1.0d0/(3.0d0-4.0d0*pc)/(7.0d0-4.0d0*pcb))/(1.0d0+yd))
    ! Loop over all wavenumbers.
    do iWavenumber=1,transferFunctionNumberPoints
       wavenumber=dexp(transferFunctionLogWavenumber(iWavenumber))
       ! Compute effective q.
       qEH=wavenumber*(Theta27**2)/Omega_Matter()/(Little_H_0()**2)
       ! Compute rescaled shape parameter (eqn. 16)
       Gammaeff=Omega_Matter()*(Little_H_0()**2)*(dsqrt(alphav)+(1.0d0-dsqrt(alphav))/(1.0d0+((0.43d0*wavenumber*s)**4)))
       qeff=wavenumber*(Theta27**2)/Gammaeff
       betac=1.0d0/(1.0d0-0.949d0*fvb)                     ! Eqn. 21.
       L=dlog(dexp(1.0d0)+1.84d0*betac*dsqrt(alphav)*qeff) ! Eqn. 19.
       C=14.4d0+325.0d0/(1.0d0+60.5d0*(qeff**1.11d0))      ! Eqn. 20.
       Tsup=L/(L+C*(qeff**2))                              ! Zero baryon form of the transfer function (eqn. 18).
       ! Apply correction for scales close to horizon.
       if (fv > 0.0d0.and.fv <= 0.3d0) then
          qv=3.92d0*qEH*dsqrt(Nv)/fv
          Bk=1.0d0+(1.2d0*(fv**0.64d0)*(Nv**(0.3d0+0.6d0*fv)))/((qv**(-1.6d0))+(qv**0.8d0))
       else
          qv=0.0d0
          Bk=1.0d0
       end if
       transferFunctionLogT(iWavenumber)=dlog(Tsup*Bk)
    end do
    return
  end subroutine Transfer_Function_Eisenstein_Hu_Make
  
  !# <galacticusStateStoreTask>
  !#  <unitName>Transfer_Function_Eisenstein_Hu_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Transfer_Function_Eisenstein_Hu_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use FGSL
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) logWavenumberMinimum,logWavenumberMaximum
    return
  end subroutine Transfer_Function_Eisenstein_Hu_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Transfer_Function_Eisenstein_Hu_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Transfer_Function_Eisenstein_Hu_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use FGSL
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) logWavenumberMinimum,logWavenumberMaximum
    return
  end subroutine Transfer_Function_Eisenstein_Hu_State_Retrieve
    
end module Transfer_Function_Eisenstein_Hu
