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

!% Contains a module which generates a tabulated transfer function using the Eisenstein \& Hu fitting formula.

module Transfer_Function_Eisenstein_Hu
  !% Implements generation of a tabulated transfer function using the Eisenstein \& Hu fitting formula.
  use ISO_Varying_String
  implicit none
  private
  public :: Transfer_Function_Eisenstein_Hu_Initialize, Transfer_Function_Eisenstein_Hu_State_Store,&
       & Transfer_Function_Eisenstein_Hu_State_Retrieve

  ! Flag to indicate if this module has been initialized.
  logical                     :: transferFunctionInitialized   =.false.

  ! Wavenumber range and fineness of gridding.
  double precision            :: logWavenumberMaximum          =log(10.0d0)
  double precision            :: logWavenumberMinimum          =log(1.0d-5)
  integer         , parameter :: numberPointsPerDecade         =1000

  ! Neutrino properties.
  double precision            :: effectiveNumberNeutrinos                  , summedNeutrinoMasses

  ! Warm dark matter cut-off scale.
  double precision            :: transferFunctionWdmCutOffScale
  double precision            :: transferFunctionWdmEpsilon
  double precision            :: transferFunctionWdmEta
  double precision            :: transferFunctionWdmNu

contains

  !# <transferFunctionMethod>
  !#  <unitName>Transfer_Function_Eisenstein_Hu_Initialize</unitName>
  !# </transferFunctionMethod>
  subroutine Transfer_Function_Eisenstein_Hu_Initialize(transferFunctionMethod,Transfer_Function_Tabulate)
    !% Initializes the ``transfer function from Eisenstein \& Hu'' module.
    use Input_Parameters
    implicit none
    type     (varying_string                      ), intent(in   )          :: transferFunctionMethod
    procedure(Transfer_Function_Eisenstein_Hu_Make), intent(inout), pointer :: Transfer_Function_Tabulate

    if (transferFunctionMethod == 'Eisenstein-Hu1999') then
       Transfer_Function_Tabulate => Transfer_Function_Eisenstein_Hu_Make
       !@ <inputParameter>
       !@   <name>effectiveNumberNeutrinos</name>
       !@   <defaultValue>3.046 \citep{mangano_relic_2005}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The effective number of neutrino species as used in the \cite{eisenstein_power_1999} transfer function.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('effectiveNumberNeutrinos',effectiveNumberNeutrinos,defaultValue=3.046d0)
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
       !@ <inputParameter>
       !@   <name>transferFunctionWdmCutOffScale</name>
       !@   <defaultValue>0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The cut-off scale in the transfer function due to warm dark matter.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('transferFunctionWdmCutOffScale',transferFunctionWdmCutOffScale,defaultValue=0.00d0)
       !@ <inputParameter>
       !@   <name>transferFunctionWdmEpsilon</name>
       !@   <defaultValue>0.361 \citep{barkana_constraints_2001}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\epsilon$ appearing in the warm dark matter transfer function \citep{barkana_constraints_2001}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('transferFunctionWdmEpsilon',transferFunctionWdmEpsilon,defaultValue=0.361d0)
       !@ <inputParameter>
       !@   <name>transferFunctionWdmEta</name>
       !@   <defaultValue>5.0 \citep{barkana_constraints_2001}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\epsilon$ appearing in the warm dark matter transfer function \citep{barkana_constraints_2001}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('transferFunctionWdmEta',transferFunctionWdmEta,defaultValue=5.0d0)
       !@ <inputParameter>
       !@   <name>transferFunctionWdmNu</name>
       !@   <defaultValue>1.2 \citep{barkana_constraints_2001}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\epsilon$ appearing in the warm dark matter transfer function \citep{barkana_constraints_2001}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('transferFunctionWdmNu',transferFunctionWdmNu,defaultValue=1.200d0)
    end if
    return
  end subroutine Transfer_Function_Eisenstein_Hu_Initialize

  subroutine Transfer_Function_Eisenstein_Hu_Make(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber&
       &,transferFunctionLogT)
    !% Build a transfer function using the \cite{eisenstein_power_1999} fitting formula. Includes a modification for warm dark
    !% matter using the fitting function of \citeauthor{bode_halo_2001}~(\citeyear{bode_halo_2001}; as re-expressed by
    !% \citealt{barkana_constraints_2001}) to impose a cut-off below a specified {\tt [transferFunctionWdmCutOffScale]}.
    use Memory_Management
    use Cosmological_Parameters
    use Numerical_Ranges
    use Numerical_Constants_Math
    implicit none
    double precision                           , intent(in   ) :: logWavenumber
    double precision, allocatable, dimension(:), intent(inout) :: transferFunctionLogT        , transferFunctionLogWavenumber
    integer                                    , intent(  out) :: transferFunctionNumberPoints
    integer                                                    :: iWavenumber
    double precision                                           :: Bk                          , C                             , &
         &                                                        Gammaeff                    , L                             , &
         &                                                        Nv                          , Theta27                       , &
         &                                                        Tsup                        , alphav                        , &
         &                                                        b1                          , b2                            , &
         &                                                        betac                       , fb                            , &
         &                                                        fc                          , fcb                           , &
         &                                                        fv                          , fvb                           , &
         &                                                        pc                          , pcb                           , &
         &                                                        qEH                         , qeff                          , &
         &                                                        qv                          , s                             , &
         &                                                        transferFunctionWdmFactor   , wavenumber                    , &
         &                                                        yd                          , zd                            , &
         &                                                        zeq

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
    s=44.5d0*log(9.83d0/Omega_Matter()/(Little_H_0()**2))/sqrt(1.0d0+10.0d0*((Omega_b()*(Little_H_0()**2))**0.75d0))
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
    pc =0.25d0*(5.0d0-sqrt(1.0d0+24.0d0*fc ))
    pcb=0.25d0*(5.0d0-sqrt(1.0d0+24.0d0*fcb))
    alphav=(fc/fcb)*((5.0d0-2.0d0*(pc+pcb))/(5.0d0-4.0d0*pcb))*((1.0d0-0.533d0*fvb+0.126d0*(fvb**3))*((1.0d0+yd)**(pcb-pc))/(1.0d0-0.193d0&
         &*sqrt(fv*Nv)+0.169d0*fv*(Nv**0.2d0)))*(1.0d0+0.5d0*(pc-pcb)*(1.0d0+1.0d0/(3.0d0-4.0d0*pc)/(7.0d0-4.0d0*pcb))/(1.0d0+yd))
    ! Loop over all wavenumbers.
    do iWavenumber=1,transferFunctionNumberPoints
       wavenumber=exp(transferFunctionLogWavenumber(iWavenumber))
       ! Compute effective q.
       qEH=wavenumber*(Theta27**2)/Omega_Matter()/(Little_H_0()**2)
       ! Compute rescaled shape parameter (eqn. 16)
       Gammaeff=Omega_Matter()*(Little_H_0()**2)*(sqrt(alphav)+(1.0d0-sqrt(alphav))/(1.0d0+((0.43d0*wavenumber*s)**4)))
       qeff=wavenumber*(Theta27**2)/Gammaeff
       betac=1.0d0/(1.0d0-0.949d0*fvb)                     ! Eqn. 21.
       L=log(exp(1.0d0)+1.84d0*betac*sqrt(alphav)*qeff) ! Eqn. 19.
       C=14.4d0+325.0d0/(1.0d0+60.5d0*(qeff**1.11d0))      ! Eqn. 20.
       Tsup=L/(L+C*(qeff**2))                              ! Zero baryon form of the transfer function (eqn. 18).
       ! Apply correction for scales close to horizon.
       if (fv > 0.0d0.and.fv <= 0.3d0) then
          qv=3.92d0*qEH*sqrt(Nv)/fv
          Bk=1.0d0+(1.2d0*(fv**0.64d0)*(Nv**(0.3d0+0.6d0*fv)))/((qv**(-1.6d0))+(qv**0.8d0))
       else
          qv=0.0d0
          Bk=1.0d0
       end if
       transferFunctionLogT(iWavenumber)=log(Tsup*Bk)
       if (transferFunctionWdmCutOffScale > 0.0d0) then
          transferFunctionWdmFactor=                                &
               & (                                                  &
               &   1.0d0                                            &
               &  +                                                 &
               &   (                                                &
               &     transferFunctionWdmEpsilon                     &
               &    *wavenumber                                     &
               &    *transferFunctionWdmCutOffScale                 &
               &   )**(2.0d0*transferFunctionWdmNu)                 &
               & )**(-transferFunctionWdmEta/transferFunctionWdmNu)
          if (transferFunctionWdmFactor > 0.0d0) then
             transferFunctionLogT(iWavenumber)=transferFunctionLogT(iWavenumber)+log(transferFunctionWdmFactor)
          else
             transferFunctionLogT(iWavenumber)=transferFunctionLogT(iWavenumber)-100.0d0
          end if
       end if
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
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

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
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) logWavenumberMinimum,logWavenumberMaximum
    return
  end subroutine Transfer_Function_Eisenstein_Hu_State_Retrieve

end module Transfer_Function_Eisenstein_Hu
