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

!% Contains a module which implements the \cite{behroozi_comprehensive_2010} fitting function descriptions of the conditional stellar mass function.

module Conditional_Stellar_Mass_Functions_Behroozi2010
  !% Implements the \cite{behroozi_comprehensive_2010} fitting function descriptions of the conditional stellar mass function.
  use FGSL
  private
  public :: Conditional_Stellar_Mass_Functions_Behroozi2010_Initialize

  ! Parameters of the Behroozi et al. fitting function.
  double precision                                               :: conditionalStellarMassFunctionBehrooziAlphaSatellite                                         
  double precision                                               :: conditionalStellarMassFunctionBehrooziLog10M1                                                
  double precision                                               :: conditionalStellarMassFunctionBehrooziLog10Mstar0                                            
  double precision                                               :: conditionalStellarMassFunctionBehrooziBeta                                                   
  double precision                                               :: conditionalStellarMassFunctionBehrooziDelta                                                  
  double precision                                               :: conditionalStellarMassFunctionBehrooziGamma                                                  
  double precision                                               :: conditionalStellarMassFunctionBehrooziSigmaLogMstar                                          
  double precision                                               :: conditionalStellarMassFunctionBehrooziBCut                                                   
  double precision                                               :: conditionalStellarMassFunctionBehrooziBSatellite                                             
  double precision                                               :: conditionalStellarMassFunctionBehrooziBetaCut                                                
  double precision                                               :: conditionalStellarMassFunctionBehrooziBetaSatellite                                          
  double precision                                               :: Mstar0                                                                                       
  
  ! Tablulation of stellar-halo mass relation.
  integer                            , parameter                 :: fMassStellarTablePointsPerDecade                    =10                                      
  integer                                                        :: fMassStellarTableCount                                                                       
  double precision                                               :: fMassStellarTableMaximum                            =1.0d13 , fMassStellarTableMinimum=1.0d8 
  double precision                                               :: fMassHaloTableMaximum                                       , fMassHaloTableMinimum          
  double precision                   , allocatable, dimension(:) :: fMassHaloTable                                              , fMassStellarTable              
  logical                                                        :: interpolationReset                                  =.false.                                 
  type            (fgsl_interp      )                            :: interpolationObject                                                                          
  type            (fgsl_interp_accel)                            :: interpolationAccelerator                                                                     
  
contains

  !# <conditionalStellarMassFunctionMethod>
  !#  <unitName>Conditional_Stellar_Mass_Functions_Behroozi2010_Initialize</unitName>
  !# </conditionalStellarMassFunctionMethod>
  subroutine Conditional_Stellar_Mass_Functions_Behroozi2010_Initialize(conditionalStellarMassFunctionMethod&
       &,Cumulative_Conditional_Stellar_Mass_Function_Get,Cumulative_Conditional_Stellar_Mass_Function_Var_Get)
    !% Initializes the ``Behroozi2010'' conditional stellar mass function method.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string  ), intent(in   )          :: conditionalStellarMassFunctionMethod                                                                   
    procedure(double precision), intent(inout), pointer :: Cumulative_Conditional_Stellar_Mass_Function_Get, Cumulative_Conditional_Stellar_Mass_Function_Var_Get 
    
    if (conditionalStellarMassFunctionMethod == 'Behroozi2010') then
       Cumulative_Conditional_Stellar_Mass_Function_Get     => Cumulative_Conditional_Stellar_Mass_Function_Behroozi2010
       Cumulative_Conditional_Stellar_Mass_Function_Var_Get => Cumulative_Conditional_Stellar_Mass_Function_Var_Behroozi2010

       ! Get parameters of the fitting function.
       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziAlphaSatellite</name>
       !@   <defaultValue>1.0 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\alpha_{\rm sat}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziAlphaSatellite',conditionalStellarMassFunctionBehrooziAlphaSatellite,defaultValue=1.0d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziLog10M1</name>
       !@   <defaultValue>12.520 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\log_{10}M_1$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziLog10M1',conditionalStellarMassFunctionBehrooziLog10M1,defaultValue=12.520d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziLog10Mstar0</name>
       !@   <defaultValue>10.916 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\log_{10}M_{\star,0}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziLog10Mstar0',conditionalStellarMassFunctionBehrooziLog10Mstar0,defaultValue=10.916d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziBeta</name>
       !@   <defaultValue>0.457 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\beta$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziBeta',conditionalStellarMassFunctionBehrooziBeta,defaultValue=0.457d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziDelta</name>
       !@   <defaultValue>0.5666 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\delta$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziDelta',conditionalStellarMassFunctionBehrooziDelta,defaultValue=0.5666d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziGamma</name>
       !@   <defaultValue>1.53 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\gamma$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziGamma',conditionalStellarMassFunctionBehrooziGamma,defaultValue=1.53d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziSigmaLogMstar</name>
       !@   <defaultValue>0.206 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\sigma_{\log M_\star}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziSigmaLogMstar',conditionalStellarMassFunctionBehrooziSigmaLogMstar,defaultValue=0.206d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziBCut</name>
       !@   <defaultValue>1.47 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $B_{\rm cut}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziBCut',conditionalStellarMassFunctionBehrooziBCut,defaultValue=1.47d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziBSatellite</name>
       !@   <defaultValue>10.62 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $B_{\rm sat}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziBSatellite',conditionalStellarMassFunctionBehrooziBSatellite,defaultValue=10.62d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziBetaCut</name>
       !@   <defaultValue>$-0.13$ (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\beta_{\rm cut}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziBetaCut',conditionalStellarMassFunctionBehrooziBetaCut,defaultValue=-0.13d0)

       !@ <inputParameter>
       !@   <name>conditionalStellarMassFunctionBehrooziBetaSatellite</name>
       !@   <defaultValue>0.859 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\beta_{\rm sat}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalStellarMassFunctionBehrooziBetaSatellite',conditionalStellarMassFunctionBehrooziBetaSatellite,defaultValue=0.859d0)

       ! Compute the actual value of the Mstar0 parameter.
       Mstar0=10.0d0**conditionalStellarMassFunctionBehrooziLog10Mstar0
    end if
    return
  end subroutine Conditional_Stellar_Mass_Functions_Behroozi2010_Initialize

  double precision function Cumulative_Conditional_Stellar_Mass_Function_Behroozi2010(massHalo,massStellar)
    !% Computes the cumulative conditional mass function, $\langle N(M_\star|M_{\rm halo}) \rangle \equiv \phi(M_\star|M_{\rm
    !% halo})$ using the fitting formula of \cite{behroozi_comprehensive_2010}.
    implicit none
    double precision, intent(in   ) :: massHalo      , massStellar      
    double precision                :: numberCentrals, numberSatellites 
    
    ! Get the number of satellites and centrals.
    call Cumulative_Conditional_Stellar_Mass_Function_Compute(massHalo,massStellar,numberCentrals,numberSatellites)

    ! Sum central and satellite contributions.
    Cumulative_Conditional_Stellar_Mass_Function_Behroozi2010=numberCentrals+numberSatellites
    return
  end function Cumulative_Conditional_Stellar_Mass_Function_Behroozi2010

  double precision function Cumulative_Conditional_Stellar_Mass_Function_Var_Behroozi2010(massHalo,massStellarLow,massStellarHigh)
    !% Computes the variance in the cumulative conditional mass function, $\langle N(M_\star|M_{\rm halo}) \rangle \equiv
    !% \phi(M_\star|M_{\rm halo})$ using the fitting formula of \cite{behroozi_comprehensive_2010}. Assumes that the number of
    !% satellite galaxies is Poisson distributed, while the number of central galaxies follows a Bernoulli distribution, and that
    !% the numbers of satellites and centrals are uncorrelated.
    implicit none
    double precision, intent(in   ) :: massHalo        , massStellarHigh     , massStellarLow         
    double precision                :: numberCentrals  , numberCentralsHigh  , numberCentralsLow  , & 
         &                             numberSatellites, numberSatellitesHigh, numberSatellitesLow    
    
    ! Get the number of satellites and centrals.
    call Cumulative_Conditional_Stellar_Mass_Function_Compute(massHalo,massStellarLow ,numberCentralsLow ,numberSatellitesLow )
    call Cumulative_Conditional_Stellar_Mass_Function_Compute(massHalo,massStellarHigh,numberCentralsHigh,numberSatellitesHigh)
    numberSatellites=max(numberSatellitesLow-numberSatellitesHigh,0.0d0)
    numberCentrals  =max(numberCentralsLow  -numberCentralsHigh  ,0.0d0)

    ! Compute the variance (using the Bienaymé formula for uncorrelated variables).
    Cumulative_Conditional_Stellar_Mass_Function_Var_Behroozi2010= &
         &  numberSatellites                                       & ! Satellites are Poisson distributed, so the variance is just their number.
         & +numberCentrals*(1.0d0-numberCentrals)                    ! Centrals are Bernoulli distributed.
    return
  end function Cumulative_Conditional_Stellar_Mass_Function_Var_Behroozi2010

  subroutine Cumulative_Conditional_Stellar_Mass_Function_Compute(massHalo,massStellar,numberCentrals,numberSatellites)
    !% Computes the cumulative conditional mass function, $\langle N(M_\star|M_{\rm halo}) \rangle \equiv \phi(M_\star|M_{\rm
    !% halo})$ using the fitting formula of \cite{behroozi_comprehensive_2010}.
    use Memory_Management
    use Numerical_Ranges
    use Numerical_Interpolation
    implicit none
    double precision                    , intent(in   ) :: massHalo                  , massStellar                   
    double precision                    , intent(  out) :: numberCentrals            , numberSatellites              
    double precision, parameter                         :: massNormalization  =1.0d12                                
    double precision                                    :: fMassHalo                 , massCut                   , & 
         &                                                 massSatellite                                             
    double precision              , save                :: fMassStellar              , massHaloPrevious   =-1.0d0    
    !$omp threadprivate(massHaloPrevious,fMassStellar)
    double precision, dimension(2), save                :: fMassHaloStored           , massCutStored             , & 
         &                                                 massSatelliteStored       , massStellarPrevious=-1.0d0    
    !$omp threadprivate(massStellarPrevious,fMassHaloStored,massSatelliteStored,massCutStored)
    ! Ensure that the stellar-halo mass relation is tabulated over a sufficient range.    !$omp critical(CSMF_Behroozi2010_Tabulate)
    do while (massHalo < fMassHaloTableMinimum .or. massHalo > fMassHaloTableMaximum)
       if (massHalo < fMassHaloTableMinimum) fMassStellarTableMinimum=0.5d0*fMassStellarTableMinimum
       if (massHalo > fMassHaloTableMaximum) fMassStellarTableMaximum=2.0d0*fMassStellarTableMaximum
       fMassStellarTableCount=int(log10(fMassStellarTableMaximum/fMassStellarTableMinimum)*fMassStellarTablePointsPerDecade)+1
       if (allocated(fMassStellarTable)) call Dealloc_Array(fMassStellarTable)
       if (allocated(fMassHaloTable   )) call Dealloc_Array(fMassHaloTable   )
       call Alloc_Array(fMassStellarTable,[fMassStellarTableCount])
       call Alloc_Array(fMassHaloTable   ,[fMassStellarTableCount])
       fMassStellarTable    =Make_Range(fMassStellarTableMinimum,fMassStellarTableMaximum,fMassStellarTableCount,rangeType=rangeTypeLogarithmic)
       fMassHaloTable       =fSHMRInverse(fMassStellarTable)
       fMassHaloTableMinimum=fMassHaloTable(                     1)
       fMassHaloTableMaximum=fMassHaloTable(fMassStellarTableCount)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,interpolationReset)
       interpolationReset=.true.
    end do
    !$omp end critical(CSMF_Behroozi2010_Tabulate)

    ! Compute the inverse stellar-halo mass relation if stellar mass has changed.
    if      (massStellar == massStellarPrevious(1)) then
       fMassHalo    =fMassHaloStored    (1)
       massSatellite=massSatelliteStored(1)
       massCut      =massCutStored      (1)
    else if (massStellar == massStellarPrevious(2)) then
       fMassHalo    =fMassHaloStored    (2)
       massSatellite=massSatelliteStored(2)
       massCut      =massCutStored      (2)
    else
       fMassHaloStored    (2)=fMassHaloStored    (1)
       massSatelliteStored(2)=massSatelliteStored(1)
       massCutStored      (2)=massCutStored      (1)
       massStellarPrevious(2)=massStellarPrevious(1)
       fMassHalo             =fSHMRInverse(massStellar)
       massSatellite         =massNormalization*conditionalStellarMassFunctionBehrooziBSatellite*(fMassHalo/massNormalization)&
            &                  **conditionalStellarMassFunctionBehrooziBetaSatellite
       massCut               =massNormalization*conditionalStellarMassFunctionBehrooziBCut      *(fMassHalo/massNormalization)&
            &                  **conditionalStellarMassFunctionBehrooziBetaCut       
       fMassHaloStored    (1)=fMassHalo
       massStellarPrevious(1)=massStellar
       massSatelliteStored(1)=massSatellite
       massCutStored      (1)=massCut
    end if

    ! Computed the forward stellar-halo mass relation is halo mass has changed.
    if (massHalo /= massHaloPrevious) then
       massHaloPrevious=massHalo
       fMassStellar=Interpolate(fMassStellarTableCount,fMassHaloTable,fMassStellarTable,interpolationObject&
            &,interpolationAccelerator,massHalo ,extrapolationType=extrapolationTypeNone,reset=interpolationReset)
    end if

    ! Compute the number of central galaxies.
    numberCentrals=                                                            &
         &          0.5d0                                                      &
         &         *(                                                          &
         &           1.0d0                                                     &
         &           -erf(                                                     &
         &                 log10(massStellar/fMassStellar)                     &
         &                /sqrt(2.0d0)                                        &
         &                /conditionalStellarMassFunctionBehrooziSigmaLogMstar &
         &               )                                                     &
         &          )

    ! Compute the number of satellites.
    numberSatellites=                                                                                &
         &            numberCentrals                                                                 &
         &           *(massHalo/massSatellite)**conditionalStellarMassFunctionBehrooziAlphaSatellite &
         &           *exp(-massCut/massHalo)

    return
  end subroutine Cumulative_Conditional_Stellar_Mass_Function_Compute

  elemental double precision function fSHMRInverse(massStellar)
    !% The median stellar mass vs. halo mass relation functional form from \cite{behroozi_comprehensive_2010}.
    implicit none
    double precision, intent(in   ) :: massStellar                  
    double precision, parameter     :: logHaloMassTransition=20.0d0 
    double precision                :: argument                     
    
    ! Compute the logarithmic halo mass for the given stellar mass.
    argument=                                                                                                                   &
         &    conditionalStellarMassFunctionBehrooziLog10M1                                                                     &
         &   +conditionalStellarMassFunctionBehrooziBeta*(log10(massStellar)-conditionalStellarMassFunctionBehrooziLog10Mstar0) &
         &   +(massStellar/Mstar0)**conditionalStellarMassFunctionBehrooziDelta                                                 &
         &   /(1.0d0+1.0d0/(massStellar/Mstar0)**conditionalStellarMassFunctionBehrooziGamma)                                   &
         &   -0.5d0
    ! For some parameter choices, the halo mass can grow unreasonably large. Therefore, above a transition value, allow the
    ! logarithmic halo mass to grow only logarithmically. Halo masses this high are irrelevant anyway (since the halo mass
    ! function will be so suppressed, while allowing the mass to continue to slowly grow allows for any tabulation to remain
    ! monotonically growing.
    if (argument > logHaloMassTransition) argument=logHaloMassTransition+log(argument/logHaloMassTransition)
    ! Compute the halo mass.
    fSHMRInverse=10.0d0**min(argument,100.0d0)
    return
  end function fSHMRInverse

end module Conditional_Stellar_Mass_Functions_Behroozi2010
