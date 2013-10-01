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

!% Contains a module which implements the \cite{behroozi_comprehensive_2010} fitting function descriptions of the conditional mass function.

module Conditional_Mass_Functions_Behroozi2010
  !% Implements the \cite{behroozi_comprehensive_2010} fitting function descriptions of the conditional mass function.
  use Tables
  private
  public :: Conditional_Mass_Functions_Behroozi2010_Initialize

  ! Parameters of the Behroozi et al. fitting function.
  double precision                                               :: conditionalMassFunctionBehrooziAlphaSatellite
  double precision                                               :: conditionalMassFunctionBehrooziLog10M1
  double precision                                               :: conditionalMassFunctionBehrooziLog10Mstar0
  double precision                                               :: conditionalMassFunctionBehrooziBeta
  double precision                                               :: conditionalMassFunctionBehrooziDelta
  double precision                                               :: conditionalMassFunctionBehrooziGamma
  double precision                                               :: conditionalMassFunctionBehrooziSigmaLogMstar
  double precision                                               :: conditionalMassFunctionBehrooziBCut
  double precision                                               :: conditionalMassFunctionBehrooziBSatellite
  double precision                                               :: conditionalMassFunctionBehrooziBetaCut
  double precision                                               :: conditionalMassFunctionBehrooziBetaSatellite
  double precision                                        :: Mstar0

  ! Tablulation of-halo mass relation.
  integer                                   , parameter   :: fMassTablePointsPerDecade                    =10
  integer                                                 :: fMassTableCount
  double precision                                        :: fMassTableMaximum                            =1.0d13 , fMassTableMinimum=1.0d8
  double precision                                        :: fMassHaloTableMaximum                                , fMassHaloTableMinimum
  type            (table1DLogarithmicLinear)              :: fMassTable
  class           (table1D                 ), allocatable :: fMassHaloTable

contains

  !# <conditionalMassFunctionMethod>
  !#  <unitName>Conditional_Mass_Functions_Behroozi2010_Initialize</unitName>
  !# </conditionalMassFunctionMethod>
  subroutine Conditional_Mass_Functions_Behroozi2010_Initialize(conditionalMassFunctionMethod&
       &,Cumulative_Conditional_Mass_Function_Get,Cumulative_Conditional_Mass_Function_Var_Get)
    !% Initializes the ``Behroozi2010'' conditional mass function method.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string  ), intent(in   )          :: conditionalMassFunctionMethod
    procedure(Cumulative_Conditional_Mass_Function_Behroozi2010), intent(inout), pointer :: Cumulative_Conditional_Mass_Function_Get
    procedure(Cumulative_Conditional_Mass_Function_Var_Behroozi2010), intent(inout), pointer :: Cumulative_Conditional_Mass_Function_Var_Get

    if (conditionalMassFunctionMethod == 'Behroozi2010') then
       Cumulative_Conditional_Mass_Function_Get     => Cumulative_Conditional_Mass_Function_Behroozi2010
       Cumulative_Conditional_Mass_Function_Var_Get => Cumulative_Conditional_Mass_Function_Var_Behroozi2010

       ! Get parameters of the fitting function.
       !@ <inputParameter>
       !@   <name>conditionalMassFunctionBehrooziAlphaSatellite</name>
       !@   <defaultValue>1.0 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\alpha_{\rm sat}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalMassFunctionBehrooziAlphaSatellite',conditionalMassFunctionBehrooziAlphaSatellite,defaultValue=1.0d0)

       !@ <inputParameter>
       !@   <name>conditionalMassFunctionBehrooziLog10M1</name>
       !@   <defaultValue>12.520 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\log_{10}M_1$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalMassFunctionBehrooziLog10M1',conditionalMassFunctionBehrooziLog10M1,defaultValue=12.520d0)

       !@ <inputParameter>
       !@   <name>conditionalMassFunctionBehrooziLog10Mstar0</name>
       !@   <defaultValue>10.916 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\log_{10}M_{\star,0}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalMassFunctionBehrooziLog10Mstar0',conditionalMassFunctionBehrooziLog10Mstar0,defaultValue=10.916d0)

       !@ <inputParameter>
       !@   <name>conditionalMassFunctionBehrooziBeta</name>
       !@   <defaultValue>0.457 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\beta$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalMassFunctionBehrooziBeta',conditionalMassFunctionBehrooziBeta,defaultValue=0.457d0)

       !@ <inputParameter>
       !@   <name>conditionalMassFunctionBehrooziDelta</name>
       !@   <defaultValue>0.5666 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\delta$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalMassFunctionBehrooziDelta',conditionalMassFunctionBehrooziDelta,defaultValue=0.5666d0)

       !@ <inputParameter>
       !@   <name>conditionalMassFunctionBehrooziGamma</name>
       !@   <defaultValue>1.53 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\gamma$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalMassFunctionBehrooziGamma',conditionalMassFunctionBehrooziGamma,defaultValue=1.53d0)

       !@ <inputParameter>
       !@   <name>conditionalMassFunctionBehrooziSigmaLogMstar</name>
       !@   <defaultValue>0.206 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\sigma_{\log M_\star}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalMassFunctionBehrooziSigmaLogMstar',conditionalMassFunctionBehrooziSigmaLogMstar,defaultValue=0.206d0)

       !@ <inputParameter>
       !@   <name>conditionalMassFunctionBehrooziBCut</name>
       !@   <defaultValue>1.47 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $B_{\rm cut}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalMassFunctionBehrooziBCut',conditionalMassFunctionBehrooziBCut,defaultValue=1.47d0)

       !@ <inputParameter>
       !@   <name>conditionalMassFunctionBehrooziBSatellite</name>
       !@   <defaultValue>10.62 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $B_{\rm sat}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalMassFunctionBehrooziBSatellite',conditionalMassFunctionBehrooziBSatellite,defaultValue=10.62d0)

       !@ <inputParameter>
       !@   <name>conditionalMassFunctionBehrooziBetaCut</name>
       !@   <defaultValue>$-0.13$ (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\beta_{\rm cut}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalMassFunctionBehrooziBetaCut',conditionalMassFunctionBehrooziBetaCut,defaultValue=-0.13d0)

       !@ <inputParameter>
       !@   <name>conditionalMassFunctionBehrooziBetaSatellite</name>
       !@   <defaultValue>0.859 (\citealt{leauthaud_new_2011}; $z_1$ sample using their {\tt SIG\_MOD1} method)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\beta_{\rm sat}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>haloModel</group>
       !@ </inputParameter>
       call Get_Input_Parameter('conditionalMassFunctionBehrooziBetaSatellite',conditionalMassFunctionBehrooziBetaSatellite,defaultValue=0.859d0)

       ! Compute the actual value of the Mstar0 parameter.
       Mstar0=10.0d0**conditionalMassFunctionBehrooziLog10Mstar0
    end if
    return
  end subroutine Conditional_Mass_Functions_Behroozi2010_Initialize

  double precision function Cumulative_Conditional_Mass_Function_Behroozi2010(massHalo,mass)
    !% Computes the cumulative conditional mass function, $\langle N(M_\star|M_{\rm halo}) \rangle \equiv \phi(M_\star|M_{\rm
    !% halo})$ using the fitting formula of \cite{behroozi_comprehensive_2010}.
    implicit none
    double precision, intent(in   ) :: massHalo      , mass
    double precision                :: numberCentrals, numberSatellites

    ! Get the number of satellites and centrals.
    call Cumulative_Conditional_Mass_Function_Compute(massHalo,mass,numberCentrals,numberSatellites)

    ! Sum central and satellite contributions.
    Cumulative_Conditional_Mass_Function_Behroozi2010=numberCentrals+numberSatellites
    return
  end function Cumulative_Conditional_Mass_Function_Behroozi2010

  double precision function Cumulative_Conditional_Mass_Function_Var_Behroozi2010(massHalo,massLow,massHigh)
    !% Computes the variance in the cumulative conditional mass function, $\langle N(M_\star|M_{\rm halo}) \rangle \equiv
    !% \phi(M_\star|M_{\rm halo})$ using the fitting formula of \cite{behroozi_comprehensive_2010}. Assumes that the number of
    !% satellite galaxies is Poisson distributed, while the number of central galaxies follows a Bernoulli distribution, and that
    !% the numbers of satellites and centrals are uncorrelated.
    implicit none
    double precision, intent(in   ) :: massHalo        , massHigh     , massLow
    double precision                :: numberCentrals  , numberCentralsHigh  , numberCentralsLow  , &
         &                             numberSatellites, numberSatellitesHigh, numberSatellitesLow

    ! Get the number of satellites and centrals.
    call Cumulative_Conditional_Mass_Function_Compute(massHalo,massLow ,numberCentralsLow ,numberSatellitesLow )
    call Cumulative_Conditional_Mass_Function_Compute(massHalo,massHigh,numberCentralsHigh,numberSatellitesHigh)
    numberSatellites=max(numberSatellitesLow-numberSatellitesHigh,0.0d0)
    numberCentrals  =max(numberCentralsLow  -numberCentralsHigh  ,0.0d0)

    ! Compute the variance (using the Bienaymé formula for uncorrelated variables).
    Cumulative_Conditional_Mass_Function_Var_Behroozi2010= &
         &  numberSatellites                                       & ! Satellites are Poisson distributed, so the variance is just their number.
         & +numberCentrals*(1.0d0-numberCentrals)                    ! Centrals are Bernoulli distributed.
    return
  end function Cumulative_Conditional_Mass_Function_Var_Behroozi2010

  subroutine Cumulative_Conditional_Mass_Function_Compute(massHalo,mass,numberCentrals,numberSatellites)
    !% Computes the cumulative conditional mass function, $\langle N(M_\star|M_{\rm halo}) \rangle \equiv \phi(M_\star|M_{\rm
    !% halo})$ using the fitting formula of \cite{behroozi_comprehensive_2010}.
    implicit none
    double precision                    , intent(in   ) :: massHalo                  , mass
    double precision                    , intent(  out) :: numberCentrals            , numberSatellites
    double precision, parameter                         :: massNormalization  =1.0d12
    double precision                                    :: fMassHalo                 , massCut                   , &
         &                                                 massSatellite
    double precision              , save                :: fMass              , massHaloPrevious   =-1.0d0
    !$omp threadprivate(massHaloPrevious,fMass)
    double precision, dimension(2), save                :: fMassHaloStored           , massCutStored             , &
         &                                                 massSatelliteStored       , massPrevious=-1.0d0
    !$omp threadprivate(massPrevious,fMassHaloStored,massSatelliteStored,massCutStored)

    !$omp critical(CSMF_Behroozi2010_Tabulate)
    do while (massHalo < fMassHaloTableMinimum .or. massHalo > fMassHaloTableMaximum)
       if (massHalo < fMassHaloTableMinimum) fMassTableMinimum=0.5d0*fMassTableMinimum
       if (massHalo > fMassHaloTableMaximum) fMassTableMaximum=2.0d0*fMassTableMaximum
       fMassTableCount=int(log10(fMassTableMaximum/fMassTableMinimum)*fMassTablePointsPerDecade)+1
       call fMassTable%destroy ()
       call fMassTable%create  (                                          &
            &                   fMassTableMinimum                       , &
            &                   fMassTableMaximum                       , &
            &                   fMassTableCount                         , &
            &                   extrapolationType=extrapolationTypeAbort  &
            &                  )
       call fMassTable%populate(                                          &
            &                   fSHMRInverse(fMassTable%xs())             &
            &                  )
       call fMassTable%reverse (                                          &
            &                   fMassHaloTable                            &
            &                  )
       fMassHaloTableMinimum=fMassTable%y(+1)
       fMassHaloTableMaximum=fMassTable%y(-1)
    end do
    !$omp end critical(CSMF_Behroozi2010_Tabulate)

    ! Compute the inverse-halo mass relation if mass has changed.
    if      (mass == massPrevious(1)) then
       fMassHalo    =fMassHaloStored    (1)
       massSatellite=massSatelliteStored(1)
       massCut      =massCutStored      (1)
    else if (mass == massPrevious(2)) then
       fMassHalo    =fMassHaloStored    (2)
       massSatellite=massSatelliteStored(2)
       massCut      =massCutStored      (2)
    else
       fMassHaloStored    (2)=fMassHaloStored    (1)
       massSatelliteStored(2)=massSatelliteStored(1)
       massCutStored      (2)=massCutStored      (1)
       massPrevious(2)=massPrevious(1)
       fMassHalo             =fSHMRInverse(mass)
       massSatellite         =massNormalization*conditionalMassFunctionBehrooziBSatellite*(fMassHalo/massNormalization)&
            &                  **conditionalMassFunctionBehrooziBetaSatellite
       massCut               =massNormalization*conditionalMassFunctionBehrooziBCut      *(fMassHalo/massNormalization)&
            &                  **conditionalMassFunctionBehrooziBetaCut
       fMassHaloStored    (1)=fMassHalo
       massPrevious(1)=mass
       massSatelliteStored(1)=massSatellite
       massCutStored      (1)=massCut
    end if

    ! Computed the forward-halo mass relation is halo mass has changed.
    if (massHalo /= massHaloPrevious) then
       massHaloPrevious=massHalo
       fMass=fMassHaloTable%interpolate(massHalo)
    end if

    ! Compute the number of central galaxies.
    numberCentrals=                                                     &
         &          0.5d0                                               &
         &         *(                                                   &
         &           1.0d0                                              &
         &           -erf(                                              &
         &                 log10(mass/fMass)                            &
         &                /sqrt(2.0d0)                                  &
         &                /conditionalMassFunctionBehrooziSigmaLogMstar &
         &               )                                              &
         &          )

    ! Compute the number of satellites.
    numberSatellites=                                                                         &
         &            numberCentrals                                                          &
         &           *(massHalo/massSatellite)**conditionalMassFunctionBehrooziAlphaSatellite &
         &           *exp(-massCut/massHalo)

    return
  end subroutine Cumulative_Conditional_Mass_Function_Compute

  elemental double precision function fSHMRInverse(mass)
    !% The median mass vs. halo mass relation functional form from \cite{behroozi_comprehensive_2010}.
    implicit none
    double precision, intent(in   ) :: mass
    double precision, parameter     :: logHaloMassTransition=20.0d0
    double precision                :: argument

    ! Compute the logarithmic halo mass for the given mass.
    argument=                                                                                              &
         &    conditionalMassFunctionBehrooziLog10M1                                                       &
         &   +conditionalMassFunctionBehrooziBeta*(log10(mass)-conditionalMassFunctionBehrooziLog10Mstar0) &
         &   +(mass/Mstar0)**conditionalMassFunctionBehrooziDelta                                          &
         &   /(1.0d0+1.0d0/(mass/Mstar0)**conditionalMassFunctionBehrooziGamma)                            &
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

end module Conditional_Mass_Functions_Behroozi2010
