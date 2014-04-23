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

!% Implements calculation of conditional mass functions using the \cite{behroozi_comprehensive_2010} fitting function.

  !# <conditionalMassFunction name="conditionalMassFunctionBehroozi2010">
  !#  <description>Computes the conditional mass function using the fiting functions of \cite{behroozi_comprehensive_2010}.</description>
  !# </conditionalMassFunction>

  use Tables

  type, extends(conditionalMassFunctionClass) :: conditionalMassFunctionBehroozi2010
     ! Parameters of the fitting function.
     double precision                                         :: alphaSatellite
     double precision                                         :: log10M1
     double precision                                         :: Mstar0
     double precision                                         :: beta
     double precision                                         :: delta
     double precision                                         :: gamma
     double precision                                         :: sigmaLogMstar
     double precision                                         :: BCut
     double precision                                         :: BSatellite
     double precision                                         :: betaCut
     double precision                                         :: betaSatellite
     ! Accelerators.
     double precision                                         :: fMass                , massHaloPrevious
     double precision                          , dimension(2) :: massSatelliteStored  , massPrevious    , &
         &                                                       fMassHaloStored      , massCutStored
     ! Tablulation of-halo mass relation.
     integer                                                  :: fMassTableCount
     double precision                                         :: fMassTableMaximum    , fMassTableMinimum
     double precision                                         :: fMassHaloTableMaximum, fMassHaloTableMinimum
     type            (table1DLogarithmicLinear)               :: fMassTable
     class           (table1D                 ), allocatable  :: fMassHaloTable
   contains
     final     ::                         behroozi2010Destructor
     procedure :: massFunction         => behroozi2010MassFunction
     procedure :: massFunctionVariance => behroozi2010MassFunctionVariance
     procedure :: compute              => behroozi2010Compute
  end type conditionalMassFunctionBehroozi2010
  
  interface conditionalMassFunctionBehroozi2010
     !% Constructors for the \cite{behroozi_comprehensive_2010} merging timescale class.
     module procedure behroozi2010DefaultConstructor
     module procedure behroozi2010Constructor
  end interface conditionalMassFunctionBehroozi2010

  ! Table resolution.
  integer                                              , parameter :: behroozi2010MassTablePointsPerDecade=10

  ! Module scope pointer to the current object.
  class           (conditionalMassFunctionBehroozi2010), pointer   :: behroozi2010SelfGlobal
  !$omp threadprivate(behroozi2010SelfGlobal)

  ! Default parameter values.
  logical                                                          :: behroozi2010Initialized             =.false.
  double precision                                                 :: behrooziAlphaSatellite                      , behrooziLog10M1, &
       &                                                              behrooziLog10Mstar0                         , behrooziBeta   , &
       &                                                              behrooziDelta                               , behrooziGamma  , &
       &                                                              behrooziSigmaLogMstar                       , behrooziBCut   , &
       &                                                              behrooziBSatellite                          , behrooziBetaCut, &
       &                                                              behrooziBetaSatellite

contains

  function behroozi2010DefaultConstructor()
    !% Default constructor for the \cite{behroozi_comprehensive_2010} conditional mass function class.
    use Input_Parameters
    implicit none
    type(conditionalMassFunctionBehroozi2010) :: behroozi2010DefaultConstructor

    if (.not.behroozi2010Initialized) then
       !$omp critical (conditionalMassFunctionBehroozi2010Initialize)
       if (.not.behroozi2010Initialized) then
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
          call Get_Input_Parameter('conditionalMassFunctionBehrooziAlphaSatellite',behrooziAlphaSatellite,defaultValue=1.0d0)
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
          call Get_Input_Parameter('conditionalMassFunctionBehrooziLog10M1',behrooziLog10M1,defaultValue=12.520d0)          
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
          call Get_Input_Parameter('conditionalMassFunctionBehrooziLog10Mstar0',behrooziLog10Mstar0,defaultValue=10.916d0)
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
          call Get_Input_Parameter('conditionalMassFunctionBehrooziBeta',behrooziBeta,defaultValue=0.457d0)
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
          call Get_Input_Parameter('conditionalMassFunctionBehrooziDelta',behrooziDelta,defaultValue=0.5666d0)
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
          call Get_Input_Parameter('conditionalMassFunctionBehrooziGamma',behrooziGamma,defaultValue=1.53d0)
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
          call Get_Input_Parameter('conditionalMassFunctionBehrooziSigmaLogMstar',behrooziSigmaLogMstar,defaultValue=0.206d0)
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
          call Get_Input_Parameter('conditionalMassFunctionBehrooziBCut',behrooziBCut,defaultValue=1.47d0)
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
          call Get_Input_Parameter('conditionalMassFunctionBehrooziBSatellite',behrooziBSatellite,defaultValue=10.62d0)
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
          call Get_Input_Parameter('conditionalMassFunctionBehrooziBetaCut',behrooziBetaCut,defaultValue=-0.13d0)
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
          call Get_Input_Parameter('conditionalMassFunctionBehrooziBetaSatellite',behrooziBetaSatellite,defaultValue=0.859d0)       
           ! Record that we are now initialized.
          behroozi2010Initialized=.true.
       end if
       !$omp end critical (conditionalMassFunctionBehroozi2010Initialize)
    end if
    behroozi2010DefaultConstructor=behroozi2010Constructor(behrooziAlphaSatellite,behrooziLog10M1,behrooziLog10Mstar0,behrooziBeta,behrooziDelta,behrooziGamma,behrooziSigmaLogMstar,behrooziBCut,behrooziBSatellite,behrooziBetaCut,behrooziBetaSatellite)
   return
  end function behroozi2010DefaultConstructor

  function behroozi2010Constructor(behrooziAlphaSatellite,behrooziLog10M1,behrooziLog10Mstar0,behrooziBeta,behrooziDelta,behrooziGamma,behrooziSigmaLogMstar,behrooziBCut,behrooziBSatellite,behrooziBetaCut,behrooziBetaSatellite)
    !% Generic constructor for the \cite{behroozi_comprehensive_2010} conditional mass function class.
    implicit none
    type            (conditionalMassFunctionBehroozi2010)                :: behroozi2010Constructor
    double precision                                     , intent(in   ) :: behrooziAlphaSatellite , behrooziLog10M1, &
         &                                                                  behrooziLog10Mstar0    , behrooziBeta   , & 
         &                                                                  behrooziDelta          , behrooziGamma  , &
         &                                                                  behrooziSigmaLogMstar  , behrooziBCut   , &
         &                                                                  behrooziBSatellite     , behrooziBetaCut, &
         &                                                                  behrooziBetaSatellite
    
    behroozi2010Constructor%alphaSatellite=        behrooziAlphaSatellite
    behroozi2010Constructor%log10M1       =        behrooziLog10M1
    behroozi2010Constructor%Mstar0        =10.0d0**behrooziLog10Mstar0
    behroozi2010Constructor%beta          =        behrooziBeta
    behroozi2010Constructor%delta         =        behrooziDelta
    behroozi2010Constructor%gamma         =        behrooziGamma
    behroozi2010Constructor%sigmaLogMstar =        behrooziSigmaLogMstar
    behroozi2010Constructor%BCut          =        behrooziBCut
    behroozi2010Constructor%BSatellite    =        behrooziBSatellite
    behroozi2010Constructor%betaCut       =        behrooziBetaCut
    behroozi2010Constructor%betaSatellite =        behrooziBetaSatellite
    ! Initialize tables and accelerators.
    behroozi2010Constructor%massPrevious     =-1.0d00
    behroozi2010Constructor%massHaloPrevious =-1.0d00
    behroozi2010Constructor%fMassTableMinimum=+1.0d08
    behroozi2010Constructor%fMassTableMaximum=+1.0d13
    return
  end function behroozi2010Constructor

  subroutine behroozi2010Destructor(self)
    !% Destructor for the \cite{{behroozi_comprehensive_2010} conditional mass function class.
     use Gaussian_Random
     implicit none
     type(conditionalMassFunctionBehroozi2010), intent(inout) :: self

     call                                     self%fMassTable    %destroy()
     if (allocated(self%fMassHaloTable)) call self%fMassHaloTable%destroy()
    return
  end subroutine behroozi2010Destructor

  double precision function behroozi2010MassFunction(self,massHalo,mass)
    !% Computes the cumulative conditional mass function, $\langle N(M_\star|M_{\rm halo}) \rangle \equiv \phi(M_\star|M_{\rm
    !% halo})$ using the fitting formula of \cite{behroozi_comprehensive_2010}.
    implicit none
    class           (conditionalMassFunctionBehroozi2010), intent(inout) :: self
    double precision                                     , intent(in   ) :: massHalo      , mass
    double precision                                                     :: numberCentrals, numberSatellites

    ! Get the number of satellites and centrals.
    call self%compute(massHalo,mass,numberCentrals,numberSatellites)
    ! Sum central and satellite contributions.
    behroozi2010MassFunction=numberCentrals+numberSatellites
    return
  end function behroozi2010MassFunction

  double precision function behroozi2010MassFunctionVariance(self,massHalo,massLow,massHigh)
    !% Computes the variance in the cumulative conditional mass function, $\langle N(M_\star|M_{\rm halo}) \rangle \equiv
    !% \phi(M_\star|M_{\rm halo})$ using the fitting formula of \cite{behroozi_comprehensive_2010}. Assumes that the number of
    !% satellite galaxies is Poisson distributed, while the number of central galaxies follows a Bernoulli distribution, and that
    !% the numbers of satellites and centrals are uncorrelated.
    implicit none
    class           (conditionalMassFunctionBehroozi2010), intent(inout) :: self
    double precision                                     , intent(in   ) :: massHalo        , massHigh     , massLow
    double precision                                                     :: numberCentrals  , numberCentralsHigh  , numberCentralsLow  , &
         &                                                                  numberSatellites, numberSatellitesHigh, numberSatellitesLow

    ! Get the number of satellites and centrals.
    call self%compute(massHalo,massLow ,numberCentralsLow ,numberSatellitesLow )
    call self%compute(massHalo,massHigh,numberCentralsHigh,numberSatellitesHigh)
    numberSatellites=max(numberSatellitesLow-numberSatellitesHigh,0.0d0)
    numberCentrals  =max(numberCentralsLow  -numberCentralsHigh  ,0.0d0)
    ! Compute the variance (using the Bienaymé formula for uncorrelated variables).
    behroozi2010MassFunctionVariance=             &
         &  numberSatellites                      & ! Satellites are Poisson distributed, so the variance is just their number.
         & +numberCentrals*(1.0d0-numberCentrals)   ! Centrals are Bernoulli distributed.
    return
  end function behroozi2010MassFunctionVariance

  subroutine behroozi2010Compute(self,massHalo,mass,numberCentrals,numberSatellites)
    !% Computes the cumulative conditional mass function, $\langle N(M_\star|M_{\rm halo}) \rangle \equiv \phi(M_\star|M_{\rm
    !% halo})$ using the fitting formula of \cite{behroozi_comprehensive_2010}.
    implicit none
    class           (conditionalMassFunctionBehroozi2010), intent(inout), target :: self
    double precision                                     , intent(in   )         :: massHalo                , mass
    double precision                                     , intent(  out)         :: numberCentrals          , numberSatellites
    double precision                                     , parameter             :: massNormalization=1.0d12
    double precision                                     , parameter             :: massMinimum      =1.0d-12
    double precision                                                             :: fMassHalo               , massCut         , &
         &                                                                          massSatellite

    behroozi2010SelfGlobal => self
    do while (                                            &
         &     .not.allocated(self%fMassHaloTable)        &
         &    .or.                                        &
         &     (                                          &
         &       massHalo    < self%fMassHaloTableMinimum &
         &      .and.                                     &
         &       massMinimum < self%fMassTableMinimum     &
         &     )                                          &
         &    .or.                                        &
         &       massHalo    > self%fMassHaloTableMaximum &
         &   )
       if (allocated(self%fMassHaloTable)) then
          if (massHalo < self%fMassHaloTableMinimum) self%fMassTableMinimum=0.5d0*self%fMassTableMinimum
          if (massHalo > self%fMassHaloTableMaximum) self%fMassTableMaximum=2.0d0*self%fMassTableMaximum
       end if
       self%fMassTableCount=int(log10(self%fMassTableMaximum/self%fMassTableMinimum)*behroozi2010MassTablePointsPerDecade)+1
       call self%fMassTable%destroy ()
       call self%fMassTable%create  (                                               &
            &                        self%fMassTableMinimum                       , &
            &                        self%fMassTableMaximum                       , &
            &                        self%fMassTableCount                         , &
            &                        extrapolationType=extrapolationTypeFix       &
            &                       )
       call self%fMassTable%populate(                                               &
            &                        behroozi2010fSHMRInverse(self%fMassTable%xs()) &
            &                       )
       call self%fMassTable%reverse (                                               &
            &                        self%fMassHaloTable                            &
            &                       )
       self%fMassHaloTableMinimum=self%fMassTable%y(+1)
       self%fMassHaloTableMaximum=self%fMassTable%y(-1)
    end do
    ! Compute the inverse-halo mass relation if mass has changed.
    if      (mass == self%massPrevious(1)) then
       fMassHalo    =self%fMassHaloStored    (1)
       massSatellite=self%massSatelliteStored(1)
       massCut      =self%massCutStored      (1)
    else if (mass == self%massPrevious(2)) then
       fMassHalo    =self%fMassHaloStored    (2)
       massSatellite=self%massSatelliteStored(2)
       massCut      =self%massCutStored      (2)
    else
       self%fMassHaloStored    (2)=self%fMassHaloStored    (1)
       self%massSatelliteStored(2)=self%massSatelliteStored(1)
       self%massCutStored      (2)=self%massCutStored      (1)
       self%massPrevious       (2)=self%massPrevious       (1)
       fMassHalo                  =behroozi2010fSHMRInverse(mass)
       massSatellite              =massNormalization*self%BSatellite*(fMassHalo/massNormalization)**self%betaSatellite
       massCut                    =massNormalization*self%BCut      *(fMassHalo/massNormalization)**self%betaCut
       self%fMassHaloStored    (1)=fMassHalo
       self%massPrevious       (1)=mass
       self%massSatelliteStored(1)=massSatellite
       self%massCutStored      (1)=massCut
    end if
    ! Computed the forward-halo mass relation is halo mass has changed.
    if (massHalo /= self%massHaloPrevious) then
       self%massHaloPrevious=                                massHalo
       self%fMass           =self%fMassHaloTable%interpolate(massHalo)
    end if
    ! Compute the number of central galaxies.
    numberCentrals=                               &
         &          0.5d0                         &
         &         *(                             &
         &           1.0d0                        &
         &           -erf(                        &
         &                 log10(mass/self%fMass) &
         &                /sqrt(2.0d0)            &
         &                /self%sigmaLogMstar     &
         &               )                        &
         &          )
    ! Compute the number of satellites.
    numberSatellites=                                               &
         &            numberCentrals                                &
         &           *(massHalo/massSatellite)**self%alphaSatellite &
         &           *exp(-massCut/massHalo)
    return
  end subroutine behroozi2010Compute

  elemental double precision function behroozi2010fSHMRInverse(mass)
    !% The median mass vs. halo mass relation functional form from \cite{behroozi_comprehensive_2010}.
    implicit none
    double precision, intent(in   ) :: mass
    double precision, parameter     :: logHaloMassTransition=20.0d0
    double precision                :: argument

    ! Compute the logarithmic halo mass for the given mass.
    argument=                                                                                  &
         &    behroozi2010SelfGlobal%log10M1                                                   &
         &   +behroozi2010SelfGlobal%beta   *log10(mass/behroozi2010SelfGlobal%Mstar0)         &
         &   +             (mass/behroozi2010SelfGlobal%Mstar0)**behroozi2010SelfGlobal%delta  &
         &   /(1.0d0+1.0d0/(mass/behroozi2010SelfGlobal%Mstar0)**behroozi2010SelfGlobal%gamma) &
         &   -0.5d0
    ! For some parameter choices, the halo mass can grow unreasonably large. Therefore, above a transition value, allow the
    ! logarithmic halo mass to grow only logarithmically. Halo masses this high are irrelevant anyway (since the halo mass
    ! function will be so suppressed, while allowing the mass to continue to slowly grow allows for any tabulation to remain
    ! monotonically growing.
    if (argument > logHaloMassTransition) argument=logHaloMassTransition+log(argument/logHaloMassTransition)
    ! Compute the halo mass.
    behroozi2010fSHMRInverse=10.0d0**min(argument,100.0d0)
    return
  end function behroozi2010fSHMRInverse
