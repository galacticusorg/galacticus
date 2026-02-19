!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!!{
Implements a class for the conditional mass functions using the \cite{behroozi_comprehensive_2010} fitting function.
!!}

  use :: Tables, only : table1D, table1DLogarithmicLinear

  !![
  <conditionalMassFunction name="conditionalMassFunctionBehroozi2010">
   <description>
    A conditional mass function class which implements the fitting functions of \cite{behroozi_comprehensive_2010}:
    \begin{equation}
     \langle N_\mathrm{c}(M_\star|M)\rangle \equiv \int_{M_\star}^\infty \phi_\mathrm{c}(M_\star^\prime) \d \ln M_\star^\prime
     = {1 \over 2} \left[ 1 - \hbox{erf}\left( {\log_{10}M_\star - \log_{10} f_\mathrm{SHMR}(M) \over \sqrt{2}\sigma_{\log
     M_\star}} \right) \right].
    \end{equation}
    Here, the function $f_\mathrm{SHMR}(M)$ is the solution of
    \begin{equation}
     \log_{10}M = \log_{10}M_1 + \beta \log_{10}\left({M_\star \over M_{\star,0}}\right) + {(M_\star/M_{\star,0})^\delta \over
     1 + (M_\star/M_{\star,0})^{-\gamma}} - {1/2}.
    \end{equation}
    For satellites,
    \begin{equation}
     \langle N_\mathrm{s}(M_\star|M)\rangle \equiv \int_{M_\star}^\infty \phi_\mathrm{s}(M_\star^\prime) \d \ln M_\star^\prime
     = \langle N_\mathrm{c}(M_\star|M)\rangle \left({f^{-1}_\mathrm{SHMR}(M_\star) \over
     M_\mathrm{sat}}\right)^{\alpha_\mathrm{sat}} \exp\left(- {M_\mathrm{cut} \over f^{-1}_\mathrm{SHMR}(M_\star)} \right),
    \end{equation}
    where
    \begin{equation}
     {M_\mathrm{sat} \over 10^{12} M_\odot} = B_\mathrm{sat} \left({f^{-1}_\mathrm{SHMR}(M_\star) \over 10^{12}
     M_\odot}\right)^{\beta_\mathrm{sat}},
    \end{equation}
    and
    \begin{equation}
     {M_\mathrm{cut} \over 10^{12} M_\odot} = B_\mathrm{cut} \left({f^{-1}_\mathrm{SHMR}(M_\star) \over 10^{12}
     M_\odot}\right)^{\beta_\mathrm{cut}}.
    \end{equation}
    By default, parameter values are taken from the fit of \cite{leauthaud_new_2011}, specifically their {\normalfont \ttfamily
    SIG\_MOD1} method for their $z_1$ sample. These default values, and the \glc\ input parameters which can be used to adjust
    them are shown in Table~\ref{table:Behroozi2010FitParameters}. This method assumes that $P_\mathrm{s}(N|M_\star,M;\delta
    \ln M_\star)$ is a Poisson distribution while $P_\mathrm{c}(N|M_\star,M;\delta \ln M_\star)$ has a Bernoulli distribution,
    with each distribution's free parameter fixed by requiring
    \begin{equation}
     \phi(M_\star;M) \delta \ln M_\star = \sum_{N=0}^\infty N P(N|M_\star,M;\delta \ln M_\star)
    \end{equation}
    \begin{table}
    \caption{Parameters of the \cite{behroozi_comprehensive_2010} conditional stellar mass function model, along with their
    default values and the corresponding \glc\ input parameters.}
    \label{table:Behroozi2010FitParameters}
    \begin{center}
    \begin{tabular}{lr@{.}ll}
    \hline
    {\normalfont \bfseries Parameter} &amp; \multicolumn{2}{c}{{\normalfont \bfseries Default}} &amp; {\normalfont \bfseries \glc\ name} \\
    \hline
    $\alpha_\mathrm{sat}$&amp; 1&amp;0&amp; {\normalfont \ttfamily [conditionalStellarMassFunctionBehrooziAlphaSatellite]} \\
    $\log_{10} M_1$&amp; 12&amp;520&amp; {\normalfont \ttfamily [conditionalStellarMassFunctionBehrooziLog10M1]} \\
    $\log_{10} M_{\star,0}$&amp; 10&amp;916&amp; {\normalfont \ttfamily [conditionalStellarMassFunctionBehrooziLog10Mstar0]} \\
    $\beta$&amp; 0&amp;457&amp; {\normalfont \ttfamily [conditionalStellarMassFunctionBehrooziBeta]} \\
    $\delta$&amp; 0&amp;5666&amp; {\normalfont \ttfamily [conditionalStellarMassFunctionBehrooziDelta]} \\
    $\gamma$&amp; 1&amp;53&amp; {\normalfont \ttfamily [conditionalStellarMassFunctionBehrooziGamma]} \\
    $\sigma_{\log M_\star}$&amp; 0&amp;206&amp; {\normalfont \ttfamily [conditionalStellarMassFunctionBehrooziSigmaLogMstar]} \\
    $B_\mathrm{cut}$&amp; 1&amp;47&amp; {\normalfont \ttfamily [conditionalStellarMassFunctionBehrooziBCut]} \\
    $B_\mathrm{sat}$&amp; 10&amp;62&amp; {\normalfont \ttfamily [conditionalStellarMassFunctionBehrooziBSatellite]} \\
    $\beta_\mathrm{cut}$&amp; $-$0&amp;13&amp; {\normalfont \ttfamily [conditionalStellarMassFunctionBehrooziBetaCut]} \\
    $\beta_\mathrm{sat}$&amp; 0&amp;859&amp; {\normalfont \ttfamily [conditionalStellarMassFunctionBehrooziBetaCut]} \\
    \hline
    \end{tabular}
    \end{center}
    \end{table}
   </description>
  </conditionalMassFunction>
  !!]
  type, extends(conditionalMassFunctionClass) :: conditionalMassFunctionBehroozi2010
     !!{
     Implements the conditional mass function using the fitting functions of \cite{behroozi_comprehensive_2010}.
     !!}
     ! Parameters of the fitting function.
     private
     double precision                                         :: alphaSatellite
     double precision                                         :: log10M1
     double precision                                         :: log10Mstar0   , Mstar0
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
     ! Tabulation of-halo mass relation.
     integer                                                  :: fMassTableCount
     double precision                                         :: fMassTableMaximum    , fMassTableMinimum
     double precision                                         :: fMassHaloTableMaximum, fMassHaloTableMinimum
     type            (table1DLogarithmicLinear)               :: fMassTable
     class           (table1D                 ), allocatable  :: fMassHaloTable
   contains
     !![
     <methods>
       <method description="Compute the cumulative conditional mass function, $\langle N(M_\star|M_\mathrm{halo}) \rangle \equiv \phi(M_\star|M_\mathrm{halo})$." method="compute" />
     </methods>
     !!]
     final     ::                         behroozi2010Destructor
     procedure :: massFunction         => behroozi2010MassFunction
     procedure :: massFunctionVariance => behroozi2010MassFunctionVariance
     procedure :: compute              => behroozi2010Compute
  end type conditionalMassFunctionBehroozi2010

  interface conditionalMassFunctionBehroozi2010
     !!{
     Constructors for the \cite{behroozi_comprehensive_2010} merging timescale class.
     !!}
     module procedure behroozi2010ConstructorParameters
     module procedure behroozi2010ConstructorInternal
  end interface conditionalMassFunctionBehroozi2010

  ! Table resolution.
  integer                                              , parameter :: massTablePointsPerDecade=10

  ! Module scope pointer to the current object.
  class           (conditionalMassFunctionBehroozi2010), pointer   :: self_
  !$omp threadprivate(self_)

contains


  function behroozi2010ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{conditionalMassFunctionBehroozi2010} conditional mass function class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (conditionalMassFunctionBehroozi2010)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    double precision                                                     :: alphaSatellite, log10M1, &
         &                                                                  log10Mstar0   , beta   , &
         &                                                                  delta         , gamma  , &
         &                                                                  sigmaLogMstar , BCut   , &
         &                                                                  BSatellite    , betaCut, &
         &                                                                  betaSatellite

    !![
    <inputParameter>
      <name>alphaSatellite</name>
      <defaultSource>(\citealt{leauthaud_new_2011}; $z_1$ sample using their {\normalfont \ttfamily SIG\_MOD1} method)</defaultSource>
      <defaultValue>1.0d0</defaultValue>
      <description>The parameter $\alpha_\mathrm{sat}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>log10M1</name>
      <defaultSource>(\citealt{leauthaud_new_2011}; $z_1$ sample using their {\normalfont \ttfamily SIG\_MOD1} method)</defaultSource>
      <defaultValue>12.520d0</defaultValue>
      <description>The parameter $\log_{10}M_1$ from the fitting functions of \cite{behroozi_comprehensive_2010}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>log10Mstar0</name>
      <defaultSource>(\citealt{leauthaud_new_2011}; $z_1$ sample using their {\normalfont \ttfamily SIG\_MOD1} method)</defaultSource>
      <defaultValue>10.916d0</defaultValue>
      <description>The parameter $\log_{10}M_{\star,0}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <defaultSource>(\citealt{leauthaud_new_2011}; $z_1$ sample using their {\normalfont \ttfamily SIG\_MOD1} method)</defaultSource>
      <defaultValue>0.457d0</defaultValue>
      <description>The parameter $\beta$ from the fitting functions of \cite{behroozi_comprehensive_2010}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>delta</name>
      <defaultSource>(\citealt{leauthaud_new_2011}; $z_1$ sample using their {\normalfont \ttfamily SIG\_MOD1} method)</defaultSource>
      <defaultValue>0.5666d0</defaultValue>
      <description>The parameter $\delta$ from the fitting functions of \cite{behroozi_comprehensive_2010}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>gamma</name>
      <defaultSource>(\citealt{leauthaud_new_2011}; $z_1$ sample using their {\normalfont \ttfamily SIG\_MOD1} method)</defaultSource>
      <defaultValue>1.53d0</defaultValue>
      <description>The parameter $\gamma$ from the fitting functions of \cite{behroozi_comprehensive_2010}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>sigmaLogMstar</name>
      <defaultSource>(\citealt{leauthaud_new_2011}; $z_1$ sample using their {\normalfont \ttfamily SIG\_MOD1} method)</defaultSource>
      <defaultValue>0.206d0</defaultValue>
      <description>The parameter $\sigma_{\log M_\star}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>BCut</name>
      <defaultSource>(\citealt{leauthaud_new_2011}; $z_1$ sample using their {\normalfont \ttfamily SIG\_MOD1} method)</defaultSource>
      <defaultValue>1.47d0</defaultValue>
      <description>The parameter $B_\mathrm{cut}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>BSatellite</name>
      <defaultSource>(\citealt{leauthaud_new_2011}; $z_1$ sample using their {\normalfont \ttfamily SIG\_MOD1} method)</defaultSource>
      <defaultValue>10.62d0</defaultValue>
      <description>The parameter $B_\mathrm{sat}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>betaCut</name>
      <defaultSource>(\citealt{leauthaud_new_2011}; $z_1$ sample using their {\normalfont \ttfamily SIG\_MOD1} method)</defaultSource>
      <defaultValue>-0.13d0</defaultValue>
      <description>The parameter $\beta_\mathrm{cut}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>betaSatellite</name>
      <defaultSource>(\citealt{leauthaud_new_2011}; $z_1$ sample using their {\normalfont \ttfamily SIG\_MOD1} method)</defaultSource>
      <defaultValue>0.859d0</defaultValue>
      <description>The parameter $\beta_\mathrm{sat}$ from the fitting functions of \cite{behroozi_comprehensive_2010}.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=conditionalMassFunctionBehroozi2010(alphaSatellite,log10M1,log10Mstar0,beta,delta,gamma,sigmaLogMstar,BCut,BSatellite,betaCut,betaSatellite)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function behroozi2010ConstructorParameters

  function behroozi2010ConstructorInternal(alphaSatellite,log10M1,log10Mstar0,beta,delta,gamma,sigmaLogMstar,BCut,BSatellite,betaCut,betaSatellite) result(self)
    !!{
    Internal constructor for the \cite{behroozi_comprehensive_2010} conditional mass function class.
    !!}
    implicit none
    type            (conditionalMassFunctionBehroozi2010)                :: self
    double precision                                     , intent(in   ) :: alphaSatellite, log10M1, &
         &                                                                  log10Mstar0   , beta   , &
         &                                                                  delta         , gamma  , &
         &                                                                  sigmaLogMstar , BCut   , &
         &                                                                  BSatellite    , betaCut, &
         &                                                                  betaSatellite
    !![
    <constructorAssign variables="alphaSatellite, log10M1, log10Mstar0, beta, delta, gamma, sigmaLogMstar, BCut, BSatellite, betaCut, betaSatellite"/>
    !!]

    self%Mstar0           =+10.0**log10Mstar0
    ! Initialize tables and accelerators.
    self%massPrevious     =-1.0d00
    self%massHaloPrevious =-1.0d00
    self%fMassTableMinimum=+1.0d08
    self%fMassTableMaximum=+1.0d13
    return
  end function behroozi2010ConstructorInternal

  subroutine behroozi2010Destructor(self)
    !!{
    Destructor for the \cite{behroozi_comprehensive_2010} conditional mass function class.
    !!}
    implicit none
    type(conditionalMassFunctionBehroozi2010), intent(inout) :: self

    call                                     self%fMassTable    %destroy()
    if (allocated(self%fMassHaloTable)) call self%fMassHaloTable%destroy()
    return
  end subroutine behroozi2010Destructor

  double precision function behroozi2010MassFunction(self,massHalo,mass,galaxyType)
    !!{
    Compute the cumulative conditional mass function, $\langle N(M_\star|M_\mathrm{halo}) \rangle \equiv
    \phi(M_\star|M_\mathrm{halo})$ using the fitting formula of \cite{behroozi_comprehensive_2010}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (conditionalMassFunctionBehroozi2010), intent(inout)           :: self
    double precision                                     , intent(in   )           :: massHalo        , mass
    type            (enumerationHaloModelGalaxyTypeType ), intent(in   ), optional :: galaxyType
    double precision                                                               :: numberCentrals  , numberSatellites
    type            (enumerationHaloModelGalaxyTypeType )                          :: galaxyTypeActual

    ! Get the number of satellites and centrals.
    call self%compute(massHalo,mass,numberCentrals,numberSatellites)
    ! Determine galaxy type required.
    if (present(galaxyType)) then
       galaxyTypeActual=galaxyType
    else
       galaxyTypeActual=haloModelGalaxyTypeAll
    end if
    ! Return required number.
    select case (galaxyTypeActual%ID)
    case (haloModelGalaxyTypeAll      %ID)
       ! Sum central and satellite contributions.
       behroozi2010MassFunction=numberCentrals+numberSatellites
    case (haloModelGalaxyTypeCentral  %ID)
       behroozi2010MassFunction=numberCentrals
    case (haloModelGalaxyTypeSatellite%ID)
       behroozi2010MassFunction=               numberSatellites
    case default
       behroozi2010MassFunction=0.0d0
       call Error_Report('unknown galaxy type'//{introspection:location})
    end select
    return
  end function behroozi2010MassFunction

  double precision function behroozi2010MassFunctionVariance(self,massHalo,massLow,massHigh)
    !!{
    Computes the variance in the cumulative conditional mass function, $\langle N(M_\star|M_\mathrm{halo}) \rangle \equiv
    \phi(M_\star|M_\mathrm{halo})$ using the fitting formula of \cite{behroozi_comprehensive_2010}. Assumes that the number of
    satellite galaxies is Poisson distributed, while the number of central galaxies follows a Bernoulli distribution, and that
    the numbers of satellites and centrals are uncorrelated.
    !!}
    implicit none
    class           (conditionalMassFunctionBehroozi2010), intent(inout) :: self
    double precision                                     , intent(in   ) :: massHalo        , massHigh            , massLow
    double precision                                                     :: numberCentrals  , numberCentralsHigh  , numberCentralsLow  , &
         &                                                                  numberSatellites, numberSatellitesHigh, numberSatellitesLow

    ! Get the number of satellites and centrals.
    call self%compute(massHalo,massLow ,numberCentralsLow ,numberSatellitesLow )
    call self%compute(massHalo,massHigh,numberCentralsHigh,numberSatellitesHigh)
    numberSatellites=max(numberSatellitesLow-numberSatellitesHigh,0.0d0)
    numberCentrals  =max(numberCentralsLow  -numberCentralsHigh  ,0.0d0)
    ! Compute the variance (using the Bienaymé formula for uncorrelated variables).
    behroozi2010MassFunctionVariance=+numberSatellites                      & ! Satellites are Poisson distributed, so the variance is just their number.
         &                           +numberCentrals*(1.0d0-numberCentrals)   ! Centrals are Bernoulli distributed.
    return
  end function behroozi2010MassFunctionVariance

  subroutine behroozi2010Compute(self,massHalo,mass,numberCentrals,numberSatellites)
    !!{
    Computes the cumulative conditional mass function, $\langle N(M_\star|M_\mathrm{halo}) \rangle \equiv \phi(M_\star|M_\mathrm{
    halo})$ using the fitting formula of \cite{behroozi_comprehensive_2010}.
    !!}
    use :: Table_Labels, only : extrapolationTypeFix
    implicit none
    class           (conditionalMassFunctionBehroozi2010), intent(inout), target :: self
    double precision                                     , intent(in   )         :: massHalo                , mass
    double precision                                     , intent(  out)         :: numberCentrals          , numberSatellites
    double precision                                     , parameter             :: massNormalization=1.0d12
    double precision                                     , parameter             :: massMinimum      =1.0d-12
    double precision                                     , parameter             :: massHaloMaximum  =1.0d+17
    double precision                                                             :: fMassHalo               , massCut         , &
         &                                                                          massSatellite

    self_ => self
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
       self%fMassTableCount=int(log10(self%fMassTableMaximum/self%fMassTableMinimum)*massTablePointsPerDecade)+1
       call          self%fMassTable    %destroy()
       if (allocated(self%fMassHaloTable)) then
          call       self%fMassHaloTable%destroy()
          deallocate(self%fMassHaloTable          )
       end if
       call self%fMassTable%create  (                                               &
            &                        self%fMassTableMinimum                       , &
            &                        self%fMassTableMaximum                       , &
            &                        self%fMassTableCount                         , &
            &                   extrapolationType=spread(extrapolationTypeFix,1,2)  &
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
       fMassHalo                  =min(behroozi2010fSHMRInverse(mass),massHaloMaximum)
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
    !!{
    The median mass vs. halo mass relation functional form from \cite{behroozi_comprehensive_2010}.
    !!}
    implicit none
    double precision, intent(in   ) :: mass
    double precision, parameter     :: logHaloMassTransition=20.0d0
    double precision                :: argument

    ! Compute the logarithmic halo mass for the given mass.
    argument=                                                      &
         &    self_%log10M1                                        &
         &   +self_%beta   *log10(mass/self_%Mstar0)               &
         &   +                   (mass/self_%Mstar0)**self_%delta  &
         &   /(                                                    &
         &     +1.0d0                                              &
         &     +1.0d0                                              &
         &     /                 (mass/self_%Mstar0)**self_%gamma  &
         &    )                                                    &
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
