!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{RST
Implements a class for the conditional mass functions using the :cite:t:`behroozi_comprehensive_2010` fitting function.
!!}

  use :: Numerical_Ranges, only : rangeLattice
  use :: Tables          , only : table1D     , table1DLogarithmicLinear

  !![
  <conditionalMassFunction name="conditionalMassFunctionBehroozi2010" docformat="rst">
   <description>
   A conditional mass function class which implements the fitting functions of :cite:t:`behroozi_comprehensive_2010`:

   .. math::

      \langle N_\mathrm{c}(M_\star|M)\rangle \equiv \int_{M_\star}^\infty \phi_\mathrm{c}(M_\star^\prime) \d \ln M_\star^\prime
      = {1 \over 2} \left[ 1 - \hbox{erf}\left( {\log_{10}M_\star - \log_{10} f_\mathrm{SHMR}(M) \over \sqrt{2}\sigma_{\log
      M_\star}} \right) \right].

   Here, the function :math:`f_\mathrm{SHMR}(M)` is the solution of

   .. math::

      \log_{10}M = \log_{10}M_1 + \beta \log_{10}\left({M_\star \over M_{\star,0}}\right) + {(M_\star/M_{\star,0})^\delta \over
      1 + (M_\star/M_{\star,0})^{-\gamma}} - {1/2}.

   For satellites,

   .. math::

      \langle N_\mathrm{s}(M_\star|M)\rangle \equiv \int_{M_\star}^\infty \phi_\mathrm{s}(M_\star^\prime) \d \ln M_\star^\prime
      = \langle N_\mathrm{c}(M_\star|M)\rangle \left({f^{-1}_\mathrm{SHMR}(M_\star) \over
      M_\mathrm{sat}}\right)^{\alpha_\mathrm{sat}} \exp\left(- {M_\mathrm{cut} \over f^{-1}_\mathrm{SHMR}(M_\star)} \right),

   where

   .. math::

      {M_\mathrm{sat} \over 10^{12} \mathrm{M}_\odot} = B_\mathrm{sat} \left({f^{-1}_\mathrm{SHMR}(M_\star) \over 10^{12}
      \mathrm{M}_\odot}\right)^{\beta_\mathrm{sat}},

   and

   .. math::

      {M_\mathrm{cut} \over 10^{12} \mathrm{M}_\odot} = B_\mathrm{cut} \left({f^{-1}_\mathrm{SHMR}(M_\star) \over 10^{12}
      \mathrm{M}_\odot}\right)^{\beta_\mathrm{cut}}.

   By default, parameter values are taken from the fit of :cite:t:`leauthaud_new_2011`, specifically their ``SIG_MOD1`` method for their :math:`z_1` sample. These default values, and the Galacticus input parameters which can be used to adjust them are shown in Table :numref:`{number} &lt;table-Behroozi2010FitParameters&gt;`. This method assumes that :math:`P_\mathrm{s}(N|M_\star,M;\delta \ln M_\star)` is a Poisson distribution while :math:`P_\mathrm{c}(N|M_\star,M;\delta \ln M_\star)` has a Bernoulli distribution, with each distribution's free parameter fixed by requiring

   .. math::

      \phi(M_\star;M) \delta \ln M_\star = \sum_{N=0}^\infty N P(N|M_\star,M;\delta \ln M_\star)

   .. list-table:: Parameters of the :cite:t:`behroozi_comprehensive_2010` conditional stellar mass function model, along with their default values and the corresponding Galacticus input parameters.
      :name: table-Behroozi2010FitParameters
      :header-rows: 1

      * - Parameter
        - Default
        - Galacticus name
      * - :math:`\alpha_\mathrm{sat}`
        - 1.0
        - ``[conditionalStellarMassFunctionBehrooziAlphaSatellite]``
      * - :math:`\log_{10} M_1`
        - 12.520
        - ``[conditionalStellarMassFunctionBehrooziLog10M1]``
      * - :math:`\log_{10} M_{\star,0}`
        - 10.916
        - ``[conditionalStellarMassFunctionBehrooziLog10Mstar0]``
      * - :math:`\beta`
        - 0.457
        - ``[conditionalStellarMassFunctionBehrooziBeta]``
      * - :math:`\delta`
        - 0.5666
        - ``[conditionalStellarMassFunctionBehrooziDelta]``
      * - :math:`\gamma`
        - 1.53
        - ``[conditionalStellarMassFunctionBehrooziGamma]``
      * - :math:`\sigma_{\log M_\star}`
        - 0.206
        - ``[conditionalStellarMassFunctionBehrooziSigmaLogMstar]``
      * - :math:`B_\mathrm{cut}`
        - 1.47
        - ``[conditionalStellarMassFunctionBehrooziBCut]``
      * - :math:`B_\mathrm{sat}`
        - 10.62
        - ``[conditionalStellarMassFunctionBehrooziBSatellite]``
      * - :math:`\beta_\mathrm{cut}`
        - :math:`-`\ 0.13
        - ``[conditionalStellarMassFunctionBehrooziBetaCut]``
      * - :math:`\beta_\mathrm{sat}`
        - 0.859
        - ``[conditionalStellarMassFunctionBehrooziBetaSatellite]``
   </description>
  </conditionalMassFunction>
  !!]
  type, extends(conditionalMassFunctionClass) :: conditionalMassFunctionBehroozi2010
     !!{RST
     Implements the conditional mass function using the fitting functions of :cite:t:`behroozi_comprehensive_2010`.
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
     ! Tabulation of-halo mass relation. The lattice to which the tabulated masses are pinned is the source of truth for the
     ! extent of the tabulation; the mass bounds are derived from it, and are retained because they are what the test for a
     ! sufficient tabulation reads, and what the extension of the range is measured from.
     type            (rangeLattice            )               :: latticeMass
     double precision                                         :: fMassTableMaximum    , fMassTableMinimum
     double precision                                         :: fMassHaloTableMaximum, fMassHaloTableMinimum
     type            (table1DLogarithmicLinear)               :: fMassTable
     class           (table1D                 ), allocatable  :: fMassHaloTable
   contains
     !![
     <methods docformat="rst">
       <method description="Compute the cumulative conditional mass function, :math:`\langle N(M_\star|M_\mathrm{halo}) \rangle \equiv \phi(M_\star|M_\mathrm{halo})`." method="compute" />
     </methods>
     !!]
     final     ::                         behroozi2010Destructor
     procedure :: massFunction         => behroozi2010MassFunction
     procedure :: massFunctionVariance => behroozi2010MassFunctionVariance
     procedure :: compute              => behroozi2010Compute
  end type conditionalMassFunctionBehroozi2010

  interface conditionalMassFunctionBehroozi2010
     !!{RST
     Constructors for the :cite:t:`behroozi_comprehensive_2010` merging timescale class.
     !!}
     module procedure behroozi2010ConstructorParameters
     module procedure behroozi2010ConstructorInternal
  end interface conditionalMassFunctionBehroozi2010

  ! Table resolution.
  integer                                              , parameter :: massTablePointsPerDecade=10

  ! Seed range of masses to tabulate. Any two tabulations therefore contain this range, and so always overlap. Note that these
  ! are exact powers of ten, so on a lattice anchored to whole decades they are themselves anchor points: a tabulation which is
  ! never asked for a halo mass outside those which they span covers exactly this range, as it did before being pinned.
  double precision                                     , parameter :: massTableSeedMinimum    =1.0d08, massTableSeedMaximum=1.0d13

  ! Module scope pointer to the current object.
  class           (conditionalMassFunctionBehroozi2010), pointer   :: self_
  !$omp threadprivate(self_)

contains


  function behroozi2010ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`conditionalMassFunctionBehroozi2010` conditional mass function class which builds the object from a parameter set.
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
    <inputParameter docformat="rst">
      <name>alphaSatellite</name>
      <defaultSource>
      (:cite:author:`leauthaud_new_2011` :cite:year:`leauthaud_new_2011`; :math:`z_1` sample using their ``SIG_MOD1`` method)
      </defaultSource>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The parameter :math:`\alpha_\mathrm{sat}` from the fitting functions of :cite:t:`behroozi_comprehensive_2010`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>log10M1</name>
      <defaultSource>
      (:cite:author:`leauthaud_new_2011` :cite:year:`leauthaud_new_2011`; :math:`z_1` sample using their ``SIG_MOD1`` method)
      </defaultSource>
      <defaultValue>12.520d0</defaultValue>
      <description>
      The parameter :math:`\log_{10}M_1` from the fitting functions of :cite:t:`behroozi_comprehensive_2010`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>log10Mstar0</name>
      <defaultSource>
      (:cite:author:`leauthaud_new_2011` :cite:year:`leauthaud_new_2011`; :math:`z_1` sample using their ``SIG_MOD1`` method)
      </defaultSource>
      <defaultValue>10.916d0</defaultValue>
      <description>
      The parameter :math:`\log_{10}M_{\star,0}` from the fitting functions of :cite:t:`behroozi_comprehensive_2010`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>beta</name>
      <defaultSource>
      (:cite:author:`leauthaud_new_2011` :cite:year:`leauthaud_new_2011`; :math:`z_1` sample using their ``SIG_MOD1`` method)
      </defaultSource>
      <defaultValue>0.457d0</defaultValue>
      <description>
      The parameter :math:`\beta` from the fitting functions of :cite:t:`behroozi_comprehensive_2010`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>delta</name>
      <defaultSource>
      (:cite:author:`leauthaud_new_2011` :cite:year:`leauthaud_new_2011`; :math:`z_1` sample using their ``SIG_MOD1`` method)
      </defaultSource>
      <defaultValue>0.5666d0</defaultValue>
      <description>
      The parameter :math:`\delta` from the fitting functions of :cite:t:`behroozi_comprehensive_2010`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>gamma</name>
      <defaultSource>
      (:cite:author:`leauthaud_new_2011` :cite:year:`leauthaud_new_2011`; :math:`z_1` sample using their ``SIG_MOD1`` method)
      </defaultSource>
      <defaultValue>1.53d0</defaultValue>
      <description>
      The parameter :math:`\gamma` from the fitting functions of :cite:t:`behroozi_comprehensive_2010`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>sigmaLogMstar</name>
      <defaultSource>
      (:cite:author:`leauthaud_new_2011` :cite:year:`leauthaud_new_2011`; :math:`z_1` sample using their ``SIG_MOD1`` method)
      </defaultSource>
      <defaultValue>0.206d0</defaultValue>
      <description>
      The parameter :math:`\sigma_{\log M_\star}` from the fitting functions of :cite:t:`behroozi_comprehensive_2010`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>BCut</name>
      <defaultSource>
      (:cite:author:`leauthaud_new_2011` :cite:year:`leauthaud_new_2011`; :math:`z_1` sample using their ``SIG_MOD1`` method)
      </defaultSource>
      <defaultValue>1.47d0</defaultValue>
      <description>
      The parameter :math:`B_\mathrm{cut}` from the fitting functions of :cite:t:`behroozi_comprehensive_2010`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>BSatellite</name>
      <defaultSource>
      (:cite:author:`leauthaud_new_2011` :cite:year:`leauthaud_new_2011`; :math:`z_1` sample using their ``SIG_MOD1`` method)
      </defaultSource>
      <defaultValue>10.62d0</defaultValue>
      <description>
      The parameter :math:`B_\mathrm{sat}` from the fitting functions of :cite:t:`behroozi_comprehensive_2010`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>betaCut</name>
      <defaultSource>
      (:cite:author:`leauthaud_new_2011` :cite:year:`leauthaud_new_2011`; :math:`z_1` sample using their ``SIG_MOD1`` method)
      </defaultSource>
      <defaultValue>-0.13d0</defaultValue>
      <description>
      The parameter :math:`\beta_\mathrm{cut}` from the fitting functions of :cite:t:`behroozi_comprehensive_2010`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>betaSatellite</name>
      <defaultSource>
      (:cite:author:`leauthaud_new_2011` :cite:year:`leauthaud_new_2011`; :math:`z_1` sample using their ``SIG_MOD1`` method)
      </defaultSource>
      <defaultValue>0.859d0</defaultValue>
      <description>
      The parameter :math:`\beta_\mathrm{sat}` from the fitting functions of :cite:t:`behroozi_comprehensive_2010`.
      </description>
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
    !!{RST
    Internal constructor for the :cite:t:`behroozi_comprehensive_2010` conditional mass function class.
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
    self%fMassTableMinimum=massTableSeedMinimum
    self%fMassTableMaximum=massTableSeedMaximum
    return
  end function behroozi2010ConstructorInternal

  subroutine behroozi2010Destructor(self)
    !!{RST
    Destructor for the :cite:t:`behroozi_comprehensive_2010` conditional mass function class.
    !!}
    implicit none
    type(conditionalMassFunctionBehroozi2010), intent(inout) :: self

    call                                     self%fMassTable    %destroy()
    if (allocated(self%fMassHaloTable)) call self%fMassHaloTable%destroy()
    return
  end subroutine behroozi2010Destructor

  double precision function behroozi2010MassFunction(self,massHalo,mass,galaxyType)
    !!{RST
    Compute the cumulative conditional mass function, :math:`\langle N(M_\star|M_\mathrm{halo}) \rangle \equiv \phi(M_\star|M_\mathrm{halo})` using the fitting formula of :cite:t:`behroozi_comprehensive_2010`.
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
    !!{RST
    Computes the variance in the cumulative conditional mass function, :math:`\langle N(M_\star|M_\mathrm{halo}) \rangle \equiv \phi(M_\star|M_\mathrm{halo})` using the fitting formula of :cite:t:`behroozi_comprehensive_2010`. Assumes that the number of satellite galaxies is Poisson distributed, while the number of central galaxies follows a Bernoulli distribution, and that the numbers of satellites and centrals are uncorrelated.
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
    !!{RST
    Computes the cumulative conditional mass function, :math:`\langle N(M_\star|M_\mathrm{halo}) \rangle \equiv \phi(M_\star|M_\mathrm{ halo})` using the fitting formula of :cite:t:`behroozi_comprehensive_2010`.
    !!}
    use :: Numerical_Ranges, only : Range_Pinned        , gridSchemePerDecade
    use :: Table_Labels    , only : extrapolationTypeFix
    implicit none
    class           (conditionalMassFunctionBehroozi2010), intent(inout), target    :: self
    double precision                                     , intent(in   )            :: massHalo                 , mass
    double precision                                     , intent(  out)            :: numberCentrals           , numberSatellites
    double precision                                     , parameter                :: massNormalization=1.0d+12
    double precision                                     , parameter                :: massMinimum      =1.0d-12
    double precision                                     , parameter                :: massHaloMaximum  =1.0d+17
    type            (rangeLattice                       )                           :: latticeMass
    logical                                              , allocatable, dimension(:):: isComputed
    double precision                                                                :: fMassHalo                , massCut         , &
         &                                                                             massSatellite            , massTableMinimum, &
         &                                                                             massTableMaximum

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
       ! Find the range of masses to tabulate, pinning it to an absolute lattice so that the masses evaluated - and therefore
       ! every value interpolated between them - depend only on which lattice points are spanned, and not on the sequence of
       ! halo masses which happened to be requested. The range already tabulated is unioned in through `latticeCurrent` rather
       ! than being folded into the target, and no safety margin is applied: the loop condition already asks for exactly the
       ! coverage which is needed.
       !
       ! The range is extended by taking a target beyond the current bound, rather than by the factors of two which were applied
       ! to those bounds before this tabulation was pinned. A factor of two need not reach the next anchor point, and an
       ! iteration which failed to enlarge the pinned range would rebuild an identical tabulation and leave the loop condition
       ! unchanged - repeating forever. Halving a bound which is itself an anchor point always falls short of it, and so always
       ! reaches the next anchor below.
       massTableMinimum=self%fMassTableMinimum
       massTableMaximum=self%fMassTableMaximum
       if (allocated(self%fMassHaloTable)) then
          if (massHalo < self%fMassHaloTableMinimum) massTableMinimum=0.5d0*massTableMinimum
          if (massHalo > self%fMassHaloTableMaximum) massTableMaximum=2.0d0*massTableMaximum
       end if
       ! Note that the lower limit on the mass tabulated is imposed through the loop condition above, and not through
       ! `limitMinimum`. A clamped bound is held *at* the limit, so the loop condition - which asks whether there is any room
       ! left below - would have to compare a lattice point with that limit for equality. The two need not agree in their final
       ! bits, and if the lattice point fell an ulp above the limit the loop would never end, rebuilding an identical tabulation
       ! forever. Letting the range step past the limit by at most one anchor interval, as it did before being pinned, avoids
       ! resting the termination of a loop on the last bit of a power of ten.
       latticeMass=Range_Pinned(                                                    &
            &                                  [massTableMinimum,massTableMaximum], &
            &                                   massTablePointsPerDecade          , &
            &                                   gridSchemePerDecade               , &
            &                   marginFactor  = 1.0d0                             , &
            &                   latticeCurrent= self%latticeMass                    &
            &                  )
       ! Extend the tabulation onto the new lattice. The tabulated masses are lattice points, so the halo mass tabulated for a
       ! given mass is identical between one tabulation and another however many each spans.
       call self%fMassTable%extend  (                                               &
            &                        latticeMass                                  , &
            &                        isComputed                                   , &
            &                   extrapolationType=spread(extrapolationTypeFix,1,2)  &
            &                       )
       call self%fMassTable%populate(                                               &
            &                        behroozi2010fSHMRInverse(self%fMassTable%xs()) &
            &                       )
       ! Rebuild the reversed table, which is tabulated at halo masses and so lies on no lattice.
       if (allocated(self%fMassHaloTable)) then
          call       self%fMassHaloTable%destroy()
          deallocate(self%fMassHaloTable          )
       end if
       call self%fMassTable%reverse (self%fMassHaloTable)
       self%latticeMass          =latticeMass
       self%fMassTableMinimum    =latticeMass%minimum()
       self%fMassTableMaximum    =latticeMass%maximum()
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
    !!{RST
    The median mass vs. halo mass relation functional form from :cite:t:`behroozi_comprehensive_2010`.
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
