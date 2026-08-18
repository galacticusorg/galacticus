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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which implements useful shared utilities for calculations of decaying dark matter.
!!}

module Decaying_Dark_Matter
  !!{RST
  Implements useful shared utilities for calculations of decaying dark matter.
  !!}
  use :: Numerical_Interpolation, only : interpolator
  use :: Numerical_Ranges       , only : rangeLattice
  private
  public :: decayingDarkMatterFractionRetained, decayingDarkMatterEnergyRetained, decayingDarkMatterFractionRetainedDerivatives, decayingDarkMatterEnergyRetainedDerivatives

  ! Tables of retained fractions and energies, together with the absolute lattices on which they are tabulated.
  type            (rangeLattice)                              :: latticeVelocityEscape                      , latticeVelocityKick
  double precision              , dimension(:  ), allocatable :: velocitiesEscapeScaleFree                  , velocitiesKickScaleFree
  double precision              , dimension(:,:), allocatable :: energyRetained                             , fractionRetained
  type            (interpolator)                , allocatable :: interpolatorVelocityEscape                 , interpolatorVelocityKick
  !$omp threadprivate(latticeVelocityEscape,latticeVelocityKick,velocitiesEscapeScaleFree,velocitiesKickScaleFree,energyRetained,fractionRetained,interpolatorVelocityEscape,interpolatorVelocityKick)

  ! Number of points per decade at which each axis of the solutions is tabulated.
  integer                                       , parameter   :: countVelocityEscapePerDecade  =    300     , countVelocityKickPerDecade    =    300

  ! Scale free velocity above which the velocity dispersion is negligible and we assume deterministic behavior.
  double precision                              , parameter   :: velocityScaleFreeLarge        =+100.0d0

  ! Minimum scale free escape velocity for which we can find a self-consistent velocity distribution function. Smaller values only
  ! occur at extreme distances outside of halos, so should not matter.
  double precision                              , parameter   :: velocityEscapeScaleFreeLimit  =+  2.3d0

  ! Seed range for the kick velocity axis. Unlike the escape velocity axis - which is bounded below by the limit above, and for
  ! which the deterministic limit provides a natural upper end - the kick velocity has no physical scale of its own, so a seed
  ! must be stated. Unity is the natural reference of the scale free variable: below it the kick is smaller than the velocity
  ! dispersion and essentially all particles are retained. Seeding matters because it guarantees that any two tabulations - built
  ! in different runs, at different requested velocities - overlap, and so can always be merged onto the common lattice.
  double precision                              , parameter   :: velocityKickScaleFreeSeed     =+  1.0d0

  ! Generate a source digest. The tabulation below is a property of the truncated Maxwell-Boltzmann distribution alone, so no
  ! object's parameters enter it and none belong in its file name. It does depend on the code which generates it - the
  ! integrands, the integration and root-finding tolerances, and the constants above - so the digest is what its cache file must
  ! be keyed on. Without it a tabulation computed by an earlier version would be silently reused after the code was changed.
  !![
  <sourceDigest name="decayingDarkMatterSourceDigest"/>
  !!]

contains

  double precision function decayingDarkMatterFractionRetained(velocityDispersion,velocityEscape,velocityKick) result(fraction)
    !!{RST
    Compute the fraction of decaying dark matter particles retained in a halo.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    double precision          , intent(in   )  :: velocityDispersion     , velocityEscape       , &
         &                                        velocityKick
    double precision                           :: velocityEscapeScaleFree, velocityKickScaleFree
    double precision          , dimension(0:1) :: hVelocityEscape        , hVelocityKick
    integer         (c_size_t), dimension(0:1) :: iVelocityEscape        , iVelocityKick
    integer                                    :: jVelocityEscape        , jVelocityKick

    if (velocityKick >= 2.0d0*velocityEscape) then
       ! For kicks above twice the escape velocity, no particles are retained.
       fraction=0.0d0
    else
       ! Compute the scale free velocities, using the velocity dispersion as our reference velocity.
       velocityEscapeScaleFree=max(velocityEscape/velocityDispersion,velocityEscapeScaleFreeLimit)
       velocityKickScaleFree  =    velocityKick  /velocityDispersion
       ! Check for deterministic limit.
       if     (                                                  &
            &   velocityEscapeScaleFree > velocityScaleFreeLarge &
            &  .and.                                             &
            &   velocityKickScaleFree   > velocityScaleFreeLarge &
            & ) then
          if  (velocityKickScaleFree > velocityEscapeScaleFree) then
             fraction=0.0d0
          else
             fraction=1.0d0
          end if
          return
       end if
       ! Ensure the tabulated solutions cover a sufficient range.
       call decayingDarkMatterRetainedTabulate(velocityEscapeScaleFree,velocityKickScaleFree)
       ! Interpolate in the tabulated solutions.
       call interpolatorVelocityEscape%linearFactors(velocityEscapeScaleFree,iVelocityEscape(0),hVelocityEscape)
       call interpolatorVelocityKick  %linearFactors(velocityKickScaleFree  ,iVelocityKick  (0),hVelocityKick  )
       iVelocityEscape(1)=iVelocityEscape(0)+1
       iVelocityKick  (1)=iVelocityKick  (0)+1
       fraction=0.0d0
       do jVelocityEscape =0,1
          do jVelocityKick=0,1
             fraction=+fraction                                                                        &
                  &   +fractionRetained(iVelocityEscape(jVelocityEscape),iVelocityKick(jVelocityKick)) &
                  &   *                 hVelocityEscape(jVelocityEscape)                               &
                  &   *                                                  hVelocityKick(jVelocityKick)
          end do
       end do
    end if
    return
  end function decayingDarkMatterFractionRetained
  
  double precision function decayingDarkMatterEnergyRetained(velocityDispersion,velocityEscape,velocityKick) result(energy)
    !!{RST
    Compute the fraction of decaying dark matter particle energy retained in a halo.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    double precision          , intent(in   )  :: velocityDispersion     , velocityEscape       , &
         &                                        velocityKick
    double precision                           :: velocityEscapeScaleFree, velocityKickScaleFree
    double precision          , dimension(0:1) :: hVelocityEscape        , hVelocityKick
    integer         (c_size_t), dimension(0:1) :: iVelocityEscape        , iVelocityKick
    integer                                    :: jVelocityEscape        , jVelocityKick

    if (velocityKick >= 2.0d0*velocityEscape) then
       ! For kicks above twice the escape velocity, no particles (and, therefore, no energy) are retained.
       energy=0.0d0
    else
       ! Compute the scale free velocities, using the velocity dispersion as our reference velocity.
       velocityEscapeScaleFree=max(velocityEscape/velocityDispersion,velocityEscapeScaleFreeLimit)
       velocityKickScaleFree  =    velocityKick  /velocityDispersion
       ! Check for deterministic limit.
       if     (                                                  &
            &   velocityEscapeScaleFree > velocityScaleFreeLarge &
            &  .and.                                             &
            &   velocityKickScaleFree   > velocityScaleFreeLarge &
            & ) then
          if  (velocityKickScaleFree > velocityEscapeScaleFree) then
             energy=+0.0d0
          else
             energy=+0.5d0*velocityKick**2
          end if
          return
       end if
       ! Ensure the tabulated solutions cover a sufficient range.
       call decayingDarkMatterRetainedTabulate(velocityEscapeScaleFree,velocityKickScaleFree)
       ! Interpolate in the tabulated solutions.
       call interpolatorVelocityEscape%linearFactors(velocityEscapeScaleFree,iVelocityEscape(0),hVelocityEscape)
       call interpolatorVelocityKick  %linearFactors(velocityKickScaleFree  ,iVelocityKick  (0),hVelocityKick  )
       iVelocityEscape(1)=iVelocityEscape(0)+1
       iVelocityKick  (1)=iVelocityKick  (0)+1
       energy=0.0d0
       do jVelocityEscape=0,1
          do jVelocityKick=0,1
             energy=+energy                                                                        &
                  & +energyRetained(iVelocityEscape(jVelocityEscape),iVelocityKick(jVelocityKick)) &
                  & *               hVelocityEscape(jVelocityEscape)                               &
                  & *                                                hVelocityKick(jVelocityKick)
          end do
       end do
    end if
    ! Scale back into physical units.
    energy=+energy                &
         & *velocityDispersion**2
    return
  end function decayingDarkMatterEnergyRetained

  subroutine decayingDarkMatterFractionRetainedDerivatives(velocityDispersion,velocityEscape,velocityKick,fractionDerivativeVelocityEscapeScaleFree,fractionDerivativeVelocityKickScaleFree)
    !!{RST
    Compute the partial derivatives of the fraction of decaying dark matter particles retained in a halo.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    double precision          , intent(in   )  :: velocityDispersion                       , velocityEscape                         , &
         &                                        velocityKick
    double precision          , intent(  out)  :: fractionDerivativeVelocityEscapeScaleFree, fractionDerivativeVelocityKickScaleFree
    double precision                           :: velocityEscapeScaleFree                  , velocityKickScaleFree
    double precision          , dimension(0:1) :: hVelocityEscape                          , hVelocityKick
    integer         (c_size_t), dimension(0:1) :: iVelocityEscape                          , iVelocityKick
    integer                                    :: jVelocityEscape                          , jVelocityKick

    if (velocityKick >= 2.0d0*velocityEscape) then
       ! For kicks above twice the escape velocity, no particles are retained.
       fractionDerivativeVelocityEscapeScaleFree=0.0d0
       fractionDerivativeVelocityKickScaleFree  =0.0d0
    else
       ! Compute the scale free velocities, using the velocity dispersion as our reference velocity.
       velocityEscapeScaleFree=max(velocityEscape/velocityDispersion,velocityEscapeScaleFreeLimit)
       velocityKickScaleFree  =    velocityKick  /velocityDispersion
       ! Check for deterministic limit.
       if     (                                                  &
            &   velocityEscapeScaleFree > velocityScaleFreeLarge &
            &  .and.                                             &
            &   velocityKickScaleFree   > velocityScaleFreeLarge &
            & ) then
          fractionDerivativeVelocityEscapeScaleFree=0.0d0
          fractionDerivativeVelocityKickScaleFree  =0.0d0
          return
       end if
       ! Ensure the tabulated solutions cover a sufficient range.
       call decayingDarkMatterRetainedTabulate(velocityEscapeScaleFree,velocityKickScaleFree)
       ! Interpolate in the tabulated solutions.
       call interpolatorVelocityEscape%linearFactors(velocityEscapeScaleFree,iVelocityEscape(0),hVelocityEscape)
       call interpolatorVelocityKick  %linearFactors(velocityKickScaleFree  ,iVelocityKick  (0),hVelocityKick  )
       iVelocityEscape(1)=iVelocityEscape(0)+1
       iVelocityKick  (1)=iVelocityKick  (0)+1
       fractionDerivativeVelocityEscapeScaleFree=0.0d0
       do jVelocityKick=0,1
          fractionDerivativeVelocityEscapeScaleFree=+fractionDerivativeVelocityEscapeScaleFree                                                                                                                                         &
               &                                    +(fractionRetained         (iVelocityEscape(1             ),iVelocityKick(jVelocityKick))-fractionRetained         (iVelocityEscape(             0),iVelocityKick(jVelocityKick))) &
               &                                    /(velocitiesEscapeScaleFree(iVelocityEscape(1             )                             )-velocitiesEscapeScaleFree(iVelocityEscape(             0)                             )) &
               &                                    *                                                           hVelocityKick(jVelocityKick)
       end do
       fractionDerivativeVelocityKickScaleFree  =0.0d0
       do jVelocityEscape=0,1
          fractionDerivativeVelocityKickScaleFree  =+fractionDerivativeVelocityKickScaleFree                                                                                                                                           &
               &                                    +(fractionRetained         (iVelocityEscape(jVelocityEscape),iVelocityKick(           1))-fractionRetained         (iVelocityEscape(jVelocityEscape),iVelocityKick(           0))) &
               &                                    /(velocitiesKickScaleFree  (                                 iVelocityKick(           1))-velocitiesKickScaleFree  (                                 iVelocityKick(           0))) &
               &                                    *                           hVelocityEscape(jVelocityEscape)
       end do
    end if
    return
  end subroutine decayingDarkMatterFractionRetainedDerivatives

  subroutine decayingDarkMatterEnergyRetainedDerivatives(velocityDispersion,velocityEscape,velocityKick,energyDerivativeVelocityEscapeScaleFree,energyDerivativeVelocityKickScaleFree)
    !!{RST
    Compute the partial derivatives of the energy of decaying dark matter particles retained in a halo.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    double precision          , intent(in   )  :: velocityDispersion                     , velocityEscape                       , &
         &                                        velocityKick
    double precision          , intent(  out)  :: energyDerivativeVelocityEscapeScaleFree, energyDerivativeVelocityKickScaleFree
    double precision                           :: velocityEscapeScaleFree                , velocityKickScaleFree
    double precision          , dimension(0:1) :: hVelocityEscape                        , hVelocityKick
    integer         (c_size_t), dimension(0:1) :: iVelocityEscape                        , iVelocityKick
    integer                                    :: jVelocityEscape                        , jVelocityKick

    if (velocityKick >= 2.0d0*velocityEscape) then
       ! For kicks above twice the escape velocity, no particles are retained.
       energyDerivativeVelocityEscapeScaleFree=0.0d0
       energyDerivativeVelocityKickScaleFree  =0.0d0
    else
       ! Compute the scale free velocities, using the velocity dispersion as our reference velocity.
       velocityEscapeScaleFree=max(velocityEscape/velocityDispersion,velocityEscapeScaleFreeLimit)
       velocityKickScaleFree  =    velocityKick  /velocityDispersion
       ! Check for deterministic limit.
       if     (                                                  &
            &   velocityEscapeScaleFree > velocityScaleFreeLarge &
            &  .and.                                             &
            &   velocityKickScaleFree   > velocityScaleFreeLarge &
            & ) then
          energyDerivativeVelocityEscapeScaleFree=0.0d0
          energyDerivativeVelocityKickScaleFree  =0.0d0
          return
       end if
       ! Ensure the tabulated solutions cover a sufficient range.
       call decayingDarkMatterRetainedTabulate(velocityEscapeScaleFree,velocityKickScaleFree)
       ! Interpolate in the tabulated solutions.
       call interpolatorVelocityEscape%linearFactors(velocityEscapeScaleFree,iVelocityEscape(0),hVelocityEscape)
       call interpolatorVelocityKick  %linearFactors(velocityKickScaleFree  ,iVelocityKick  (0),hVelocityKick  )
       iVelocityEscape(1)=iVelocityEscape(0)+1
       iVelocityKick  (1)=iVelocityKick  (0)+1
       energyDerivativeVelocityEscapeScaleFree=0.0d0
       do jVelocityKick =0,1
          energyDerivativeVelocityEscapeScaleFree=+energyDerivativeVelocityEscapeScaleFree                                                                                                                                             &
               &                                  +(energyRetained           (iVelocityEscape(1              ),iVelocityKick(jVelocityKick))-energyRetained           (iVelocityEscape(              0),iVelocityKick(jVelocityKick))) &
               &                                  /(velocitiesEscapeScaleFree(iVelocityEscape(1              )                             )-velocitiesEscapeScaleFree(iVelocityEscape(              0)                             )) &
               &                                  *                                                            hVelocityKick(jVelocityKick)
       end do
       energyDerivativeVelocityEscapeScaleFree=energyDerivativeVelocityEscapeScaleFree*velocityDispersion**2
       energyDerivativeVelocityKickScaleFree  =0.0d0
       do jVelocityEscape =0,1        
          ! Note that, since we expect ε/σ² ~ xₖ² when all energy is retained, we evaluate this derivative as:
          !
          !   d(ε/σ²)/dxₖ =  d(ε/σ²)/d(xₖ²) d(xₖ²)/dx = d(ε/σ²)/d(xₖ²) 2 xₖ ≅ Δ(ε/σ²)/Δ(xₖ²) 2 xₖ
          !
          ! thereby ensuring that the finite difference derivative is precise in the limit of all energy being retained.
          energyDerivativeVelocityKickScaleFree  =+energyDerivativeVelocityKickScaleFree                                                                                                                                                     &
               &                                  +(energyRetained           (iVelocityEscape(jVelocityEscape),iVelocityKick(            1))   -energyRetained           (iVelocityEscape(jVelocityEscape),iVelocityKick(            0))   ) &
               &                                  /(velocitiesKickScaleFree  (                                 iVelocityKick(            1))**2-velocitiesKickScaleFree  (                                 iVelocityKick(            0))**2) &
               &                                  *                           hVelocityEscape(jVelocityEscape)                                                                                                                               &
               &                                  *2.0d0                                                                                                                                                                                     &
               &                                  *velocityKickScaleFree
       end do
       energyDerivativeVelocityKickScaleFree  =energyDerivativeVelocityKickScaleFree  *velocityDispersion**2
    end if
    return
  end subroutine decayingDarkMatterEnergyRetainedDerivatives

  subroutine decayingDarkMatterRetainedTabulate(velocityEscapeScaleFree,velocityKickScaleFree)
    !!{RST
    Compute the fraction of the decaying dark matter particles and kick energy that is retained. Assumes that the initial distribution of particle velocities is a Maxwell-Boltzmann distribution, truncated at the escape velocity. The mean energy of retained particles, minus their mean energy pre-kick is computed.

    The truncated Maxwell-Boltzmann distribution is

    .. math::

       p(v,\theta|s) = \left\{ \begin{array}{ll} A^{-1} v^2 \exp\left(-\frac{1}{2}\left[\frac{v}{s}\right]^2\right) & \hbox{if } v < v_\mathrm{e} \\ 0 &  \hbox{if } v \ge v_\mathrm{e}, \end{array} \right.

    where :math:`s` is the velocity width, :math:`v_\mathrm{e}` is the escape velocity, :math:`v` is particle speed, :math:`\theta` is the direction of the particle velocity relative to the :math:`z`-axis, and

    .. math::

       A = \sqrt{2 \pi} s^3 \hbox{erf}\left( \frac{v_\mathrm{e}}{\sqrt{2}s} \right) - 2 v_\mathrm{e} s^2 \exp \left( -\frac{1}{2}\left[\frac{v_\mathrm{e}}{s}\right]^2\right)

    is a normalization factor.

    Assume, without loss of generality, that the kick is along the :math:`z`-axis. The specific kinetic energy of the retained particles, in excess of their original energy is:

    .. math::

       \epsilon =  \int_{-1}^{+1} \mathrm{d}\cos\theta \int_0^{v_\mathrm{e}} \mathrm{d}v \frac{1}{2} \left( v^2 + v_\mathrm{k}^2 + 2 v v_\mathrm{k} \cos\theta - v^2 \right) p(v,\theta|s) H\left( v^2 + v_\mathrm{k}^2 + 2 v v_\mathrm{k} \cos\theta < v_\mathrm{e}^2 \right)

    where :math:`v_\mathrm{k}` is the scale-free kick velocity, and :math:`H(x) = 1` if :math:`x` is true, and 0 otherwise. Solving the inequality for the velocity, :math:`v`, that will remain bound as a function of :math:`\theta`, gives :math:`v_\mathrm{max|min}(\theta) = \pm \left( v_\mathrm{e}^2 - v_\mathrm{k}^2 \sin^2\theta \right)^{1/2} - v_\mathrm{k} \cos \theta`, so:

    .. math::
       :label: eq-decayingDMRetainedEnergy

       \epsilon =  \int_{-1}^{+_1} \mathrm{d}\cos\theta \int_{v_\mathrm{min}(\theta)}^{v_\mathrm{max}(\theta)} \mathrm{d}v \frac{1}{2} \left( v_\mathrm{k}^2 + 2 v v_\mathrm{k} \cos\theta \right) p(v,\theta|s).

    Similarly, the retained fraction of particles is simply

    .. math::
       :label: eq-decayingDMRetainedFraction

       f =  \int_{-1}^{+_1} \mathrm{d}\cos\theta \int_{v_\mathrm{min}(\theta)}^{v_\mathrm{max}(\theta)} \mathrm{d}v p(v,\theta|s).
    !!}
    use :: Display                 , only : displayCounter, displayCounterClear          , displayIndent                , displayUnindent          , &
         &                                  displayMessage, verbosityLevelStandard
    use :: ISO_Varying_String      , only : varying_string, operator(//)                 , assignment(=)
    use :: Numerical_Constants_Math, only : Pi
    use :: Root_Finder             , only : rootFinder    , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rangeExpandMultiplicative
    use :: Numerical_Integration   , only : integrator
    use :: Numerical_Ranges        , only : Range_Pinned  , rangeLattice                 , gridSchemePerDecade          , Range_Lattice_Offset
    use :: File_Utilities          , only : File_Exists   , Directory_Make               , File_Path, File_Lock         , File_Unlock, lockDescriptor
    use :: HDF5_Access             , only : hdf5Access
    use :: IO_HDF5                 , only : hdf5File
    use :: Input_Paths             , only : inputPath     , pathTypeDataDynamic
    use :: String_Handling         , only : String_C_To_Fortran
    implicit none
    double precision              , intent(inout)               :: velocityEscapeScaleFree
    double precision              , intent(in   )               :: velocityKickScaleFree
    double precision              , parameter                   :: toleranceAbsolute         =0.0d+0, toleranceRelative       =1.0d-6
    ! Note that the abscissae and solutions are accumulated in these local arrays rather than directly in the module-scope arrays:
    ! the latter are threadprivate, and so are not shared with the threads of the parallel region below.
    double precision              , dimension(:  ), allocatable :: velocitiesEscapeScaleFree_       , velocitiesKickScaleFree_
    double precision              , dimension(:,:), allocatable :: energyRetained_                  , fractionRetained_
    double precision              , save                        :: velocityWidthScaleFree           , normalization                  , &
         &                                                         velocityEscapeScaleFree_         , velocityKickScaleFree_         , &
         &                                                         cosThetaMaximum
    type            (rootFinder  ), save          , allocatable :: finder
    type            (integrator  ), save          , allocatable :: integratorEnergy                 , integratorFraction
    integer                       , save                        :: iVelocityKick
    integer                                                     :: countVelocityEscape              , countVelocityKick              , &
         &                                                         iVelocityEscape                  , offsetVelocityEscape           , &
         &                                                         offsetVelocityKick
    type            (rangeLattice)                              :: latticeEscape                    , latticeKick
    logical                       , dimension(:,:), allocatable :: isComputed
    !$omp threadprivate(iVelocityKick,velocityEscapeScaleFree_,velocityKickScaleFree_,velocityWidthScaleFree,normalization,cosThetaMaximum,finder,integratorFraction,integratorEnergy)

    ! Find the ranges of escape and kick velocity required, pinned to absolute lattices so that the tabulation - and hence every
    ! value interpolated from it - is independent of the velocities at which it happened to be first requested, and so that it
    ! can be extended without recomputing any solution already found.
    call rangesRequired()
    if     (                                             &
         &   latticeVelocityEscape%covers(latticeEscape) &
         &  .and.                                        &
         &   latticeVelocityKick  %covers(latticeKick  ) &
         & ) then
       call velocityEscapeScaleFreeClamp()
       return
    end if
    block
      type     (varying_string) :: message     , fileName
      character(len=8         ) :: labelMinimum, labelMaximum
      type     (lockDescriptor) :: fileLock

      fileName=inputPath(pathTypeDataDynamic)                     // &
           &   'darkMatter/decayingDarkMatterRetention_'          // &
           &   String_C_To_Fortran(decayingDarkMatterSourceDigest)// &
           &   '.hdf5'
      call Directory_Make(File_Path(fileName))
      ! Adopt any tabulation already cached on disk. Since the tabulation lies on absolute lattices the cached and in-memory
      ! solutions share abscissae wherever they overlap, so the cache is adopted whenever it spans everything we already hold -
      ! rather than, as previously, being read over the top of whatever is in memory, which could discard a wider range. Note
      ! that the file name deliberately carries no descriptor of any object's parameters: the retained fractions and energies are
      ! a property of the truncated Maxwell-Boltzmann distribution alone, and so are common to every model. It is keyed instead
      ! on this module's source digest, since the code which generates the tabulation is the only thing its values depend on.
      if (File_Exists(fileName)) then
         ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
         call File_Lock(fileName,fileLock,lockIsShared=.true.)
         call retainedRestore(fileName)
         call File_Unlock(fileLock)
      end if
      ! Re-evaluate the ranges required in the light of anything restored.
      call rangesRequired()
      if     (                                                &
           &  .not.latticeVelocityEscape%covers(latticeEscape) &
           &  .or.                                            &
           &  .not.latticeVelocityKick  %covers(latticeKick  ) &
           & ) then
         ! Extend the tabulation onto the new lattices, preserving every solution already found. Both axes may grow, so the
         ! surviving block is a rectangle offset within the new arrays.
         countVelocityEscape=latticeEscape%count
         countVelocityKick  =latticeKick  %count
         ! Take the abscissae from the lattices rather than by subdividing the range, so that they are bit-identical to those of
         ! any other tabulation built on the same lattice.
         velocitiesEscapeScaleFree_=latticeEscape%values()
         velocitiesKickScaleFree_  =latticeKick  %values()
         allocate(fractionRetained_(countVelocityEscape,countVelocityKick))
         allocate(energyRetained_  (countVelocityEscape,countVelocityKick))
         allocate(isComputed       (countVelocityEscape,countVelocityKick))
         fractionRetained_=0.0d0
         energyRetained_  =0.0d0
         isComputed       =.false.
         if     (                                   &
              &   allocated(fractionRetained)       &
              &  .and.                              &
              &   allocated(energyRetained  )       &
              &  .and.                              &
              &   latticeVelocityEscape%isDefined() &
              &  .and.                              &
              &   latticeVelocityKick  %isDefined() &
              & ) then
            offsetVelocityEscape=Range_Lattice_Offset(latticeVelocityEscape,latticeEscape)
            offsetVelocityKick  =Range_Lattice_Offset(latticeVelocityKick  ,latticeKick  )
            fractionRetained_(offsetVelocityEscape+1:offsetVelocityEscape+size(fractionRetained,dim=1),offsetVelocityKick+1:offsetVelocityKick+size(fractionRetained,dim=2))=fractionRetained
            energyRetained_  (offsetVelocityEscape+1:offsetVelocityEscape+size(energyRetained  ,dim=1),offsetVelocityKick+1:offsetVelocityKick+size(energyRetained  ,dim=2))=energyRetained
            isComputed       (offsetVelocityEscape+1:offsetVelocityEscape+size(fractionRetained,dim=1),offsetVelocityKick+1:offsetVelocityKick+size(fractionRetained,dim=2))=.true.
         end if
         call displayIndent     ('Tabulating decaying dark matter retained fractions',verbosity=verbosityLevelStandard)
         write (labelMinimum,'(f8.2)') latticeEscape%minimum()
         write (labelMaximum,'(f8.2)') latticeEscape%maximum()
         message=labelMinimum//' ≤ vₑ/σ  ≤ '//labelMaximum
         call displayMessage(message,verbosity=verbosityLevelStandard)
         write (labelMinimum,'(f8.2)') latticeKick  %minimum()
         write (labelMaximum,'(f8.2)') latticeKick  %maximum()
         message=labelMinimum//' ≤ vₖ/σ  ≤ '//labelMaximum
         call displayMessage(message,verbosity=verbosityLevelStandard)
         ! Begin parallel execution.
         call displayCounter(                                  &
              &                        0                     , &
              &              isNew    =.true.                , &
              &              verbosity=verbosityLevelStandard  &
              &             )
         !$omp parallel
         ! Build root finders and integrators.
         allocate(integratorFraction)
         allocate(integratorEnergy  )
         allocate(finder            )
         integratorFraction=integrator(fractionRetainedIntegrand,toleranceRelative=1.0d-6,toleranceAbsolute=1.0d-9)
         integratorEnergy  =integrator(energyRetainedIntegrand  ,toleranceRelative=1.0d-6,toleranceAbsolute=1.0d-9)
         finder            =rootFinder(                                     &
              &                        rootFunction     =velocityWidthRoot, &
              &                        toleranceAbsolute=toleranceAbsolute, &
              &                        toleranceRelative=toleranceRelative  &
              &                       )
         call finder%rangeExpand(                                                             &
              &                  rangeExpandUpward            =2.0d0                        , &
              &                  rangeExpandDownward          =0.50d0                       , &
              &                  rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
              &                  rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
              &                  rangeExpandType              =rangeExpandMultiplicative      &
              &                 )
         !$omp do schedule(dynamic)
         ! Iterate over escape and kick velocities, skipping any solution which is already known.
         do iVelocityEscape=1,countVelocityEscape
            call displayCounter(                                                                          &
                 &                        int(100.0d0*dble(iVelocityEscape-1)/dble(countVelocityEscape)), &
                 &              isNew    =.false.                                                       , &
                 &              verbosity=verbosityLevelStandard                                          &
                 &             )
            ! Skip the whole row - including the solution for the velocity width - if every solution in it is already known.
            if (all(isComputed(iVelocityEscape,:))) cycle
            velocityEscapeScaleFree_=velocitiesEscapeScaleFree_(iVelocityEscape)
            ! Solve for the velocity width of the Maxwell-Boltzmann distribution that gives the required root-mean-squared velocity.
            velocityWidthScaleFree=finder%find(rootGuess=1.0d0)
            ! Evaluate the normalization of the truncated Maxwell-Boltzmann distribution.
            normalization=+sqrt(2.0d0*Pi)*velocityWidthScaleFree**3                         *erf(        velocityEscapeScaleFree_/sqrt(2.0d0)/velocityWidthScaleFree    ) &
                 &        -     2.0d0    *velocityWidthScaleFree**2*velocityEscapeScaleFree_*exp(-0.5d0*(velocityEscapeScaleFree_            /velocityWidthScaleFree)**2)
            do iVelocityKick=1,countVelocityKick
               if (isComputed(iVelocityEscape,iVelocityKick)) cycle
               velocityKickScaleFree_=velocitiesKickScaleFree_(iVelocityKick)
               ! Compute the kick energy retained. If the kick velocity exceeds twice the escape velocity then no particles remain bound.
               if (velocityKickScaleFree_ >= 2.0d0*velocityEscapeScaleFree_ .or. velocityKickScaleFree_ <= 0.0d0) then
                  ! Return zero retained fraction and energy in this case.
                  energyRetained_  (iVelocityEscape,iVelocityKick)=0.0d0
                  fractionRetained_(iVelocityEscape,iVelocityKick)=0.0d0
               else
                  ! If the kick velocity exceeds the escape velocity, there is a maximum angle θ beyond which no velocity remains bound.
                  if (velocityKickScaleFree_ >= velocityEscapeScaleFree_) then
                     ! Compute the maximum angle.
                     cosThetaMaximum=cos(Pi-asin(velocityEscapeScaleFree_/velocityKickScaleFree_))
                  else
                     ! For smaller kick velocities, all angles have some velocity for which particles remain bound.
                     cosThetaMaximum=+1.0d0
                  end if
                  ! Integrate to find the retained fraction and energy.
                  energyRetained_  (iVelocityEscape,iVelocityKick)=integratorEnergy  %integrate(-1.0d0,cosThetaMaximum)/normalization
                  fractionRetained_(iVelocityEscape,iVelocityKick)=integratorFraction%integrate(-1.0d0,cosThetaMaximum)/normalization
               end if
            end do
         end do
         !$omp end do
         ! Transfer solutions to module-scope. This is done inside the parallel region because the module-scope variables are
         ! threadprivate, so each thread must be given its own copy.
         if (allocated(velocitiesEscapeScaleFree)) deallocate(velocitiesEscapeScaleFree)
         if (allocated(velocitiesKickScaleFree  )) deallocate(velocitiesKickScaleFree  )
         if (allocated(energyRetained           )) deallocate(energyRetained           )
         if (allocated(fractionRetained         )) deallocate(fractionRetained         )
         velocitiesEscapeScaleFree=velocitiesEscapeScaleFree_
         velocitiesKickScaleFree  =velocitiesKickScaleFree_
         energyRetained           =energyRetained_
         fractionRetained         =fractionRetained_
         latticeVelocityEscape    =latticeEscape
         latticeVelocityKick      =latticeKick
         ! Clean up objects.
         deallocate(integratorFraction)
         deallocate(integratorEnergy  )
         deallocate(finder            )
         !$omp end parallel
         call displayCounterClear(        verbosity=verbosityLevelStandard)
         call displayUnindent     ('done',verbosity=verbosityLevelStandard)
         ! Store the results to file, recording the lattices so that the tabulation can be restored onto exactly these abscissae.
         ! The lock is taken only for the write - the solutions above were computed without holding it, so that concurrent
         ! processes are not serialized behind a tabulation which may take minutes.
         call File_Lock(fileName,fileLock,lockIsShared=.false.)
         block
           type(hdf5File  ) :: file
           !$ call hdf5Access%set()
           file=hdf5File(fileName,overWrite=.true.,readOnly=.false.)
           call file%writeDataset  (                      velocitiesEscapeScaleFree,'velocitiesEscape'          )
           call file%writeDataset  (                      velocitiesKickScaleFree  ,'velocitiesKick'            )
           call file%writeDataset  (                      energyRetained           ,'energyRetained'            )
           call file%writeDataset  (                      fractionRetained         ,'fractionRetained'          )
           call file%writeAttribute(latticeVelocityEscape%pointsPer                ,'pointsPerVelocityEscape'   )
           call file%writeAttribute(latticeVelocityEscape%indexMinimum             ,'indexMinimumVelocityEscape')
           call file%writeAttribute(latticeVelocityEscape%count                    ,'countVelocityEscape'       )
           call file%writeAttribute(latticeVelocityKick  %pointsPer                ,'pointsPerVelocityKick'     )
           call file%writeAttribute(latticeVelocityKick  %indexMinimum             ,'indexMinimumVelocityKick'  )
           call file%writeAttribute(latticeVelocityKick  %count                    ,'countVelocityKick'         )
           !$ call hdf5Access%unset()
         end block
         call File_Unlock(fileLock)
      end if
      ! Build the interpolators.
      if (allocated(interpolatorVelocityEscape)) deallocate(interpolatorVelocityEscape)
      if (allocated(interpolatorVelocityKick  )) deallocate(interpolatorVelocityKick  )
      allocate(interpolatorVelocityEscape)
      allocate(interpolatorVelocityKick  )
      interpolatorVelocityEscape=interpolator(velocitiesEscapeScaleFree)
      interpolatorVelocityKick  =interpolator(velocitiesKickScaleFree  )
    end block
    call velocityEscapeScaleFreeClamp()
    return

  contains

    subroutine rangesRequired()
      !!{RST
      Find the ranges of scale-free escape and kick velocity over which solutions must be tabulated, as sets of points on
      absolute lattices. Both axes are pinned to whole decades, which gives the coarsest - and so most reproducible - set of
      possible ranges. Each is seeded with a default range, which guarantees that any two tabulations overlap and hence can
      always be merged.
      !!}
      implicit none

      latticeEscape=Range_Pinned(                                                                      &
           &                                    velocityEscapeScaleFree                              , &
           &                                    countVelocityEscapePerDecade                         , &
           &                                    gridSchemePerDecade                                  , &
           &                     rangeCurrent  =[velocityEscapeScaleFreeLimit,velocityScaleFreeLarge], &
           &                     latticeCurrent=latticeVelocityEscape                                , &
           &                     limitMinimum  =velocityEscapeScaleFreeLimit                           &
           &                    )
      latticeKick  =Range_Pinned(                                                                      &
           &                                    velocityKickScaleFree                                , &
           &                                    countVelocityKickPerDecade                           , &
           &                                    gridSchemePerDecade                                  , &
           &                     rangeCurrent  =[velocityKickScaleFreeSeed   ,velocityScaleFreeLarge], &
           &                     latticeCurrent=latticeVelocityKick                                    &
           &                    )
      return
    end subroutine rangesRequired

    subroutine velocityEscapeScaleFreeClamp()
      !!{RST
      Clamp the requested scale-free escape velocity to the lower bound of the tabulation. That bound is the limiting value,
      ``velocityEscapeScaleFreeLimit``, snapped *inward* to a lattice point, and so lies marginally above the limit itself -
      while callers clamp their requests to precisely the limit. Without this the interpolators, which abort on out-of-range
      arguments, would be asked to extrapolate below the first tabulated point.
      !!}
      implicit none

      velocityEscapeScaleFree=max(velocityEscapeScaleFree,latticeVelocityEscape%minimum())
      return
    end subroutine velocityEscapeScaleFreeClamp

    subroutine retainedRestore(fileName)
      !!{RST
      Adopt the tabulation cached in ``fileName``, if there is one and if it spans everything already held in memory. Because
      the tabulation lies on absolute lattices, the cached solutions share abscissae with those in memory wherever the two
      overlap, so adopting the cache in this case discards nothing.
      !!}
      use :: ISO_Varying_String, only : varying_string
      implicit none
      type   (varying_string), intent(in   ) :: fileName
      type   (rangeLattice   )               :: latticeEscapeCached   , latticeKickCached
      integer                                :: pointsPerEscapeCached , indexMinimumEscapeCached, &
           &                                    countEscapeCached     , pointsPerKickCached     , &
           &                                    indexMinimumKickCached, countKickCached

      !$ call hdf5Access%set()
      hdf5FileScope: block
        type            (hdf5File)                              :: file
        double precision          , dimension(:  ), allocatable :: velocitiesEscapeCached, velocitiesKickCached
        double precision          , dimension(:,:), allocatable :: energyCached          , fractionCached
        ! Open read-only: opening read-write would have HDF5 take an exclusive lock on the file, so that a second thread reading
        ! it concurrently would fail to open it at all.
        file=hdf5File(fileName,overWrite=.false.,readOnly=.true.)
        ! A file written before the tabulation was pinned carries no record of its lattices, and cannot be placed on them, so it
        ! is simply ignored and overwritten.
        if (file%hasAttribute('pointsPerVelocityEscape')) then
           call file%readAttribute('pointsPerVelocityEscape'   ,pointsPerEscapeCached   )
           call file%readAttribute('indexMinimumVelocityEscape',indexMinimumEscapeCached)
           call file%readAttribute('countVelocityEscape'       ,countEscapeCached       )
           call file%readAttribute('pointsPerVelocityKick'     ,pointsPerKickCached     )
           call file%readAttribute('indexMinimumVelocityKick'  ,indexMinimumKickCached  )
           call file%readAttribute('countVelocityKick'         ,countKickCached         )
           latticeEscapeCached=rangeLattice(gridSchemePerDecade,pointsPerEscapeCached,indexMinimumEscapeCached,countEscapeCached)
           latticeKickCached  =rangeLattice(gridSchemePerDecade,pointsPerKickCached  ,indexMinimumKickCached  ,countKickCached  )
           ! Adopt the cache if we hold nothing yet, or if it spans everything we do hold - in both dimensions. Note that
           ! `covers` is false whenever either lattice is undefined, so the "nothing held yet" case must be tested explicitly;
           ! it is the usual one, since this is how a freshly started thread acquires the tabulation.
           if     (                                                           &
                &   latticeEscapeCached%isDefined()                           &
                &  .and.                                                      &
                &   latticeKickCached  %isDefined()                           &
                &  .and.                                                      &
                &   (                                                         &
                &     .not.                                                   &
                &      latticeVelocityEscape%isDefined(                     ) &
                &    .or.                                                     &
                &      latticeEscapeCached  %covers   (latticeVelocityEscape) &
                &   )                                                         &
                &  .and.                                                      &
                &   (                                                         &
                &     .not.                                                   &
                &      latticeVelocityKick  %isDefined(                     ) &
                &    .or.                                                     &
                &      latticeKickCached    %covers   (latticeVelocityKick  ) &
                &   )                                                         &
                & ) then
              call file%readDataset('velocitiesEscape',velocitiesEscapeCached)
              call file%readDataset('velocitiesKick'  ,velocitiesKickCached  )
              call file%readDataset('energyRetained'  ,energyCached          )
              call file%readDataset('fractionRetained',fractionCached        )
              ! Adopt only if the datasets match the lattices recorded alongside them; otherwise the file is not
              ! self-consistent, and everything read from it is discarded rather than leaving a partially-restored tabulation.
              if     (                                                                 &
                   &   size(velocitiesEscapeCached      ) == latticeEscapeCached%count &
                   &  .and.                                                            &
                   &   size(velocitiesKickCached        ) == latticeKickCached  %count &
                   &  .and.                                                            &
                   &   size(fractionCached        ,dim=1) == latticeEscapeCached%count &
                   &  .and.                                                            &
                   &   size(fractionCached        ,dim=2) == latticeKickCached  %count &
                   &  .and.                                                            &
                   &   size(energyCached          ,dim=1) == latticeEscapeCached%count &
                   &  .and.                                                            &
                   &   size(energyCached          ,dim=2) == latticeKickCached  %count &
                   & ) then
                 call Move_Alloc(fractionCached,fractionRetained)
                 call Move_Alloc(energyCached  ,energyRetained  )
                 ! Take the abscissae from the lattices rather than from the file, so that they are bit-identical to those of
                 ! any other tabulation built on the same lattice.
                 if (allocated(velocitiesEscapeScaleFree)) deallocate(velocitiesEscapeScaleFree)
                 if (allocated(velocitiesKickScaleFree  )) deallocate(velocitiesKickScaleFree  )
                 allocate(velocitiesEscapeScaleFree(latticeEscapeCached%count))
                 allocate(velocitiesKickScaleFree  (latticeKickCached  %count))
                 velocitiesEscapeScaleFree=latticeEscapeCached%values()
                 velocitiesKickScaleFree  =latticeKickCached  %values()
                 latticeVelocityEscape    =latticeEscapeCached
                 latticeVelocityKick      =latticeKickCached
              end if
           end if
        end if
      end block hdf5FileScope
      !$ call hdf5Access%unset()
      return
    end subroutine retainedRestore

    double precision function fractionRetainedIntegrand(cosTheta) result(integrand)
      !!{RST
      The integrand used to find the retained fraction of particles.
      !!}
      implicit none
      double precision, intent(in   ) :: cosTheta
      double precision                :: velocityMinimumScaleFree, velocityMaximumScaleFree
      
      call velocityLimitsScaleFree(cosTheta,velocityMinimumScaleFree,velocityMaximumScaleFree)
      if (velocityMaximumScaleFree > velocityMinimumScaleFree) then
         integrand=+fractionRetainedIntegrandIndefinite(velocityMaximumScaleFree) &
              &    -fractionRetainedIntegrandIndefinite(velocityMinimumScaleFree)
      else
         integrand=+0.0d0
      end if
      return
    end function fractionRetainedIntegrand

    double precision function energyRetainedIntegrand(cosTheta) result(integrand)
      !!{RST
      The integrand used to find the retained energy of particles.
      !!}
      implicit none
      double precision, intent(in   ) :: cosTheta
      double precision                :: velocityMinimumScaleFree, velocityMaximumScaleFree
      
      call velocityLimitsScaleFree(cosTheta,velocityMinimumScaleFree,velocityMaximumScaleFree)
      if (velocityMaximumScaleFree > velocityMinimumScaleFree) then
         integrand=+energyRetainedIntegrandIndefinite(cosTheta,velocityMaximumScaleFree) &
              &    -energyRetainedIntegrandIndefinite(cosTheta,velocityMinimumScaleFree)
      else
         integrand=+0.0d0
      end if
      return
    end function energyRetainedIntegrand

    double precision function fractionRetainedIntegrandIndefinite(velocityScaleFree) result(integrand)
      !!{RST
      The indefinite integral over velocity of the retained fraction of particles (note that the normalization factor is not included here; see eqn. :eq:`eq-decayingDMRetainedFraction`):

      .. math::

         \sqrt{2 \pi} s^3 \hbox{erf}\left( \frac{v}{\sqrt{2}s} \right) - 2 v s^2 \exp \left( -\frac{1}{2}\left[\frac{v}{s}\right]^2\right).
      !!}
      implicit none
      double precision, intent(in   ) :: velocityScaleFree

      integrand=+sqrt(Pi/2.0d0)*                  velocityWidthScaleFree**3*erf(        velocityScaleFree/sqrt(2.0d0)/velocityWidthScaleFree    ) &
           &    -               velocityScaleFree*velocityWidthScaleFree**2*exp(-0.5d0*(velocityScaleFree            /velocityWidthScaleFree)**2)
      return
    end function fractionRetainedIntegrandIndefinite
    
    double precision function energyRetainedIntegrandIndefinite(cosTheta,velocityScaleFree) result(integrand)
      !!{RST
      The indefinite integral over velocity of the retained energy of particles (note that the normalization factor is not included here; see eqn. :eq:`eq-decayingDMRetainedEnergy`):

      .. math::

         \frac{1}{4} v_\mathrm{k} s^2 \left\{ 8 s^2 \cos\theta - 2 (v v_\mathrm{k} + 2 [v^2+2s^2]\cos\theta) \exp\left(-\frac{1}{2}\left[\frac{v}{s}\right]^2\right) + \sqrt{2 \pi} v_\mathrm{k} s \, \hbox{erf}\left(\frac{v}{\sqrt{2}s}\right) \right\}.
      !!}
      implicit none
      double precision, intent(in   ) :: cosTheta,velocityScaleFree


      integrand=+       0.25d0    * velocityKickScaleFree_*velocityWidthScaleFree**2                                                                                                                                              &
           &    *(                                                                                                                                                                                                                &
           &      +     8.00d0    *                                                                                    velocityWidthScaleFree**2 *cosTheta                                                                        &
           &      -     2.00d0    *(velocityKickScaleFree_*velocityScaleFree        +2.0d0*(velocityScaleFree**2+2.0d0*velocityWidthScaleFree**2)*cosTheta)*exp(-0.5d0*(velocityScaleFree            /velocityWidthScaleFree)**2) &
           &      +sqrt(2.00d0*Pi)* velocityKickScaleFree_*velocityWidthScaleFree                                                                          *erf(        velocityScaleFree/sqrt(2.0d0)/velocityWidthScaleFree    ) &
           &     )
      return
    end function energyRetainedIntegrandIndefinite
    
    subroutine velocityLimitsScaleFree(cosTheta,velocityMinimumScaleFree,velocityMaximumScaleFree)
      !!{RST
      Find the minimum and maximum velocity, as a function of :math:`\cos\theta`, for which particles remain bound: :math:`v_\mathrm{max|min}(\theta) = \pm \left( v_\mathrm{e}^2 - v_\mathrm{k}^2 \sin^2\theta \right)^{1/2} - v_\mathrm{k} \cos \theta`.
      !!}
      implicit none
      double precision, intent(in   ) :: cosTheta
      double precision, intent(  out) :: velocityMinimumScaleFree, velocityMaximumScaleFree
      double precision                :: argument

      argument=velocityEscapeScaleFree_**2-velocityKickScaleFree_**2*(1.0d0-cosTheta**2)
      if (argument < 0.0d0) then
         ! No solution exists - particles do not remain bound for any velocity.
         velocityMinimumScaleFree=                                0.0d0
         velocityMaximumScaleFree=                                0.0d0
      else
         ! Evaluate the minimum and maximum velocities.
         velocityMinimumScaleFree=min(velocityEscapeScaleFree_,max(0.0d0,-sqrt(argument)-velocityKickScaleFree_*cosTheta))
         velocityMaximumScaleFree=min(velocityEscapeScaleFree_,max(0.0d0,+sqrt(argument)-velocityKickScaleFree_*cosTheta))
      end if
      return
    end subroutine velocityLimitsScaleFree
    
    double precision function velocityWidthRoot(velocityWidthScaleFree)
      !!{RST
      The root function used to find the velocity width parameter of the truncated Maxwell-Boltzmann distribution such that the 3D, scale-free (i.e. all velocities are expressed in units of the 1D root-mean-squared velocity, :math:`\sigma`) root-mean-squared speed is 3 (i.e. a non-scale-free root-mean-squared velocity of :math:`3\sigma` as expected for a non-truncated Maxwell-Boltzmann distribution). This ensures that the kinetic energy density in the distribution is as predicted by a Jeans analysis.

      This requires that:

      .. math::

         \int_{-1}^{+_1} \mathrm{d}\cos\theta \int_0^{v_\mathrm{e}} v^2 p(v,\theta|s) \mathrm{d}v = 3,

      which implies that:

      .. math::

         3 \sqrt{2 \pi} s^3 \, \hbox{erf}\left(\frac{v_\mathrm{e}}{\sqrt{2}s}\right) - 2 (3 s^2 + v_\mathrm{e}^2) v_\mathrm{e} \exp\left(-\frac{1}{2}\left[\frac{v_\mathrm{e}}{s}\right]^2\right) = 3 \left\{ \sqrt{2\pi} s \, \hbox{erf}\left(\frac{v_\mathrm{e}}{\sqrt{2}s}\right) -2 v_\mathrm{e} \exp\left(-\frac{1}{2}\left[\frac{v_\mathrm{e}}{s}\right]^2\right)  \right\}.
      !!}
      implicit none
      double precision, intent(in   ) :: velocityWidthScaleFree

      velocityWidthRoot=+  3.0d0*sqrt(2.0d0*Pi)*velocityWidthScaleFree **3                                                              *erf(        velocityEscapeScaleFree_/sqrt(2.0d0)/velocityWidthScaleFree    ) &
           &            -  2.0d0               *velocityEscapeScaleFree_  *(3.0d0*velocityWidthScaleFree**2+velocityEscapeScaleFree_**2)*exp(-0.5d0*(velocityEscapeScaleFree_            /velocityWidthScaleFree)**2) &
           &            -  3.0d0                                                                                                                                                                                      &
           &            *(                                                                                                                                                                                            &
           &              +      sqrt(2.0d0*Pi)*velocityWidthScaleFree                                                                  *erf(        velocityEscapeScaleFree_/sqrt(2.0d0)/velocityWidthScaleFree    ) &
           &              -2.0d0               *velocityEscapeScaleFree_                                                                *exp(-0.5d0*(velocityEscapeScaleFree_            /velocityWidthScaleFree)**2) &
           &             )
      return
    end function velocityWidthRoot

  end subroutine decayingDarkMatterRetainedTabulate

end module Decaying_Dark_Matter
