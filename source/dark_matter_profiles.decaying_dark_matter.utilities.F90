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

!!{
Contains a module which implements useful shared utilities for calculations of decaying dark matter.
!!}

module Decaying_Dark_Matter
  !!{
  Implements useful shared utilities for calculations of decaying dark matter.
  !!}
  use :: Numerical_Interpolation, only : interpolator
  private
  public :: decayingDarkMatterFractionRetained, decayingDarkMatterEnergyRetained, decayingDarkMatterFractionRetainedDerivatives, decayingDarkMatterEnergyRetainedDerivatives

  ! Tables of retained fractions and energies.
  double precision                                            :: velocityEscapeScaleFreeMinimum=+huge(0.0d0), velocityEscapeScaleFreeMaximum=-huge(0.0d0), &
       &                                                         velocityKickScaleFreeMinimum  =+huge(0.0d0), velocityKickScaleFreeMaximum  =-huge(0.0d0)
  double precision              , dimension(:  ), allocatable :: velocitiesEscapeScaleFree                  , velocitiesKickScaleFree
  double precision              , dimension(:,:), allocatable :: energyRetained                             , fractionRetained
  type            (interpolator)                , allocatable :: interpolatorVelocityEscape                 , interpolatorVelocityKick
  !$omp threadprivate(velocityEscapeScaleFreeMinimum,velocityEscapeScaleFreeMaximum,velocityKickScaleFreeMinimum,velocityKickScaleFreeMaximum,velocitiesEscapeScaleFree,velocitiesKickScaleFree,energyRetained,fractionRetained,interpolatorVelocityEscape,interpolatorVelocityKick)

  ! Scale free velocity above which the velocity dispersion is negligible and we assume deterministic behavior.
  double precision                              , parameter   :: velocityScaleFreeLarge        =+100.0d0

  ! Minimum scale free escape velocity for which we can find a self-consistent velocity distribution function. Smaller values only
  ! occur at extreme distances outside of halos, so should not matter.
  double precision                              , parameter   :: velocityEscapeScaleFreeLimit  =+  2.3d0

contains

  double precision function decayingDarkMatterFractionRetained(velocityDispersion,velocityEscape,velocityKick) result(fraction)
    !!{
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
    !!{
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
    !!{
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
    !!{
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
    !!{
    Compute the fraction of the decaying dark matter particles and kick energy that is retained. Assumes that the initial
    distribution of particle velocities is a Maxwell-Boltzmann distribution, truncated at the escape velocity. The mean energy of
    retained particles, minus their mean energy pre-kick is computed.

    The truncated Maxwell-Boltzmann distribution is
    \begin{equation}
      p(v,\theta|s) = \left\{ \begin{array}{ll} A^{-1} v^2 \exp\left(-\frac{1}{2}\left[\frac{v}{s}\right]^2\right) & \hbox{if } v < v_\mathrm{e} \\ 0 &  \hbox{if } v \ge v_\mathrm{e}, \end{array} \right.
    \end{equation}
    where $s$ is the velocity width, $v_\mathrm{e}$ is the escape velocity, $v$ is particle speed, $\theta$ is the direction of
    the particle velocity relative to the $z$-axis, and
    \begin{equation}
      A = \sqrt{2 \pi} s^3 \hbox{erf}\left( \frac{v_\mathrm{e}}{\sqrt{2}s} \right) - 2 v_\mathrm{e} s^2 \exp \left( -\frac{1}{2}\left[\frac{v_\mathrm{e}}{s}\right]^2\right)
    \end{equation}
    is a normalization factor.

    Assume, without loss of generality, that the kick is along the $z$-axis. The specific kinetic energy of the retained particles,
    in excess of their original energy is:
    \begin{equation}
      \epsilon =  \int_{-1}^{+1} \mathrm{d}\cos\theta \int_0^{v_\mathrm{e}} \mathrm{d}v \frac{1}{2} \left( v^2 + v_\mathrm{k}^2 + 2 v v_\mathrm{k} \cos\theta - v^2 \right) p(v,\theta|s) H\left( v^2 + v_\mathrm{k}^2 + 2 v v_\mathrm{k} \cos\theta < v_\mathrm{e}^2 \right)
    \end{equation}

    where $v_\mathrm{k}$ is the scale-free kick velocity, and $H(x) = 1$ if $x$ is true, and 0 otherwise. Solving the inequality
    for the velocity, $v$, that will remain bound as a function of $\theta$, gives
    $v_\mathrm{max|min}(\theta) = \pm \left( v_\mathrm{e}^2 - v_\mathrm{k}^2 \sin^2\theta \right)^{1/2} - v_\mathrm{k} \cos
    \theta$, so:
    \begin{equation}
      \epsilon =  \int_{-1}^{+_1} \mathrm{d}\cos\theta \int_{v_\mathrm{min}(\theta)}^{v_\mathrm{max}(\theta)} \mathrm{d}v \frac{1}{2} \left( v_\mathrm{k}^2 + 2 v v_\mathrm{k} \cos\theta \right) p(v,\theta|s).
      \label{eq:decayingDMRetainedEnergy}
    \end{equation}
    Similarly, the retained fraction of particles is simply
    \begin{equation}
      f =  \int_{-1}^{+_1} \mathrm{d}\cos\theta \int_{v_\mathrm{min}(\theta)}^{v_\mathrm{max}(\theta)} \mathrm{d}v p(v,\theta|s).
      \label{eq:decayingDMRetainedFraction}
    \end{equation}
    !!}
    use :: Display                 , only : displayCounter, displayCounterClear          , displayIndent                , displayUnindent          , &
         &                                  displayMessage, verbosityLevelStandard
    use :: ISO_Varying_String      , only : varying_string, operator(//)                 , assignment(=)                , char
    use :: Numerical_Constants_Math, only : Pi
    use :: Root_Finder             , only : rootFinder    , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rangeExpandMultiplicative
    use :: Numerical_Integration   , only : integrator
    use :: Numerical_Ranges        , only : Make_Range    , rangeTypeLogarithmic
    use :: File_Utilities          , only : File_Exists   , Directory_Make               , File_Path, File_Lock         , File_Unlock, lockDescriptor
    use :: HDF5_Access             , only : hdf5Access
    use :: IO_HDF5                 , only : hdf5Object
    use :: Input_Paths             , only : inputPath     , pathTypeDataDynamic
    implicit none
    double precision            , intent(in   )               :: velocityEscapeScaleFree              , velocityKickScaleFree
    double precision            , parameter                   :: toleranceAbsolute           =  0.0d+0, toleranceRelative         =  1.0d-6
    double precision            , parameter                   :: countVelocityEscapePerDecade=300.0d+0, countVelocityKickPerDecade=300.0d+0
    double precision            , dimension(:  ), allocatable :: velocitiesEscapeScaleFree_           , velocitiesKickScaleFree_
    double precision            , dimension(:,:), allocatable :: energyRetained_                      , fractionRetained_
    double precision            , save                        :: velocityWidthScaleFree               , normalization                      , &
         &                                                       velocityEscapeScaleFree_             , velocityKickScaleFree_             , &
         &                                                       cosThetaMaximum
    type            (rootFinder), save          , allocatable :: finder
    type            (integrator), save          , allocatable :: integratorEnergy                     , integratorFraction
    integer                     , save                        :: iVelocityKick
    integer                                                   :: countVelocityEscape                  , countVelocityKick                  , &
         &                                                       iVelocityEscape
    logical                                                   :: remakeTable
    !$omp threadprivate(iVelocityKick,velocityEscapeScaleFree_,velocityKickScaleFree_,velocityWidthScaleFree,normalization,cosThetaMaximum,finder,integratorFraction,integratorEnergy)
    
    ! Determine if the tables must be remade.
    remakeTable= velocityEscapeScaleFree < velocityEscapeScaleFreeMinimum .or. velocityEscapeScaleFree > velocityEscapeScaleFreeMaximum &
         &      .or.                                                                                                                    &
         &       velocityKickScaleFree   < velocityKickScaleFreeMinimum   .or. velocityKickScaleFree   > velocityKickScaleFreeMaximum
    if (.not.remakeTable) return
    block
      type     (varying_string) :: message     , fileName
      character(len=8         ) :: labelMinimum, labelMaximum
      type     (hdf5Object    ) :: file
      type     (lockDescriptor) :: fileLock

      ! Get a lock on the file.
      fileName=inputPath(pathTypeDataDynamic)//'darkMatter/decayingDarkMatterRetention.hdf5'
      call Directory_Make(File_Path(fileName)                              )
      call File_Lock     (          fileName ,fileLock,lockIsShared=.false.)
      ! Attempt to read existing data from file.
      if (File_Exists(fileName)) then
         if (allocated(velocitiesEscapeScaleFree)) deallocate(velocitiesEscapeScaleFree)
         if (allocated(velocitiesKickScaleFree  )) deallocate(velocitiesKickScaleFree  )
         if (allocated(fractionRetained         )) deallocate(fractionRetained         )
         if (allocated(energyRetained           )) deallocate(energyRetained           )
         !$ call hdf5Access%set()
         call file%openFile   (char(fileName),overWrite=.false.,readOnly=.true.)
         call file%readDataset('velocitiesEscape',velocitiesEscapeScaleFree)
         call file%readDataset('velocitiesKick'  ,velocitiesKickScaleFree  )
         call file%readDataset('energyRetained'  ,energyRetained           )
         call file%readDataset('fractionRetained',fractionRetained         )
         call file%close      (                                            )
         !$ call hdf5Access%unset()
         velocityEscapeScaleFreeMinimum=velocitiesEscapeScaleFree(                             1 )
         velocityEscapeScaleFreeMaximum=velocitiesEscapeScaleFree(size(velocitiesEscapeScaleFree))
         velocityKickScaleFreeMinimum  =velocitiesKickScaleFree  (                             1 )
         velocityKickScaleFreeMaximum  =velocitiesKickScaleFree  (size(velocitiesKickScaleFree  ))
         remakeTable= velocityEscapeScaleFree < velocityEscapeScaleFreeMinimum .or. velocityEscapeScaleFree > velocityEscapeScaleFreeMaximum &
              &      .or.                                                                                                                    &
              &       velocityKickScaleFree   < velocityKickScaleFreeMinimum   .or. velocityKickScaleFree   > velocityKickScaleFreeMaximum
      end if
      if (remakeTable) then
         ! Find suitable ranges of escape velocity and kick velocity over which to tabulate. For the escape velocity we limit the
         ! minimum value to a fixed multiple of the velocity dispersion. Below this it becomes difficult to find a velocity width
         ! for the truncated Maxwell-Boltzmann distribution.
         velocityEscapeScaleFreeMinimum=min(velocityEscapeScaleFreeMinimum,max(velocityEscapeScaleFreeLimit,0.5d0*velocityEscapeScaleFree))
         velocityEscapeScaleFreeMaximum=max(velocityEscapeScaleFreeMaximum,                                 2.0d0*velocityEscapeScaleFree )
         velocityKickScaleFreeMinimum  =min(velocityKickScaleFreeMinimum  ,                                 0.5d0*velocityKickScaleFree   )
         velocityKickScaleFreeMaximum  =max(velocityKickScaleFreeMaximum  ,                                 2.0d0*velocityKickScaleFree   )
         call displayIndent     ('Tabulating decaying dark matter retained fractions',verbosity=verbosityLevelStandard)
         write (labelMinimum,'(f8.2)') velocityEscapeScaleFreeMinimum
         write (labelMaximum,'(f8.2)') velocityEscapeScaleFreeMaximum
         message=labelMinimum//' ≤ vₑ/σ  ≤ '//labelMaximum
         call displayMessage(message,verbosity=verbosityLevelStandard)
         write (labelMinimum,'(f8.2)') velocityKickScaleFreeMinimum
         write (labelMaximum,'(f8.2)') velocityKickScaleFreeMaximum
         message=labelMinimum//' ≤ vₖ/σ  ≤ '//labelMaximum
         call displayMessage(message,verbosity=verbosityLevelStandard)
         ! Construct grids.
         if (allocated(velocitiesEscapeScaleFree)) deallocate(velocitiesEscapeScaleFree)
         if (allocated(velocitiesKickScaleFree  )) deallocate(velocitiesKickScaleFree  )
         if (allocated(fractionRetained         )) deallocate(fractionRetained         )
         if (allocated(energyRetained           )) deallocate(energyRetained           )
         countVelocityEscape       =int(log10(velocityEscapeScaleFreeMaximum/velocityEscapeScaleFreeMinimum)*countVelocityEscapePerDecade+1.0d0)
         countVelocityKick         =int(log10(velocityKickScaleFreeMaximum  /velocityKickScaleFreeMinimum  )*countVelocityKickPerDecade  +1.0d0)
         velocitiesEscapeScaleFree_=Make_Range(velocityEscapeScaleFreeMinimum,velocityEscapeScaleFreeMaximum,countVelocityEscape,rangeTypeLogarithmic)
         velocitiesKickScaleFree_  =Make_Range(velocityKickScaleFreeMinimum  ,velocityKickScaleFreeMaximum  ,countVelocityKick  ,rangeTypeLogarithmic)
         allocate(fractionRetained_(countVelocityEscape,countVelocityKick))
         allocate(energyRetained_  (countVelocityEscape,countVelocityKick))
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
         ! Iterate over escape and kick velocities.
         do iVelocityEscape=1,countVelocityEscape
            call displayCounter(                                                                          &
                 &                        int(100.0d0*dble(iVelocityEscape-1)/dble(countVelocityEscape)), &
                 &              isNew    =.false.                                                       , &
                 &              verbosity=verbosityLevelStandard                                          &
                 &             )
            velocityEscapeScaleFree_=velocitiesEscapeScaleFree_(iVelocityEscape)
            ! Solve for the velocity width of the Maxwell-Boltzmann distribution that gives the required root-mean-squared velocity.
            velocityWidthScaleFree=finder%find(rootGuess=1.0d0)
            ! Evaluate the normalization of the truncated Maxwell-Boltzmann distribution.
            normalization=+sqrt(2.0d0*Pi)*velocityWidthScaleFree**3                         *erf(        velocityEscapeScaleFree_/sqrt(2.0d0)/velocityWidthScaleFree    ) &
                 &        -     2.0d0    *velocityWidthScaleFree**2*velocityEscapeScaleFree_*exp(-0.5d0*(velocityEscapeScaleFree_            /velocityWidthScaleFree)**2)
            do iVelocityKick=1,countVelocityKick
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
         ! Transfer solutions to module-scope.
         if (allocated(velocitiesEscapeScaleFree)) deallocate(velocitiesEscapeScaleFree)
         if (allocated(velocitiesKickScaleFree  )) deallocate(velocitiesKickScaleFree  )
         if (allocated(energyRetained           )) deallocate(energyRetained           )
         if (allocated(fractionRetained         )) deallocate(fractionRetained         )
         velocitiesEscapeScaleFree=velocitiesEscapeScaleFree_
         velocitiesKickScaleFree  =velocitiesKickScaleFree_
         energyRetained           =energyRetained_
         fractionRetained         =fractionRetained_
         ! Clean up objects.
         deallocate(integratorFraction)
         deallocate(integratorEnergy  )
         deallocate(finder            )
         !$omp end parallel
         call displayCounterClear(        verbosity=verbosityLevelStandard)
         call displayUnindent     ('done',verbosity=verbosityLevelStandard)
         ! Store results to file.
         !$ call hdf5Access%set()
         call file%openFile    (char(fileName),overWrite=.true.,readOnly=.false.)
         call file%writeDataset(velocitiesEscapeScaleFree,'velocitiesEscape')
         call file%writeDataset(velocitiesKickScaleFree  ,'velocitiesKick'  )
         call file%writeDataset(energyRetained           ,'energyRetained'  )
         call file%writeDataset(fractionRetained         ,'fractionRetained')
         call file%close       (                                            )
         !$ call hdf5Access%unset()
      end if
      call File_Unlock(fileLock)
      ! Build the interpolators.
      if (allocated(interpolatorVelocityEscape)) deallocate(interpolatorVelocityEscape)
      if (allocated(interpolatorVelocityKick  )) deallocate(interpolatorVelocityKick  )
      interpolatorVelocityEscape=interpolator(velocitiesEscapeScaleFree)
      interpolatorVelocityKick  =interpolator(velocitiesKickScaleFree  )
    end block
    return

  contains

    double precision function fractionRetainedIntegrand(cosTheta) result(integrand)
      !!{
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
      !!{
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
      !!{      
      The indefinite integral over velocity of the retained fraction of particles (note that the normalization factor is not
      included here; see eqn.~\ref{eq:decayingDMRetainedFraction}):
      \begin{equation}
        \sqrt{2 \pi} s^3 \hbox{erf}\left( \frac{v}{\sqrt{2}s} \right) - 2 v s^2 \exp \left( -\frac{1}{2}\left[\frac{v}{s}\right]^2\right).
      \end{equation}
      !!}
      implicit none
      double precision, intent(in   ) :: velocityScaleFree

      integrand=+sqrt(Pi/2.0d0)*                  velocityWidthScaleFree**3*erf(        velocityScaleFree/sqrt(2.0d0)/velocityWidthScaleFree    ) &
           &    -               velocityScaleFree*velocityWidthScaleFree**2*exp(-0.5d0*(velocityScaleFree            /velocityWidthScaleFree)**2)
      return
    end function fractionRetainedIntegrandIndefinite
    
    double precision function energyRetainedIntegrandIndefinite(cosTheta,velocityScaleFree) result(integrand)
      !!{      
      The indefinite integral over velocity of the retained energy of particles (note that the normalization factor is not
      included here; see eqn.~\ref{eq:decayingDMRetainedEnergy}):
      \begin{equation}
        \frac{1}{4} v_\mathrm{k} s^2 \left\{ 8 s^2 \cos\theta - 2 (v v_\mathrm{k} + 2 [v^2+2s^2]\cos\theta) \exp\left(-\frac{1}{2}\left[\frac{v}{s}\right]^2\right) + \sqrt{2 \pi} v_\mathrm{k} s \, \hbox{erf}\left(\frac{v}{\sqrt{2}s}\right) \right\}.
      \end{equation}
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
      !!{
      Find the minimum and maximum velocity, as a function of $\cos\theta$, for which particles remain bound:
      $v_\mathrm{max|min}(\theta) = \pm \left( v_\mathrm{e}^2 - v_\mathrm{k}^2 \sin^2\theta \right)^{1/2} - v_\mathrm{k} \cos
      \theta$.      
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
      !!{      
      The root function used to find the velocity width parameter of the truncated Maxwell-Boltzmann distribution such that the
      3D, scale-free (i.e. all velocities are expressed in units of the 1D root-mean-squared velocity, $\sigma$) root-mean-squared
      speed is 3 (i.e. a non-scale-free root-mean-squared velocity of $3\sigma$ as expected for a non-truncated Maxwell-Boltzmann
      distribution). This ensures that the kinetic energy density in the distribution is as predicted by a Jeans analysis.
      
      This requires that:
      \begin{equation}
        \int_{-1}^{+_1} \mathrm{d}\cos\theta \int_0^{v_\mathrm{e}} v^2 p(v,\theta|s) \mathrm{d}v = 3,
      \end{equation}
      which implies that:
      \begin{equation}
        3 \sqrt{2 \pi} s^3 \, \hbox{erf}\left(\frac{v_\mathrm{e}}{\sqrt{2}s}\right) - 2 (3 s^2 + v_\mathrm{e}^2) v_\mathrm{e} \exp\left(-\frac{1}{2}\left[\frac{v_\mathrm{e}}{s}\right]^2\right) = 3 \left\{ \sqrt{2\pi} s \, \hbox{erf}\left(\frac{v_\mathrm{e}}{\sqrt{2}s}\right) -2 v_\mathrm{e} \exp\left(-\frac{1}{2}\left[\frac{v_\mathrm{e}}{s}\right]^2\right)  \right\}.
      \end{equation}
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
