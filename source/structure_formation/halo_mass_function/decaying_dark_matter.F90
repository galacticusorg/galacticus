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
Contains a module which implements a dark matter halo mass function class for decaying dark matter (DDM)
models, following the revised spherical collapse model of :cite:t:`montandon_decaying_2026`.
!!}

  use :: Dark_Matter_Particles  , only : darkMatterParticleClass
  use :: Numerical_Interpolation, only : interpolator

  !![
  <haloMassFunction name="haloMassFunctionDecayingDarkMatter" docformat="rst">
   <description>
   A dark matter halo mass function class for decaying dark matter (:term:`DDM`) models, following
   the revised spherical collapse model of :cite:t:`montandon_decaying_2026`. The decay of dark matter
   causes each halo to lose mass, so that the observed collapsed mass :math:`M_\mathrm{coll}` is smaller
   than the initial Lagrangian mass :math:`M_0` from which the perturbation formed. This class wraps
   another halo mass function (which should be evaluated in terms of the Lagrangian mass :math:`M_0`)
   and maps it to the observed collapsed mass via the change of variables (their eq. 8),

   .. math::

      \frac{\mathrm{d}n}{\mathrm{d}M_\mathrm{coll}} = \left.\frac{\mathrm{d}n}{\mathrm{d}M_0}\right|_{M_0(M_\mathrm{coll})} \frac{\mathrm{d}M_0}{\mathrm{d}M_\mathrm{coll}} ,

   where the mapping :math:`M_\mathrm{coll}(M_0)` is given by their eq. 46. The mapping is tabulated over
   a grid of Lagrangian masses at fixed resolution ``[pointsPerDecade]`` and inverted by interpolation.
   The grid initially spans ``[massMinimum]`` to ``[massMaximum]``, but is expanded dynamically (by whole
   grid steps, so the resolution is preserved) whenever a collapsed mass outside the currently-tabulated
   range is requested. The decay lifetime and velocity kick are taken from a
   :cite:t:`montandon_decaying_2026` decaying dark matter particle
   (:galacticus-class:`darkMatterParticleDecayingDarkMatter`).

   .. warning::

      In the model of :cite:t:`montandon_decaying_2026` the DDM suppression is carried by a
      mass-dependent critical overdensity for collapse (see
      :galacticus-class:`criticalOverdensityDecayingDarkMatter`), with the variance :math:`\sigma(M)` computed from
      the *unmodified* :math:`\Lambda`CDM linear power spectrum. The wrapped halo mass function must
      therefore (i) be an :math:`f(\nu)`-type mass function that consumes a
      :galacticus-class:`criticalOverdensityDecayingDarkMatter` critical overdensity (this class emits a warning if
      it does not; see ``isCriticalOverdensityDependent``), and (ii) use a standard :math:`\Lambda`CDM
      ``cosmologicalMassVariance``/``transferFunction`` (combining with a suppressed transfer function
      would double-count the suppression).

   .. note::

      **Halo mass definition.** The collapsed mass :math:`M_\mathrm{coll}` returned here is the mass
      associated with the spherical-collapse (:math:`M_\mathrm{virial}`) definition inherited from the
      wrapped mass function, whereas :cite:t:`montandon_decaying_2026` calibrate and measure their halo
      mass function using the :math:`M_\mathrm{200m}` definition (the mass within a sphere enclosing a
      mean density of 200 times the cosmic mean matter density). Converting between these definitions
      requires a model for the halo density profile, which is not currently applied here, so no attempt
      is made to match the :math:`M_\mathrm{200m}` definition. This introduces a systematic offset in the
      mass (and hence in the abundance at fixed mass) when comparing directly to their results.
   </description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionDecayingDarkMatter
     !!{RST
     A halo mass function class for decaying dark matter cosmologies :cite:t:`montandon_decaying_2026`.
     !!}
     private
     class           (haloMassFunctionClass  ), pointer                   :: massFunction_       => null()
     class           (darkMatterParticleClass), pointer                   :: darkMatterParticle_ => null()
     double precision                                                     :: lifetime                           , velocityKick         , &
          &                                                                  massMinimum                        , massMaximum          , &
          &                                                                  pointsPerDecade
     ! Cache of the (epoch-dependent) M_coll(M₀) mapping. The grid of Lagrangian masses spans
     ! [massMinimum,massMaximum], which are expanded dynamically as needed; massCollapsedMinimum and
     ! massCollapsedMaximum are the corresponding limits of the (monotonic) mapped collapsed mass, used to
     ! detect when a requested collapsed mass falls outside the tabulated range.
     double precision                         , allocatable, dimension(:) :: logMass0Table                      , logMassCollapsedTable
     type            (interpolator           )                            :: interpolator_
     double precision                                                     :: timeTabulated       =  -huge(0.0d0)
     double precision                                                     :: massCollapsedMinimum               , massCollapsedMaximum
     logical                                                              :: tabulated           =  .false.
   contains
     final     ::                                   decayingDarkMatterDestructor
     procedure :: differential                   => decayingDarkMatterDifferential
     procedure :: isCriticalOverdensityDependent => decayingDarkMatterIsCriticalOverdensityDependent
  end type haloMassFunctionDecayingDarkMatter

  interface haloMassFunctionDecayingDarkMatter
     !!{RST
     Constructors for the :galacticus-class:`haloMassFunctionDecayingDarkMatter` halo mass function class.
     !!}
     module procedure decayingDarkMatterConstructorParameters
     module procedure decayingDarkMatterConstructorInternal
  end interface haloMassFunctionDecayingDarkMatter

contains

  function decayingDarkMatterConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`haloMassFunctionDecayingDarkMatter` halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloMassFunctionDecayingDarkMatter)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (haloMassFunctionClass             ), pointer       :: massFunction_
    class           (cosmologyParametersClass          ), pointer       :: cosmologyParameters_
    class           (darkMatterParticleClass           ), pointer       :: darkMatterParticle_
    double precision                                                    :: massMinimum         , massMaximum, &
         &                                                                 pointsPerDecade

    !![
    <inputParameter docformat="rst">
      <name>massMinimum</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The initial minimum Lagrangian mass (in :math:`\mathrm{M}_\odot`) used in tabulating the mapping between Lagrangian and collapsed mass. The tabulated range is expanded automatically if a smaller mass is required.</description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>massMaximum</name>
      <source>parameters</source>
      <defaultValue>1.0d18</defaultValue>
      <description>The initial maximum Lagrangian mass (in :math:`\mathrm{M}_\odot`) used in tabulating the mapping between Lagrangian and collapsed mass. The tabulated range is expanded automatically if a larger mass is required.</description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>pointsPerDecade</name>
      <source>parameters</source>
      <defaultValue>100.0d0</defaultValue>
      <description>The number of points per decade of mass used in tabulating the mapping between Lagrangian and collapsed mass.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="haloMassFunction"    name="massFunction_"        source="parameters"/>
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    !!]
    self=haloMassFunctionDecayingDarkMatter(massFunction_,cosmologyParameters_,darkMatterParticle_,massMinimum,massMaximum,pointsPerDecade)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="massFunction_"       />
    <objectDestructor name="darkMatterParticle_" />
    !!]
    return
  end function decayingDarkMatterConstructorParameters

  function decayingDarkMatterConstructorInternal(massFunction_,cosmologyParameters_,darkMatterParticle_,massMinimum,massMaximum,pointsPerDecade) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`haloMassFunctionDecayingDarkMatter` halo mass function class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleDecayingDarkMatter
    use :: Error                , only : Error_Report                        , Warn
    implicit none
    type            (haloMassFunctionDecayingDarkMatter)                        :: self
    class           (haloMassFunctionClass             ), target, intent(in   ) :: massFunction_
    class           (cosmologyParametersClass          ), target, intent(in   ) :: cosmologyParameters_
    class           (darkMatterParticleClass           ), target, intent(in   ) :: darkMatterParticle_
    double precision                                            , intent(in   ) :: massMinimum         , massMaximum, &
         &                                                                         pointsPerDecade
    !![
    <constructorAssign variables="*massFunction_, *cosmologyParameters_, *darkMatterParticle_, massMinimum, massMaximum, pointsPerDecade"/>
    !!]

    ! Extract the decay lifetime and velocity kick from the (decaying dark matter) particle. Note that we
    ! must select on our own (assigned) pointer to the particle, as the accessor methods require an
    ! `intent(inout)` object.
    select type (particle_ => self%darkMatterParticle_)
    class is (darkMatterParticleDecayingDarkMatter)
       self%lifetime    =particle_%lifetime    ()
       self%velocityKick=particle_%velocityKick()
    class default
       call Error_Report('a decaying dark matter particle ([darkMatterParticleDecayingDarkMatter]) is required'//{introspection:location})
    end select
    ! Warn if the wrapped mass function will not respond to the (mass-dependent) DDM critical overdensity,
    ! in which case the DDM suppression---which is carried entirely by δ_c(M₀)---would be lost and
    ! only the mass remapping applied.
    if (.not.self%massFunction_%isCriticalOverdensityDependent())                                                                            &
         & call Warn(                                                                                                                        &
         &           'haloMassFunctionDecayingDarkMatter: the wrapped halo mass function does not depend on the critical overdensity for' // &
         &           ' collapse, so the decaying dark matter suppression (which is carried by a mass-dependent critical overdensity) will'// &
         &           ' be lost and only the mass remapping applied. Use an f(ν)-type mass function (e.g. shethTormen) configured with a'  // &
         &           ' [criticalOverdensityDecayingDarkMatter] critical overdensity.'                                                        &
         &          )
    ! The mapping between Lagrangian and collapsed mass is tabulated lazily (and re-tabulated as the epoch
    ! changes or the required mass range grows) in decayingDarkMatterTabulate.
    return
  end function decayingDarkMatterConstructorInternal

  subroutine decayingDarkMatterDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`haloMassFunctionDecayingDarkMatter` halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionDecayingDarkMatter), intent(inout) :: self

    !![
    <objectDestructor name="self%massFunction_"       />
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%darkMatterParticle_" />
    !!]
    return
  end subroutine decayingDarkMatterDestructor

  subroutine decayingDarkMatterTabulate(self,time,mass)
    !!{RST
    (Re)build the tabulated mapping between Lagrangian mass :math:`M_0` and collapsed mass :math:`M_\mathrm{coll}`
    at the given cosmic ``time``. The Lagrangian-mass grid spans ``self%massMinimum`` to ``self%massMaximum``
    at fixed resolution ``self%pointsPerDecade``; if the (optional) requested collapsed ``mass`` falls
    outside the resulting range these limits are first expanded, by whole grid steps so that the
    resolution is preserved. The mapping is monotonic, so the tabulated collapsed masses are strictly
    increasing and can be used to build an interpolator for the inverse mapping.
    !!}
    use :: Decaying_Dark_Matter_Spherical_Collapse, only : decayingDarkMatterMassCollapsed
    use :: Numerical_Interpolation                , only : GSL_Interp_Linear
    use :: Numerical_Ranges                       , only : Make_Range                     , rangeTypeLogarithmic
    use :: Table_Labels                           , only : extrapolationTypeExtrapolate
    implicit none
    class           (haloMassFunctionDecayingDarkMatter), intent(inout)           :: self
    double precision                                    , intent(in   )           :: time
    double precision                                    , intent(in   ), optional :: mass
    integer                                                                       :: i          , countTable
    double precision                                                              :: massMinimum, massMaximum, &
         &                                                                           stepFactor

    massMinimum=self%massMinimum
    massMaximum=self%massMaximum
    ! Expand the Lagrangian-mass range (by whole grid steps, preserving the resolution) until the
    ! corresponding collapsed-mass range brackets any requested collapsed mass. The mapping is monotonic
    ! increasing with M_coll ≤ M₀, so the low and high limits can be tested independently and expansion is
    ! guaranteed to terminate.
    if (present(mass)) then
       stepFactor=10.0d0**(1.0d0/self%pointsPerDecade)
       do while (decayingDarkMatterMassCollapsed(massMinimum,time,self%lifetime,self%velocityKick) > mass)
          massMinimum=massMinimum/stepFactor
       end do
       do while (decayingDarkMatterMassCollapsed(massMaximum,time,self%lifetime,self%velocityKick) < mass)
          massMaximum=massMaximum*stepFactor
       end do
    end if
    self%massMinimum=massMinimum
    self%massMaximum=massMaximum
    ! Build the grid of Lagrangian masses at the fixed resolution.
    countTable=int(log10(massMaximum/massMinimum)*self%pointsPerDecade)+1
    if (allocated(self%logMass0Table        )) deallocate(self%logMass0Table        )
    if (allocated(self%logMassCollapsedTable)) deallocate(self%logMassCollapsedTable)
    allocate(self%logMass0Table        (countTable))
    allocate(self%logMassCollapsedTable(countTable))
    self%logMass0Table=log(Make_Range(massMinimum,massMaximum,countTable,rangeTypeLogarithmic))
    do i=1,countTable
       self%logMassCollapsedTable(i)=+log(decayingDarkMatterMassCollapsed(exp(self%logMass0Table(i)),time,self%lifetime,self%velocityKick))
    end do
    ! Build an interpolator for the inverse mapping, log(M₀) as a function of log(M_coll).
    self%interpolator_       =interpolator(self%logMassCollapsedTable,self%logMass0Table,interpolationType=GSL_Interp_Linear,extrapolationType=extrapolationTypeExtrapolate)
    self%massCollapsedMinimum=exp(self%logMassCollapsedTable(         1))
    self%massCollapsedMaximum=exp(self%logMassCollapsedTable(countTable))
    self%timeTabulated       =time
    self%tabulated           =.true.
    return
  end subroutine decayingDarkMatterTabulate

  double precision function decayingDarkMatterDifferential(self,time,mass,node)
    !!{RST
    Return the differential halo mass function at the given time and (collapsed) mass. The wrapped mass
    function is evaluated at the corresponding Lagrangian mass and multiplied by the Jacobian of the
    mass mapping.
    !!}
    implicit none
    class           (haloMassFunctionDecayingDarkMatter), intent(inout), target   :: self
    double precision                                    , intent(in   )           :: time            , mass
    type            (treeNode                          ), intent(inout), optional :: node
    double precision                                                              :: logMassCollapsed, mass0   , &
         &                                                                           gradientLog     , jacobian

    ! (Re)tabulate the mass mapping on demand: when the epoch changes, or when the requested collapsed
    ! mass falls outside the currently-tabulated range (which extends the range and rebuilds the table).
    ! No locking is required as the halo mass function task deep-copies this object per OpenMP thread, so
    ! the tabulated state mutated here is thread-private.
    if (.not.self%tabulated .or. time /= self%timeTabulated .or. mass < self%massCollapsedMinimum .or. mass > self%massCollapsedMaximum) &
         & call decayingDarkMatterTabulate(self,time,mass)
    ! Invert the mapping: find the Lagrangian mass M₀ corresponding to the collapsed mass, and the
    ! logarithmic gradient d ln M₀ / d ln M_coll.
    logMassCollapsed=+log(mass)
    mass0           =+exp(self%interpolator_%interpolate(logMassCollapsed))
    gradientLog     =+    self%interpolator_%derivative (logMassCollapsed)
    ! Jacobian of the change of variables, dM₀/dM_coll = (M₀/M_coll) (d ln M₀ / d ln M_coll)
    ! [Montandon et al. (2026), their eq. 8].
    jacobian        =+mass0       &
         &           /mass        &
         &           *gradientLog
    decayingDarkMatterDifferential=+self%massFunction_%differential(time,mass0,node=node) &
         &                         *jacobian
    return
  end function decayingDarkMatterDifferential

  logical function decayingDarkMatterIsCriticalOverdensityDependent(self)
    !!{RST
    Return whether the differential halo mass function depends on the critical overdensity for
    collapse, by forwarding the query to the wrapped halo mass function.
    !!}
    implicit none
    class(haloMassFunctionDecayingDarkMatter), intent(inout) :: self

    decayingDarkMatterIsCriticalOverdensityDependent=self%massFunction_%isCriticalOverdensityDependent()
    return
  end function decayingDarkMatterIsCriticalOverdensityDependent
