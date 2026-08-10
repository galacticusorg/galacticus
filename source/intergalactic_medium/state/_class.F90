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
Contains a module which provides a class for calculations of the intergalactic medium thermal and ionization state.
!!}

module Intergalactic_Medium_State
  !!{RST
  Provides a class for calculations of the intergalactic medium thermal and ionization state.
  !!}
  use :: Cosmology_Functions , only : cosmologyFunctions , cosmologyFunctionsClass
  use :: Cosmology_Parameters, only : cosmologyParameters, cosmologyParametersClass
  use :: Numerical_Ranges    , only : rangeLattice
  use :: Tables              , only : table              , table1DGeneric          , table1DLogarithmicLinear
  private

  !![
  <functionClass docformat="rst">
   <name>intergalacticMediumState</name>
   <descriptiveName>Intergalactic Medium State</descriptiveName>
   <description>
   Class providing the thermal and ionization state of the :term:`IGM`---the hydrogen and helium neutral, singly-ionized, and doubly-ionized fractions, the electron fraction, the temperature, and the electron-scattering optical depth as functions of cosmic time. These quantities evolve through the epoch of reionization and affect the cooling rates, photo-ionization suppression, and the UV background modelled elsewhere in Galacticus. The instantaneous Jeans mass computed from the IGM temperature governs the filtering mass scale for baryon accretion onto low-mass halos.
   </description>
   <default>recFast</default>
   <method name="electronFraction" >
    <description>
    Return the electron fraction (relative to hydrogen) in the :term:`IGM` at the given time.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="neutralHydrogenFraction" >
    <description>
    Return the neutral fraction of hydrogen in the :term:`IGM` at the given time.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="singlyIonizedHydrogenFraction" >
    <description>
    Return the singly-ionized fraction of hydrogen in the :term:`IGM` at the given time.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: time</argument>
    <code>intergalacticMediumStateSinglyIonizedHydrogenFraction=1.0d0-self%neutralHydrogenFraction(time)</code>
   </method>
   <method name="neutralHeliumFraction" >
    <description>
    Return the neutral fraction of helium in the :term:`IGM` at the given time.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="singlyIonizedHeliumFraction" >
    <description>
    Return the singly-ionized fraction of helium in the :term:`IGM` at the given time.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="doublyIonizedHeliumFraction" >
    <description>
    Return the doubly-ionized fraction of helium in the :term:`IGM` at the given time.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: time</argument>
    <code>intergalacticMediumStateDoublyIonizedHeliumFraction=1.0d0-self%singlyIonizedHeliumFraction(time)-self%neutralHeliumFraction(time)</code>
   </method>
   <method name="metallicity" >
    <description>
    Return the metallicity (mass fraction) in the :term:`IGM` at the given time.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: time</argument>
    <code>intergalacticMediumStateMetallicity=0.0d0</code>
   </method>
   <method name="temperature" >
    <description>
    Return the temperature (in Kelvin) of the :term:`IGM` at the given time.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="massJeans" >
    <description>
    Return the instantaneous Jeans mass (in :math:`\mathrm{M}_\odot`) at the given time.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <modules>Numerical_Constants_Physical Numerical_Constants_Astronomical Numerical_Constants_Atomic Numerical_Constants_Prefixes</modules>
    <argument>double precision, intent(in   ) :: time</argument>
    <code>
     double precision :: massParticleMean, speedSound
     massParticleMean                 =+(hydrogenByMassPrimordial*(1.0d0+self%electronFraction(time)*electronMass/massHydrogenAtom)                 +heliumByMassPrimordial               ) &amp;
          &amp;                        /(hydrogenByMassPrimordial*(1.0d0+self%electronFraction(time)                              )/massHydrogenAtom+heliumByMassPrimordial/massHeliumAtom)
     speedSound                       =+sqrt(                          &amp;
          &amp;                              +boltzmannsConstant       &amp;
          &amp;                              *self%temperature  (time) &amp;
          &amp;                              /massParticleMean         &amp;
          &amp;                             )                          &amp;
          &amp;                        /kilo
     intergalacticMediumStateMassJeans=+4.0d0                                                       &amp;
          &amp;                        *Pi                                                          &amp;
          &amp;                        /3.0d0                                                       &amp;
          &amp;                        *        self%cosmologyFunctions_%matterDensityEpochal(time) &amp;
          &amp;                        *(                                                           &amp;
          &amp;                          +2.0d0                                                     &amp;
          &amp;                          *Pi                                                        &amp;
          &amp;                          *speedSound                                                &amp;
          &amp;                          /sqrt(                                                     &amp;
          &amp;                                +4.0d0                                               &amp;
          &amp;                                *Pi                                                  &amp;
          &amp;                                *gravitationalConstant_internal                      &amp;
          &amp;                                *self%cosmologyFunctions_%matterDensityEpochal(time) &amp;
          &amp;                               )                                                     &amp;
          &amp;                         )**3
    </code>
   </method>
   <method name="electronScatteringOpticalDepth" >
    <description>
    Return the electron scattering optical depth from the present day back to the given ``time`` in the :term:`IGM`.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <modules>Error</modules>
    <argument>double precision, intent(in   )           :: time</argument>
    <argument>logical         , intent(in   ), optional :: assumeFullyIonized</argument>
    <code>
     logical :: assumeFullyIonizedActual
     ! Ensure that the table is initialized.
     call intergalacticMediumStateElectronScatteringTabulate(self,time)
     ! Check for invalid input.
     if (time &gt; self%electronScatteringTableTimeMaximum)                  &amp;
        &amp; call Error_Report(                                             &amp;
        &amp;                   'time exceeds present age of the universe'// &amp;
        &amp;                   {introspection:location}                     &amp;
        &amp;                  )
     assumeFullyIonizedActual=.false.
     if (present(assumeFullyIonized)) assumeFullyIonizedActual=assumeFullyIonized
     ! The tabulations are made against time in units of the present age of the universe - see
     ! `intergalacticMediumStateElectronScatteringTabulate`.
     if (assumeFullyIonizedActual) then
        intergalacticMediumStateElectronScatteringOpticalDepth=-self%electronScatteringFullyIonized%interpolate(time/self%electronScatteringTableTimeMaximum)
     else
        intergalacticMediumStateElectronScatteringOpticalDepth=-self%electronScattering            %interpolate(time/self%electronScatteringTableTimeMaximum)
     end if
    </code>
   </method>
   <method name="electronScatteringTime" >
    <description>
    Return the cosmological time at which the given electron scattering ``opticalDepth`` is reached (integrating from the present day) in the :term:`IGM`.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <modules>Error</modules>
    <argument>double precision, intent(in   )           :: opticalDepth</argument>
    <argument>logical         , intent(in   ), optional :: assumeFullyIonized</argument>
    <code>
     logical                                            :: assumeFullyIonizedActual
     double precision                                   :: time
     ! Check for invalid input.
     if (opticalDepth &lt; 0.0d0) call Error_Report('optical depth must be non-negative'//{introspection:location})
     ! Determine which optical depth to use.
     assumeFullyIonizedActual=.false.
     if (present(assumeFullyIonized)) assumeFullyIonizedActual=assumeFullyIonized
     ! Ensure that the table is initialized.
     time=self%cosmologyFunctions_%cosmicTime(1.0d0)
     call intergalacticMediumStateElectronScatteringTabulate(self,time)
     do while (                                                                                                       &amp;
          &amp;     (.not.assumeFullyIonizedActual .and. self%electronScattering            %y(1) &gt; -opticalDepth) &amp;
          &amp;    .or.                                                                                               &amp;
          &amp;     (     assumeFullyIonizedActual .and. self%electronScatteringFullyIonized%y(1) &gt; -opticalDepth) &amp;
          &amp;   )
        time=time/2.0d0
        call intergalacticMediumStateElectronScatteringTabulate(self,time)
     end do
     if (assumeFullyIonizedActual) then
        intergalacticMediumStateElectronScatteringTime=self%electronScatteringFullyIonizedInverse%interpolate(-opticalDepth)
     else
        intergalacticMediumStateElectronScatteringTime=self%electronScatteringInverse            %interpolate(-opticalDepth)
     end if
    </code>
   </method>
   <!-- Objects. -->
   <data scope="self"  >class           (cosmologyParametersClass), pointer   :: cosmologyParameters_                   => null()                                  </data>
   <data scope="self"  >class           (cosmologyFunctionsClass ), pointer   :: cosmologyFunctions_                    => null()                                  </data>
   <!-- Electron scattering optical depth tables. -->
   <data scope="module">integer                                   , parameter :: electronScatteringTablePointsPerDecade =  100                                 </data>
   <data scope="self"  >logical                                               :: electronScatteringTableInitialized     =  .false.                             </data>
   <data scope="self"  >double precision                                      :: electronScatteringTableTimeMaximum                                            </data>
   <!-- Lattice to which the tabulated times are pinned. This is the source of truth for the extent of the tabulation. -->
   <data scope="self"  >type            (rangeLattice            )            :: electronScatteringTableLattice                                                </data>
   <data scope="self"  >type            (table1DLogarithmicLinear)            :: electronScattering                        , electronScatteringFullyIonized    </data>
   <data scope="self"  >type            (table1DGeneric          )            :: electronScatteringFullyIonizedInverse     , electronScatteringInverse         </data>
  </functionClass>
  !!]

contains

  subroutine intergalacticMediumStateElectronScatteringTabulate(self,time)
    !!{RST
    Construct a table of electron scattering optical depth as a function of cosmological time.

    The optical depths are tabulated against time *in units of the present age of the universe*, and that axis is pinned to an
    absolute lattice, so that the times at which the optical depth is evaluated---and therefore every value interpolated between
    them---depend only on which lattice points are spanned, and not on the sequence of times which happened to be requested. In
    those units the present day is exactly the lattice point of index zero, which is what allows the axis to be pinned while
    still ending exactly at the present day: that is the upper limit of every integral performed here, the epoch at which the
    optical depth is exactly zero, and the largest time for which this tabulation is ever asked. On an axis of absolute time no
    lattice point coincides with it, and pinning would necessarily either fall short of it or run beyond it.
    !!}
    use :: Error                , only : Error_Report
    use :: Numerical_Integration, only : integrator
    use :: Numerical_Ranges     , only : Range_Pinned                    , gridSchemePerDecade
    implicit none
    class           (intergalacticMediumStateClass), intent(inout), target      :: self
    double precision                               , intent(in   )              :: time
    double precision                               , dimension(:) , allocatable :: time_                     , opticalDepth              , &
         &                                                                         opticalDepthFullyIonized
    logical                                        , dimension(:) , allocatable :: isComputed                , isComputedFullyIonized
    type            (integrator                   )                             :: integrator_
    type            (rangeLattice                 )                             :: latticeTime
    integer                                                                     :: iTime                     , iTimeMonotonic            , &
         &                                                                         iTimeMonotonicFullyIonized
    double precision                                                            :: timeRelative
    logical                                                                     :: fullyIonized              , makeTable

    ! Discard any previous tabulation which has been marked as invalid - the reionization history on which it was built may have
    ! changed, in which case no value in it may be carried over into the new tabulation.
    if (.not.self%electronScatteringTableInitialized) then
       ! Validate cosmological parameters.
       if (self%cosmologyParameters_%OmegaBaryon() <= 0.0d0) call Error_Report('can not compute electron scattering optical depths in a universe with no baryons'//{introspection:location})
       call self%electronScattering            %destroy()
       call self%electronScatteringFullyIonized%destroy()
       self%electronScatteringTableLattice    =rangeLattice()
       self%electronScatteringTableTimeMaximum=self%cosmologyFunctions_%cosmicTime(1.0d0)
    end if
    ! Find the range of times to tabulate. The request is passed as the target and the range already tabulated is unioned in
    ! through `latticeCurrent`; folding the latter into the target instead would apply the safety margin to an already margined
    ! bound and so ratchet the range outward on every retabulation. The margin is a factor of two below the request, which is the
    ! margin applied before this tabulation was pinned. The bounds are anchored to half decades rather than to whole decades
    ! because this is a cosmic time axis on which a whole decade is a large fraction of the history of the universe, and each
    ! additional point costs an integration over the whole span to the present day.
    !
    ! Every tabulation runs to the present day, whatever was requested of it: that is the upper limit of every integral performed
    ! here, so a tabulation which stopped short of it would be a table of integrals to some other epoch. It is imposed through
    ! `rangeCurrent`, which is applied after the safety margin and so raises the upper bound without also lowering the lower one,
    ! while `limitMaximum` stops the margin carrying the range beyond the present day when the time requested is close to it.
    timeRelative=min(time,self%electronScatteringTableTimeMaximum)/self%electronScatteringTableTimeMaximum
    latticeTime =Range_Pinned(                                                             &
         &                                  [timeRelative]                               , &
         &                                   electronScatteringTablePointsPerDecade       , &
         &                                   gridSchemePerDecade                          , &
         &                   marginFactor  = 2.0d0                                        , &
         &                   rangeCurrent  =[1.0d0,1.0d0]                                 , &
         &                   limitMaximum  = 1.0d0                                        , &
         &                   anchorEvery   = electronScatteringTablePointsPerDecade/2     , &
         &                   latticeCurrent= self%electronScatteringTableLattice            &
         &                  )
    ! Decide whether the tabulation must be built or extended. The decision is taken from the pinned lattice rather than from the
    ! bounds directly: the safety margin is applied to the request, so testing the request against the bound would apply that
    ! margin only when the request happened to arrive while it was the earliest seen. The range reached would then depend on the
    ! order in which times were asked for - which is precisely the dependence that pinning the lattice exists to remove.
    makeTable=.not.self%electronScatteringTableInitialized
    if (.not.makeTable) makeTable=.not.self%electronScatteringTableLattice%covers(latticeTime)
    if (makeTable) then
       ! Extend the tabulations onto the new lattice, preserving the optical depths already computed. Each is the integral from
       ! its own time to the present day, neither limit of which moves when the range is extended, so a preserved value is
       ! precisely the value which would have been computed afresh.
       call self%electronScattering            %extend(latticeTime,isComputed            )
       call self%electronScatteringFullyIonized%extend(latticeTime,isComputedFullyIonized)
       ! Take the times to tabulate from the abscissae of the table itself, so that the limits of the integrals performed here
       ! are precisely the times at which interpolation will later be made.
       time_=self%electronScatteringTableTimeMaximum*self%electronScattering%xs()
       allocate(opticalDepth            (latticeTime%count))
       allocate(opticalDepthFullyIonized(latticeTime%count))
       ! Loop over times, carrying over the optical depths already computed and evaluating only those which are new. Also find
       ! where the optical depth last decreases: only the monotonic tail beyond that point can be inverted.
       integrator_               =integrator(intergalacticMediumStateElectronScatteringIntegrand,toleranceRelative=1.0d-6)
       iTimeMonotonic            =1
       iTimeMonotonicFullyIonized=1
       do iTime=1,latticeTime%count
          if      (isComputed            (iTime)        ) then
             opticalDepth            (iTime)=self%electronScattering            %y(iTime)
          else if (latticeTime%indexMinimum+iTime-1 == 0) then
             ! This point is the present day - the lattice point of index zero - at which the optical depth is zero by
             ! definition. Testing the lattice index rather than the position in the table means that a range which did not reach
             ! the present day would simply integrate every one of its points, rather than silently setting its last to zero.
             opticalDepth            (iTime)=0.0d0
          else
             fullyIonized                   =.false.
             opticalDepth            (iTime)=-integrator_%integrate(                                                  &
                  &                                                 time_                                    (iTime), &
                  &                                                 self%electronScatteringTableTimeMaximum           &
                  &                                                )
          end if
          if      (isComputedFullyIonized(iTime)        ) then
             opticalDepthFullyIonized(iTime)=self%electronScatteringFullyIonized%y(iTime)
          else if (latticeTime%indexMinimum+iTime-1 == 0) then
             opticalDepthFullyIonized(iTime)=0.0d0
          else
             fullyIonized                   =.true.
             opticalDepthFullyIonized(iTime)=-integrator_%integrate(                                                  &
                  &                                                 time_                                    (iTime), &
                  &                                                 self%electronScatteringTableTimeMaximum           &
                  &                                                )
          end if
          if (iTime > 1) then
             if (opticalDepth            (iTime) < opticalDepth            (iTime-1)) &
                  & iTimeMonotonic            =iTime
             if (opticalDepthFullyIonized(iTime) < opticalDepthFullyIonized(iTime-1)) &
                  & iTimeMonotonicFullyIonized=iTime
          end if
       end do
       call self%electronScattering                   %populate(opticalDepth            )
       call self%electronScatteringFullyIonized       %populate(opticalDepthFullyIonized)
       ! Build the inverse tables, which are used to find the time at which a given optical depth is reached. Their abscissae are
       ! optical depths, which lie on no lattice - but they are now a function of the lattice points spanned alone, as the forward
       ! tabulation is.
       call self%electronScatteringInverse            %destroy ()
       call self%electronScatteringFullyIonizedInverse%destroy ()
       call self%electronScatteringInverse            %create  (opticalDepth            (iTimeMonotonic            :latticeTime%count))
       call self%electronScatteringFullyIonizedInverse%create  (opticalDepthFullyIonized(iTimeMonotonicFullyIonized:latticeTime%count))
       call self%electronScatteringInverse            %populate(time_                   (iTimeMonotonic            :latticeTime%count))
       call self%electronScatteringFullyIonizedInverse%populate(time_                   (iTimeMonotonicFullyIonized:latticeTime%count))
       ! Specify that tabulation has been made.
       self%electronScatteringTableLattice    =latticeTime
       self%electronScatteringTableInitialized=.true.
    end if
    return

  contains

    double precision function intergalacticMediumStateElectronScatteringIntegrand(time)
      !!{RST
      Integrand for electron scattering optical depth calculations.
      !!}
      use :: Numerical_Constants_Astronomical, only : gigaYear        , heliumByMassPrimordial, hydrogenByMassPrimordial, massSolar, &
          &                                           megaParsec
      use :: Numerical_Constants_Atomic      , only : atomicMassHelium, atomicMassHydrogen    , atomicMassUnit
      use :: Numerical_Constants_Physical    , only : speedLight      , thomsonCrossSection
      implicit none
      double precision, intent(in   ) :: time
      double precision                :: electronFraction, expansionFactor

      expansionFactor=self%cosmologyFunctions_%expansionFactor(time)
      if (fullyIonized) then
         electronFraction=+      hydrogenByMassPrimordial                            /atomicMassHydrogen &
              &           +2.0d0*heliumByMassPrimordial                              /atomicMassHelium
      else
         electronFraction=+      hydrogenByMassPrimordial*self%electronFraction(time)/atomicMassHydrogen
      end if
      intergalacticMediumStateElectronScatteringIntegrand  &
           & =+speedLight                                  &
           &  *gigaYear                                    &
           &  *thomsonCrossSection                         &
           &  *massSolar                                   &
           &  /atomicMassUnit                              &
           &  /megaParsec         **3                      &
           &  /expansionFactor    **3                      &
           &  *self%cosmologyParameters_%OmegaBaryon    () &
           &  *self%cosmologyParameters_%densityCritical() &
           &  *electronFraction
      return
    end function intergalacticMediumStateElectronScatteringIntegrand

  end subroutine intergalacticMediumStateElectronScatteringTabulate

end module Intergalactic_Medium_State
