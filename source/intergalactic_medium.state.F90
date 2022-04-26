!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which provides a class for calculations of the intergalactic medium thermal and ionization state.
!!}

module Intergalactic_Medium_State
  !!{
  Provides a class for calculations of the intergalactic medium thermal and ionization state.
  !!}
  use :: Cosmology_Functions , only : cosmologyFunctions , cosmologyFunctionsClass
  use :: Cosmology_Parameters, only : cosmologyParameters, cosmologyParametersClass
  use :: Tables              , only : table              , table1D                 , table1DLogarithmicLinear
  private

  !![
  <functionClass>
   <name>intergalacticMediumState</name>
   <descriptiveName>Intergalactic Medium State</descriptiveName>
   <description>
    Class providing the thermal and ionization state of the intergalactic medium.
   </description>
   <default>recFast</default>
   <method name="electronFraction" >
    <description>Return the electron fraction (relative to hydrogen) in the \gls{igm} at the given time.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="neutralHydrogenFraction" >
    <description>Return the neutral fraction of hydrogen in the \gls{igm} at the given time.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="singlyIonizedHydrogenFraction" >
    <description>Return the singly-ionized fraction of hydrogen in the \gls{igm} at the given time.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: time</argument>
    <code>intergalacticMediumStateSinglyIonizedHydrogenFraction=1.0d0-self%neutralHydrogenFraction(time)</code>
   </method>
   <method name="neutralHeliumFraction" >
    <description>Return the neutral fraction of helium in the \gls{igm} at the given time.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="singlyIonizedHeliumFraction" >
    <description>Return the singly-ionized fraction of helium in the \gls{igm} at the given time.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="doublyIonizedHeliumFraction" >
    <description>Return the doubly-ionized fraction of helium in the \gls{igm} at the given time.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: time</argument>
    <code>intergalacticMediumStateDoublyIonizedHeliumFraction=1.0d0-self%singlyIonizedHeliumFraction(time)-self%neutralHeliumFraction(time)</code>
   </method>
   <method name="temperature" >
    <description>Return the temperature (in Kelvin) of the \gls{igm} at the given time.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="massJeans" >
    <description>Return the instantaneus Jeans mass (in $\mathrm{M}_\odot$) at the given time.</description>
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
          &amp;                                *gravitationalConstantGalacticus                     &amp;
          &amp;                                *self%cosmologyFunctions_%matterDensityEpochal(time) &amp;
          &amp;                               )                                                     &amp;
          &amp;                         )**3
    </code>
   </method>
   <method name="electronScatteringOpticalDepth" >
    <description>Return the electron scattering optical depth from the present day back to the given {\normalfont \ttfamily time} in the \gls{igm}.</description>
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
     if (assumeFullyIonizedActual) then
        intergalacticMediumStateElectronScatteringOpticalDepth=-self%electronScatteringFullyIonized%interpolate(time)
     else
        intergalacticMediumStateElectronScatteringOpticalDepth=-self%electronScattering            %interpolate(time)
     end if
    </code>
   </method>
   <method name="electronScatteringTime" >
    <description>Return the cosmological time at which the given electron scattering {\normalfont \ttfamily opticalDepth} is reached (integrating from the present day) in the \gls{igm}.</description>
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
   <data scope="self"  >class           (cosmologyParametersClass), pointer     :: cosmologyParameters_                   => null()                                  </data>
   <data scope="self"  >class           (cosmologyFunctionsClass ), pointer     :: cosmologyFunctions_                    => null()                                  </data>
   <!-- Electron scattering optical depth tables. -->
   <data scope="module">integer                                   , parameter   :: electronScatteringTablePointsPerDecade =  100                                     </data>
   <data scope="self"  >logical                                                 :: electronScatteringTableInitialized     =  .false.                                 </data>
   <data scope="self"  >integer                                                 :: electronScatteringTableNumberPoints                                               </data>
   <data scope="self"  >double precision                                        :: electronScatteringTableTimeMaximum            , electronScatteringTableTimeMinimum</data>
   <data scope="self"  >type            (table1DLogarithmicLinear)              :: electronScattering                            , electronScatteringFullyIonized    </data>
   <data scope="self"  >class           (table1D                 ), allocatable :: electronScatteringFullyIonizedInverse         , electronScatteringInverse         </data>
  </functionClass>
  !!]

contains

  subroutine intergalacticMediumStateElectronScatteringTabulate(self,time)
    !!{
    Construct a table of electron scattering optical depth as a function of cosmological time.
    !!}
    use :: Error                , only : Error_Report
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (intergalacticMediumStateClass), intent(inout), target :: self
    double precision                               , intent(in   )         :: time
    type            (integrator                   )                        :: integrator_
    integer                                                                :: iTime
    logical                                                                :: fullyIonized

    if (.not.self%electronScatteringTableInitialized.or.time < self%electronScatteringTableTimeMinimum) then
       ! Validate cosmological parameters.
       if (self%cosmologyParameters_%OmegaBaryon() <= 0.0d0) call Error_Report('can not compute electron scattering optical depths in a universe with no baryons'//{introspection:location})
       ! Find minimum and maximum times to tabulate.
       self%electronScatteringTableTimeMaximum=    self%cosmologyFunctions_%cosmicTime(1.0d0)
       self%electronScatteringTableTimeMinimum=min(self%cosmologyFunctions_%cosmicTime(1.0d0),time)/2.0d0
       ! Decide how many points to tabulate and allocate table arrays.
       self%electronScatteringTableNumberPoints=int(log10(self%electronScatteringTableTimeMaximum/self%electronScatteringTableTimeMinimum)&
            &*dble(electronScatteringTablePointsPerDecade))+1
       ! Create the tables.
       call self%electronScattering            %destroy()
       call self%electronScatteringFullyIonized%destroy()
       call self%electronScattering            %create (                                          &
            &                                           self%electronScatteringTableTimeMinimum , &
            &                                           self%electronScatteringTableTimeMaximum , &
            &                                           self%electronScatteringTableNumberPoints  &
            &                                          )
       call self%electronScatteringFullyIonized%create (                                          &
            &                                           self%electronScatteringTableTimeMinimum , &
            &                                           self%electronScatteringTableTimeMaximum , &
            &                                           self%electronScatteringTableNumberPoints  &
            &                                          )
       ! Loop over times and populate tables.
       integrator_=integrator(intergalacticMediumStateElectronScatteringIntegrand,toleranceRelative=1.0d-6)
       do iTime=1,self%electronScatteringTableNumberPoints-1
          fullyIonized=.false.
          call self%electronScattering            %populate(                                                                         &
               &                                            -integrator_%integrate(                                                  &
               &                                                                   self%electronScattering                %x(iTime), &
               &                                                                   self%electronScatteringTableTimeMaximum           &
               &                                                                  )                                                , &
               &                                                                                                             iTime   &
               &                                           )
          fullyIonized=.true.
          call self%electronScatteringFullyIonized%populate(                                                                         &
               &                                            -integrator_%integrate(                                                  &
               &                                                                   self%electronScatteringFullyIonized    %x(iTime), &
               &                                                                   self%electronScatteringTableTimeMaximum           &
               &                                                                  )                                                , &
               &                                                                                                             iTime   &
               &                                           )
       end do
       call self%electronScattering            %populate(0.0d0,self%electronScatteringTableNumberPoints)
       call self%electronScatteringFullyIonized%populate(0.0d0,self%electronScatteringTableNumberPoints)
       call self%electronScattering            %reverse (self%electronScatteringInverse                )
       call self%electronScatteringFullyIonized%reverse (self%electronScatteringFullyIonizedInverse    )
       ! Specify that tabulation has been made.
       self%electronScatteringTableInitialized=.true.
    end if
    return

  contains

    double precision function intergalacticMediumStateElectronScatteringIntegrand(time)
      !!{
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
