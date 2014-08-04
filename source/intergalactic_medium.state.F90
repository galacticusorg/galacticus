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

!% Contains a module which provides a class for calculations of the intergalactic medium thermal and ionization state.

module Intergalactic_Medium_State
  !% Provides a class for calculations of the intergalactic medium thermal and ionization state.
  use, intrinsic :: ISO_C_Binding
  use               ISO_Varying_String
  use               Tables
  !# <include directive="intergalacticMediumState" type="functionModules" >
  include 'intergalacticMediumState.functionModules.inc'
  !# </include>
  private

  !# <include directive="intergalacticMediumState" type="function" >
  !#  <descriptiveName>Intergalactic Medium State</descriptiveName>
  !#  <description>Class providing intergalactic medium state.</description>
  !#  <default>recFast</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <method name="electronFraction" >
  !#   <description>Return the electron fraction (relative to hydrogen) in the \gls{igm} at the given time.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: time</argument>
  !#  </method>
  !#  <method name="neutralHydrogenFraction" >
  !#   <description>Return the neutral fraction of hydrogen in the \gls{igm} at the given time.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: time</argument>
  !#  </method>
  !#  <method name="singlyIonizedHydrogenFraction" >
  !#   <description>Return the singly-ionized fraction of hydrogen in the \gls{igm} at the given time.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   )           :: time</argument>
  !#   <code>intergalacticMediumStateSinglyIonizedHydrogenFraction=1.0d0-self%neutralHydrogenFraction(time)</code>
  !#  </method>
  !#  <method name="neutralHeliumFraction" >
  !#   <description>Return the neutral fraction of helium in the \gls{igm} at the given time.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: time</argument>
  !#  </method>
  !#  <method name="singlyIonizedHeliumFraction" >
  !#   <description>Return the singly-ionized fraction of helium in the \gls{igm} at the given time.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: time</argument>
  !#  </method>
  !#  <method name="doublyIonizedHeliumFraction" >
  !#   <description>Return the doubly-ionized fraction of helium in the \gls{igm} at the given time.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   )           :: time</argument>
  !#   <code>intergalacticMediumStateDoublyIonizedHeliumFraction=1.0d0-self%singlyIonizedHeliumFraction(time)-self%neutralHeliumFraction(time)</code>
  !#  </method>
  !#  <method name="temperature" >
  !#   <description>Return the temperature (in Kelvin) of the \gls{igm} at the given time.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: time</argument>
  !#  </method>
  !#  <method name="electronScatteringOpticalDepth" >
  !#   <description>Return the electron scattering optical depth from the present day back to the given {\tt time} in the \gls{igm}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <modules>Galacticus_Error</modules>
  !#   <argument>double precision, intent(in   )           :: time</argument>
  !#   <argument>logical         , intent(in   ), optional :: assumeFullyIonized</argument>
  !#   <code>
  !#    logical :: assumeFullyIonizedActual
  !#    ! Ensure that the table is initialized.
  !#    call intergalacticMediumStateElectronScatteringTabulate(self,time)
  !#    ! Check for invalid input.
  !#    if (time > self%electronScatteringTableTimeMaximum)                               &amp;
  !#       &amp; call Galacticus_Error_Report(                                            &amp;
  !#       &amp;                              'electronScatteringOpticalDepth'          , &amp;
  !#       &amp;                              'time exceeds present age of the universe'  &amp;
  !#       &amp;                             )
  !#    !$omp critical     (igmStateElectronScatteringInterpolation)
  !#    assumeFullyIonizedActual=.false.
  !#    if (present(assumeFullyIonized)) assumeFullyIonizedActual=assumeFullyIonized
  !#    if (assumeFullyIonizedActual) then
  !#       intergalacticMediumStateElectronScatteringOpticalDepth=-self%electronScatteringFullyIonized%interpolate(time)
  !#    else
  !#       intergalacticMediumStateElectronScatteringOpticalDepth=-self%electronScattering            %interpolate(time)
  !#    end if
  !#    !$omp end critical (igmStateElectronScatteringInterpolation)
  !#   </code>
  !#  </method>
  !#  <method name="electronScatteringTime" >
  !#   <description>Return the cosmological time at which the given electron scattering {\tt opticalDepth} is reached (integrating from the present day) in the \gls{igm}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <modules>Galacticus_Error Cosmology_Functions</modules>
  !#   <argument>double precision, intent(in   )           :: opticalDepth</argument>
  !#   <argument>logical         , intent(in   ), optional :: assumeFullyIonized</argument>
  !#   <code>
  !#    class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_
  !#    logical                                            :: assumeFullyIonizedActual
  !#    double precision                                   :: time
  !#    ! Check for invalid input.
  !#    if (opticalDepth &lt; 0.0d0) call Galacticus_Error_Report('electronScatteringTime','optical depth must be non-negative')
  !#    ! Determine which optical depth to use.
  !#    assumeFullyIonizedActual=.false.
  !#    if (present(assumeFullyIonized)) assumeFullyIonizedActual=assumeFullyIonized
  !#    ! Get the default cosmology functions object.
  !#    cosmologyFunctions_ => cosmologyFunctions()
  !#    ! Ensure that the table is initialized.
  !#    time=cosmologyFunctions_%cosmicTime(1.0d0)
  !#    call intergalacticMediumStateElectronScatteringTabulate(self,time)
  !#    do while (                                                                                                    &amp;
  !#         &amp;     (.not.assumeFullyIonizedActual .and. self%electronScattering            %y(1) > -opticalDepth) &amp;
  !#         &amp;    .or.                                                                                            &amp;
  !#         &amp;     (     assumeFullyIonizedActual .and. self%electronScatteringFullyIonized%y(1) > -opticalDepth) &amp;
  !#         &amp;   )
  !#       time=time/2.0d0
  !#       call intergalacticMediumStateElectronScatteringTabulate(self,time)
  !#    end do
  !#    !$omp critical     (igmStateElectronScatteringInterpolation)
  !#    if (assumeFullyIonizedActual) then
  !#       intergalacticMediumStateElectronScatteringTime=self%electronScatteringFullyIonizedInverse%interpolate(-opticalDepth)
  !#    else
  !#       intergalacticMediumStateElectronScatteringTime=self%electronScatteringInverse            %interpolate(-opticalDepth)
  !#    end if
  !#    !$omp end critical (igmStateElectronScatteringInterpolation)
  !#   </code>
  !#  </method>
  !#  <!-- Electron scattering optical depth tables. -->
  !#  <data scope="module">integer                                        , parameter   :: electronScatteringTablePointsPerDecade=100</data>
  !#  <data scope="self"  >logical                                                      :: electronScatteringTableInitialized    =.false.</data>
  !#  <data scope="self"  >integer                                                      :: electronScatteringTableNumberPoints</data>
  !#  <data scope="self"  >double precision                                             :: electronScatteringTableTimeMaximum            , electronScatteringTableTimeMinimum</data>
  !#  <data scope="self"  >type            (table1DLogarithmicLinear)                   :: electronScattering                            , electronScatteringFullyIonized</data>
  !#  <data scope="self"  >class           (table1D                      ), allocatable :: electronScatteringFullyIonizedInverse         , electronScatteringInverse</data>
  !#  <!-- Option controlling whether electron scattering optical depth calculations should assume a fully ionized universe. -->
  !#  <data scope="module" threadprivate="yes">class (intergalacticMediumStateClass), pointer :: selfGlobal</data>
  !#  <data scope="module" threadprivate="no" >logical                                        :: fullyIonized</data>
  include 'intergalacticMediumState.type.inc'
  !# </include>

  subroutine intergalacticMediumStateElectronScatteringTabulate(self,time)
    !% Construct a table of electron scattering optical depth as a function of cosmological time.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    use Cosmology_Functions
    use FGSL
    implicit none
    class           (intergalacticMediumStateClass), intent(inout), target :: self
    double precision                               , intent(in   )         :: time
    class           (cosmologyFunctionsClass      ), pointer               :: cosmologyFunctions_
    type            (c_ptr                        )                        :: parameterPointer
    type            (fgsl_function                )                        :: integrandFunction
    type            (fgsl_integration_workspace   )                        :: integrationWorkspace
    integer                                                                :: iTime

    !$omp critical (igmStateElectronScatteringInterpolation)
    if (.not.self%electronScatteringTableInitialized.or.time < self%electronScatteringTableTimeMinimum) then
       ! Set module-scope pointer to self.
       selfGlobal => self
       ! Get the default cosmology functions object.
       cosmologyFunctions_ => cosmologyFunctions()
       ! Find minimum and maximum times to tabulate.
       self%electronScatteringTableTimeMaximum=    cosmologyFunctions_%cosmicTime(1.0d0)
       self%electronScatteringTableTimeMinimum=min(cosmologyFunctions_%cosmicTime(1.0d0),time)/2.0d0
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
       do iTime=1,self%electronScatteringTableNumberPoints-1
          fullyIonized=.false.
          call self%electronScattering%populate(                                                             &
               &                      -Integrate(                                                            &
               &                                 self%electronScattering%x(iTime)                          , &
               &                                 self%electronScatteringTableTimeMaximum                   , &
               &                                 intergalacticMediumStateElectronScatteringIntegrand       , &
               &                                 parameterPointer                                          , &
               &                                 integrandFunction                                         , &
               &                                 integrationWorkspace                                      , &
               &                                 toleranceAbsolute                                  =0.0d+0, &
               &                                 toleranceRelative                                  =1.0d-3  &
               &                                )                                                          , &
               &                       iTime                                                                 &
               &                      )
          call Integrate_Done(integrandFunction,integrationWorkspace)
          fullyIonized=.true.
          call self%electronScatteringFullyIonized%populate(                                                             &
               &                                  -Integrate(                                                            &
               &                                             self%electronScatteringFullyIonized%x(iTime)              , &
               &                                             self%electronScatteringTableTimeMaximum                   , &
               &                                             intergalacticMediumStateElectronScatteringIntegrand       , &
               &                                             parameterPointer                                          , &
               &                                             integrandFunction                                         , &
               &                                             integrationWorkspace                                      , &
               &                                             toleranceAbsolute                                  =0.0d+0, &
               &                                             toleranceRelative                                  =1.0d-3  &
               &                                            )                                                          , &
               &                       iTime                                                                             &
               &                      )
          call Integrate_Done(integrandFunction,integrationWorkspace)
       end do
       call self%electronScattering            %populate(0.0d0,self%electronScatteringTableNumberPoints)
       call self%electronScatteringFullyIonized%populate(0.0d0,self%electronScatteringTableNumberPoints)
       call self%electronScattering            %reverse (self%electronScatteringInverse                )
       call self%electronScatteringFullyIonized%reverse (self%electronScatteringFullyIonizedInverse    )
       ! Specify that tabulation has been made.
       self%electronScatteringTableInitialized=.true.
    end if
    !$omp end critical (igmStateElectronScatteringInterpolation)
    return
  end subroutine intergalacticMediumStateElectronScatteringTabulate

  function intergalacticMediumStateElectronScatteringIntegrand(time,parameterPointer) bind(c)
    !% Integrand for electron scattering optical depth calculations.
    use, intrinsic :: ISO_C_Binding
    use Cosmology_Functions
    use Cosmology_Parameters
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    implicit none
    real            (kind=c_double           )          :: intergalacticMediumStateElectronScatteringIntegrand
    real            (kind=c_double           ), value   :: time
    type            (c_ptr                   ), value   :: parameterPointer
    class           (cosmologyParametersClass), pointer :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_
    double precision                                    :: electronFraction                                    , expansionFactor

    ! Get the default cosmology.
    cosmologyParameters_ => cosmologyParameters                (    )
    ! Get the default cosmology functions object.
    cosmologyFunctions_  => cosmologyFunctions                 (    )
    expansionFactor      =  cosmologyFunctions_%expansionFactor(time)
    if (fullyIonized) then
       electronFraction=+      hydrogenByMassPrimordial                                  /atomicMassHydrogen &
            &           +2.0d0*heliumByMassPrimordial                                    /atomicMassHelium
    else
       electronFraction=       hydrogenByMassPrimordial*selfGlobal%electronFraction(time)/atomicMassHydrogen
    end if
    intergalacticMediumStateElectronScatteringIntegrand &
         & =+speedLight                                 &
         &  *gigaYear                                   &
         &  *thomsonCrossSection                        &
         &  *massSolar                                  &
         &  /atomicMassUnit                             &
         &  /megaParsec         **3                     &
         &  /expansionFactor    **3                     &
         &  *cosmologyParameters_%OmegaBaryon    ()     &
         &  *cosmologyParameters_%densityCritical()     &
         &  *electronFraction
    return
  end function intergalacticMediumStateElectronScatteringIntegrand

end module Intergalactic_Medium_State
