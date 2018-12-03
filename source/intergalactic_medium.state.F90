!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which provides a class for calculations of the intergalactic medium thermal and ionization state.

module Intergalactic_Medium_State
  !% Provides a class for calculations of the intergalactic medium thermal and ionization state.
  use Cosmology_Parameters, only : cosmologyParameters, cosmologyParametersClass
  use Cosmology_Functions , only : cosmologyFunctions , cosmologyFunctionsClass
  use Linear_Growth       , only : linearGrowth       , linearGrowthClass
  use Tables
  private
  
  !# <functionClass>
  !#  <name>intergalacticMediumState</name>
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
  !#   <description>Return the electron scattering optical depth from the present day back to the given {\normalfont \ttfamily time} in the \gls{igm}.</description>
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
  !#    if (time > self%electronScatteringTableTimeMaximum)                                &amp;
  !#       &amp; call Galacticus_Error_Report(                                             &amp;
  !#       &amp;                              'time exceeds present age of the universe'// &amp;
  !#       &amp;                              {introspection:location}                     &amp;
  !#       &amp;                             )
  !#    assumeFullyIonizedActual=.false.
  !#    if (present(assumeFullyIonized)) assumeFullyIonizedActual=assumeFullyIonized
  !#    if (assumeFullyIonizedActual) then
  !#       intergalacticMediumStateElectronScatteringOpticalDepth=-self%electronScatteringFullyIonized%interpolate(time)
  !#    else
  !#       intergalacticMediumStateElectronScatteringOpticalDepth=-self%electronScattering            %interpolate(time)
  !#    end if
  !#   </code>
  !#  </method>
  !#  <method name="electronScatteringTime" >
  !#   <description>Return the cosmological time at which the given electron scattering {\normalfont \ttfamily opticalDepth} is reached (integrating from the present day) in the \gls{igm}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <modules>Galacticus_Error</modules>
  !#   <argument>double precision, intent(in   )           :: opticalDepth</argument>
  !#   <argument>logical         , intent(in   ), optional :: assumeFullyIonized</argument>
  !#   <code>
  !#    logical                                            :: assumeFullyIonizedActual
  !#    double precision                                   :: time
  !#    ! Check for invalid input.
  !#    if (opticalDepth &lt; 0.0d0) call Galacticus_Error_Report('optical depth must be non-negative'//{introspection:location})
  !#    ! Determine which optical depth to use.
  !#    assumeFullyIonizedActual=.false.
  !#    if (present(assumeFullyIonized)) assumeFullyIonizedActual=assumeFullyIonized
  !#    ! Ensure that the table is initialized.
  !#    time=self%cosmologyFunctions_%cosmicTime(1.0d0)
  !#    call intergalacticMediumStateElectronScatteringTabulate(self,time)
  !#    do while (                                                                                                    &amp;
  !#         &amp;     (.not.assumeFullyIonizedActual .and. self%electronScattering            %y(1) > -opticalDepth) &amp;
  !#         &amp;    .or.                                                                                            &amp;
  !#         &amp;     (     assumeFullyIonizedActual .and. self%electronScatteringFullyIonized%y(1) > -opticalDepth) &amp;
  !#         &amp;   )
  !#       time=time/2.0d0
  !#       call intergalacticMediumStateElectronScatteringTabulate(self,time)
  !#    end do
  !#    if (assumeFullyIonizedActual) then
  !#       intergalacticMediumStateElectronScatteringTime=self%electronScatteringFullyIonizedInverse%interpolate(-opticalDepth)
  !#    else
  !#       intergalacticMediumStateElectronScatteringTime=self%electronScatteringInverse            %interpolate(-opticalDepth)
  !#    end if
  !#   </code>
  !#  </method>
  !#  <method name="filteringMass" >
  !#   <description>Return the filtering mass at the given {\normalfont \ttfamily time}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <modules>Galacticus_Error</modules>
  !#   <argument>double precision, intent(in   ) :: time</argument>
  !#   <code>
  !#    ! Ensure that the table is initialized.
  !#    call intergalacticMediumStateFilteringMassTabulate(self,time)
  !#    intergalacticMediumStateFilteringMass=self%filteringMassTable%interpolate(time)
  !#   </code>
  !#  </method>
  !#  <!-- Objects. -->
  !#  <data scope="self"  >class           (cosmologyParametersClass), pointer     :: cosmologyParameters_                   => null()                                  </data>
  !#  <data scope="self"  >class           (cosmologyFunctionsClass ), pointer     :: cosmologyFunctions_                    => null()                                  </data>
  !#  <data scope="self"  >class           (linearGrowthClass       ), pointer     :: linearGrowth_                          => null()                                  </data>
  !#  <!-- Filtering mass table. -->
  !#  <data scope="module">integer                                   , parameter   :: filteringMassTablePointsPerDecade      =  100                                     </data>
  !#  <data scope="self"  >logical                                                 :: filteringMassTableInitialized          =  .false.                                 </data>
  !#  <data scope="self"  >integer                                                 :: filteringMassTableNumberPoints                                                    </data>
  !#  <data scope="self"  >double precision                                        :: filteringMassTableTimeMaximum                   , filteringMassTableTimeMinimum   </data>
  !#  <data scope="self"  >type            (table1DLogarithmicLinear)              :: filteringMassTable                                                                </data>
  !#  <!-- Electron scattering optical depth tables. -->
  !#  <data scope="module">integer                                   , parameter   :: electronScatteringTablePointsPerDecade =  100                                     </data>
  !#  <data scope="self"  >logical                                                 :: electronScatteringTableInitialized     =  .false.                                 </data>
  !#  <data scope="self"  >integer                                                 :: electronScatteringTableNumberPoints                                               </data>
  !#  <data scope="self"  >double precision                                        :: electronScatteringTableTimeMaximum            , electronScatteringTableTimeMinimum</data>
  !#  <data scope="self"  >type            (table1DLogarithmicLinear)              :: electronScattering                            , electronScatteringFullyIonized    </data>
  !#  <data scope="self"  >class           (table1D                 ), allocatable :: electronScatteringFullyIonizedInverse         , electronScatteringInverse         </data>
  !# </functionClass>

contains
  
  subroutine intergalacticMediumStateElectronScatteringTabulate(self,time)
    !% Construct a table of electron scattering optical depth as a function of cosmological time.
    use Numerical_Integration
    use FGSL                 , only : fgsl_function, fgsl_integration_workspace
    implicit none
    class           (intergalacticMediumStateClass), intent(inout), target :: self
    double precision                               , intent(in   )         :: time
    type            (fgsl_function                )                        :: integrandFunction
    type            (fgsl_integration_workspace   )                        :: integrationWorkspace
    integer                                                                :: iTime
    logical                                                                :: fullyIonized

    if (.not.self%electronScatteringTableInitialized.or.time < self%electronScatteringTableTimeMinimum) then
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
       do iTime=1,self%electronScatteringTableNumberPoints-1
          fullyIonized=.false.
          call self%electronScattering%populate(                                                             &
               &                      -Integrate(                                                            &
               &                                 self%electronScattering%x(iTime)                          , &
               &                                 self%electronScatteringTableTimeMaximum                   , &
               &                                 intergalacticMediumStateElectronScatteringIntegrand       , &
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
    return
    
  contains

    double precision function intergalacticMediumStateElectronScatteringIntegrand(time)
      !% Integrand for electron scattering optical depth calculations.
      use Numerical_Constants_Physical
      use Numerical_Constants_Astronomical
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

  subroutine intergalacticMediumStateFilteringMassTabulate(self,time)
    !% Construct a table of filtering mass as a function of cosmological time.
    use FODEIV2
    use ODEIV2_Solver
    use Galacticus_Error
    use Numerical_Constants_Math
    use Intergalactic_Medium_Filtering_Masses
    implicit none
    class           (intergalacticMediumStateClass), intent(inout), target :: self
    double precision                               , intent(in   )         :: time
    double precision                               , parameter             :: redshiftMaximumNaozBarkana=150.0d0 ! Maximum redshift at which fitting function of Naoz & Barkana is valid.
    double precision                               , dimension(3)          :: massFiltering                     , massFilteringScales
    double precision                               , parameter             :: odeToleranceAbsolute      =1.0d-03, odeToleranceRelative      =1.0d-03
    type            (fodeiv2_system               )                        :: ode2System
    type            (fodeiv2_driver               )                        :: ode2Driver
    logical                                                                :: odeReset
    integer                                                                :: iTime
    double precision                                                       :: timeInitial                       , timeCurrent

    if (.not.self%filteringMassTableInitialized .or. time < self%filteringMassTableTimeMinimum) then
       ! Find minimum and maximum times to tabulate.
       self%filteringMassTableTimeMaximum=max(self%cosmologyFunctions_%cosmicTime(1.0d0),time      )
       self%filteringMassTableTimeMinimum=min(self%cosmologyFunctions_%cosmicTime(1.0d0),time/2.0d0)
       ! Decide how many points to tabulate and allocate table arrays.
       self%filteringMassTableNumberPoints=int(log10(self%filteringMassTableTimeMaximum/self%filteringMassTableTimeMinimum)&
            & *dble(filteringMassTablePointsPerDecade))+1
       ! Create the tables.
       call self%filteringMassTable%destroy()
       call self%filteringMassTable%create (                                     &
            &                               self%filteringMassTableTimeMinimum , &
            &                               self%filteringMassTableTimeMaximum , &
            &                               self%filteringMassTableNumberPoints  &
            &                              )
       ! Evaluate a suitable starting time for filtering mass calculations.
       timeInitial=self%cosmologyFunctions_ %cosmicTime                 (                            &
            &       self%cosmologyFunctions_%expansionFactorFromRedshift (                           &
            &                                                             redshiftMaximumNaozBarkana &
            &                                                            )                           &
            &                                                           )
       ! Loop over times and populate tables.
       do iTime=1,self%filteringMassTableNumberPoints
          ! Abort if time is too early.
          if (self%filteringMassTable%x(iTime) <= timeInitial) call Galacticus_Error_Report('time is too early'//{introspection:location})
          ! Set the composite variables used to solve for filtering mass.
          call Mass_Filtering_ODE_Initial_Conditions(timeInitial,self%cosmologyParameters_,self%cosmologyFunctions_,self%linearGrowth_,massFiltering,massFilteringScales)
          ! Solve the ODE system
          odeReset   =.true.
          timeCurrent=timeInitial
          call ODEIV2_Solve(                                                     &
               &            ode2Driver                                         , &
               &            ode2System                                         , &
               &            timeCurrent                                        , &
               &            self%filteringMassTable%x(iTime)                   , &
               &            3                                                  , &
               &            massFiltering                                      , &
               &            Intergalactic_Medium_State_ODEs                    , &
               &            odeToleranceAbsolute                               , &
               &            odeToleranceRelative                               , &
               &            yScale                         =massFilteringScales, &
               &            reset                          =odeReset             &
               &           )
          call self%filteringMassTable%populate(massFiltering(3),iTime)
       end do
       ! Specify that tabulation has been made.
       self%filteringMassTableInitialized=.true.
    end if
    return

  contains

    integer function Intergalactic_Medium_State_ODEs(time,properties,propertiesRateOfChange)
      !% Evaluates the ODEs controlling the evolution temperature.
      use ODE_Solver_Error_Codes
      use Numerical_Constants_Astronomical
      use Numerical_Constants_Atomic
      use Intergalactic_Medium_Filtering_Masses
      use FGSL                                 , only : FGSL_Success
      implicit none
      double precision, intent(in  )                :: time
      double precision, intent(in   ), dimension(:) :: properties
      double precision, intent(  out), dimension(:) :: propertiesRateOfChange
      double precision                              :: temperature            , massParticleMean

       ! Find mean particle mass.
      massParticleMean=+(hydrogenByMassPrimordial*(1.0d0+self%electronFraction(time)*electronMass/massHydrogenAtom)                 +heliumByMassPrimordial               ) &
           &           /(hydrogenByMassPrimordial*(1.0d0+self%electronFraction(time)                              )/massHydrogenAtom+heliumByMassPrimordial/massHeliumAtom)
      ! Get the temperature.
      temperature=self%temperature(time)
      ! Compute the rates of change for the ODE system.
      propertiesRateOfChange=Mass_Filtering_ODE_System(self%cosmologyParameters_,self%cosmologyFunctions_,self%linearGrowth_,time,massParticleMean,temperature,properties)
      ! Return success.
      Intergalactic_Medium_State_ODEs=FGSL_Success
      return
    end function Intergalactic_Medium_State_ODEs

  end subroutine intergalacticMediumStateFilteringMassTabulate

end module Intergalactic_Medium_State
