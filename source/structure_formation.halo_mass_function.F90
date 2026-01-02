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
Contains a module which provides a class that implements halo mass functions.
!!}

module Halo_Mass_Functions
  use :: Cosmology_Parameters, only : cosmologyParameters, cosmologyParametersClass
  use :: Galacticus_Nodes    , only : treeNode
  implicit none
  private

  !![
  <functionClass>
   <name>haloMassFunction</name>
   <descriptiveName>Halo Mass Function</descriptiveName>
   <description>Class providing halo mass functions.</description>
   <default>tinker2008</default>
   <data>double precision                                    :: time_                         </data>
   <data>class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()</data>
   <data>type            (treeNode                ), pointer :: node                 => null()</data>
   <data>
    <scope>module</scope>
    <threadprivate>yes</threadprivate>
    <content>class(haloMassFunctionClass), pointer :: globalSelf</content>
   </data>
   <method name="differential" >
    <description>Return the differential halo mass function for {\normalfont \ttfamily mass} [$M_\odot$] at {\normalfont \ttfamily time} [Gyr].</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision          , intent(in   )           :: time, mass</argument>
    <argument>type            (treeNode), intent(inout), optional :: node      </argument>
   </method>
   <method name="integrated" >
    <description>Return the halo mass function at {\normalfont \ttfamily time} [Gyr] integrated between {\normalfont \ttfamily massLow} and {\normalfont \ttfamily massHigh} [$M_\odot$].</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision          , intent(in   )                   :: time, massLow, massHigh</argument>
    <argument>type            (treeNode), intent(inout), target, optional :: node                   </argument>
    <modules>Numerical_Integration</modules>
    <code>
     double precision             :: logMassHigh, logMassLow
     type            (integrator) :: integrator_
     globalSelf => self
     if (present(node)) then
        self%node => node
     else
        self%node => null()
     end if
     self%time_ =  time
     logMassLow =log(massLow )
     logMassHigh=log(massHigh)
     integrator_=integrator(integratedIntegrand,toleranceRelative=1.0d-3,toleranceAbsolute=1.0d-100,integrationRule=GSL_Integ_Gauss15)
     haloMassFunctionIntegrated=integrator_%integrate(logMassLow,logMassHigh)
     return
    </code>
   </method>
   <method name="massFraction" >
    <description>Return the halo mass fraction at {\normalfont \ttfamily time} [Gyr] integrated between {\normalfont \ttfamily massLow} and {\normalfont \ttfamily massHigh} [$M_\odot$].</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision          , intent(in   )                   :: time, massLow, massHigh</argument>
    <argument>type            (treeNode), intent(inout), target, optional :: node                   </argument>
    <modules>Numerical_Integration</modules>
    <code>
     double precision             :: logMassHigh, logMassLow
     type            (integrator) :: integrator_
     globalSelf => self
     if (present(node)) then
        self%node => node
     else
        self%node => null()
     end if
     self%time_ =  time
     logMassLow =log(massLow )
     logMassHigh=log(massHigh)
     integrator_=integrator(massFractionIntegrand,toleranceRelative=1.0d-3,toleranceAbsolute=1.0d-100,integrationRule=GSL_Integ_Gauss15)
     haloMassFunctionMassFraction=integrator_%integrate(logMassLow,logMassHigh)
     ! Convert to a mass fraction.
     haloMassFunctionMassFraction=+haloMassFunctionMassFraction                &amp;
         &amp;                    /self%cosmologyParameters_%densityCritical() &amp;
         &amp;                    /self%cosmologyParameters_%OmegaMatter    ()
     return
    </code>
   </method>
  </functionClass>
  !!]

contains

  double precision function integratedIntegrand(logMass)
    !!{
    Integrand function used to integrate the dark matter halo mass function.
    !!}
    implicit none
    double precision, intent(in   ) :: logMass
    double precision                :: mass

    ! Extract integrand parameters.
    mass=exp(logMass)
    ! Return the differential mass function multiplied by mass since we are integrating over log(mass).
    if (associated(globalSelf%node)) then
       integratedIntegrand=+globalSelf%differential(globalSelf%time_,mass,node=globalSelf%node) &
            &              *                                         mass
    else
       integratedIntegrand=+globalSelf%differential(globalSelf%time_,mass                     ) &
            &              *                                         mass
    end if
    return
  end function integratedIntegrand

  double precision function massFractionIntegrand(logMass)
    !!{
    Integrand function used in computing the halo mass fraction.
    !!}
    implicit none
    double precision, intent(in   ) :: logMass
    double precision                :: mass

    ! Extract integrand parameters.
    mass=exp(logMass)
    ! Return the differential mass function multiplied by mass since we are integrating over log(mass) and by mass again to get
    ! the mass fraction.
    if (associated(globalSelf%node)) then
       massFractionIntegrand=+globalSelf%differential(globalSelf%time_,mass,node=globalSelf%node)    &
            &                *                                         mass                      **2
    else
       massFractionIntegrand=+globalSelf%differential(globalSelf%time_,mass                     )    &
            &                *                                         mass                      **2
    end if
    return
  end function massFractionIntegrand

end module Halo_Mass_Functions
