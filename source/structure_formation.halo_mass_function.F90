!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which provides a class that implements halo mass functions.

module Halo_Mass_Functions
  use, intrinsic :: ISO_C_Binding
  implicit none
  private

  !# <functionClass>
  !#  <name>haloMassFunction</name>
  !#  <descriptiveName>Halo Mass Function</descriptiveName>
  !#  <description>Class providing halo mass functions.</description>
  !#  <default>tinker2008</default>
  !#  <stateful>no</stateful>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <data>double precision                                    :: time_                         </data>
  !#  <data>class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()</data>
  !#  <data>
  !#   <scope>module</scope>
  !#   <threadprivate>yes</threadprivate>
  !#   <content>class(haloMassFunctionClass), pointer :: globalSelf</content>
  !#  </data>
  !#  <method name="differential" >
  !#   <description>Return the differential halo mass function for {\normalfont \ttfamily mass} [$M_\odot$] at {\normalfont \ttfamily time} [Gyr].</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: time, mass</argument>
  !#  </method>
  !#  <method name="integrated" >
  !#   <description>Return the halo mass function at {\normalfont \ttfamily time} [Gyr] integrated between {\normalfont \ttfamily massLow} and {\normalfont \ttfamily massHigh} [$M_\odot$].</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <selfTarget>yes</selfTarget>
  !#   <argument>double precision, intent(in   ) :: time, massLow, massHigh</argument>
  !#   <modules>Numerical_Integration</modules>
  !#   <code>
  !#    double precision                             :: logMassHigh         , logMassLow
  !#    type            (c_ptr                     ) :: parameterPointer
  !#    type            (fgsl_function             ) :: integrandFunction
  !#    type            (fgsl_integration_workspace) :: integrationWorkspace
  !#    globalSelf => self
  !#    self%time_ =  time
  !#    logMassLow =log(massLow )
  !#    logMassHigh=log(massHigh)
  !#    haloMassFunctionIntegrated=Integrate(                                         &amp;
  !#        &amp;                            logMassLow                             , &amp;
  !#        &amp;                            logMassHigh                            , &amp;
  !#        &amp;                            integratedIntegrand                    , &amp;
  !#        &amp;                            parameterPointer                       , &amp;
  !#        &amp;                            integrandFunction                      , &amp;
  !#        &amp;                            integrationWorkspace                   , &amp;
  !#        &amp;                            toleranceAbsolute   =0.0d+0            , &amp;
  !#        &amp;                            toleranceRelative   =1.0d-4            , &amp;
  !#        &amp;                            integrationRule     =FGSL_Integ_Gauss15  &amp;
  !#        &amp;                           )
  !#    call Integrate_Done(integrandFunction,integrationWorkspace)
  !#    return
  !#   </code>
  !#  </method>
  !#  <method name="massFraction" >
  !#   <description>Return the halo mass fraction at {\normalfont \ttfamily time} [Gyr] integrated between {\normalfont \ttfamily massLow} and {\normalfont \ttfamily massHigh} [$M_\odot$].</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <selfTarget>yes</selfTarget>
  !#   <argument>double precision, intent(in   ) :: time, massLow, massHigh</argument>
  !#   <modules>Numerical_Integration</modules>
  !#   <code>
  !#    double precision                             :: logMassHigh         , logMassLow
  !#    type            (c_ptr                     ) :: parameterPointer
  !#    type            (fgsl_function             ) :: integrandFunction
  !#    type            (fgsl_integration_workspace) :: integrationWorkspace
  !#    globalSelf => self
  !#    self%time_ =  time
  !#    logMassLow =log(massLow )
  !#    logMassHigh=log(massHigh)
  !#    haloMassFunctionMassFraction=Integrate(                                          &amp;
  !#         &amp;                             logMassLow                              , &amp;
  !#         &amp;                             logMassHigh                             , &amp;
  !#         &amp;                             massFractionIntegrand                   , &amp;
  !#         &amp;                             parameterPointer                        , &amp;
  !#         &amp;                             integrandFunction                       , &amp;
  !#         &amp;                             integrationWorkspace                    , &amp;
  !#         &amp;                             toleranceAbsolute    =0.0d+0            , &amp;
  !#         &amp;                             toleranceRelative    =1.0d-4            , &amp;
  !#         &amp;                             integrationRule      =FGSL_Integ_Gauss15  &amp;
  !#         &amp;                            )
  !#    call Integrate_Done(integrandFunction,integrationWorkspace)
  !#    ! Convert to a mass fraction.
  !#    haloMassFunctionMassFraction=+haloMassFunctionMassFraction                &amp;
  !#        &amp;                    /self%cosmologyParameters_%densityCritical() &amp;
  !#        &amp;                    /self%cosmologyParameters_%OmegaMatter    ()
  !#    return
  !#   </code>
  !#  </method>
  !# </functionClass>

contains

  function integratedIntegrand(logMass,parameterPointer) bind(c)
    !% Integrand function used to integrate the dark matter halo mass function.
    implicit none
    real(kind=c_double)        :: integratedIntegrand
    real(kind=c_double), value :: logMass
    type(c_ptr        ), value :: parameterPointer
    real(kind=c_double)        :: mass

    ! Extract integrand parameters.
    mass=exp(logMass)
    ! Return the differential mass function multiplied by mass since we are integrating over log(mass).
    integratedIntegrand=+globalSelf%differential(globalSelf%time_,mass) &
         &              *                                         mass
    return
  end function integratedIntegrand

  function massFractionIntegrand(logMass,parameterPointer) bind(c)
    !% Integrand function used in computing the halo mass fraction.
    implicit none
    real(kind=c_double)        :: massFractionIntegrand
    real(kind=c_double), value :: logMass
    type(c_ptr        ), value :: parameterPointer
    real(kind=c_double)        :: mass

    ! Extract integrand parameters.
    mass=exp(logMass)
    ! Return the differential mass function multiplied by mass since we are integrating over log(mass) and by mass again to get
    ! the mass fraction.
    massFractionIntegrand=+globalSelf%differential(globalSelf%time_,mass)    &
         &                *                                         mass **2
    return
  end function massFractionIntegrand

end module Halo_Mass_Functions
