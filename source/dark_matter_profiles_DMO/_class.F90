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
Contains a module which provides an object that implements dark matter halo profiles.
!!}

module Dark_Matter_Profiles_DMO
  !!{RST
  Provides an object that implements dark matter halo profiles.
  !!}
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScale              , darkMatterHaloScaleClass
  use :: Galacticus_Nodes          , only : treeNode
  use :: Mass_Distributions        , only : massDistributionClass            , massDistributionHeatingClass
  use :: Galactic_Structure_Options, only : enumerationStructureErrorCodeType, enumerationWeightByType
  private

  !![
  <functionClass docformat="rst">
   <name>darkMatterProfileDMO</name>
   <descriptiveName>Dark Matter Only Halo Profiles</descriptiveName>
   <description>
   Class providing dark matter-only halo density profiles, i.e. the profile a halo would have in the absence of baryonic effects. This returns a :galacticus-class:`massDistributionClass` object for the specified node. Common implementations include NFW and Einasto profiles parameterized by a scale radius or concentration. This class is used in calculations of dynamical friction, tidal stripping, and other processes where the unmodified dark matter profile is needed.
   </description>
   <default>NFW</default>
   <method name="get" >
    <description>
    Return the mass distribution of the dark matter-only profile.
    </description>
    <type>class(massDistributionClass)</type>
    <pass>yes</pass>
    <argument>type   (treeNode               ), intent(inout)           :: node       </argument>
    <argument>type   (enumerationWeightByType), intent(in   ), optional :: weightBy   </argument>
    <argument>integer                         , intent(in   ), optional :: weightIndex</argument>
   </method>
   <method name="scaleRadiusValidated" >
    <description>
    Return the scale radius of the dark matter profile of ``node``, reporting a fatal error if it is not positive. Profiles
    parameterized by a scale radius should obtain it through this method, so that a scale radius which was never set, or which
    was set to zero, is reported rather than propagated into the profile.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <modules>Display Error Galacticus_Nodes ISO_Varying_String String_Handling</modules>
    <argument>type(treeNode), intent(inout) :: node</argument>
    <code>
     class           (nodeComponentDarkMatterProfile), pointer   :: darkMatterProfile
     type            (varying_string                )            :: message
     character       (len=16                        )            :: label
     double precision                                , parameter :: radiusScaleUnset =-1.0d0
     darkMatterProfile                        => node             %darkMatterProfile()
     darkMatterProfileDMOScaleRadiusValidated =  darkMatterProfile%scale            ()
     if (darkMatterProfileDMOScaleRadiusValidated &gt; 0.0d0) return
     message='the ['//char(self%objectType())//'] dark matter profile requires a positive scale radius, but '
     if (darkMatterProfileDMOScaleRadiusValidated == radiusScaleUnset) then
        ! This is the class default value of the scale radius, so it was never set.
        message=message                                                                                                                // &amp;
             &amp; 'the scale radius of node '//node%index()//' has not been set'//char(10)                                            // &amp;
             &amp; displayGreen()//'HELP:'//displayReset()                                                                             // &amp;
             &amp; ' scale radii are set by [nodeOperator]=darkMatterProfileScaleSet - check that it is present, and that it is'       // &amp;
             &amp; ' applied to every halo (an operator nested inside a filtering operator, such as [nodeOperator]=filteredMainBranch,'// &amp;
             &amp; ' is applied to only some halos)'
     else
        write (label,'(e12.6)') darkMatterProfileDMOScaleRadiusValidated
        message=message                                                                                                                // &amp;
             &amp; 'the scale radius of node '//node%index()//' is '//trim(adjustl(label))//' Mpc'//char(10)                           // &amp;
             &amp; displayGreen()//'HELP:'//displayReset()                                                                             // &amp;
             &amp; ' check [darkMatterProfileScaleRadius] - the "zero" class, for example, is suitable only for profiles which have no'// &amp;
             &amp; ' scale radius, such as [darkMatterProfileDMO]=isothermal'
     end if
     call Error_Report(message//{introspection:location})
    </code>
   </method>
  </functionClass>
  !!]

  !![
  <functionClass docformat="rst">
   <name>darkMatterProfileHeating</name>
   <descriptiveName>Dark Matter Profile Heating</descriptiveName>
   <description>
   Class providing models of heating applied to dark matter-only halo profiles. Heating can modify the density profile of a dark matter halo, for example due to tidal shocks, dynamical heating from baryons, or other perturbative processes. This class returns a :galacticus-class:`massDistributionHeatingClass` object encoding the heating distribution for a given node.
   </description>
   <default>null</default>
   <method name="get" >
    <description>
    Return the dark matter profile heating in the dark matter-only profile.
    </description>
    <type>class(massDistributionHeatingClass)</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Dark_Matter_Profiles_DMO
