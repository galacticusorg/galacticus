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
Contains a module which implements a class for stellar spectra postprocessors.
!!}

module Stellar_Population_Spectra_Postprocess
  !!{RST
  Implements a class for stellar spectra postprocessors.
  !!}
  private
  public :: stellarPopulationSpectraPostprocessorList

  !![
  <functionClass docformat="rst">
   <name>stellarPopulationSpectraPostprocessor</name>
   <descriptiveName>Postprocessors for stellar population spectra</descriptiveName>
   <description>
   Class providing postprocessors for stellar population spectra---multiplicative corrections to the spectral energy distribution applied after the intrinsic stellar population spectrum is computed. Each postprocessor returns a wavelength-, age-, and redshift-dependent multiplicative factor to account for effects such as :term:`IGM` absorption of Lyman-series photons, intergalactic dust attenuation, or other propagation effects between the source and the observer. Multiple postprocessors can be chained to apply cumulative corrections for different physical processes.
   </description>
   <default>inoue2014</default>
   <method name="multiplier" >
    <description>
    Return the multiplicative correction factor (between 0 and 1) to apply to the stellar population spectrum at the given wavelength (Å), stellar population age (Gyr), and source redshift, accounting for propagation effects such as IGM absorption.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: wavelength, age, redshift</argument>
   </method>
   <method name="ageRange" >
    <description>
    Return the range of stellar population ages, in Gyr, over which this postprocessor's multiplier is non-zero. The
    default is unbounded; postprocessors which suppress emission outside a range of ages override it.

    This is what allows a consumer to discover which of a set of named postprocessing chains corresponds to which
    range of ages---so that, for example, a luminosity computed with a chain restricted to young populations can be
    subtracted from one computed with an unrestricted chain to isolate the light of old populations---rather than
    requiring the user to declare that correspondence and keep it consistent by hand.
    </description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, intent(  out) :: ageMinimum, ageMaximum</argument>
    <code>
     !$GLC attributes unused :: self
     ageMinimum=0.0d0
     ageMaximum=huge(0.0d0)
    </code>
   </method>
   <method name="ageWindowIsSharp" >
    <description>
    Return true if the age dependence of this postprocessor's multiplier is a sharp window: that is, if the multiplier
    is independent of age within the range reported by ``ageRange`` and zero outside it.

    Postprocessors which have no age dependence at all satisfy this trivially, and so the default is true. Those which
    taper---:galacticus-class:`stellarPopulationSpectraPostprocessorBirthCloudsLacey2016`, which declines linearly
    between one and two birth cloud lifetimes, and
    :galacticus-class:`stellarPopulationSpectraPostprocessorUnescaped`, which decays exponentially with age---do not,
    and must return false.

    The distinction matters because only a sharp window lets a luminosity be split cleanly into age bins by
    differencing chains. Applying that construction to a tapering postprocessor would silently mis-apportion light
    between the bins, so a consumer must check this before relying on ``ageRange``.
    </description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     stellarPopulationSpectraPostprocessorAgeWindowIsSharp=.true.
    </code>
   </method>
   <method name="isRedshiftDependent" >
    <description>
    Return true if this postprocessor's correction factor depends on the source redshift, allowing the ODE solver to determine whether the luminosity must be recomputed when the redshift changes rather than using a cached value.
    </description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
      !$GLC attributes unused :: self
      stellarPopulationSpectraPostprocessorIsRedshiftDependent=.true.
    </code>
   </method>
   <method name="isAgeDependent" >
    <description>
    Return true if this postprocessor's correction factor depends on the age of the stellar population. The default is true,
    which is always safe; postprocessors whose multiplier is independent of age override it.
    </description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
      !$GLC attributes unused :: self
      stellarPopulationSpectraPostprocessorIsAgeDependent=.true.
    </code>
   </method>
   <method name="isSeparable" >
    <description>
    Return true if this postprocessor's multiplier separates into the product of a factor independent of redshift and a
    factor independent of age,

    .. math::

       M(\lambda,\tau,z) = M_\tau(\lambda,\tau) M_z(\lambda,z),

    where :math:`\tau` is the age of the stellar population, the two factors being returned by ``multiplierAge`` and
    ``multiplierRedshift`` respectively.

    This allows the redshift-dependent factor to be applied after the stellar population has been integrated over its star
    formation history, so that a single template spectrum, built with only the redshift-independent factor, can serve
    galaxies at all redshifts. A postprocessor is separable unless it depends on both age and redshift. A
    :galacticus-class:`stellarPopulationSpectraPostprocessorSequence` is separable if each of its members is, even though it
    may depend on both.
    </description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
      stellarPopulationSpectraPostprocessorIsSeparable=.not.(self%isRedshiftDependent().and.self%isAgeDependent())
    </code>
   </method>
   <method name="multiplierAge" >
    <description>
    Return the redshift-independent factor, :math:`M_\tau(\lambda,\tau)`, of the multiplier of a separable postprocessor (see
    ``isSeparable``) at the given wavelength (Å) and stellar population age (Gyr).
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <modules>Error</modules>
    <argument>double precision, intent(in   ) :: wavelength, age</argument>
    <code>
      if (self%isRedshiftDependent()) then
         if (self%isAgeDependent()) then
            stellarPopulationSpectraPostprocessorMultiplierAge=0.0d0
            call Error_Report('postprocessor is not separable'//{introspection:location})
         else
            stellarPopulationSpectraPostprocessorMultiplierAge=1.0d0
         end if
      else
         ! The multiplier is independent of redshift, so any redshift may be passed.
         stellarPopulationSpectraPostprocessorMultiplierAge=self%multiplier(wavelength,age,redshift=0.0d0)
      end if
    </code>
   </method>
   <method name="multiplierRedshift" >
    <description>
    Return the age-independent factor, :math:`M_z(\lambda,z)`, of the multiplier of a separable postprocessor (see
    ``isSeparable``) at the given wavelength (Å, in the rest frame of the source) and source redshift.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <modules>Error</modules>
    <argument>double precision, intent(in   ) :: wavelength, redshift</argument>
    <code>
      if (self%isRedshiftDependent()) then
         if (self%isAgeDependent()) then
            stellarPopulationSpectraPostprocessorMultiplierRedshift=0.0d0
            call Error_Report('postprocessor is not separable'//{introspection:location})
         else
            ! The multiplier is independent of age, so any age may be passed.
            stellarPopulationSpectraPostprocessorMultiplierRedshift=self%multiplier(wavelength,age=0.0d0,redshift=redshift)
         end if
      else
         stellarPopulationSpectraPostprocessorMultiplierRedshift=1.0d0
      end if
    </code>
   </method>
  </functionClass>
  !!]

  type :: stellarPopulationSpectraPostprocessorList
     !!{RST
     Type used to build linked list of stellar population spectra postprocessors.
     !!}
     class(stellarPopulationSpectraPostprocessorClass), pointer :: stellarPopulationSpectraPostprocessor_ => null()
   contains
     !![
     <methods docformat="rst">
       <method method="assignment(=)" description="Assign postprocessor list objects."/>
     </methods>
     !!]
     final     ::                  stellarPopulationSpectraPostprocessorListDestructor
     procedure ::                  stellarPopulationSpectraPostprocessorListAssign
     generic   :: assignment(=) => stellarPopulationSpectraPostprocessorListAssign
  end type stellarPopulationSpectraPostprocessorList

  !![
  <functionClass docformat="rst">
   <name>stellarPopulationSpectraPostprocessorBuilder</name>
   <descriptiveName>Stellar Population Spectra Postprocessor Builders</descriptiveName>
   <description>
   Class providing builders for postprocessors for stellar population spectra. These act as a factory for ``stellarPopulationSpectraPostprocessor`` objects. Different postprocessors can be applied to different filters. The ``[luminosityPostprocessSet]`` parameter specifies, for each filter, a descriptor which is passed to the builder object, which then uses that descriptor to build a postprocessor. (If this parameter is not present then "``default``" is assumed for all filters.)
   </description>
   <default>lookup</default>
   <method name="build" >
    <description>
    Build and return a :galacticus-class:`stellarPopulationSpectraPostprocessorClass` object corresponding to the given string descriptor, acting as a factory that constructs the appropriate postprocessor (e.g., IGM absorption, dust) for each photometric filter.
    </description>
    <type>class(stellarPopulationSpectraPostprocessorClass)</type>
    <pass>yes</pass>
    <argument>type(varying_string), intent(in   ) :: descriptor</argument>
   </method>
  </functionClass>
  !!]

contains

  subroutine stellarPopulationSpectraPostprocessorListDestructor(self)
    !!{RST
    Destructor for elements of stellar population spectra postprocessor lists.
    !!}
    implicit none
    type(stellarPopulationSpectraPostprocessorList), intent(inout) :: self

    !![
    <objectDestructor name="self%stellarPopulationSpectraPostprocessor_"/>
    !!]
    return
  end subroutine stellarPopulationSpectraPostprocessorListDestructor

  recursive subroutine stellarPopulationSpectraPostprocessorListAssign(self,from)
    !!{RST
    Perform assignment for the :galacticus-class:`stellarPopulationSpectraPostprocessorList` class.
    !!}
    implicit none
    class(stellarPopulationSpectraPostprocessorList), intent(  out) :: self
    class(stellarPopulationSpectraPostprocessorList), intent(in   ) :: from

    self%stellarPopulationSpectraPostprocessor_ => from%stellarPopulationSpectraPostprocessor_
    !![
    <referenceCountIncrement owner="self" object="stellarPopulationSpectraPostprocessor_"/>
    !!]
    return
  end subroutine stellarPopulationSpectraPostprocessorListAssign

end module Stellar_Population_Spectra_Postprocess
