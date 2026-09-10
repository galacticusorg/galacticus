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
Contains a module which provides a class implementing dust extinction curves.
!!}

module Dust_Extinction_Curves
  !!{RST
  Provides a class implementing dust extinction curves---the wavelength dependence of the opacity of interstellar
  dust, normalized to its value in the :math:`V`-band.

  An extinction curve carries no information about *how much* dust lies along a line of sight, only about the relative
  efficiency with which dust absorbs and scatters light as a function of wavelength. The normalization is supplied
  separately by a ``dustAttenuation`` object, which multiplies the curve by a :math:`V`-band optical depth computed
  from the properties of a galaxy. Keeping the two apart means any curve can be combined with any normalization.

  Curves are normalized to their value in the :math:`V`-band, so that for a simple screen the transmission is

  .. math::

     T(\lambda) = \exp\left[-\tau_\mathrm{V} \frac{k(\lambda)}{k_\mathrm{V}}\right].

  Note that the normalization is only approximate for the empirically fit curves. Each published fit anchors its
  :math:`V`-band normalization at its own nominal wavelength---:cite:t:`cardelli_relationship_1989` at
  :math:`x=1.82\,\mu\mathrm{m}^{-1}`, i.e. :math:`0.5495\,\mu\mathrm{m}`---which differs slightly from the
  :math:`5504.6\,`Å effective wavelength of the Buser :math:`V` filter that Galacticus uses, and the fitting functions
  do not evaluate to exactly unity even there. The residual is a few parts in a thousand (0.14% for
  :cite:t:`calzetti_dust_2000`, 0.21% for :cite:t:`cardelli_relationship_1989`). The curves are deliberately *not*
  renormalized to remove it, both because doing so would depart from the published fits and because it would change
  the results of existing models.

  Many published curves are fit only over a limited range of wavelengths. The ``wavelengthRange`` method reports that
  range, and the description of each implementation states what it does outside it---several return zero, which is to
  say *no* attenuation rather than clamped attenuation, and that distinction matters when a curve is evaluated over a
  broad spectrum rather than through an optical filter.
  !!}
  implicit none
  private
  public :: wavelengthVBand

  ! Effective wavelength of the Buser V-band filter, in Angstroms. Extinction curves are normalized to unity here.
  double precision, parameter :: wavelengthVBand=5504.61227375652d0

  !![
  <functionClass docformat="rst">
   <name>dustExtinctionCurve</name>
   <descriptiveName>Dust Extinction Curves</descriptiveName>
   <description>
   Class providing the wavelength dependence of dust opacity, normalized to unity in the :math:`V`-band.
   </description>
   <default>calzetti2000</default>
   <method name="attenuationRelative" >
    <description>
    Return the opacity of dust at the given rest-frame ``wavelength`` (in Å), relative to its value in the
    :math:`V`-band; that is, :math:`k(\lambda)/k_\mathrm{V}`. The result is dimensionless, and is approximately (but,
    for the empirically fit curves, not exactly) unity at the :math:`V`-band effective wavelength---see the discussion
    of normalization in the class description.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: wavelength</argument>
   </method>
   <method name="wavelengthRange" >
    <description>
    Return the range of rest-frame wavelengths (in Å) over which this curve is defined. The default is unbounded;
    implementations fit over a limited range override it. Callers which evaluate a curve across a broad spectrum---the
    dust emission calculations in particular---should consult this rather than assume the curve is meaningful
    everywhere.
    </description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, intent(  out) :: wavelengthMinimum, wavelengthMaximum</argument>
    <code>
     !$GLC attributes unused :: self
     wavelengthMinimum=0.0d0
     wavelengthMaximum=huge(0.0d0)
    </code>
   </method>
  </functionClass>
  !!]

end module Dust_Extinction_Curves
