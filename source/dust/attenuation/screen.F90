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
  Implements a dust attenuation class describing a uniform screen of dust in front of the emitting stars.
  !!}

  use :: Dust_Extinction_Curves, only : dustExtinctionCurveClass

  !![
  <dustAttenuation name="dustAttenuationScreen" abstract="yes" docformat="rst">
   <description>
   An abstract dust attenuation class describing a uniform screen of dust lying between the emitting stars and the
   observer, so that

   .. math::

      T(\lambda) = \exp\left[ -\tau_\mathrm{V} \frac{k(\lambda)}{k_\mathrm{V}} \right],

   where the wavelength dependence :math:`k(\lambda)/k_\mathrm{V}` is supplied by a
   :galacticus-class:`dustExtinctionCurveClass` object and the normalization :math:`\tau_\mathrm{V}` by the
   ``depthOpticalV`` method of a concrete subclass. Separating the two means any extinction curve may be combined with
   any prescription for how much dust a galaxy has.

   A screen is the simplest possible geometry. It attenuates all of the light it is applied to equally, and---unlike a
   mixed distribution of dust and stars---has no upper limit to the attenuation it can produce.
   </description>
  </dustAttenuation>
  !!]

  ! Note: declared as a concrete type despite `abstract="yes"` above, so that it may carry a finalizer to release the
  ! extinction curve object; a final subroutine's `self` cannot be declared `type(...)` for an abstract type. The
  ! `depthOpticalV` method therefore cannot be `deferred`, and instead reports an error if a subclass fails to
  ! override it -- the same idiom used by `nodePropertyExtractor%decompose`.
  type, extends(dustAttenuationClass) :: dustAttenuationScreen
     !!{RST
     A dust attenuation class describing a uniform screen of dust.
     !!}
     private
     class(dustExtinctionCurveClass), pointer :: dustExtinctionCurve_ => null()
   contains
     !![
     <methods docformat="rst">
       <method method="depthOpticalV" description="Return the :math:`V`-band optical depth of the screen for one component of one node."/>
     </methods>
     !!]
     final     ::                  screenDestructor
     procedure :: transmission  => screenTransmission
     procedure :: depthOpticalV => screenDepthOpticalV
  end type dustAttenuationScreen

contains

  subroutine screenDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`dustAttenuationScreen` dust attenuation class.
    !!}
    implicit none
    type(dustAttenuationScreen), intent(inout) :: self

    !![
    <objectDestructor name="self%dustExtinctionCurve_"/>
    !!]
    return
  end subroutine screenDestructor

  double precision function screenDepthOpticalV(self,node,componentType) result(depthOpticalV)
    !!{RST
    Return the :math:`V`-band optical depth of the screen. The default implementation reports an error: every
    concrete subclass of :galacticus-class:`dustAttenuationScreen` must override it.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(dustAttenuationScreen         ), intent(inout)         :: self
    type (treeNode                      ), intent(inout), target :: node
    type (enumerationComponentTypeType  ), intent(in   )         :: componentType
    !$GLC attributes unused :: self, node, componentType

    depthOpticalV=0.0d0
    call Error_Report('this screen does not implement depthOpticalV'//{introspection:location})
    return
  end function screenDepthOpticalV

  function screenTransmission(self,node,descriptors) result(transmission)
    !!{RST
    Return the transmission through a uniform screen of dust.

    The :math:`V`-band optical depth is evaluated once per distinct component appearing in ``descriptors``, rather
    than once per parcel: a spectrum decomposed at hundreds of wavelengths would otherwise repeat the same structural
    calculation at each one.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeMax, componentTypeMin
    implicit none
    class           (dustAttenuationScreen), intent(inout)                                               :: self
    type            (treeNode             ), intent(inout), target                                       :: node
    type            (emissionDescriptor   ), intent(in   ), dimension(:                                ) :: descriptors
    double precision                                      , dimension(size(descriptors)                ) :: transmission
    double precision                                      , dimension(componentTypeMin:componentTypeMax) :: depthOpticalV
    logical                                               , dimension(componentTypeMin:componentTypeMax) :: computed
    integer                                                                                              :: i            , indexComponent

    computed=.false.
    do i=1,size(descriptors)
       indexComponent=descriptors(i)%componentType%ID
       if (.not.computed(indexComponent)) then
          depthOpticalV(indexComponent)=self%depthOpticalV(node,descriptors(i)%componentType)
          computed     (indexComponent)=.true.
       end if
       transmission(i)=exp(                                                                          &
            &              -depthOpticalV(indexComponent)                                            &
            &              *self%dustExtinctionCurve_%attenuationRelative(descriptors(i)%wavelength) &
            &             )
    end do
    return
  end function screenTransmission
