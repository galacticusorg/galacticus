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

!+    Contributions to this file made by: Andrew Benson
!+    Implemented with assistance from Claude (Anthropic).

!!{RST
Contains a module which provides the data structures used to describe emission to be attenuated by dust.
!!}

module Dust_Attenuation_Descriptors
  !!{RST
  Provides the data structures which mediate between objects that *produce* luminosity (typically
  :galacticus-class:`nodePropertyExtractorClass` objects) and objects that *attenuate* it (``dustAttenuation``
  objects).

  Dust attenuation is not, in general, a single multiplicative factor applied to a galaxy's total luminosity: it
  depends on the wavelength of the emission, on which component the emission comes from, and---where birth clouds are
  modelled---on the age of the emitting stellar population. A luminosity-producing object must therefore be able to
  hand over its luminosity *split* along whichever of those axes matter, have each piece attenuated, and then
  recombine the pieces.

  Three types cooperate to make that possible:

  ``emissionDescriptor``
     Describes a single, homogeneous parcel of emission---one wavelength, one component, one source type, one age
     range. This is what a ``dustAttenuation`` object is asked to attenuate.

  ``luminosityDecomposition``
     A set of such parcels together with their luminosities, plus a map recording which output element each parcel
     contributes to. This is what a luminosity-producing object hands over.

  ``decompositionRequest``
     Describes how finely the decomposition must be split. A ``dustAttenuation`` object builds this to declare which
     axes it actually depends upon, so that a producer can collapse the axes which are irrelevant. This matters for
     performance: an unnegotiated decomposition of a spectrum would have one term per (wavelength, age, metallicity)
     triple, while an age-independent screen applied to a broad-band luminosity needs only one term per component.
  !!}
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType, componentTypeUnknown
  implicit none
  private
  public :: emissionDescriptor, luminosityDecomposition, decompositionRequest

  !![
  <enumeration docformat="rst">
   <name>emissionSource</name>
   <description>
   Specifies the physical origin of a parcel of emission which is to be attenuated by dust. Dust attenuation models may
   treat these differently---for example, nebular emission arises from H II regions which are embedded in their birth
   clouds and so may be more heavily obscured than the stellar continuum of the same age, while emission from an
   accretion disk originates at the very centre of the galaxy and so is seen through the full column of the galaxy's
   dust.
   </description>
   <encodeFunction>yes</encodeFunction>
   <decodeFunction>yes</decodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="stellar"       description="Stellar continuum emission"               />
   <entry label="nebular"       description="Nebular continuum and line emission"      />
   <entry label="accretionDisk" description="Emission from a black hole accretion disk"/>
   <entry label="dust"          description="Thermal re-emission from dust grains"     />
   <entry label="unknown"       description="Emission of unknown or unspecified origin"/>
  </enumeration>
  !!]

  type :: emissionDescriptor
     !!{RST
     Describes a single, homogeneous parcel of emission to be attenuated by dust.

     The fields are:

     ``wavelength``
        Rest-frame wavelength of the emission, in Å.
     ``componentType``
        The galactic component from which the emission originates.
     ``sourceType``
        The physical origin of the emission.
     ``ageMinimum``, ``ageMaximum``
        The range of stellar population ages contributing to this parcel, in Gyr. A parcel whose age is unresolved
        ---because the attenuation model does not depend on age, or because the producer cannot resolve it---spans
        the full range, from zero to ``huge(0.0d0)``. Attenuators which distinguish young from old populations should
        therefore test ``ageMaximum`` against their birth cloud lifetime, so that an unresolved parcel is treated as
        old rather than young.
     ``metallicity``
        Metallicity of the emitting material (linear, by mass); negative if unresolved.
     ``radius``
        Galactocentric radius at which the emission originates, in Mpc; negative if unresolved.
     !!}
     double precision                                :: wavelength   =-1.0d0
     type            (enumerationComponentTypeType ) :: componentType= componentTypeUnknown
     type            (enumerationEmissionSourceType) :: sourceType   = emissionSourceUnknown
     double precision                                :: ageMinimum   = 0.0d0
     double precision                                :: ageMaximum   = huge(0.0d0)
     double precision                                :: metallicity  =-1.0d0
     double precision                                :: radius       =-1.0d0
  end type emissionDescriptor

  type :: luminosityDecomposition
     !!{RST
     A luminosity, split into parcels which may each be attenuated differently, together with a map recording which
     output element each parcel contributes to.

     ``elementIndex(i)`` gives the (1-based) index of the output element to which parcel ``i`` contributes, and must
     lie in the range :math:`[1,` ``countElements`` :math:`]`. Several parcels may share an output element---indeed
     that is the usual case, since the whole point of the decomposition is that a single output luminosity is built
     from parcels which attenuate differently.
     !!}
     type            (emissionDescriptor), allocatable, dimension(:) :: descriptors
     double precision                    , allocatable, dimension(:) :: luminosities
     integer                             , allocatable, dimension(:) :: elementIndex
     integer                                                         :: countElements=0
   contains
     !![
     <methods docformat="rst">
       <method method="initialize" description="Allocate storage for the given numbers of parcels and output elements."/>
       <method method="countTerms" description="Return the number of parcels in the decomposition."                    />
       <method method="reduce"     description="Sum attenuated parcel luminosities into their output elements."        />
     </methods>
     !!]
     procedure :: initialize => luminosityDecompositionInitialize
     procedure :: countTerms => luminosityDecompositionCountTerms
     procedure :: reduce     => luminosityDecompositionReduce
  end type luminosityDecomposition

  type :: decompositionRequest
     !!{RST
     Declares how finely a ``luminosityDecomposition`` must be split for a given dust attenuation model.

     ``ageBoundaries`` holds the *internal* boundaries of the age bins required, in Gyr, in increasing order. An empty
     (or unallocated) array requests no age resolution at all; a single boundary at the birth cloud lifetime requests
     a young/old split; and so on. A producer should return parcels which do not straddle any of these boundaries.

     The ``resolve*`` flags request splitting along the remaining axes. Producers are permitted to return a *finer*
     decomposition than requested---the request is a lower bound on resolution, not an upper one---but returning a
     coarser one will give wrong answers.
     !!}
     double precision, allocatable, dimension(:) :: ageBoundaries
     logical                                     :: resolveComponents =.true.
     logical                                     :: resolveMetallicity=.false.
     logical                                     :: resolveRadius     =.false.
   contains
     !![
     <methods docformat="rst">
       <method method="countAgeBins" description="Return the number of age bins requested."                       />
       <method method="ageBinIndex"  description="Return the index of the age bin containing the given age."       />
       <method method="ageBinRange"  description="Return the age range spanned by the age bin of the given index."/>
     </methods>
     !!]
     procedure :: countAgeBins => decompositionRequestCountAgeBins
     procedure :: ageBinIndex  => decompositionRequestAgeBinIndex
     procedure :: ageBinRange  => decompositionRequestAgeBinRange
  end type decompositionRequest

contains

  subroutine luminosityDecompositionInitialize(self,countTerms,countElements)
    !!{RST
    Allocate storage in a ``luminosityDecomposition`` for ``countTerms`` parcels contributing to ``countElements``
    output elements. Luminosities are initialized to zero and element indices to one, so that a producer need only
    fill in the parcels which are actually non-zero.
    !!}
    implicit none
    class  (luminosityDecomposition), intent(inout) :: self
    integer                         , intent(in   ) :: countTerms, countElements

    if (allocated(self%descriptors )) deallocate(self%descriptors )
    if (allocated(self%luminosities)) deallocate(self%luminosities)
    if (allocated(self%elementIndex)) deallocate(self%elementIndex)
    allocate(self%descriptors (countTerms))
    allocate(self%luminosities(countTerms))
    allocate(self%elementIndex(countTerms))
    self%luminosities =0.0d0
    self%elementIndex =1
    self%countElements=countElements
    return
  end subroutine luminosityDecompositionInitialize

  integer function luminosityDecompositionCountTerms(self) result(countTerms)
    !!{RST
    Return the number of parcels in a ``luminosityDecomposition``.
    !!}
    implicit none
    class(luminosityDecomposition), intent(in   ) :: self

    if (allocated(self%luminosities)) then
       countTerms=size(self%luminosities)
    else
       countTerms=0
    end if
    return
  end function luminosityDecompositionCountTerms

  subroutine luminosityDecompositionReduce(self,transmission,values)
    !!{RST
    Sum the parcel luminosities of a ``luminosityDecomposition``, each multiplied by the corresponding
    ``transmission`` factor, into the output elements to which they contribute. This is the default reduction, and is
    correct for any output which is linear in luminosity. Producers whose output is *not* linear in luminosity---
    magnitudes, colours, or ratios, for example---must reduce their decomposition themselves.

    ``values`` is allocated to ``countElements`` and zeroed before accumulation.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (luminosityDecomposition), intent(in   )                            :: self
    double precision                         , intent(in   ), dimension(:)              :: transmission
    double precision                         , intent(inout), dimension(:), allocatable :: values
    integer                                                                             :: i

    if (size(transmission) /= self%countTerms())                                                                      &
         & call Error_Report('number of transmission factors does not match number of terms'//{introspection:location})
    if (allocated(values)) deallocate(values)
    allocate(values(self%countElements))
    values=0.0d0
    do i=1,self%countTerms()
       if (self%elementIndex(i) < 1 .or. self%elementIndex(i) > self%countElements)                                   &
            & call Error_Report('element index out of range'                             //{introspection:location})
       values(self%elementIndex(i))=+values      (self%elementIndex(i)) &
            &                       +self        %luminosities    (i)   &
            &                       *transmission                 (i)
    end do
    return
  end subroutine luminosityDecompositionReduce

  integer function decompositionRequestCountAgeBins(self) result(countAgeBins)
    !!{RST
    Return the number of age bins requested by a ``decompositionRequest``. This is one more than the number of
    internal boundaries, and is always at least one.
    !!}
    implicit none
    class(decompositionRequest), intent(in   ) :: self

    if (allocated(self%ageBoundaries)) then
       countAgeBins=size(self%ageBoundaries)+1
    else
       countAgeBins=1
    end if
    return
  end function decompositionRequestCountAgeBins

  integer function decompositionRequestAgeBinIndex(self,age) result(indexBin)
    !!{RST
    Return the index of the age bin of a ``decompositionRequest`` which contains the given ``age`` (in Gyr). Bins are
    closed at their lower boundary and open at their upper boundary, so an age exactly equal to a boundary falls in
    the *older* bin.
    !!}
    implicit none
    class           (decompositionRequest), intent(in   ) :: self
    double precision                      , intent(in   ) :: age

    indexBin=1
    if (.not.allocated(self%ageBoundaries)) return
    do while (indexBin <= size(self%ageBoundaries))
       if (age < self%ageBoundaries(indexBin)) exit
       indexBin=indexBin+1
    end do
    return
  end function decompositionRequestAgeBinIndex

  subroutine decompositionRequestAgeBinRange(self,indexBin,ageMinimum,ageMaximum)
    !!{RST
    Return the range of ages (in Gyr) spanned by the age bin of a ``decompositionRequest`` with the given index. The
    lowest bin begins at zero age, and the highest extends to ``huge(0.0d0)``.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (decompositionRequest), intent(in   ) :: self
    integer                               , intent(in   ) :: indexBin
    double precision                      , intent(  out) :: ageMinimum, ageMaximum

    if (indexBin < 1 .or. indexBin > self%countAgeBins()) then
       ageMinimum=0.0d0
       ageMaximum=0.0d0
       call Error_Report('age bin index out of range'//{introspection:location})
       return
    end if
    if (indexBin == 1                   ) then
       ageMinimum=0.0d0
    else
       ageMinimum=self%ageBoundaries(indexBin-1)
    end if
    if (indexBin == self%countAgeBins()) then
       ageMaximum=huge(0.0d0)
    else
       ageMaximum=self%ageBoundaries(indexBin  )
    end if
    return
  end subroutine decompositionRequestAgeBinRange

end module Dust_Attenuation_Descriptors
