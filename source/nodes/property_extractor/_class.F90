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

!!{RST
Contains a module which provides a class that implements extraction of properties from nodes.
!!}

module Node_Property_Extractors
  !!{RST
  Provides a class that implements extraction of properties from nodes.
  !!}
  use :: Dust_Attenuation_Descriptors, only : decompositionRequest                         , luminosityDecomposition
  use :: Dust_Attenuations           , only : dustAttenuationClass
  use :: Galacticus_Nodes            , only : treeNode
  use :: Multi_Counters              , only : multiCounter
  use :: Output_Analyses_Options     , only : enumerationOutputAnalysisPropertyQuantityType, enumerationOutputAnalysisPropertyTypeType, outputAnalysisPropertyQuantityUnknown, outputAnalysisPropertyTypeLinear
  use :: ISO_Varying_String          , only : varying_string
  private
  public :: dustAbsorbedLuminosities, dustPhaseLabels

  !![
  <functionClass docformat="rst">
   <name>nodePropertyExtractor</name>
   <descriptiveName>Node Property Extractor</descriptiveName>
   <description>
   Class providing extraction of scalar, 1D, or multi-D properties from merger tree nodes for output. Property extractors are used by output analysis classes to retrieve galaxy and halo properties (e.g. stellar mass, dark matter halo mass, star formation rate, positions, velocities) and convert them to formats suitable for comparison with observational data or for writing to the Galacticus output file.
   </description>
   <default>nodeIndices</default>
   <method name="type" >
    <description>
    Return the type of the extracted property.
    </description>
    <type>type(enumerationOutputAnalysisPropertyTypeType)</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     nodePropertyExtractorType=outputAnalysisPropertyTypeLinear
    </code>
   </method>
   <method name="quantity" >
    <description>
    Return the class of the extracted property.
    </description>
    <type>type(enumerationOutputAnalysisPropertyQuantityType)</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     nodePropertyExtractorQuantity=outputAnalysisPropertyQuantityUnknown
    </code>
   </method>
   <method name="addInstances" >
    <description>
    Add multiple instances of this property to a ``multiCounter`` object.
    </description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type(treeNode    ), intent(inout) :: node</argument>
    <argument>type(multiCounter), intent(inout) :: instance</argument>
    <code>
     !$GLC attributes unused :: self, node, instance
     ! Nothing to do.
    </code>
   </method>
   <method name="extractScalar" >
    <description>
    Extract a scalar property from the given ``node``. This is a convenience method for callers (such as the
    ``libgalacticus`` library interface) that hold an extractor through the base class: it dispatches to the
    ``extract`` method of the scalar extractor subclass, and reports an error for extractors of any other class
    (whose ``extract`` methods have different signatures).
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <modules>Error</modules>
    <argument>type(treeNode), intent(inout), target :: node</argument>
    <code>
     select type (self)
     class is (nodePropertyExtractorScalar)
        nodePropertyExtractorExtractScalar=self%extract(node)
     class default
        nodePropertyExtractorExtractScalar=0.0d0
        call Error_Report('extractScalar requires an extractor of the scalar class'//{introspection:location})
     end select
    </code>
   </method>
   <method name="supportsAttenuation" >
    <description>
    Return true if this extractor is able to decompose its output into parcels of emission which can be attenuated by
    dust---that is, if it implements the ``decompose`` method. Extractors which do not produce a luminosity, and
    luminosity-producing extractors for which no decomposition has yet been implemented, return false. The
    ``dustAttenuation`` property extractor uses this to reject unusable children at construction time, rather than
    failing part way through a run.
    </description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     nodePropertyExtractorSupportsAttenuation=.false.
    </code>
   </method>
   <method name="decompose" >
    <description>
    Return the luminosity extracted from the given ``node`` at the given ``time``, split into parcels of emission which
    may each be attenuated differently by dust, at a resolution at least as fine as that specified by ``request``. Each
    parcel records the output element to which it contributes, so that the parcels can be recombined by the
    ``recompose`` method once they have been attenuated.

    Output elements are numbered as follows. A scalar extractor has a single element. A tuple extractor's elements are
    numbered in the order its ``extract`` method returns them. An array extractor, whose ``extract`` returns an array
    of shape ``(size,elementCount)``, numbers its elements in Fortran's column-major order, so that row :math:`i` of
    column :math:`j` is element :math:`i + (j-1) \times` ``size``---the position that element occupies in memory.

    The default implementation reports an error: only extractors whose ``supportsAttenuation`` method returns true
    override it. This mirrors the treatment of ``extractScalar``, and is used because Fortran provides no way to mix a
    decomposition interface into the several rank-specific extractor classes independently.
    </description>
    <type>type(luminosityDecomposition)</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <modules>Error</modules>
    <argument>type            (treeNode            ), intent(inout), target :: node</argument>
    <argument>double precision                      , intent(in   )         :: time</argument>
    <argument>type            (decompositionRequest), intent(in   )         :: request</argument>
    <code>
     !$GLC attributes unused :: self, node, time, request
     call nodePropertyExtractorDecompose%initialize(0,0)
     call Error_Report('this property extractor does not support dust attenuation'//{introspection:location})
    </code>
   </method>
   <method name="recompose" >
    <description>
    Recombine the parcels of a ``luminosityDecomposition``, each multiplied by the corresponding ``transmission``
    factor, into the output element values of this extractor. On return ``values`` has one element per output element
    of the decomposition, in the same order as the extractor's own ``extract`` method would return them.

    The default implementation simply sums each parcel into the output element which it contributes to, which is
    correct for any output that is linear in luminosity. Extractors whose output is not linear in luminosity---
    magnitudes, colors, or ratios, for example---must override this method.
    </description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type            (luminosityDecomposition), intent(in   )                            :: decomposition</argument>
    <argument>double precision                         , intent(in   ), dimension(:)              :: transmission</argument>
    <argument>double precision                         , intent(inout), dimension(:), allocatable :: values</argument>
    <code>
     !$GLC attributes unused :: self
     call decomposition%reduce(transmission,values)
    </code>
   </method>
  </functionClass>
  !!]

  ! Enumerations for galactic components.
  !![
  <enumeration docformat="rst">
   <name>galacticComponent</name>
   <description>
   Specifies the galactic component for various node property extractors.
   </description>
   <encodeFunction>yes</encodeFunction>
   <visibility>public</visibility>
   <entry label="disk"              />
   <entry label="spheroid"          />
   <entry label="nuclearStarCluster"/>
   <entry label="total"             />
  </enumeration>
  !!]

contains

  subroutine dustAbsorbedLuminosities(dustAttenuation_,extractor_,node,time,attenuated,absorbed,cosineInclination,weight)
    !!{RST
    Decompose the luminosity of ``extractor_`` into parcels of emission, and recompose from those parcels the luminosity
    transmitted through the dust of ``dustAttenuation_`` (``attenuated``, with one value per output element of the
    extractor), and the luminosity absorbed by each phase of that dust (``absorbed``, whose first index runs over output
    elements and second over phases). Absorption is averaged over orientation, where the attenuator depends on it, with
    the quadrature rule given by ``cosineInclination`` and ``weight``, which must be supplied if ``absorbed`` is
    requested.

    Both the dust attenuation and the dust emission property extractors form their luminosities through this routine,
    so that the luminosity one reports as absorbed is exactly the luminosity the other re-emits. It is a module
    procedure, rather than a method of either, because implementations of this class are generated into separate
    submodules which can see this module but not one another.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (dustAttenuationClass      ), intent(inout)                                        :: dustAttenuation_
    class           (nodePropertyExtractorClass), intent(inout), target                                :: extractor_
    type            (treeNode                  ), intent(inout), target                                :: node
    double precision                            , intent(in   )                                        :: time
    double precision                            , intent(inout), optional, allocatable, dimension(:  ) :: attenuated
    double precision                            , intent(inout), optional, allocatable, dimension(:,:) :: absorbed
    double precision                            , intent(in   ), optional             , dimension(:  ) :: cosineInclination, weight
    type            (luminosityDecomposition   )                                                       :: decomposition
    double precision                            , allocatable                         , dimension(:  ) :: transmission     , valuesPhase
    double precision                            , allocatable                         , dimension(:,:) :: fractions
    integer                                                                                            :: i                , k          , &
         &                                                                                                countPhases

    decomposition=extractor_%decompose(node,time,dustAttenuation_%request())
    ! The attenuator may refuse the component a parcel came from -- most refuse a luminosity summed over components,
    ! which can not be attenuated meaningfully. The child's component is not visible until it decomposes, so this is
    ! checked here rather than at construction.
    do i=1,decomposition%countTerms()
       if (.not.dustAttenuation_%supportsComponent(decomposition%descriptors(i)%componentType))                     &
            & call Error_Report(                                                                                    &
            &                   'the dust attenuation model refuses the component of a parcel of emission from "'// &
            &                   extractor_%objectType()                                                          // &
            &                   '" - a luminosity summed over components can not be attenuated'                  // &
            &                   {introspection:location}                                                            &
            &                  )
    end do
    if (present(attenuated)) then
       allocate(transmission(decomposition%countTerms()))
       if (decomposition%countTerms() > 0) transmission=dustAttenuation_%transmission(node,decomposition%descriptors)
       call extractor_%recompose(decomposition,transmission,attenuated)
    end if
    ! The luminosity absorbed by each phase of dust is recomposed from the same parcels, each weighted by the fraction its
    ! phase absorbs, averaged over orientation where the attenuator depends on it.
    if (present(absorbed)) then
       if (.not.(present(cosineInclination).and.present(weight))) &
            & call Error_Report('a quadrature rule over orientation is required to compute absorbed luminosities'//{introspection:location})
       countPhases=dustAttenuation_%countPhases()
       if (decomposition%countTerms() > 0) then
          fractions=dustAttenuation_%absorbedFractions(node,decomposition%descriptors,cosineInclination,weight)
       else
          allocate(fractions(0,countPhases))
       end if
       if (allocated(absorbed)) deallocate(absorbed)
       do k=1,countPhases
          call extractor_%recompose(decomposition,fractions(:,k),valuesPhase)
          if (.not.allocated(absorbed)) allocate(absorbed(size(valuesPhase),countPhases))
          absorbed(:,k)=valuesPhase
       end do
       if (.not.allocated(absorbed)) allocate(absorbed(0,0))
    end if
    return
  end subroutine dustAbsorbedLuminosities

  function dustPhaseLabels(dustAttenuation_) result(labels)
    !!{RST
    Return labels for the phases of dust of ``dustAttenuation_``, made unique: where two phases share a label---two
    members of a sequence of the same class, say---each such label has the index of its phase appended, so that
    properties named for them do not collide.
    !!}
    use :: ISO_Varying_String, only : operator(//), operator(==)
    implicit none
    type     (varying_string      ), allocatable  , dimension(:) :: labels
    class    (dustAttenuationClass), intent(inout)               :: dustAttenuation_
    logical                        , allocatable  , dimension(:) :: duplicated
    integer                                                      :: j               , k
    character(len=16              )                              :: labelIndex

    allocate(labels    (dustAttenuation_%countPhases()))
    allocate(duplicated(size(labels)                  ))
    do k=1,size(labels)
       labels(k)=dustAttenuation_%labelPhase(k)
    end do
    duplicated=.false.
    do k=1,size(labels)
       do j=1,size(labels)
          if (j /= k .and. labels(j) == labels(k)) duplicated(k)=.true.
       end do
    end do
    do k=1,size(labels)
       if (duplicated(k)) then
          write (labelIndex,'(i0)') k
          labels(k)=labels(k)//trim(labelIndex)
       end if
    end do
    return
  end function dustPhaseLabels

end module Node_Property_Extractors
