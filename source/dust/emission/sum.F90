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
  Implements a dust emission spectrum which sums other spectra.
  !!}

  type, public :: dustEmissionSpectrumList
     class(dustEmissionSpectrumClass), pointer :: dustEmissionSpectrum_ => null()
     type (dustEmissionSpectrumList ), pointer :: next                  => null()
  end type dustEmissionSpectrumList

  !![
  <dustEmissionSpectrum name="dustEmissionSpectrumSum" docformat="rst">
   <description>
   A dust emission spectrum which is the sum of other spectra, each emitting a fixed fraction of the absorbed luminosity
   from a fixed fraction of the dust mass. Member :math:`i` is given a luminosity :math:`f_{L,i} L_\mathrm{abs}` and a
   mass :math:`f_{M,i} M_\mathrm{dust}`, where :math:`f_{L,i}` are the ``fractionsLuminosity`` and :math:`f_{M,i}` the
   ``fractionsMass``, each of which must sum to one. ``fractionsMass`` defaults to ``fractionsLuminosity``.

   The canonical use is dust at two temperatures---a warm component heated strongly, and a cool one heated weakly---as
   in :cite:t:`da_cunha_simple_2008`. Given fixed temperatures, only the luminosity fractions matter; given energy
   balance, the mass fractions set each component's temperature. A member which is given absorbed luminosity but no
   mass can not find its temperature by energy balance, and reports an error.
   </description>
   <linkedList type="dustEmissionSpectrumList" variable="dustEmissionSpectra" next="next" object="dustEmissionSpectrum_" objectType="dustEmissionSpectrumClass"/>
  </dustEmissionSpectrum>
  !!]
  type, extends(dustEmissionSpectrumClass) :: dustEmissionSpectrumSum
     !!{RST
     A dust emission spectrum which sums other spectra.
     !!}
     private
     type            (dustEmissionSpectrumList), pointer                   :: dustEmissionSpectra => null()
     double precision                          , allocatable, dimension(:) :: fractionsLuminosity          , fractionsMass
   contains
     final     ::               sumDestructor
     procedure :: luminosity => sumLuminosity
  end type dustEmissionSpectrumSum

  interface dustEmissionSpectrumSum
     !!{RST
     Constructors for the :galacticus-class:`dustEmissionSpectrumSum` dust emission spectrum class.
     !!}
     module procedure sumConstructorParameters
     module procedure sumConstructorInternal
  end interface dustEmissionSpectrumSum

contains

  function sumConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustEmissionSpectrumSum` dust emission spectrum class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (dustEmissionSpectrumSum )                              :: self
    type            (inputParameters         ), intent(inout)               :: parameters
    type            (dustEmissionSpectrumList), pointer                     :: dustEmissionSpectrum_
    double precision                          , allocatable  , dimension(:) :: fractionsLuminosity  , fractionsMass
    integer                                                                 :: i

    self                 %dustEmissionSpectra => null()
    dustEmissionSpectrum_                     => null()
    do i=1,parameters%copiesCount('dustEmissionSpectrum',zeroIfNotPresent=.true.)
       if (associated(dustEmissionSpectrum_)) then
          allocate(dustEmissionSpectrum_%next)
          dustEmissionSpectrum_ => dustEmissionSpectrum_%next
       else
          allocate(self%dustEmissionSpectra)
          dustEmissionSpectrum_ => self%dustEmissionSpectra
       end if
       !![
       <objectBuilder class="dustEmissionSpectrum" name="dustEmissionSpectrum_%dustEmissionSpectrum_" source="parameters" copy="i" />
       !!]
    end do
    allocate(fractionsLuminosity(parameters%count('fractionsLuminosity')))
    !![
    <inputParameter docformat="rst">
      <name>fractionsLuminosity</name>
      <description>
      The fraction of the absorbed luminosity emitted by each member spectrum, in the order the members are listed.
      These must sum to one.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    if (parameters%isPresent('fractionsMass')) then
       allocate(fractionsMass(parameters%count('fractionsMass')))
       !![
       <inputParameter docformat="rst">
         <name>fractionsMass</name>
         <description>
         The fraction of the dust mass given to each member spectrum, in the order the members are listed. These must sum
         to one. If absent, the luminosity fractions are used.
         </description>
         <source>parameters</source>
       </inputParameter>
       !!]
    else
       fractionsMass=fractionsLuminosity
    end if
    self%fractionsLuminosity=fractionsLuminosity
    self%fractionsMass      =fractionsMass
    call sumValidate(self)
    !![
    <inputParametersValidate source="parameters" multiParameters="dustEmissionSpectrum"/>
    !!]
    return
  end function sumConstructorParameters

  function sumConstructorInternal(dustEmissionSpectra,fractionsLuminosity,fractionsMass) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustEmissionSpectrumSum` dust emission spectrum class.
    !!}
    implicit none
    type            (dustEmissionSpectrumSum )                                        :: self
    type            (dustEmissionSpectrumList), intent(in   ), target                 :: dustEmissionSpectra
    double precision                          , intent(in   ), dimension(:)           :: fractionsLuminosity
    double precision                          , intent(in   ), dimension(:), optional :: fractionsMass
    type            (dustEmissionSpectrumList), pointer                               :: dustEmissionSpectrum_

    self                 %dustEmissionSpectra => dustEmissionSpectra
    dustEmissionSpectrum_                     => dustEmissionSpectra
    do while (associated(dustEmissionSpectrum_))
       !![
       <referenceCountIncrement owner="dustEmissionSpectrum_" object="dustEmissionSpectrum_"/>
       !!]
       dustEmissionSpectrum_ => dustEmissionSpectrum_%next
    end do
    self%fractionsLuminosity=fractionsLuminosity
    if (present(fractionsMass)) then
       self%fractionsMass=fractionsMass
    else
       self%fractionsMass=fractionsLuminosity
    end if
    call sumValidate(self)
    return
  end function sumConstructorInternal

  subroutine sumValidate(self)
    !!{RST
    Check that there is one luminosity fraction and one mass fraction for each member spectrum, that none is negative,
    and that each set sums to one.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (dustEmissionSpectrumSum ), intent(inout) :: self
    type            (dustEmissionSpectrumList), pointer       :: dustEmissionSpectrum_
    ! Tolerance on the fractions summing to one, loose enough to admit values written to a few significant figures.
    double precision                          , parameter     :: toleranceSum         =1.0d-6
    integer                                                   :: countMembers

    countMembers          =  0
    dustEmissionSpectrum_ => self%dustEmissionSpectra
    do while (associated(dustEmissionSpectrum_))
       countMembers          =  countMembers+1
       dustEmissionSpectrum_ => dustEmissionSpectrum_%next
    end do
    if (countMembers == 0                                                                 ) call Error_Report('no member spectra to sum'                                                   //{introspection:location})
    if (size(self%fractionsLuminosity) /= countMembers .or. size(self%fractionsMass) /= countMembers) call Error_Report('there must be one luminosity fraction and one mass fraction for each member spectrum'//{introspection:location})
    if (any(self%fractionsLuminosity < 0.0d0) .or. any(self%fractionsMass < 0.0d0)       ) call Error_Report('fractions must be non-negative'                                             //{introspection:location})
    if (abs(sum(self%fractionsLuminosity)-1.0d0) > toleranceSum                           ) call Error_Report('`fractionsLuminosity` must sum to one'                                      //{introspection:location})
    if (abs(sum(self%fractionsMass      )-1.0d0) > toleranceSum                           ) call Error_Report('`fractionsMass` must sum to one'                                            //{introspection:location})
    return
  end subroutine sumValidate

  subroutine sumDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`dustEmissionSpectrumSum` dust emission spectrum class.
    !!}
    implicit none
    type(dustEmissionSpectrumSum ), intent(inout) :: self
    type(dustEmissionSpectrumList), pointer       :: dustEmissionSpectrum_, dustEmissionSpectrumNext

    if (associated(self%dustEmissionSpectra)) then
       dustEmissionSpectrum_ => self%dustEmissionSpectra
       do while (associated(dustEmissionSpectrum_))
          dustEmissionSpectrumNext => dustEmissionSpectrum_%next
          !![
          <objectDestructor name="dustEmissionSpectrum_%dustEmissionSpectrum_"/>
          !!]
          deallocate(dustEmissionSpectrum_)
          dustEmissionSpectrum_ => dustEmissionSpectrumNext
       end do
    end if
    return
  end subroutine sumDestructor

  function sumLuminosity(self,wavelengths,luminosityAbsorbed,massDust,time) result(luminosity)
    !!{RST
    Return the sum of the luminosities of the member spectra, each given its fractions of the absorbed luminosity and of
    the dust mass.
    !!}
    implicit none
    class           (dustEmissionSpectrumSum ), intent(inout)                  :: self
    double precision                          , intent(in   ), dimension(:   ) :: wavelengths
    double precision                          , intent(in   )                  :: luminosityAbsorbed, massDust, &
         &                                                                        time
    double precision                          , dimension(size(wavelengths))   :: luminosity
    type            (dustEmissionSpectrumList), pointer                        :: dustEmissionSpectrum_
    integer                                                                    :: i

    luminosity            =  0.0d0
    i                     =  0
    dustEmissionSpectrum_ => self%dustEmissionSpectra
    do while (associated(dustEmissionSpectrum_))
       i                     =  i+1
       luminosity            =  +luminosity                                                                                             &
            &                   +dustEmissionSpectrum_%dustEmissionSpectrum_%luminosity(                                                &
            &                                                                           wavelengths                                   , &
            &                                                                           self%fractionsLuminosity(i)*luminosityAbsorbed, &
            &                                                                           self%fractionsMass      (i)*massDust          , &
            &                                                                           time                                            &
            &                                                                          )
       dustEmissionSpectrum_ => dustEmissionSpectrum_%next
    end do
    return
  end function sumLuminosity
