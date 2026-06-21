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
  Implements a stellar initial mass function class for piecewise power-law :term:`IMF`\ s.
  !!}

  !![
  <initialMassFunction name="initialMassFunctionPiecewisePowerLaw" docformat="rst">
   <description>
   A stellar initial mass function class for piecewise power-law :term:`IMF`\ s. Arbitrary piecewise power-law :term:`IMF`\ s can be defined using the ``PiecewisePowerLaw`` method. The :term:`IMF` will be constructed such that:

   .. math::

      \phi(M) \propto M^{\alpha_i} \hbox{ if } M_i \le M &lt; M_{i+1},

   where :math:`i=1`\ …\ :math:`N`, the :math:`M_i` are given by ``[mass]`` and the :math:`\alpha_i` are given by ``[exponent]``. (Note that ``[mass]`` must contain :math:`N+1` elements, while ``[exponent]`` contains only :math:`N` elements.) The normalization of each power-law piece is chosen to ensure a continuous :term:`IMF` that is normalized to unit mass overall.
   </description>
  </initialMassFunction>
  !!]
  type, extends(initialMassFunctionClass) :: initialMassFunctionPiecewisePowerLaw
     !!{RST
     A stellar initial mass function class for piecewise power-law :term:`IMF`\ s.
     !!}
     private
     integer                                     :: countPieces
     double precision, allocatable, dimension(:) :: normalization, exponent , &
          &                                         massLower    , massUpper, &
          &                                         mass
   contains
     procedure :: massMinimum      => piecewisePowerLawMassMinimum
     procedure :: massMaximum      => piecewisePowerLawMassMaximum
     procedure :: phi              => piecewisePowerLawPhi
     procedure :: numberCumulative => piecewisePowerLawNumberCumulative
     procedure :: tabulate         => piecewisePowerLawTabulate
     procedure :: label            => piecewisePowerLawLabel
  end type initialMassFunctionPiecewisePowerLaw

  interface initialMassFunctionPiecewisePowerLaw
     !!{RST
     Constructors for the :galacticus-class:`initialMassFunctionPiecewisePowerLaw` initial mass function class.
     !!}
     module procedure piecewisePowerLawConstructorParameters
     module procedure piecewisePowerLawConstructorInternal
  end interface initialMassFunctionPiecewisePowerLaw

contains

  function piecewisePowerLawConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`initialMassFunctionPiecewisePowerLaw` initial mass function class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (initialMassFunctionPiecewisePowerLaw)                            :: self
    type            (inputParameters                     ), intent(inout)             :: parameters
    double precision                                      , allocatable, dimension(:) :: mass      , exponent
    integer                                                                           :: countMass , countExponent

    countMass    =parameters%count('mass'    ,zeroIfNotPresent=.true.)
    countExponent=parameters%count('exponent',zeroIfNotPresent=.true.)
    allocate(mass    (countMass    ))
    allocate(exponent(countExponent))
    !![
    <inputParameter docformat="rst">
      <name>mass</name>
      <defaultValue>[0.1d0,125.0d0]</defaultValue>
      <source>parameters</source>
      <description>
      The mass points used to define a piecewise power-law initial mass function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponent</name>
      <defaultValue>[-2.35d0]</defaultValue>
      <source>parameters</source>
      <description>
      The exponents used to define a piecewise power-law initial mass function.
      </description>
    </inputParameter>
    !!]
    self=initialMassFunctionPiecewisePowerLaw(mass,exponent)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function piecewisePowerLawConstructorParameters

  function piecewisePowerLawConstructorInternal(mass,exponent) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`initialMassFunctionPiecewisePowerLaw` initial mass function.
    !!}
    use :: Array_Utilities, only : Array_Is_Monotonic, directionIncreasing
    use :: Error          , only : Error_Report
    implicit none
    type            (initialMassFunctionPiecewisePowerLaw)                              :: self
    double precision                                      , intent(in   ), dimension(:) :: mass     , exponent
    integer                                                                             :: i
    double precision                                                                    :: massTotal

    ! Validate.
    if (size(mass) <                 2                             ) call Error_Report('at least 2 mass points are required to define the IMF'//{introspection:location})
    if (size(mass) /= size(exponent)+1                             ) call Error_Report('number of exponents must match number of intervals'   //{introspection:location})
    if (.not.Array_Is_Monotonic(mass,direction=directionIncreasing)) call Error_Report('masses must be monotonically increasing'              //{introspection:location})
    ! Set upper and lower masses.
    allocate(self%massLower    (size(exponent)  ))
    allocate(self%massUpper    (size(exponent)  ))
    allocate(self%mass         (size(exponent)+1))
    allocate(self%exponent     (size(exponent)  ))
    allocate(self%normalization(size(exponent)  ))
    self%massLower  =mass    (1:size(exponent)  )
    self%massUpper  =mass    (2:size(exponent)+1)
    self%mass       =mass
    self%exponent   =exponent
    self%countPieces=size(exponent)
    ! Find total mass in the IMF.
    massTotal=0.0d0
    do i=1,size(exponent)
       ! Ensure continuity in the IMF.
       if (i == 1) then
          ! Arbitrary normalization for first piece.
          self%normalization(i)=1.0d0
       else
          ! Normalize to get continuity at the lower boundary.
          self%normalization(i)=+self%normalization(i-1)                &
               &                *self%massUpper    (i-1)**exponent(i-1) &
               &                /self%massLower    (i  )**exponent(i  )
       end if
       ! Sum the mass contributed by this range.
       if (self%exponent(i) == -2.0d0) then
          ! Integral is logarithmic in this case.
          massTotal=+massTotal                  &
               &    +     self%normalization(i) &
               &    *log(                       &
               &         +self%massUpper    (i) &
               &         /self%massLower    (i) &
               &        )
       else
          ! Integral is a power-law in all other cases.
          massTotal=+massTotal                                         &
               &    +  self%normalization(i)                           &
               &    *(                                                 &
               &      +self%massUpper    (i)**(2.0d0+self%exponent(i)) &
               &      -self%massLower    (i)**(2.0d0+self%exponent(i)) &
               &     )                                                 &
               &    /                         (2.0d0+self%exponent(i))
       end if
    end do
    ! Divide through by the total mass to get the correct normalization.
    self%normalization=+self%normalization &
         &             /massTotal
    return
  end function piecewisePowerLawConstructorInternal

  double precision function piecewisePowerLawMassMinimum(self)
    !!{RST
    Return the minimum mass of stars in the :cite:t:`chabrier_galactic_2001` :term:`IMF`.
    !!}
    implicit none
    class(initialMassFunctionPiecewisePowerLaw), intent(inout) :: self

    piecewisePowerLawMassMinimum=self%massLower(1)
    return
  end function piecewisePowerLawMassMinimum

  double precision function piecewisePowerLawMassMaximum(self)
    !!{RST
    Return the maximum mass of stars in a piecewise power-law :term:`IMF`.
    !!}
    implicit none
    class(initialMassFunctionPiecewisePowerLaw), intent(inout) :: self

    piecewisePowerLawMassMaximum=self%massUpper(self%countPieces)
    return
  end function piecewisePowerLawMassMaximum

  double precision function piecewisePowerLawPhi(self,massInitial)
    !!{RST
    Evaluate a piecewise power-law stellar initial mass function.
    !!}
    implicit none
    class           (initialMassFunctionPiecewisePowerLaw), intent(inout) :: self
    double precision                                      , intent(in   ) :: massInitial
    integer                                                               :: i

    piecewisePowerLawPhi=0.0d0
    do i=1,self%countPieces
       if     (                                        &
            &       massInitial >= self%massLower  (i) &
            &  .and.                                   &
            &   (                                      &
            &       massInitial <  self%massUpper  (i) &
            &    .or.                                  &
            &     (                                    &
            &       massInitial == self%massUpper  (i) &
            &      .and.                               &
            &       i           == self%countPieces    &
            &     )                                    &
            &   )                                      &
            & ) then
          piecewisePowerLawPhi=+             self%normalization(i) &
               &               *massInitial**self%exponent     (i)
          exit
       end if
    end do
    return
  end function piecewisePowerLawPhi

  double precision function piecewisePowerLawNumberCumulative(self,massLower,massUpper) result(number)
    !!{RST
    Evaluate a piecewise power-law stellar initial mass function.
    !!}
    implicit none
    class           (initialMassFunctionPiecewisePowerLaw), intent(inout) :: self
    double precision                                      , intent(in   ) :: massLower , massUpper
    double precision                                                      :: massLower_, massUpper_
    integer                                                               :: i

    number=0.0d0
    do i=1,self%countPieces
       massLower_=max(massLower,self%massLower(i))
       massUpper_=min(massUpper,self%massUpper(i))
       if (massUpper_ > massLower_)                               &
            & number=+     number                                 &
            &        +                     self%normalization(i)  &
            &        /              (1.0d0+self%exponent     (i)) &
            &        *(                                           &
            &          +massUpper_**(1.0d0+self%exponent     (i)) &
            &          -massLower_**(1.0d0+self%exponent     (i)) &
            &        )
    end do
    return
  end function piecewisePowerLawNumberCumulative

  subroutine piecewisePowerLawTabulate(self,imfTable)
    !!{RST
    Construct and return a tabulation of a piecewise power-law :term:`IMF`.
    !!}
    use :: Tables, only : table1DLogarithmicLinear
    implicit none
    class  (initialMassFunctionPiecewisePowerLaw)             , intent(inout) :: self
    class  (table1D                             ), allocatable, intent(inout) :: imfTable
    integer                                      , parameter                  :: countPoints=100
    integer                                                                   :: i

    allocate(table1DLogarithmicLinear :: imfTable)
    select type (imfTable)
    type is (table1DLogarithmicLinear)
       call imfTable%create(                                  &
            &               self%massLower(1               ), &
            &               self%massUpper(self%countPieces), &
            &               countPoints                       &
            &           )
       do i=1,countPoints
          call imfTable%populate(self%phi(imfTable%x(i)),i)
       end do
    end select
    return
  end subroutine piecewisePowerLawTabulate

  function piecewisePowerLawLabel(self)
    !!{RST
    Return a label for this :term:`IMF`.
    !!}
    implicit none
    class(initialMassFunctionPiecewisePowerLaw), intent(inout) :: self
    type (varying_string                      )                :: piecewisePowerLawLabel
    !$GLC attributes unused :: self

    piecewisePowerLawLabel="PiecewisePowerLaw"
    return
  end function piecewisePowerLawLabel
