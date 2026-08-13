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
Contains a module providing tabulations of monotonic functions, built for inversion.
!!}

module Tabulations_Inverse
  !!{RST
  Provides a tabulation of a monotonic function of a single, positive variable, built so that the function may be
  *inverted*---given a value of the function, the abscissa at which it takes that value is returned.

  Such tabulations occur throughout the analytic mass distributions, which invert the enclosed mass, the mean enclosed
  density, the specific angular momentum, and the freefall time to recover a radius. Each was built by bracketing the
  requested value---halving the lower bound, or doubling the upper, until the tabulated range spanned it---and then
  redistributing a fixed number of points per decade across the resulting range. That redistribution moves every abscissa
  whenever the range grows, so the tabulation depended on the sequence of values it happened to be asked for, and every
  previously computed point had to be discarded and recomputed.

  Here the abscissae are instead points of an absolute lattice, ``pointsPerOctave`` to the octave. The bounds reached by
  halving and doubling are themselves exact powers of two, so they are already points of that lattice: pinning to it costs
  nothing in range, and the abscissae become a function of which octaves have been visited rather than of the order in which
  values were requested. Extending the tabulation then preserves every point already computed, and only the newly exposed
  points are evaluated.

  The function is evaluated by the caller rather than through a procedure pointer, which keeps each evaluator where it is
  written---with access to whatever host state it needs---and avoids imposing a common interface on evaluators which range
  from a single ``elemental`` expression to a numerical quadrature per point. The pattern of use is:

  .. code-block:: fortran

     do while (.not.tabulation%brackets(value))
        call tabulation%expand(value)
        radii=tabulation%abscissae()
        do i=1,size(radii)
           if (tabulation%isComputed(i)) cycle
           call tabulation%set(i,f(radii(i)))
        end do
        call tabulation%build()
     end do
     abscissa=tabulation%invert(value)
  !!}
  use :: Numerical_Interpolation, only : interpolator
  use :: Numerical_Ranges       , only : rangeLattice

  private

  type, public :: tabulationInverse
     !!{RST
     A tabulation of a monotonic function, built for inversion.
     !!}
     type            (rangeLattice)                            :: lattice
     double precision              , allocatable, dimension(:) :: values
     logical                       , allocatable, dimension(:) :: isComputed
     type            (interpolator), allocatable               :: interpolator_
     double precision                                          :: valueMinimum   =+huge(0.0d0), &
          &                                                       valueMaximum   =-huge(0.0d0)
     integer                                                   :: pointsPerOctave=0
     logical                                                   :: increasing     =.true.      , &
          &                                                       initialized    =.false.
   contains
     !![
     <methods docformat="rst">
       <method description="Initialize the tabulation, discarding anything previously tabulated."                               method="reset"       />
       <method description="Return true if the tabulated range of the function spans the given value."                          method="brackets"    />
       <method description="Extend the tabulation by one octave, in whichever direction is needed to approach the given value." method="expand"      />
       <method description="Return the abscissae of the tabulation."                                                            method="abscissae"   />
       <method description="Set the value of the function at the given point of the tabulation."                                method="set"         />
       <method description="Complete the tabulation, recording its range and constructing the interpolator used to invert it."  method="build"       />
       <method description="Return the abscissa at which the function takes the given value."                                   method="invert"      />
       <method description="Return the derivative of the abscissa with respect to the value of the function."                   method="derivative"  />
       <method description="Write the tabulation to an open state file."                                                        method="stateStore"  />
       <method description="Read the tabulation back from an open state file."                                                  method="stateRestore"/>
     </methods>
     !!]
     procedure :: reset        => tabulationInverseReset
     procedure :: brackets     => tabulationInverseBrackets
     procedure :: expand       => tabulationInverseExpand
     procedure :: abscissae    => tabulationInverseAbscissae
     procedure :: set          => tabulationInverseSet
     procedure :: build        => tabulationInverseBuild
     procedure :: invert       => tabulationInverseInvert
     procedure :: derivative   => tabulationInverseDerivative
     procedure :: stateStore   => tabulationInverseStateStore
     procedure :: stateRestore => tabulationInverseStateRestore
  end type tabulationInverse

contains

  subroutine tabulationInverseReset(self,pointsPerOctave,increasing)
    !!{RST
    Initialize the tabulation, discarding anything previously tabulated. ``increasing`` states whether the function increases
    with its abscissa, which fixes the direction in which the tabulation must be extended to reach a given value.
    !!}
    implicit none
    class  (tabulationInverse), intent(inout) :: self
    integer                   , intent(in   ) :: pointsPerOctave
    logical                   , intent(in   ) :: increasing

    if (allocated(self%values       )) deallocate(self%values       )
    if (allocated(self%isComputed   )) deallocate(self%isComputed   )
    if (allocated(self%interpolator_)) deallocate(self%interpolator_)
    self%lattice        =rangeLattice()
    self%pointsPerOctave=pointsPerOctave
    self%increasing     =increasing
    self%valueMinimum   =+huge(0.0d0)
    self%valueMaximum   =-huge(0.0d0)
    self%initialized    =.false.
    return
  end subroutine tabulationInverseReset

  logical function tabulationInverseBrackets(self,value) result(brackets)
    !!{RST
    Return true if the tabulated range of the function spans ``value``, so that it can be inverted without extending the
    tabulation.
    !!}
    implicit none
    class           (tabulationInverse), intent(in   ) :: self
    double precision                   , intent(in   ) :: value

    brackets=               self%initialized  &
         &   .and. value >= self%valueMinimum &
         &   .and. value <= self%valueMaximum
    return
  end function tabulationInverseBrackets

  subroutine tabulationInverseExpand(self,value)
    !!{RST
    Extend the tabulation by one octave, in whichever direction brings the tabulated range closer to spanning ``value``,
    preserving every point already computed. An unbuilt tabulation is seeded with the single octave about unity, which is
    where the bracketing loops these tabulations replace also began.
    !!}
    use :: Error           , only : Error_Report
    use :: Numerical_Ranges, only : Range_Lattice_Extend, gridSchemePerOctave
    implicit none
    class           (tabulationInverse), intent(inout) :: self
    double precision                   , intent(in   ) :: value
    type            (rangeLattice     )                :: latticeNew
    logical                                            :: extendUpward

    if (self%pointsPerOctave <= 0) call Error_Report('tabulation has not been initialized'//{introspection:location})
    if (.not.self%lattice%isDefined()) then
       ! Seed with the octave spanning [1/2,2].
       latticeNew=rangeLattice(gridSchemePerOctave,self%pointsPerOctave,-self%pointsPerOctave,2*self%pointsPerOctave+1)
    else
       ! Extend towards the requested value. A value above the tabulated maximum is reached by extending to larger abscissae
       ! for an increasing function, and to smaller for a decreasing one.
       extendUpward=(value > self%valueMaximum) .eqv. self%increasing
       if (extendUpward) then
          latticeNew=rangeLattice(gridSchemePerOctave,self%pointsPerOctave,self%lattice%indexMinimum                     ,self%lattice%count+self%pointsPerOctave)
       else
          latticeNew=rangeLattice(gridSchemePerOctave,self%pointsPerOctave,self%lattice%indexMinimum-self%pointsPerOctave,self%lattice%count+self%pointsPerOctave)
       end if
    end if
    call Range_Lattice_Extend(self%lattice,latticeNew,self%values,self%isComputed)
    self%lattice=latticeNew
    return
  end subroutine tabulationInverseExpand

  function tabulationInverseAbscissae(self) result(abscissae)
    !!{RST
    Return the abscissae of the tabulation.
    !!}
    implicit none
    class           (tabulationInverse), intent(in   )               :: self
    double precision                   , allocatable  , dimension(:) :: abscissae

    abscissae=self%lattice%values()
    return
  end function tabulationInverseAbscissae

  subroutine tabulationInverseSet(self,index,value)
    !!{RST
    Set the value of the function at the given point of the tabulation.
    !!}
    implicit none
    class           (tabulationInverse), intent(inout) :: self
    integer                            , intent(in   ) :: index
    double precision                   , intent(in   ) :: value

    self%values    (index)=value
    self%isComputed(index)=.true.
    return
  end subroutine tabulationInverseSet

  subroutine tabulationInverseBuild(self)
    !!{RST
    Complete the tabulation, recording the range of the function over it and constructing the interpolator through which it
    is inverted. The interpolator requires ascending abscissae, so for a decreasing function the tabulated values are negated
    ---as is the value presented to it when the function is inverted---rather than reversing the arrays.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(tabulationInverse), intent(inout) :: self

    if (.not.self%lattice%isDefined()) call Error_Report('tabulation has no abscissae'//{introspection:location})
    if (.not.all(self%isComputed)    ) call Error_Report('tabulation is incomplete'   //{introspection:location})
    if (allocated(self%interpolator_)) deallocate(self%interpolator_)
    allocate(self%interpolator_)
    if (self%increasing) then
       self%valueMinimum =self%values(                 1)
       self%valueMaximum =self%values(self%lattice%count)
       self%interpolator_=interpolator(+self%values,self%lattice%values())
    else
       self%valueMinimum =self%values(self%lattice%count)
       self%valueMaximum =self%values(                 1)
       self%interpolator_=interpolator(-self%values,self%lattice%values())
    end if
    self%initialized=.true.
    return
  end subroutine tabulationInverseBuild

  double precision function tabulationInverseInvert(self,value) result(abscissa)
    !!{RST
    Return the abscissa at which the function takes the given value.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (tabulationInverse), intent(inout) :: self
    double precision                   , intent(in   ) :: value

    if (.not.self%initialized) call Error_Report('tabulation has not been built'//{introspection:location})
    if (self%increasing) then
       abscissa=self%interpolator_%interpolate(+value)
    else
       abscissa=self%interpolator_%interpolate(-value)
    end if
    return
  end function tabulationInverseInvert

  double precision function tabulationInverseDerivative(self,value) result(derivative)
    !!{RST
    Return the derivative of the abscissa with respect to the value of the function. For a decreasing function the tabulation
    is held against the negated value, so the derivative with respect to the value itself carries the corresponding sign.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (tabulationInverse), intent(inout) :: self
    double precision                   , intent(in   ) :: value

    if (.not.self%initialized) call Error_Report('tabulation has not been built'//{introspection:location})
    if (self%increasing) then
       derivative=+self%interpolator_%derivative(+value)
    else
       derivative=-self%interpolator_%derivative(-value)
    end if
    return
  end function tabulationInverseDerivative

  subroutine tabulationInverseStateStore(self,stateFile)
    !!{RST
    Write the tabulation to an open state file. The tabulated values are written along with the lattice on which they were
    computed, so that the tabulation is restored in full. Note that recording the range alone would not be sufficient: it
    would leave a restored object which reports that it spans a value while holding nothing with which to invert it.
    !!}
    implicit none
    class  (tabulationInverse), intent(in   ) :: self
    integer                   , intent(in   ) :: stateFile

    write (stateFile) self%pointsPerOctave,self%increasing,self%initialized
    write (stateFile) self%lattice%scheme%ID,self%lattice%pointsPer,self%lattice%indexMinimum,self%lattice%count
    write (stateFile) allocated(self%values)
    if (allocated(self%values)) write (stateFile) self%values
    return
  end subroutine tabulationInverseStateStore

  subroutine tabulationInverseStateRestore(self,stateFile)
    !!{RST
    Read the tabulation back from an open state file, and reconstruct the interpolator through which it is inverted.
    !!}
    use :: Numerical_Ranges, only : enumerationGridSchemeType
    implicit none
    class  (tabulationInverse), intent(inout) :: self
    integer                   , intent(in   ) :: stateFile
    integer                                   :: schemeID    , pointsPer, &
         &                                       indexMinimum, count
    logical                                   :: hasValues

    if (allocated(self%values       )) deallocate(self%values       )
    if (allocated(self%isComputed   )) deallocate(self%isComputed   )
    if (allocated(self%interpolator_)) deallocate(self%interpolator_)
    read (stateFile) self%pointsPerOctave,self%increasing,self%initialized
    read (stateFile) schemeID,pointsPer,indexMinimum,count
    self%lattice=rangeLattice(enumerationGridSchemeType(schemeID),pointsPer,indexMinimum,count)
    read (stateFile) hasValues
    if (hasValues) then
       allocate(self%values    (count))
       allocate(self%isComputed(count))
       read (stateFile) self%values
       self%isComputed=.true.
    end if
    ! Rebuild the interpolator from the restored values, rather than attempting to serialize it.
    if (self%initialized .and. hasValues) then
       self%initialized=.false.
       call self%build()
    end if
    return
  end subroutine tabulationInverseStateRestore

end module Tabulations_Inverse
