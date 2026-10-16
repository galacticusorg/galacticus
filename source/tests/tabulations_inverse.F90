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
Contains a program to test tabulations of monotonic functions built for inversion.
!!}

program Test_Tabulations_Inverse
  !!{RST
  Tests tabulations of monotonic functions built for inversion.

  The property these tabulations exist to provide is that they do *not* depend on the sequence of values they happen to be
  asked for. That is asserted directly here, by building two tabulations of the same function from the same set of requests
  presented in opposite orders and requiring that they agree exactly - which is the check that exposed an order-dependent
  rebuild trigger when the transfer function tabulations were converted, a defect which pinning the abscissae alone did not
  cure.
  !!}
  use :: Display            , only : displayVerbositySet, verbosityLevelStandard
  use :: Tabulations_Inverse, only : tabulationInverse
  use :: Unit_Tests         , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (tabulationInverse)                              :: increasingAscending, increasingDescending, &
       &                                                              decreasingAscending, decreasingDescending, &
       &                                                              extended           , restored
  integer                            , parameter                   :: pointsPerOctave      =12
  double precision                   , dimension(7)                :: resultAscending    , resultDescending
  double precision                   , allocatable, dimension(:)   :: valuesPreserved    , abscissaePreserved  , &
       &                                                              abscissaeExtended
  integer                                                          :: i                  , fileUnit            , &
       &                                                              countPreserved
  logical                                                          :: preserved
  ! Values at which each tabulation is queried. They span several octaves in each direction, so that both must be extended
  ! repeatedly, and they are deliberately not points of the lattice.
  double precision                   , dimension(7), parameter     :: valuesIncreasing     =[3.7d-3,4.1d-2,0.83d0,7.3d0,51.0d0,410.0d0,3300.0d0], &
       &                                                              valuesDecreasing     =[3.1d-3,2.7d-2,0.31d0,2.9d0,44.0d0,370.0d0,2900.0d0]

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Inverse tabulations")

  ! An increasing function, f(x)=x², requested in ascending and in descending order.
  call tabulate(increasingAscending ,valuesIncreasing         ,increasing=.true.)
  call tabulate(increasingDescending,reverse(valuesIncreasing),increasing=.true.)
  do i=1,size(valuesIncreasing)
     resultAscending (i)=increasingAscending %invert(valuesIncreasing(i))
     resultDescending(i)=increasingDescending%invert(valuesIncreasing(i))
  end do
  call Assert("increasing: inverse is independent of the order in which values are requested",resultAscending,resultDescending,absTol=0.0d0)
  call Assert("increasing: the two tabulations span the same lattice"                       , &
       &      [increasingAscending %lattice%indexMinimum,increasingAscending %lattice%count], &
       &      [increasingDescending%lattice%indexMinimum,increasingDescending%lattice%count]  )
  ! The inverse should recover the abscissa: f(x)=x² so x=√f.
  call Assert("increasing: inverse recovers the abscissa at a tabulated point",increasingAscending%invert(4.0d0),2.0d0,relTol=1.0d-12)

  ! A decreasing function, f(x)=1/x, likewise.
  call tabulate(decreasingAscending ,valuesDecreasing              ,increasing=.false.)
  call tabulate(decreasingDescending,reverse(valuesDecreasing)     ,increasing=.false.)
  do i=1,size(valuesDecreasing)
     resultAscending (i)=decreasingAscending %invert(valuesDecreasing(i))
     resultDescending(i)=decreasingDescending%invert(valuesDecreasing(i))
  end do
  call Assert("decreasing: inverse is independent of the order in which values are requested",resultAscending,resultDescending,absTol=0.0d0)
  call Assert("decreasing: the two tabulations span the same lattice"                       , &
       &      [decreasingAscending %lattice%indexMinimum,decreasingAscending %lattice%count], &
       &      [decreasingDescending%lattice%indexMinimum,decreasingDescending%lattice%count]  )
  call Assert("decreasing: inverse recovers the abscissa at a tabulated point",decreasingAscending%invert(0.25d0),4.0d0,relTol=1.0d-12)

  ! Extension must preserve every previously computed value bit-for-bit, and must place them at the abscissae they were
  ! computed for.
  call tabulate(extended,valuesIncreasing(4:4),increasing=.true.)
  countPreserved    =extended%lattice%count
  valuesPreserved   =extended%values
  abscissaePreserved=extended%abscissae()
  call tabulate(extended,valuesIncreasing,increasing=.true.,reset=.false.)
  abscissaeExtended =extended%abscissae()
  preserved         =.false.
  do i=1,size(abscissaeExtended)-countPreserved+1
     if (all(abscissaeExtended(i:i+countPreserved-1) == abscissaePreserved)) then
        preserved=all(extended%values(i:i+countPreserved-1) == valuesPreserved)
        exit
     end if
  end do
  call Assert("extension preserves previously computed values, at their own abscissae",preserved,.true.)

  ! A tabulation must survive a store and restore unchanged - including its ability to invert, which requires that the
  ! values, and not merely the range, be stored.
  open(newUnit=fileUnit,file='tabulationsInverse.state',status='replace',form='unformatted')
  call increasingAscending%stateStore  (fileUnit)
  rewind(fileUnit)
  call restored           %stateRestore(fileUnit)
  close(fileUnit,status='delete')
  do i=1,size(valuesIncreasing)
     resultDescending(i)=restored%invert(valuesIncreasing(i))
  end do
  do i=1,size(valuesIncreasing)
     resultAscending (i)=increasingAscending%invert(valuesIncreasing(i))
  end do
  call Assert("a restored tabulation inverts identically to the one it was stored from",resultAscending,resultDescending,absTol=0.0d0)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

contains

  function reverse(values)
    !!{RST
    Return ``values`` in reverse order.
    !!}
    implicit none
    double precision, intent(in   ), dimension(           :) :: values
    double precision               , dimension(size(values)) :: reverse
    integer                                                  :: i

    do i=1,size(values)
       reverse(i)=values(size(values)-i+1)
    end do
    return
  end function reverse

  subroutine tabulate(tabulation,values,increasing,reset)
    !!{RST
    Build ``tabulation`` by requesting each of ``values`` in turn, in the order given.
    !!}
    implicit none
    type            (tabulationInverse), intent(inout)               :: tabulation
    double precision                   , intent(in   ), dimension(:) :: values
    logical                            , intent(in   )               :: increasing
    logical                            , intent(in   ), optional     :: reset
    double precision                   , allocatable  , dimension(:) :: abscissae
    integer                                                          :: i         , j
    logical                                                          :: reset_

    reset_=.true.
    if (present(reset)) reset_=reset
    if (reset_) call tabulation%reset(pointsPerOctave,increasing)
    do i=1,size(values)
       do while (.not.tabulation%brackets(values(i)))
          call tabulation%expand(values(i))
          abscissae=tabulation%abscissae()
          do j=1,size(abscissae)
             if (tabulation%isComputed(j)) cycle
             if (increasing) then
                call tabulation%set(j,      abscissae(j)**2)
             else
                call tabulation%set(j,1.0d0/abscissae(j)   )
             end if
          end do
          call tabulation%build()
       end do
    end do
    return
  end subroutine tabulate

end program Test_Tabulations_Inverse
