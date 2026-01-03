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

!!{
Contains a module which provides poly-ranked types (i.e. types which can store data in arrays of different ranks).
!!}

module Poly_Ranks
  !!{
  Provides poly-ranked types (i.e. types which can store data in arrays of different ranks).
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_size_t
  use            :: Kind_Numbers , only : kind_int8
  implicit none
  public

  !![
  <generic identifier="Type">
   <instance label="integer" intrinsic="integer(kind_int8)"/>
   <instance label="double"  intrinsic="double precision"/>
  </generic>
  !!]

  type :: polyRank{Type¦label}
     !!{
     A type which provides poly-ranked double precision data.
     !!}
     integer         (c_size_t), allocatable, dimension(:) :: shape_
     {Type¦intrinsic}          , allocatable, dimension(:) :: data
   contains
     !![
     <methods>
       <method description="Return the rank of the data."  method="rank"  />
       <method description="Return the shape of the data." method="shape" />
     </methods>
     !!]
     procedure :: rank   => {Type¦label}Rank
     procedure :: shape  => {Type¦label}Shape
  end type polyRank{Type¦label}

  interface polyRank{Type¦label}
     module procedure {Type¦label}Constructor
  end interface polyRank{Type¦label}

  interface assignment(=)
     module procedure {Type¦label}Assign
  end interface assignment(=)

contains

  function {Type¦label}Constructor(array) result(self)
    !!{
    Constructor for poly-ranked arrays.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (polyRank{Type¦label})                               :: self
    {Type¦intrinsic}                      , dimension(..), intent(in   ) :: array

    allocate   (self%shape_(rank (array)))
    if (rank(array) > 0) then
       allocate(self%data  (size (array)))
       self%shape_=         shape(array)
    else
       allocate(self%data(            1 ))
    end if
    select rank (array)
    rank (0)
       self%data(1)=        array
    rank (1)
       self%data   =reshape(array,shape(self%data))
    rank (2)
       self%data   =reshape(array,shape(self%data))
    rank (3)
       self%data   =reshape(array,shape(self%data))
    rank (4)
       self%data   =reshape(array,shape(self%data))
    rank (5)
       self%data   =reshape(array,shape(self%data))
    rank (6)
       self%data   =reshape(array,shape(self%data))
    rank (7)
       self%data   =reshape(array,shape(self%data))
    rank default
       call Error_Report('unsupported rank'//{introspection:location})
    end select
    return
  end function {Type¦label}Constructor
  
  subroutine {Type¦label}Assign(array,self)
    !!{
    Assign to an array.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (polyRank{Type¦label}), intent(in   )                :: self
    {Type¦intrinsic}                      , dimension(..), intent(  out) :: array
    integer         (c_size_t            ), dimension( 1)                :: shape1
    integer         (c_size_t            ), dimension( 2)                :: shape2
    integer         (c_size_t            ), dimension( 3)                :: shape3
    integer         (c_size_t            ), dimension( 4)                :: shape4
    integer         (c_size_t            ), dimension( 5)                :: shape5
    integer         (c_size_t            ), dimension( 6)                :: shape6
    integer         (c_size_t            ), dimension( 7)                :: shape7

    if (rank(array) /= size(self%shape_)) call Error_Report('rank mismatch'//{introspection:location})
    select rank (array)
    rank (0)
       array=        self%data(1)
    rank (1)
       shape1=self%shape_
       array=reshape(self%data   ,shape1)
    rank (2)
       shape2=self%shape_
       array=reshape(self%data   ,shape2)
    rank (3)
       shape3=self%shape_
       array=reshape(self%data   ,shape3)
    rank (4)
       shape4=self%shape_
       array=reshape(self%data   ,shape4)
    rank (5)
       shape5=self%shape_
       array=reshape(self%data   ,shape5)
    rank (6)
       shape6=self%shape_
       array=reshape(self%data   ,shape6)
    rank (7)
       shape7=self%shape_
       array=reshape(self%data   ,shape7)
    rank default
       call Error_Report('unsupported rank'//{introspection:location})
    end select
    return
  end subroutine {Type¦label}Assign
  
  integer function {Type¦label}Rank(self)
    !!{
    Return the rank of the object.
    !!}
    implicit none
    class(polyRank{Type¦label}), intent(in   ) :: self

    {Type¦label}Rank=size(self%shape_)
    return
  end function {Type¦label}Rank
  
  function {Type¦label}Shape(self) result(shape_)
    !!{
    Return the shape of the object.
    !!}
    implicit none
    integer(c_size_t            ), allocatable  , dimension(:) :: shape_
    class  (polyRank{Type¦label}), intent(in   )               :: self

    allocate(shape_(size(self%shape_)))
    shape_=self%shape_
    return
  end function {Type¦label}Shape
  
end module Poly_Ranks
