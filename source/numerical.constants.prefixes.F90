!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a module of useful numerical prefixes.
!!}

module Numerical_Constants_Prefixes
  !!{
  Contains useful numerical prefixes.
  !!}
  implicit none
  public

  double precision       , parameter                    :: quetta  =1.0d+30
  double precision       , parameter                    :: ronna   =1.0d+27
  double precision       , parameter                    :: yotta   =1.0d+24
  double precision       , parameter                    :: zetta   =1.0d+21
  double precision       , parameter                    :: exa     =1.0d+18
  double precision       , parameter                    :: peta    =1.0d+15
  double precision       , parameter                    :: tera    =1.0d+12
  double precision       , parameter                    :: giga    =1.0d+09
  double precision       , parameter                    :: mega    =1.0d+06
  double precision       , parameter                    :: kilo    =1.0d+03
  double precision       , parameter                    :: hecto   =1.0d+02
  double precision       , parameter                    :: deca    =1.0d+01
  double precision       , parameter                    :: deci    =1.0d-01
  double precision       , parameter                    :: centi   =1.0d-02
  double precision       , parameter                    :: milli   =1.0d-03
  double precision       , parameter                    :: micro   =1.0d-06
  double precision       , parameter                    :: nano    =1.0d-09
  double precision       , parameter                    :: pico    =1.0d-12
  double precision       , parameter                    :: femto   =1.0d-15
  double precision       , parameter                    :: atto    =1.0d-18
  double precision       , parameter                    :: zepto   =1.0d-21
  double precision       , parameter                    :: yocto   =1.0d-24
  double precision       , parameter                    :: ronto   =1.0d-27
  double precision       , parameter                    :: quecto  =1.0d-30
  character       (len=2), parameter, dimension(-10:10) :: siPrefix=['q ','r ','y ','z ','a ','f ','p ','n ','μ','m ',' ','k ','M ','G ','T ','P ','E ','Z ','Y ','R ','Q ']

contains

  function siFormat(value_,format_)
    !!{
    Format a value using SI prefixes.
    !!}   
    implicit none
    character       (len=64            )                :: siFormat
    double precision                    , intent(in   ) :: value_
    character       (len=*             ), intent(in   ) :: format_
    integer                                             :: iPrefix
    character       (len=len(format_)+5)                :: formatExtended
    
    if (value_ == 0.0d0) then
       iPrefix=0
    else
       iPrefix=min(max(floor(dble(int(log10(abs(value_)))/3.0d0)),lbound(siPrefix,dim=1)),ubound(siPrefix,dim=1))
    end if
    formatExtended="("//trim(format_)//",a2)"
    write (siFormat,formatExtended) value_/10.0d0**(3*iPrefix),siPrefix(iPrefix)
    return
  end function siFormat
  
end module Numerical_Constants_Prefixes
