!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which imports the FFTW3 library Fortran interface.

module FFTW3
  !% Imports the FFTW3 library Fortran interface.
  use, intrinsic :: ISO_C_Binding
  public
  include 'fftw3.f03'

contains

  double precision function FFTW_Wavenumber(k,n)
    !% Return the wavenumber (in units of $1/L$ where $L$ is the box length) corresponding to element {\tt k} out of {\tt n} of a
    !% 1-D FFT using the FFTW convention.
    implicit none
    integer, intent(in   ) :: k , n
    integer                :: kk

    kk=k-1
    if (kk < n/2) then
       FFTW_Wavenumber=dble(kk  )
    else
       FFTW_Wavenumber=dble(kk-n)
    end if
    return
  end function FFTW_Wavenumber

end module FFTW3
