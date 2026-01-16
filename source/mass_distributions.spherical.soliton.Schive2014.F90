!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
  Contains a module which provides utility variables for the \cite{schive_understanding_2014} soliton density profile.
  !!}

  module Mass_Distribution_Soliton_Schive2014
    !!{
    Provides utility variables for the \cite{schive_understanding_2014} soliton density profile.
    !!}
    implicit none
    private
    
    ! Coefficient of the dimensionless radius in the soliton profile.
    double precision, parameter, public :: coefficientCore=0.091d0 ! Schive et al. (2014; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S; equation 3).

  end module Mass_Distribution_Soliton_Schive2014
