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

  !% Implements calculations of incompleteness assuming a complete sample.
  
  !# <massFunctionIncompleteness name="massFunctionIncompletenessComplete">
  !#  <description>Computes incompleteness for a complete survey.</description>
  !# </massFunctionIncompleteness>

  type, extends(massFunctionIncompletenessClass) :: massFunctionIncompletenessComplete
     !% A class implementing incompleteness calculations for a complete survey.
     private
   contains
     final     ::                 completeDestructor
     procedure :: completeness => completeCompleteness
  end type massFunctionIncompletenessComplete

  interface massFunctionIncompletenessComplete
     !% Constructors for the ``complete'' incompleteness class.
     module procedure completeDefaultConstructor
  end interface massFunctionIncompletenessComplete

contains

  function completeDefaultConstructor()
    !% Default constructor for the ``complete'' incompleteness class.
    implicit none
    type(massFunctionIncompletenessComplete) :: completeDefaultConstructor

   return
  end function completeDefaultConstructor

  subroutine completeDestructor(self)
    !% Destructor for the ``complete'' incompleteness class.
     use Gaussian_Random
     implicit none
     type(massFunctionIncompletenessComplete), intent(inout) :: self
    return
  end subroutine completeDestructor

  double precision function completeCompleteness(self,mass)
    !% Return the completeness.
    implicit none
    class           (massFunctionIncompletenessComplete), intent(inout) :: self
    double precision                                    , intent(in   ) :: mass

    completeCompleteness=1.0d0
    return
  end function completeCompleteness
